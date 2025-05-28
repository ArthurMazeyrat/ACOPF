import gurobipy as gp
from gurobipy import GRB
from gurobi_optimods import opf
from gurobi_optimods import datasets
import numpy as np
import conversion
from tqdm import tqdm
import time
import itertools
import math
import matplotlib.pyplot as plt
import networkx as nx
import hyperplanes
import SCIP
import ObjBound
import Theta

class Node:
    def __init__(self,number:int, v_i: int, parent: int, bounds_v: list[float],theta_i:int,bounds_theta:list[float], val: float,next_v:list,next_theta:list,v_or_theta:bool):
        self.number = number
        self.v_i = v_i        
        self.parent = parent
        self.bounds_v = bounds_v
        self.theta_i = theta_i
        self.bounds_theta = bounds_theta
        self.val = val 
        self.next_v = next_v
        self.next_theta = next_theta
        self.v_or_theta = v_or_theta

    def __str__(self):
        parent_str = f"Parent: {self.parent}"
        if self.bounds_v!=None:
            bounds_str_v = f"Bounds_v: [{self.bounds_v[0]}, {self.bounds_v[1]}]"
        else:
            bounds_str_v = f"Bounds_theta: [None]"
            
        if self.bounds_theta!=None:
            bounds_str_theta = f"Bounds_theta: [{self.bounds_theta[0]}, {self.bounds_theta[1]}]"
        else:
            bounds_str_theta = f"Bounds_theta: [None]"
        val_str = f"Value: {self.val}"
        next_str_v = f"Next nodes_v: {self.next_v}"
        next_str_theta = f"Next nodes_theta: {self.next_theta}"
        str_vort = f"V/Theta : {self.v_or_theta}"
        return f"Node(n°={self.number}, v_i={self.v_i}, {parent_str}, {bounds_str_v}, {bounds_str_theta}, {val_str}, {next_str_v}, {next_str_theta}, {str_vort})"


def Solve(case,plot,theta_lim_,node : Node,nodes:dict):
    duration=time.time()

    Bus_Gen={}
    for i in range(len(case["gen"])):
        bus = case["gen"][i]["bus"]
        if Bus_Gen.get(bus)==None:
            Bus_Gen[bus] = [i]
        else:
            if i in Bus_Gen[bus]:
                raise ValueError
            else:
                Bus_Gen[bus].append(i)

    N=len(case['branch'])
    # print(f"Number of branches : {N}")

    Branch_Bus={}
    for i in range(len(case["branch"])):
        fbus_ = case["branch"][i]["fbus"]
        tbus_ = case["branch"][i]["tbus"]
        if fbus_ not in Branch_Bus:
            Branch_Bus[fbus_] = [i] # Branch i starts from fbus_
        else:
            if i not in Branch_Bus[fbus_]:
                Branch_Bus[fbus_].append(i) # Branch i starts from fbus_
            else:
                # print(f"{fbus_,tbus_} double")
                raise ValueError
        if tbus_ not in Branch_Bus:
            Branch_Bus[tbus_] = [i+N] # Branch i+N starts from tbus_
        else:
            if i+N not in Branch_Bus[tbus_]:
                Branch_Bus[tbus_].append(i+N) # Branch i+N starts from tbus_
            else:
                # print(f"{tbus_,fbus_} double")
                raise ValueError

    Bus_SH={} #Dict that makes a bus_id refers to its gs and bs
    for bus in case["bus"]:
        Bus_SH[bus["bus_i"]]=[bus["Gs"],bus["Bs"]]

    def create_admittance_matrix(branch):
        G={}
        B={}
        #for branch in case["branch"]:
        fbus = branch["fbus"]
        tbus = branch["tbus"]
        # if branch==case['branch'][1]:
        #     G[(fbus,fbus)] = 0
        #     G[(fbus,tbus)] =0
        #     G[(tbus,fbus)] = 0
        #     G[(tbus,tbus)] = 0

        #     B[(fbus,fbus)] = 0
        #     B[(fbus,tbus)] = 0
        #     B[(tbus,fbus)] = 0
        #     B[(tbus,tbus)] = 0
        #     return G,B
        r = branch["r"] # real part of impendence
        x = branch["x"] # complex part of impendence
        b = branch["b"] # line charging susceptance
        t = branch["ratio"]
        v = branch["angle"]
        z = r + x*1j # impendence
        y = 1/z
        A = (y + (b/2)*1j) #no trasformer admittance

        if t == 0:
            t = 1 # optimod lo fa quindi io lo faccio

        Y22 = A
        Y11 = A/(t**2)
        Y12 = - y / (t * np.exp(-v*1j))
        Y21 = - y / (t * np.exp(v*1j))

        G[(fbus,fbus)] = Y11.real
        G[(fbus,tbus)] = Y12.real
        G[(tbus,fbus)] = Y21.real
        G[(tbus,tbus)] = Y22.real

        B[(fbus,fbus)] = Y11.imag
        B[(fbus,tbus)] = Y12.imag
        B[(tbus,fbus)] = Y21.imag
        B[(tbus,tbus)] = Y22.imag
        if fbus==tbus:
            raise ValueError
        return G,B

    baseMVA=case["baseMVA"]

    # Create a new model
    model = gp.Model("Jabr-ACOPF")
    model.setParam(GRB.Param.OutputFlag, 0)

    ### Variables definition & intervals

    Pg = []  # Active power generation
    Qg = []  # Reactive power generation

    for i in range(len(case['gen'])):
        Pg.append(model.addVar(lb=case['gen'][i]['Pmin'], ub=case['gen'][i]['Pmax'], name=f"Pg{i}"))
        Qg.append(model.addVar(lb=case['gen'][i]['Qmin'], ub=case['gen'][i]['Qmax'], name=f"Qg{i}"))

    # print("Power Gen var generated")

    v_lim={}
    for i in range(len(case["bus"])):
        bus_i=case["bus"][i]["bus_i"]
        v_lim[bus_i]=[case["bus"][i]["Vmin"],case["bus"][i]["Vmax"]]
    current_node=node
    while current_node.parent!=-1:
        if current_node.bounds_v!=None:
            if current_node.bounds_v[1]==0:
                if current_node.bounds_v[0]>v_lim[current_node.v_i][0]:
                    v_lim[current_node.v_i][0]=current_node.bounds_v[0]
                    # print(f"V {current_node.v_i} : {v_lim[current_node.v_i]}")
            else:
                if current_node.bounds_v[0]<v_lim[current_node.v_i][1]:
                    v_lim[current_node.v_i][1]=current_node.bounds_v[0]
                    # print(f"V {current_node.v_i} : {v_lim[current_node.v_i]}")
        current_node=nodes[current_node.parent]

    theta={}
    theta_lim={}
    for i in range(len(case["branch"])):
        theta_lim[i]=[-theta_lim_,theta_lim_]
    current_node=node
    while current_node.parent!=-1:
        if current_node.bounds_theta!=None:
            if current_node.bounds_theta[1]==0:
                if current_node.bounds_theta[0]>theta_lim[current_node.theta_i][0]:
                    theta_lim[current_node.theta_i][0]=current_node.bounds_theta[0]
                    # print(f"T {current_node.theta_i} : {theta_lim[current_node.theta_i]}")
            else:
                if current_node.bounds_theta[0]<theta_lim[current_node.theta_i][1]:
                    theta_lim[current_node.theta_i][1]=current_node.bounds_theta[0]
                    # print(f"T {current_node.theta_i} : {theta_lim[current_node.theta_i]}")
        current_node=nodes[current_node.parent]

    score_t_lim=0
    for i in range(len(case["branch"])):
        score_t_lim+=(theta_lim[i][1]-theta_lim[i][0])
    print(f"Score T = {score_t_lim}   {node.number}")
    
    score_v_lim=0
    for i in range(len(case["bus"])):
        bus_i=case["bus"][i]["bus_i"]
        score_v_lim+=(v_lim[bus_i][1]-v_lim[bus_i][0])
    print(f"Score V = {score_v_lim}   {node.number}")

    v = {}
    for i in range(len(case["bus"])):
        bus_i = int(case["bus"][i]["bus_i"])
        v[bus_i] = model.addVar(lb=v_lim[bus_i][0], ub=v_lim[bus_i][1], name=f"v {bus_i}")
        if v_lim[bus_i][1]-v_lim[bus_i][0]<1e-3:
            model.addConstr(v[bus_i]==(v_lim[bus_i][1]-v_lim[bus_i][0])/2)
        
    v2 = {}
    for i in range(len(case["bus"])):
        bus_i = int(case["bus"][i]["bus_i"])
        v2[bus_i] = model.addVar(lb=v_lim[bus_i][0]**2, ub=v_lim[bus_i][1]**2, name=f"v2 {bus_i}")

    theta_real={}
    for i in range(len(case["bus"])):
        bus_i = int(case["bus"][i]["bus_i"])
        theta_real[bus_i]=model.addVar(lb=-2*math.pi,ub=2*math.pi)

    # print("v2 var generated")

    VMAX={}
    for i in range(len(case["bus"])):
        VMAX[case["bus"][i]["bus_i"]]=case["bus"][i]["Vmax"]
    VMIN={}
    for i in range(len(case["bus"])):
        VMIN[case["bus"][i]["bus_i"]]=case["bus"][i]["Vmin"]

    P={}
    Q={}
    c={}
    s={}
    c_mult={}
    s_mult={}
    c_lim={}
    s_lim={}

    for i in range(len(case['branch'])):
        fbus=int(case['branch'][i]['fbus'])
        tbus=int(case['branch'][i]['tbus'])
        theta[i] = model.addVar(lb=theta_lim[i][0], ub=theta_lim[i][1], name=f"theta {i}")   
        if theta_lim[i][1]-theta_lim[i][0]<1e-3:
            model.addConstr(theta[i]==(theta_lim[i][1]-theta_lim[i][0])/2)
        model.addConstr(theta[i]==theta_real[fbus]-theta_real[tbus],name=f"Linking Thetakm {i}, {fbus}-{tbus}")
        vk_max=v_lim[fbus][1]
        vm_max=v_lim[tbus][1]
        vk_min=v_lim[fbus][0]
        vm_min=v_lim[tbus][0]
        if vm_max*vk_max<=0:
            raise ValueError
        P[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"P{(fbus,tbus,i)}") #P_km
        P[i+N]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"P{(tbus,fbus,i+N)}")
        Q[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"Q{(fbus,tbus,i)}") #Q_km
        Q[i+N]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"Q{(tbus,fbus,i+N)}")
        if theta_lim[i][0]*theta_lim[i][1]<0:
            c[i]=model.addVar(lb=min(math.cos(theta_lim[i][0]),math.cos(theta_lim[i][1])),ub=1,name=f"c{(fbus,tbus,i)}")
            c_lim[i]=[min(math.cos(theta_lim[i][0]),math.cos(theta_lim[i][1])),1]
        else:
            c[i]=model.addVar(lb=min(math.cos(theta_lim[i][0]),math.cos(theta_lim[i][1])),ub=max(math.cos(theta_lim[i][0]),math.cos(theta_lim[i][1])),name=f"c{(fbus,tbus,i)}")
            c_lim[i]=[min(math.cos(theta_lim[i][0]),math.cos(theta_lim[i][1])),max(math.cos(theta_lim[i][0]),math.cos(theta_lim[i][1]))]
        s[i]=model.addVar(lb=math.sin(theta_lim[i][0]),ub=math.sin(theta_lim[i][1]),name=f"s{(fbus,tbus,i)}")
        s_lim[i]=[math.sin(theta_lim[i][0]),math.sin(theta_lim[i][1])]
        c_mult[i]=model.addVar(lb=vm_min*vk_min*c_lim[i][0],ub=vm_max*vk_max*c_lim[i][1],name=f"c_mult{(fbus,tbus,i)}") #jabr c_km
        s_mul_vals=[vm_max*vk_max*s_lim[i][0],vm_max*vk_max*s_lim[i][1],vm_min*vk_min*s_lim[i][0],vm_min*vk_min*s_lim[i][1]]
        s_mult[i]=model.addVar(lb=min(s_mul_vals),ub=max(s_mul_vals),name=f"s_mult{(fbus,tbus,i)}") #jabr s_km

    nb_step=20
    for i in range(len(case['branch'])):
        fbus=int(case['branch'][i]['fbus'])
        tbus=int(case['branch'][i]['tbus'])
        # model.addConstr(s[i]==(theta[fbus]-theta[tbus]))
        gaps=[s_lim[i][0]-theta_lim[i][0],s_lim[i][1]-theta_lim[i][1]]
        model.addConstr(s[i]<=theta[i]+max(gaps))
        model.addConstr(s[i]>=theta[i]+min(gaps))
        a=(math.cos(theta_lim[i][1])-math.cos(theta_lim[i][0]))/(theta_lim[i][1]-theta_lim[i][0])
        b=math.cos(theta_lim[i][0])-theta_lim[i][0]*a
        model.addConstr(c[i]>=a*theta[i]+b)
        for k in range(0,nb_step+1):
            x_bar=theta_lim[i][0]+k/nb_step*(theta_lim[i][1]-theta_lim[i][0])
            a=-math.sin(x_bar)
            b=math.cos(x_bar)-x_bar*a
            model.addConstr(c[i]<=a*(theta[i])+b)

    nb_step=5
    for i in range(len(case["bus"])):
        fbus=case["bus"][i]["bus_i"]
        vk_max=v_lim[fbus][1]
        vk_min=v_lim[fbus][0]
        a=vk_max+vk_min
        b=-vk_min*vk_max
        model.addConstr(a*v[fbus]+b >= v2[fbus],name=f"v2 {fbus}")
        step=(v_lim[fbus][1]-v_lim[fbus][0])/nb_step
        for k in range(nb_step+1):
            x_bar=v_lim[fbus][0]+step*k
            a=2*x_bar
            b=-a*x_bar+x_bar**2
            model.addConstr(a*v[fbus]+b <= v2[fbus],name=f"v2 {fbus}")

    for i in range(len(case['branch'])): ### Linear mult
        fbus=int(case['branch'][i]['fbus'])
        tbus=int(case['branch'][i]['tbus'])
        vars=[c[i],v[fbus],v[tbus]]
        cuboid=[c_lim[i],v_lim[fbus],v_lim[tbus]]
        HP=hyperplanes.generate_all_neighbor_hp_local_MAXMIN(cuboid)
        for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
            hyperplane=hp[0][:-1]
            constant=hp[0][-1]
            constraint=np.dot(hyperplane,vars)+constant
            if hp[1]>0:
                model.addConstr(constraint >= c_mult[i],name=f"c_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]<0:
                model.addConstr(constraint <= c_mult[i],name=f"c_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]==0:
                model.addConstr(constraint == c_mult[i],name=f"c_mult {fbus} {tbus}")
        vars=[s[i],v[fbus],v[tbus]]
        cuboid=[s_lim[i],v_lim[fbus],v_lim[tbus]]
        HP=hyperplanes.generate_all_neighbor_hp_local_MAXMIN(cuboid)
        for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
            hyperplane=hp[0][:-1]
            constant=hp[0][-1]
            constraint=np.dot(hyperplane,vars)+constant
            if hp[1]>0:
                model.addConstr(constraint >= s_mult[i],name=f"s_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]<0:
                model.addConstr(constraint <= s_mult[i],name=f"s_mult {fbus} {tbus}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]==0:
                model.addConstr(constraint == s_mult[i],name=f"s_mult {fbus} {tbus}")

    ### Constraints 
    # (1b) & (1c)

    for bus in range(len(case['bus'])):
        bus_number = int(case['bus'][bus]['bus_i'])
        Pd=case['bus'][bus]['Pd']  # Active load
        Qd=case['bus'][bus]['Qd']  # Reactive load
        GenBus=Bus_Gen.get(bus_number)
        Pg_bus=0
        Qg_bus=0
        if GenBus!=None:
            Pg_bus=sum(Pg[i] for i in GenBus) # Active power gen on bus k
            Qg_bus=sum(Qg[i] for i in GenBus) # Reactive power gen on bus k
        P_bus=0
        Q_bus=0
        BranchBus=Branch_Bus.get(bus_number)
        if BranchBus!=None:
            P_bus=sum(P[i] for i in BranchBus) 
            Q_bus=sum(Q[i] for i in BranchBus) 
        # for i,branch in enumerate(case["branch"]): Conventional way
        #     if branch["fbus"]==bus_number:
        #         P_bus+=P[i]
        #         Q_bus+=Q[i]
        #     if branch["tbus"]==bus_number:
        #         P_bus+=P[i+N]
        #         Q_bus+=Q[i+N]
        model.addLConstr(P_bus == Pg_bus - Pd , name=f"ActivePowerBalance_{bus_number}")
        model.addLConstr(Q_bus == Qg_bus - Qd , name=f"ReactivePowerBalance_{bus_number}")

    # print("Const. (1b) & (1c) generated")

    # (6) Jabr const.

    for i in range(len(case['branch'])):
        
        branch_=case["branch"][i]
        G,B=create_admittance_matrix(branch_)
        k=branch_['fbus']
        m=branch_['tbus']

        model.addLConstr(P[i]/baseMVA==G[(k,k)]*v2[k]+G[(k,m)]*c_mult[i]+B[(k,m)]*s_mult[i], name=f"6a_{k},{m},{i}")
        model.addLConstr(P[i+N]/baseMVA==G[(m,m)]*v2[m]+G[(m,k)]*c_mult[i]-B[(m,k)]*s_mult[i], name=f"6b_{m},{k},{i}")
        model.addLConstr(Q[i]/baseMVA==-B[(k,k)]*v2[k]-B[(k,m)]*c_mult[i]+G[(k,m)]*s_mult[i], name=f"6c_{k},{m},{i}")
        model.addLConstr(Q[i+N]/baseMVA==-B[(m,m)]*v2[m]-B[(m,k)]*c_mult[i]-G[(m,k)]*s_mult[i], name=f"6d_{m},{k},{i}")

    # for i in range(len(case["branch"])): #SOCP Jabr
    #     fbus=case["branch"][i]["fbus"]
    #     tbus=case["branch"][i]["tbus"]
    #     model.addConstr(c_mult[i] ** 2 + s_mult[i] ** 2 <= v2[fbus] * v2[tbus]) # (2)

    # for i in range(len(case["branch"])): #SOCP Jabr
    #     fbus=case["branch"][i]["fbus"]
    #     tbus=case["branch"][i]["tbus"]
    #     model.addConstr(c[i] ** 2 + s[i] ** 2 <= 1) # (2)

    # print("Const. (6) Jabr generated")

    SumPd=0
    for i in range(len(case["bus"])):
        if case["bus"][i]["type"] in [1,2,3]:
            SumPd+=case["bus"][i]["Pd"]

    i2={}
    for i in range(len(case['branch'])):
        branch_=case["branch"][i]
        rateA=case["branch"][i]["rateA"]
        r = branch_["r"] # real part of impendence
        x = branch_["x"] # complex part of impendence
        z = r + x*1j # impendence
        y = 1/z
        g=y.real
        b=y.imag

        fbus=case['branch'][i]["fbus"]
        tbus=case['branch'][i]["tbus"]

        ratio = case["branch"][i]["ratio"]
        bshunt = case["branch"][i]["b"]
        angle = case["branch"][i]["angle"]

        if ratio==0:
            ratio=1

        # mp_cbusf=v2[fbus]
        # mp_cbust=v2[tbus]
        # mp_c=c[i]
        # mp_s=s[i]

        # i2_var  = 0
        # i2_var += (g*g + b*b)/(ratio*ratio) * ( (mp_cbusf/(ratio*ratio)) + mp_cbust - (2/ratio) * ( mp_c * math.cos(angle) + mp_s * math.sin(angle) ) )
        # i2_var += b*bshunt/(ratio**3) * ( (mp_cbusf/ratio) - (mp_c * math.cos(angle) + mp_s * math.sin(angle) ))
        # i2_var += g*bshunt/(ratio**3) * ( mp_s * math.cos(angle) - mp_c * math.sin(angle) )
        # i2_var += (bshunt*bshunt*mp_cbusf/(4*(ratio**4)) )

        tau=ratio
        sigma=angle

        gsh=0
        bsh=bshunt

        Alpha=(1/(tau**4))*((g**2+b**2)+(g*gsh+b*bsh)+((gsh**2+bsh**2)/4))
        Beta=(g**2+b**2)/(tau**2)
        Gamma=(1/(tau**3))*(math.cos(sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))+math.sin(sigma)*(b*gsh-g*bsh))
        Ki=(1/(tau**3))*(math.sin(sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))-math.cos(sigma)*(b*gsh-g*bsh))

        i2like=Alpha*v2[fbus]+Beta*v2[tbus]+Gamma*c_mult[i]+Ki*s_mult[i]

        v_min_k=VMIN[fbus]

        rho=100
        
        if Alpha<rho:
            if rateA>0:
                rateA=case["branch"][i]["rateA"]/baseMVA
                i2[i]=model.addVar(lb=0,ub=(rateA/v_min_k)**2,name=f"i2{(fbus,tbus,i)}") #i2km
                model.addConstr(i2like==i2[i],name=f"goodi2 {fbus},{tbus},{i}")
                # i2[i+N]=model.addVar(lb=0,ub=(rateA/v_min_m)**2,name=f"i2{(tbus,fbus,i)}") #i2km
                # model.addConstr(i2like_mk==i2[i+N],name=f"goodi2 {tbus},{fbus},{i}")
                # model.addConstr((P[i]**2+Q[i]**2)/(baseMVA**2)<=i2[i]*v2[fbus])
            else:
                rateA=(2*SumPd)/baseMVA
                i2[i]=model.addVar(lb=0,ub=(rateA/v_min_k)**2,name=f"i2{(fbus,tbus,i)}") #i2km
                model.addConstr(i2like==i2[i],name=f"goodi2 {fbus},{tbus},{i}")
                # i2[i+N]=model.addVar(lb=0,ub=(rateA/v_min_m)**2,name=f"i2{(fbus,tbus,i)}") #i2km
                # model.addConstr(i2like_mk==i2[i+N],name=f"goodi2 {fbus},{tbus},{i}")
            # model.addConstr((P[i]**2+Q[i]**2)/(baseMVA**2)<=i2[i]*v2[fbus])

        else:
            if rateA>0:
                rateA=case["branch"][i]["rateA"]/baseMVA
                model.addConstr(i2like/Alpha>=0,name=f"badi2-2 {fbus},{tbus},{i}")
                model.addConstr(i2like/Alpha<=((rateA/v_min_k)**2)/Alpha,name=f"badi2-1 {fbus},{tbus},{i}")
            else:
                rateA=(2*SumPd)/baseMVA
                model.addConstr(i2like/Alpha>=0,name=f"badi2-2 {fbus},{tbus},{i}")
                model.addConstr(i2like/Alpha<=((rateA/v_min_k)**2)/Alpha,name=f"badi2-1 {fbus},{tbus},{i}")

    # print("Const. i2 generated")

    # for i in range(len(case['branch'])): #SOCP limit
    #     rateA=case['branch'][i]["rateA"]/baseMVA
    #     if rateA!=0:
    #         model.addConstr(P[i] ** 2 + Q[i] ** 2 <= rateA**2*(baseMVA**2),name=f"limit {i},{rateA}")
    #         model.addConstr(P[i+N] ** 2 + Q[i+N] ** 2 <= rateA**2*(baseMVA**2),name=f"limit {i+N},{rateA}")
    #     else:
    #         rateA=2*SumPd/baseMVA
    #         model.addConstr(P[i] ** 2 + Q[i] ** 2 <= rateA**2*(baseMVA**2),name=f"limit {i},{rateA}")
    #         model.addConstr(P[i+N] ** 2 + Q[i+N] ** 2 <= rateA**2*(baseMVA**2),name=f"limit {i+N},{rateA}")
            

    # print("Const. limit generated")

    model.setParam('Threads', 32) 
    model.setParam('BarHomogeneous',1)
    model.setParam('NumericFocus',1)
    model.setParam('Crossover',0)
    model.setParam('Method', 2)

    objective = gp.QuadExpr()

    for i in range(len(case['gen'])):
        ca,cb,cc = case['gencost'][i]['costvector']
        objective += ca*Pg[i]*Pg[i]+cb*Pg[i]+cc

    model.setObjective(objective, GRB.MINIMIZE)

    model.optimize()

    # if model.status == GRB.INFEASIBLE:
        # print("Model is infeasible.")
        # You can try to find infeasibilities by calling the following:
        # model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
        # model.write("infeasible_model.ilp")

    # print(f'Temps : {time.time()-duration}')

    Cuts=[]

    def JabrCuts(eps,pT,Const_to_add): 
        violation_amount=[]
        for i in range(len(case['branch'])):
            fbus=case['branch'][i]["fbus"]
            tbus=case['branch'][i]["tbus"]
            amount=c_mult[i].X*c_mult[i].X+s_mult[i].X*s_mult[i].X-v2[fbus].X*v2[tbus].X
            if amount>eps:
                violation_amount.append([amount,i,c_mult[i].X,s_mult[i].X,v2[fbus].X,v2[tbus].X])
        violation_amount.sort(reverse=True) # Order the violation amounts
        violation_amount=violation_amount[:round(len(violation_amount)*pT)] # Keep the top p_T percent
    # for i in range(len(violation_amount)):
        for i in range(len(violation_amount)):
            branch_i=violation_amount[i][1]
            fbus=case['branch'][branch_i]["fbus"]
            tbus=case['branch'][branch_i]["tbus"]
            c_=violation_amount[i][2]
            s_=violation_amount[i][3]
            v2k=violation_amount[i][4]
            v2m=violation_amount[i][5]
            if v2k+v2m>eps:
                n0=((2*c_)**2+(2*s_)**2+(v2k-v2m)**2)**(1/2)
                Const_to_add.append([4*c_*c_mult[branch_i]+4*s_*s_mult[branch_i]+(v2k-v2m-n0)*v2[fbus]+(v2m-v2k-n0)*v2[tbus],f"CutJabr_{fbus},{tbus},{branch_i}"])
        return len(violation_amount)

    def i2Cuts(eps,pT,Const_to_add): 
            # model.addConstr((P[i]**2+Q[i]**2)/(baseMVA**2)<=i2[i]*v[fbus]**2)
        violation_amount=[]

        for i in range(len(case['branch'])):
            if i in i2:
                fbus=int(case['branch'][i]["fbus"])
                tbus=int(case['branch'][i]["tbus"])
                if i>=N:
                    tbus=int(case['branch'][i]["fbus"])
                    fbus=int(case['branch'][i]["tbus"])
                amount=(P[i].X**2+Q[i].X**2)/(baseMVA**2)-i2[i].X*v2[fbus].X
                if amount>eps:
                    violation_amount.append([amount,0,P[i].X/(baseMVA),Q[i].X/(baseMVA),v2[fbus].X,i2[i].X,i])  

        violation_amount.sort(key=lambda x: x[0],reverse=True) # Order the violation amounts
        violation_amount=violation_amount[:round(len(violation_amount)*pT)] # Keep the top p_T percent
        for i in range(len(violation_amount)):
            branch_i=violation_amount[i][-1]
            fbus=int(case['branch'][branch_i]["fbus"])
            tbus=int(case['branch'][branch_i]["tbus"])
            P_=violation_amount[i][2]
            Q_=violation_amount[i][3]
            v2k=violation_amount[i][4]
            i2kX=violation_amount[i][5]
            if v2k+i2kX>eps :
                n0=((2*P_)**2+(2*Q_)**2+(v2k-i2kX)**2)**(1/2)
                Const_to_add.append([(4*P_*P[branch_i]+4*Q_*Q[branch_i])/baseMVA+(v2k-i2kX-n0)*v2[fbus]+(i2kX-v2k-n0)*i2[branch_i],f"Cuti2_{fbus},{tbus},{branch_i}"])
        return len(violation_amount)

    def LimitCuts(eps,pT,Const_to_add): # P[i] ** 2 + Q[i] ** 2 <= rate*baseMVA**2
        violation_amount=[]
        for i in range(len(case['branch'])):
            fbus=case['branch'][i]["fbus"]
            tbus=case['branch'][i]["tbus"]
            rateA=case["branch"][i]["rateA"]/baseMVA
            if rateA==0:
                rateA=2*SumPd/baseMVA
            amount=(P[i].X*P[i].X+Q[i].X*Q[i].X)-(rateA**2)*baseMVA**2
            if amount>eps:
                violation_amount.append([amount,i,P[i].X,Q[i].X])
        violation_amount.sort(reverse=True) # Order the violation amounts
        violation_amount=violation_amount[:round(len(violation_amount)*pT)] # Keep the top p_T percent
    # for i in range(len(violation_amount)):
        for i in range(len(violation_amount)):
            branch_i=violation_amount[i][1]
            if branch_i<N:
                fbus=case['branch'][branch_i]["fbus"]
                tbus=case['branch'][branch_i]["tbus"]
                rateA=case["branch"][branch_i]["rateA"]/baseMVA
            else:
                fbus=case['branch'][branch_i-N]["fbus"]
                tbus=case['branch'][branch_i-N]["tbus"]
                rateA=case["branch"][branch_i-N]["rateA"]/baseMVA
            if rateA==0:
                rateA=2*SumPd/baseMVA
            P_=violation_amount[i][2]
            Q_=violation_amount[i][3]
            if violation_amount[i][0]>eps and rateA>0:
                norm_PQ=(P_**2+Q_**2)**(1/2)
                Const_to_add.append([(P_*P[branch_i]+Q_*Q[branch_i])-norm_PQ*rateA*baseMVA,f"CutLimit_{fbus},{tbus},{branch_i}"])
        return len(violation_amount)
    
    def SC_Cuts(eps,pT,Const_to_add): 
        violation_amount=[]
        for i in range(len(case['branch'])):
            fbus=case['branch'][i]["fbus"]
            tbus=case['branch'][i]["tbus"]
            amount=c[i].X*c[i].X+s[i].X*s[i].X-1
            if amount>eps:
                violation_amount.append([amount,i,c[i].X,s[i].X])
        violation_amount.sort(reverse=True) # Order the violation amounts
        violation_amount=violation_amount[:round(len(violation_amount)*pT)] # Keep the top p_T percent
    # for i in range(len(violation_amount)):
        for i in range(len(violation_amount)):
            branch_i=violation_amount[i][1]
            fbus=case['branch'][branch_i]["fbus"]
            tbus=case['branch'][branch_i]["tbus"]
            c_=violation_amount[i][2]
            s_=violation_amount[i][3]
            n0=((2*c_)**2+(2*s_)**2)**(1/2)
            Const_to_add.append([4*c_*c_mult[branch_i]+4*s_*s_mult[branch_i]+(-n0)*v2[fbus]+(-n0)*v2[tbus],f"CutSC_{fbus},{tbus},{branch_i}"])
        return len(violation_amount)

    def f_eps_par(new_cut, Cuts):
        dict_new_cut={}
        for term in range(new_cut.size()):
            dict_new_cut[new_cut.getVar(term)]=new_cut.getCoeff(term)
        new_cut_norm=0
        for i in dict_new_cut:
            new_cut_norm+=dict_new_cut[i]**2
        new_cut_norm=new_cut_norm**(1/2)

        for cuts in Cuts:
            dot_product=0
            cut=cuts[0]
            row = model.getRow(cut)
            for k in range(row.size()):
                if row.getVar(k) in dict_new_cut:
                    dot_product+=dict_new_cut[row.getVar(k)]*row.getCoeff(k)
            if dot_product!=0:
                cut_norm=0
                for k in range(row.size()):
                    cut_norm+=row.getCoeff(k)**2
                cut_norm=cut_norm**(1/2)

                cos_theta = dot_product / (new_cut_norm * cut_norm)

                if cos_theta > 1 - eps_par:
                    return False 
        return True

    T=0
    Tlim=1000
    pjabr=1
    pi2=1
    plim=1
    Tage=5
    Tftol=5
    eps=10E-5
    eps_par=10E-5
    eps_ftol=10E-5
    eps_slack=10E-5
    r=0


    nb_const=[[0]]
    nb_const_tot=[0]
    nb_const_deleted=[0]
    obj=model.getObjective()
    val=[obj.getValue()]
    temps=[time.time()-duration]

    initval=val[0]
    ftol=True
    init=True #I am still on the root ?
    while(T<Tlim and (r<Tftol or init)):

        nb_const.append([])
        T+=1
        # print(f"{T}--------------------CUTS-------------------")
        
        to_remove=[]
        for i,cut in enumerate(Cuts): #Cut old ones that aren't slack
            constraint = cut[0]
            slack=-model.getRow(constraint).getValue()
            if slack > eps_slack and (T - cut[1]) >= Tage and not init: # Removing condition
                to_remove.append([constraint,cut])
        remove_compt=0
        for rr in range(len(to_remove)):
            remove=to_remove[rr]
            model.remove(remove[0])  # Remove from the model
            Cuts.remove(remove[1])  # Remove from Cuts
            remove_compt+=1
        nb_const_deleted.append(remove_compt)
        
        Const_to_add=[]
        nb_const[-1].append(JabrCuts(eps,pjabr,Const_to_add))
        nb_const[-1].append(SC_Cuts(eps,pjabr,Const_to_add))
        nb_const[-1].append(i2Cuts(eps,pi2,Const_to_add))
        nb_const[-1].append(LimitCuts(eps,plim,Const_to_add))
        # print(nb_const[-1])
        nb_const_tot.append(sum(nb_const[-1]))
        if(nb_const_tot[-1]==0):
            nb_const_tot.remove(nb_const_tot[-1])
            nb_const_deleted.remove(nb_const_deleted[-1])
            break
        cuts_to_add=[True for i in range(len(Const_to_add))]
        # for i,new_cuts in enumerate(Const_to_add): # Test de parallélisme, enelever la boucle.
        #     if not f_eps_par(new_cuts[0],Cuts):
        #         cuts_to_add[i]=False
        #         nb_const_tot[-1]-=1
        for i in range(len(cuts_to_add)):
            if cuts_to_add[i]:
                Cuts.append([model.addLConstr(Const_to_add[i][0] <= 0, name=Const_to_add[i][1]),T])

        model.update()  # Synchronize the model after removal
        
        model.optimize()

        temps.append(time.time()-duration)
        if model.status == GRB.INFEASIBLE:
            # print("Model is infeasible.")
            # model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
            # model.write("infeasible_model.ilp")
            break
        obj=model.getObjective()
        val.append(obj.getValue())
        if(val[-1]-val[-2]<val[-2]*eps_ftol):
            r+=1
        else:
            r=0
        if val[-1]-initval>10**-5:
            init=False
        
        
    # print(val)
    # print(nb_const_tot)
    # print(nb_const)

    ### Finding max v_i violation
    max_v_i=-1
    max_so_far=-1
    for i in range(len(case["bus"])):
        fbus=case["bus"][i]["bus_i"]
        # qtte_v=min(abs(v2[fbus].X-v[fbus].X**2)+abs(v_lim[fbus][1]-v_lim[fbus][0]),(v_lim[fbus][1]-v_lim[fbus][0]))
        qtte_v=abs(v2[fbus].X-v[fbus].X**2)
        if qtte_v>max_so_far:
            max_so_far=qtte_v
            max_v_i=int(fbus)
    print(f"MAX VIOL V : {qtte_v}")

    ### Finding max theta_i violation
    max_theta_i=-1
    max_so_far=-1
    for i in range(len(case["branch"])):
        # qtte_theta=min(abs(math.cos(theta[i].X)-c[i].X)+abs(math.sin(theta[i].X)-s[i].X),(theta_lim[i][1]-theta_lim[i][0]))
        qtte_theta=abs(math.cos(theta[i].X)-c[i].X)+abs(math.sin(theta[i].X)-s[i].X)
        if qtte_theta>max_so_far:
            max_so_far=qtte_theta
            max_theta_i=int(i)
    print(f"MAX VIOL theta : {qtte_theta}")

    ### Computing violation quantity ###

    tot_val_viol=0
    for i in range(len(case['branch'])):
        
        branch_=case["branch"][i]
        fbus=branch_["fbus"]
        tbus=branch_["tbus"]
        G,B=create_admittance_matrix(branch_)
        k=branch_['fbus']
        m=branch_['tbus']
        aureol=c_mult[i].X**2+s_mult[i].X**2
        new_c=c_mult[i].X/math.sqrt(aureol)
        new_s=s_mult[i].X/math.sqrt(aureol)

        val_viol=0
        val_viol+=abs(P[i].X/baseMVA-(G[(k,k)]*v2[k].X+G[(k,m)]*new_c+B[(k,m)]*new_s))
        val_viol+=abs(P[i+N].X/baseMVA-(G[(m,m)]*v2[m].X+G[(m,k)]*new_c-B[(m,k)]*new_s))
        val_viol+=abs(Q[i].X/baseMVA-(-B[(k,k)]*v2[k].X-B[(k,m)]*new_c+G[(k,m)]*new_s))
        val_viol+=abs(Q[i+N].X/baseMVA-(-B[(m,m)]*v2[m].X-B[(m,k)]*new_c-G[(m,k)]*new_s))
        tot_val_viol+=val_viol
        # print(f"Branch {fbus},{tbus} = {val_viol:.3e}")
    # print(f"TOT VIOLATION ={tot_val_viol}")

    if plot:
        plt.show()
    
    if qtte_v+qtte_theta<1e-3: ##VALID SOLUTION
        # return [val[-1],tot_val_viol,max_v_i,v[max_v_i].X,max_theta_i,theta[max_theta_i].X,-1]
        return [val[-1],tot_val_viol,max_v_i,v[max_v_i].X,max_theta_i,(theta_lim[max_theta_i][0]+theta_lim[max_theta_i][1])/4+theta[max_theta_i].X/2,-2]
      
    # return [val[-1],tot_val_viol,max_v_i,(v_lim[max_v_i][0]+v_lim[max_v_i][1])/2,max_theta_i,(theta_lim[max_theta_i][0]+theta_lim[max_theta_i][1])/2,qtte_v>qtte_theta]
    return [val[-1],tot_val_viol,max_v_i,v[max_v_i].X/2+(v_lim[max_v_i][0]+v_lim[max_v_i][1])/4,max_theta_i,theta[max_theta_i].X/2+(theta_lim[max_theta_i][1]+theta_lim[max_theta_i][0])/4,10*qtte_v>qtte_theta]

def Solve_tot(case,plot,theta_lim_):

    def GenerateFigli(node,next_num):
        
        if node.v_or_theta==True:
            node1=Node(next_num,node.next_v[0],node.number,[node.next_v[1],0],node.next_theta[0],None,node.val,None,None,None)
            node2=Node(next_num+1,node.next_v[0],node.number,[node.next_v[1],1],node.next_theta[0],None,node.val,None,None,None)
        elif node.v_or_theta==False:
            node1=Node(next_num,node.next_v[0],node.number,None,node.next_theta[0],[node.next_theta[1],0],node.val,None,None,None)
            node2=Node(next_num+1,node.next_v[0],node.number,None,node.next_theta[0],[node.next_theta[1],1],node.val,None,None,None)

        nodes[next_num]=node1
        nodes[next_num+1]=node2
        return([node1,node2])
    def Chose_next_node(nodes,nodes_to_visit):
        min_val=nodes[nodes_to_visit[0]].val
        val_i=nodes[nodes_to_visit[0]].number
        for i in nodes_to_visit:
            if nodes[i].val<=min_val and nodes[i].next_v!=[-1,-1] and nodes[i].next_theta!=[-1,-1] and nodes[i].v_or_theta!=-1:
                min_val=nodes[i].val
                val_i=nodes[i].number
        print(f"CURRENT VAL : {min_val}")
        if min_val>1000000000:
            return -1
        return val_i

    nodes={}
    nodes_to_visit=[]
    nodes[0]=Node(0,-1,-1,None,-1,None,None,None,None,None)
    current_sol=Solve(case,False,theta_lim,nodes[0],nodes)
    nodes[0].val=current_sol[0]
    nodes[0].next_v=[current_sol[2],current_sol[3]]
    nodes[0].next_theta=[current_sol[4],current_sol[5]]
    nodes[0].v_or_theta=current_sol[6]
    next_node=0
    FeasibleVal=[]
    for _ in range(10000):
        compteur=len(nodes)
        # print(nodes[compteur-1])
        if nodes[next_node].v_or_theta!=-1 and current_sol[6]!=-2:
            Figli=GenerateFigli(nodes[next_node],compteur)
            nodes_to_visit.append(Figli[0].number)
            nodes_to_visit.append(Figli[1].number)
        next_node=Chose_next_node(nodes,nodes_to_visit)
        if next_node==-1:
            raise ValueError
        nodes_to_visit.remove(next_node)
        try:
            current_sol=Solve(case,False,theta_lim_,nodes[next_node],nodes)
        except:
            current_sol=[100000000000000,10000000000,-1,-1,-1,-1,-1]
        if current_sol[6]==-2:
            print(f"FEASIBLE VAL : {current_sol[0]}")
            FeasibleVal.append(current_sol[0])
            break

        nodes[next_node].val=current_sol[0]
        nodes[next_node].next_v=[current_sol[2],current_sol[3]]
        nodes[next_node].next_theta=[current_sol[4],current_sol[5]]
        nodes[next_node].v_or_theta=current_sol[6]

        
    # for i in nodes:
    #     print(nodes[i])
    print(f"Feasibles solutions : {FeasibleVal}")
    return nodes

if __name__=="__main__": 
    local_path="cases/"
    
    name="case9"
    # name="case14"
    name="case30"
    # name="case57"
    # name="case118"
    # name="caseNY"
    # name="case_ACTIVSg500"
    # name="case1354pegase"
    # name="case9241pegase"
    # name="case13659pegase"
    # name="case_ACTIVSg25k"
    path=local_path+name+".m"
    case = conversion.conversion(path)
    theta_lim=math.pi/2
    Solve_tot(case,True,theta_lim)