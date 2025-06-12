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
import ObjBound


def Solve(case,plot,tightbound):
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
    print(f"Number of branches : {N}")

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

    ### Variables definition & intervals

    Pg = []  # Active power generation
    Qg = []  # Reactive power generation

    for i in range(len(case['gen'])):
        Pg.append(model.addVar(lb=case['gen'][i]['Pmin'], ub=case['gen'][i]['Pmax'], name=f"Pg{i}"))
        Qg.append(model.addVar(lb=case['gen'][i]['Qmin'], ub=case['gen'][i]['Qmax'], name=f"Qg{i}"))

    print("Power Gen var generated")

    SumPd=0
    for i in range(len(case["bus"])):
        if case["bus"][i]["type"] in [1,2,3]:
            SumPd+=case["bus"][i]["Pd"]

    v2 = {}
    for i in range(len(case["bus"])):
        bus_i = int(case["bus"][i]["bus_i"])
        v2[bus_i] = model.addVar(lb=case["bus"][i]["Vmin"]**2, ub=case["bus"][i]["Vmax"]**2, name=f"v2_{bus_i}")
        v2[bus_i].start=1

    print("v2 var generated")

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

    for i in range(len(case['branch'])):
        fbus=int(case['branch'][i]['fbus'])
        tbus=int(case['branch'][i]['tbus'])
        vk_max=VMAX[fbus]
        vm_max=VMAX[tbus]
        vk_min=VMIN[fbus]
        vm_min=VMIN[tbus]
        if vm_max*vk_max<=0:
            raise ValueError
        P[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"P{(fbus,tbus,i)}") #P_km
        P[i+N]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"P{(tbus,fbus,i+N)}")
        Q[i]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"Q{(fbus,tbus,i)}") #Q_km
        Q[i+N]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"Q{(tbus,fbus,i+N)}")
        c[i]=model.addVar(lb=0,ub=vm_max*vk_max,name=f"c{(fbus,tbus,i)}") #jabr c_km
        s[i]=model.addVar(lb=-vm_max*vk_max,ub=vm_max*vk_max,name=f"s{(fbus,tbus,i)}") #jabr s_km
        c[i].start=1
        s[i].start=0

    print("Branch var generated")

    ### Constraints 
    # (1b) & (1c)

    for bus in tqdm(range(len(case['bus'])), desc="Const. (1b) & (1c) ...", ncols=100):
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

    print("Const. (1b) & (1c) generated")

    # (6) Jabr const.

    for i in tqdm(range(len(case['branch'])), desc="Const. (6) Jabr ..." ,ncols=100):
        
        branch_=case["branch"][i]
        G,B=create_admittance_matrix(branch_)
        k=branch_['fbus']
        m=branch_['tbus']

        model.addLConstr(P[i]/baseMVA==G[(k,k)]*v2[k]+G[(k,m)]*c[i]+B[(k,m)]*s[i], name=f"6a_{k},{m},{i}")
        model.addLConstr(P[i+N]/baseMVA==G[(m,m)]*v2[m]+G[(m,k)]*c[i]-B[(m,k)]*s[i], name=f"6b_{m},{k},{i}")
        model.addLConstr(Q[i]/baseMVA==-B[(k,k)]*v2[k]-B[(k,m)]*c[i]+G[(k,m)]*s[i], name=f"6c_{k},{m},{i}")
        model.addLConstr(Q[i+N]/baseMVA==-B[(m,m)]*v2[m]-B[(m,k)]*c[i]-G[(m,k)]*s[i], name=f"6d_{m},{k},{i}")

    # for i in range(len(case["branch"])): #SCOP Jabr
    #     fbus=case["branch"][i]["fbus"]
    #     tbus=case["branch"][i]["tbus"]
    #     model.addConstr(c[i] ** 2 + s[i] ** 2 <= v2[fbus] * v2[tbus]) # (2)

    print("Const. (6) Jabr generated")

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

        # i2like  = 0
        # i2like += (g*g + b*b)/(ratio*ratio) * ( (mp_cbusf/(ratio*ratio)) + mp_cbust - (2/ratio) * ( mp_c * math.cos(angle) + mp_s * math.sin(angle) ) )
        # i2like += b*bshunt/(ratio**3) * ( (mp_cbusf/ratio) - (mp_c * math.cos(angle) + mp_s * math.sin(angle) ))
        # i2like += g*bshunt/(ratio**3) * ( mp_s * math.cos(angle) - mp_c * math.sin(angle) )
        # i2like += (bshunt*bshunt*mp_cbusf/(4*(ratio**4)) )

        tau=ratio
        sigma=angle

        gsh=0
        bsh=bshunt

        Alpha=(1/(tau**4))*((g**2+b**2)+(g*gsh+b*bsh)+((gsh**2+bsh**2)/4))
        Beta=(g**2+b**2)/(tau**2)
        Gamma=(1/(tau**3))*(math.cos(sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))+math.sin(sigma)*(b*gsh-g*bsh))
        Ki=(1/(tau**3))*(math.sin(sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))-math.cos(sigma)*(b*gsh-g*bsh))

        # Gamma_ck=(1/(tau**3))*(math.cos(-sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))+math.sin(-sigma)*(b*gsh-g*bsh))
        # Ki_ck=(1/(tau**3))*(math.sin(-sigma)*(-2*(g**2+b**2)-(g*gsh+b*bsh))-math.cos(-sigma)*(b*gsh-g*bsh))

        i2like=Alpha*v2[fbus]+Beta*v2[tbus]+Gamma*c[i]+Ki*s[i]
        # i2like_mk=Alpha*v2[tbus]+Beta*v2[fbus]+Gamma_ck*c[i]-Ki_ck*s[i]

        v_min_k=VMIN[fbus]
        # v_min_m=VMIN[tbus]

        rho=100

        if Alpha<rho:
            if rateA>0:
                rateA=case["branch"][i]["rateA"]/baseMVA
                i2[i]=model.addVar(lb=0,ub=(rateA/v_min_k)**2,name=f"i2{(fbus,tbus,i)}") #i2km
                model.addConstr(i2like==i2[i],name=f"goodi2 {fbus},{tbus},{i}")
            else:
                rateA=(2*SumPd)/baseMVA
                i2[i]=model.addVar(lb=0,ub=(rateA/v_min_k)**2,name=f"i2{(fbus,tbus,i)}") #i2km
                model.addConstr(i2like==i2[i],name=f"goodi2 {fbus},{tbus},{i}")
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

    # print("ALPHA MEAN "+str(alpha_tot/len(case["branch"])))
    # print("BETA/ALPHA MEAN "+str(beta_tot/len(case["branch"])))
    # print("GAMMA/ALPHA MEAN "+str(gamma_tot/len(case["branch"])))
    # print("KI/ALPHA MEAN "+str(Ki_tot/len(case["branch"])))

    print("Const. i2 generated")

    # for i in range(len(case['branch'])): #SCOP limit
    #     fbus=case["branch"][i]["fbus"]
    #     tbus=case["branch"][i]["tbus"]
    #     rateA=case['branch'][i]["rateA"]
    #     if rateA!=0:
    #         model.addConstr(P[i] ** 2 + Q[i] ** 2 <= rateA**2/(baseMVA**2),name=f"limit {i},{rateA}")
    #         model.addConstr(P[i+N] ** 2 + Q[i+N] ** 2 <= rateA**2/(baseMVA**2),name=f"limit {i+N},{rateA}")

    print("Const. limit generated")

    model.setParam('Threads', 32) 
    model.setParam('BarHomogeneous',1)
    model.setParam('NumericFocus',1)
    model.setParam('Crossover',0)
    model.setParam('Method', 2)
    ### Loop Constraints


    G = nx.Graph()  # building the support graph of the network
    for i in range(len(case["branch"])):
        fbus=case["branch"][i]["fbus"]
        tbus=case["branch"][i]["tbus"]
        G.add_edge(fbus,tbus)
    loops = nx.minimum_cycle_basis(G) 
    print("LOOPS : ")
    print(loops)
    # nx.draw(G, with_labels = True)
    init_loops=loops.copy()

    fbustbus_ID_new={}
    c_new={}
    s_new={}
    id=0
    for loop in loops: # Construction des 3 loops
        if len(loop)>3:
            for i in range(2,len(loop)-1):
                id+=1
                G.add_edge(loop[0],loop[i])
                fbus=loop[0]
                tbus=loop[i]
                fbustbus_ID_new[(loop[0],loop[i])]=id
                fbustbus_ID_new[(loop[i],loop[0])]=-id
                vk_max=VMAX[fbus]
                vm_max=VMAX[tbus]
                c_new[id]=model.addVar(lb=0,ub=vm_max*vk_max,name=f"c{(fbus,tbus,id)}") #jabr c_km
                s_new[id]=model.addVar(lb=-vm_max*vk_max,ub=vm_max*vk_max,name=f"s{(fbus,tbus,id)}") #jabr s_km

    new_loops=[]
    for loop in loops: 
        if len(loop)<=3:
            new_loops.append(loop)
        else:
            for i in range(2,len(loop)-1):
                new_loops.append([loop[0],loop[i],loop[i-1]])
            new_loops.append([loop[0],loop[-1],loop[-2]])
    loops=new_loops
    # nx.draw(G, with_labels = True)
    print("LOOPS : ")
    print(loops)
    loops_3=loops.copy()

    z={}
    z_var={}

    def gen_var_on_loop(loops):
        """Generate the mult variable on all loops
        z_var(loop n°,k)=[...]"""
        for j in range(len(loops)):
            loop=loops[j]
            for k in range(len(loop)//2+1):
                z[(j,k)]=[]
                z_var[(j,k)]=[]
                combinations = list(itertools.combinations(loop, 2*k))
                # print(f'{j,k} {combinations}')
                for i in range(len(combinations)):
                    z[(j,k)].append(combinations[i])
                    z_var[(j,k)].append(model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"z{k,i}"))
        for j in range(len(loops)):
            z_var[j+len(loops)]=model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"z{k,i}")


    fbustbus_ID={}
    for i in range(len(case["branch"])):
        fbus=case["branch"][i]["fbus"]
        tbus=case["branch"][i]["tbus"]
        fbustbus_ID[(fbus,tbus)]=i
        fbustbus_ID[(tbus,fbus)]=i+len(case["branch"])


    def get_cuboid_var(loop_number,k,id_comb,cs_bounds,cs_bounds_new):
        s_list=z[(loop_number,k)][id_comb]
        s_var_list=[]
        c_var_list=[]
        for i in loops[loop_number]:
            if i in s_list:
                s_var_list.append((i,loops[loop_number][(loops[loop_number].index(i)+1)%len(loops[loop_number])]))
            else:
                c_var_list.append((i,loops[loop_number][(loops[loop_number].index(i)+1)%len(loops[loop_number])]))

        var_bounds=[]

        for i in s_var_list:
            if (i[0],i[1]) in fbustbus_ID:
                id=fbustbus_ID[i[0],i[1]]
                index=s_var_list.index(i)
                sLB=cs_bounds[0][index]
                sUB=cs_bounds[1][index]
                if id<len(case['branch']):
                    var_bounds.append([sLB,sUB])
                else:
                    id=id-len(case["branch"])
                    var_bounds.append([-sUB,-sLB])
                    # var_bounds.append([sLB,sUB])
            else:
                id=fbustbus_ID_new[i[0],i[1]]
                if id>0:
                    sLB=cs_bounds_new[0][id-1]
                    sUB=cs_bounds_new[1][id-1]
                    var_bounds.append([sLB,sUB])
                else:
                    id=-id
                    sLB=cs_bounds_new[0][id-1]
                    sUB=cs_bounds_new[1][id-1]
                    var_bounds.append([-sUB,-sLB])
                    # var_bounds.append([sLB,sUB])


        for i in c_var_list:
            if (i[0],i[1]) in fbustbus_ID:
                id=fbustbus_ID[i[0],i[1]]
                index=c_var_list.index(i)
                cLB=cs_bounds[2][index]
                cUB=cs_bounds[3][index]
                if id<len(case['branch']):
                    var_bounds.append([cLB,cUB])
                else:
                    id=id-len(case["branch"])
                    # var_bounds.append([-cUB,-cLB])
                    var_bounds.append([cLB,cUB])
            else:
                id=fbustbus_ID_new[i[0],i[1]]
                if id>0:
                    cLB=cs_bounds_new[2][id-1]
                    cUB=cs_bounds_new[3][id-1]
                    var_bounds.append([cLB,cUB])
                else:
                    id=-id
                    cLB=cs_bounds_new[2][id-1]
                    cUB=cs_bounds_new[3][id-1]
                    # var_bounds.append([-cUB,-cLB])
                    var_bounds.append([cLB,cUB])

        return var_bounds

    def gen_constraints(loop_number,id_k,id_comb,cs_bounds,cs_bounds_new):
        s_list=z[(loop_number,id_k)][id_comb]
        s_var_list=[]
        c_var_list=[]
        for i in loops[loop_number]:
            if i in s_list:
                s_var_list.append((i,loops[loop_number][(loops[loop_number].index(i)+1)%len(loops[loop_number])]))
            else:
                c_var_list.append((i,loops[loop_number][(loops[loop_number].index(i)+1)%len(loops[loop_number])]))
        cuboid=get_cuboid_var(loop_number,k,id_comb,cs_bounds,cs_bounds_new)
        HP=hyperplanes.generate_all_neighbor_hp_local_MAXMIN(cuboid)
        compteur=0
        vars=[]
        for i in s_var_list:
            if (i[0],i[1]) in fbustbus_ID:
                id=fbustbus_ID[i[0],i[1]]
                if id<len(case["branch"]):
                    vars.append(s[id])
                else:
                    id=id-len(case["branch"])
                    vars.append(-s[id])
                    compteur+=1
                    # vars.append(s[id])
            else:
                id=fbustbus_ID_new[i[0],i[1]]
                if id>0:
                    vars.append(s_new[id])
                if id<0:
                    id=-id
                    vars.append(-s_new[id])
                    compteur+=1
                    # vars.append(s_new[id])


        for i in c_var_list:
            if (i[0],i[1]) in fbustbus_ID:
                id=fbustbus_ID[i[0],i[1]]
                if id<len(case["branch"]):
                    vars.append(c[id])
                else:
                    id=id-len(case["branch"])
                    # vars.append(-c[id])
                    vars.append(c[id])
            else:
                id=fbustbus_ID_new[i[0],i[1]]
                if id>0:
                    vars.append(c_new[id])
                if id<0:
                    id=-id
                    # vars.append(-c_new[id])
                    vars.append(c_new[id])


        for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
            hyperplane=hp[0][:-1]
            constant=hp[0][-1]
            constraint=np.dot(hyperplane,vars)+constant
            if hp[1]>0:
                model.addConstr(constraint >= z_var[(loop_number,k)][id_comb],name=f"HP, {loop_number,k,id_comb}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]<0:
                model.addConstr(constraint <= z_var[(loop_number,k)][id_comb],name=f"HP, {loop_number,k,id_comb}")
                # print(f"HP, {loop_number,k,id_comb}")
            elif hp[1]==0:
                model.addConstr(constraint == z_var[(loop_number,k)][id_comb],name=f"HP, {loop_number,k,id_comb}")

    model.update()

    mem_cs=[[],[],[],[]]
    for i in range(len(case["branch"])):
        mem_cs[0].append(s[i].LB)
        mem_cs[1].append(s[i].UB)
        mem_cs[2].append(c[i].LB)
        mem_cs[3].append(c[i].UB)
    
    ecart_tot=[]
    moy=0
    for i in range(len(case["branch"])):
        moy+=s[i].UB-s[i].LB
        moy+=c[i].UB-c[i].LB
    for i in range(1,len(s_new)+1):
        moy+=s_new[i].UB-s_new[i].LB
        moy+=c_new[i].UB-c_new[i].LB
    ecart_tot.append(moy)

    if tightbound==True:
        ecart_eps=1
    else:
        ecart_eps=0
    print(f"Ecart tot = {ecart_tot[-1]}")
    # while ecart_eps>0.01:
    for k in range(1):
        print(f"Bounding n° {len(ecart_tot)}")
        model.setParam('OutputFlag', 0)
        # cs_bounds=ObjBound.GetBounds(model,s,c,v2) #Desactivate too delete variable bounding
        # cs_bounds_new=ObjBound.GetBoundsnew(model,s_new,c_new,v2)
        model.setParam('OutputFlag', 1)

        
        cs_bounds=[[],[],[],[]] #Activate to delete variable bounding
        for i in range(len(case['branch'])):
            cs_bounds[0].append(s[i].LB)
            cs_bounds[1].append(s[i].UB)
            cs_bounds[2].append(c[i].LB)
            cs_bounds[3].append(c[i].UB)
        cs_bounds_new=[[],[],[],[]] 
        for i in s_new:
            cs_bounds_new[0].append(s_new[i].LB)
            cs_bounds_new[1].append(s_new[i].UB)
            cs_bounds_new[2].append(c_new[i].LB)
            cs_bounds_new[3].append(c_new[i].UB)

        
        for i in range(len(case['branch'])):
            for j in range(0,4):
                # print(cs_bounds[j][i])
                if cs_bounds[j][i]=="Error":
                    cs_bounds[j][i]=mem_cs[j][i]

        mem_cs=cs_bounds

        moy=0
        for i in range(len(case["branch"])):
            fbus=case["branch"][i]['fbus']
            tbus=case["branch"][i]['tbus']
            model.addConstr(cs_bounds[0][i]<=s[i],name=f"LB s{fbus,tbus,i}")
            model.addConstr(cs_bounds[1][i]>=s[i],name=f"LB s{fbus,tbus,i}")
            model.addConstr(cs_bounds[2][i]<=c[i],name=f"LB c{fbus,tbus,i}")
            model.addConstr(cs_bounds[3][i]>=c[i],name=f"LB c{fbus,tbus,i}")
            moy+=cs_bounds[1][i]-cs_bounds[0][i]
            moy+=cs_bounds[3][i]-cs_bounds[2][i]
        for i in range(len(s_new)):
            model.addConstr(cs_bounds_new[0][i]<=s_new[i+1],name=f"LB s_new{i}")
            model.addConstr(cs_bounds_new[1][i]>=s_new[i+1],name=f"LB s_new{i}")
            model.addConstr(cs_bounds_new[2][i]<=c_new[i+1],name=f"LB c_new{i}")
            model.addConstr(cs_bounds_new[3][i]>=c_new[i+1],name=f"LB c_new{i}")
            moy+=cs_bounds_new[1][i]-cs_bounds_new[0][i]
            moy+=cs_bounds_new[3][i]-cs_bounds_new[2][i]

        ecart_tot.append(moy)
        print(f"Ecart tot = {ecart_tot[-1]}")
        if ecart_tot[-2]!=0:
            ecart_eps=(ecart_tot[-2]-ecart_tot[-1])/ecart_tot[-2]
        else:
            ecart_eps=0
        # for loop_number in range(len(loops)):
        #     while model.getConstrByName(f"HP_loop_kk, {loop_number}")!=None:
        #         constr_to_remove = model.getConstrByName(f"HP_loop_kk, {loop_number}")
        #         model.remove(constr_to_remove)  
        # for loop_number in range(len(loops)): 
        #     for k in range(len(loops[loop_number])//2+1):
        #         combinations = math.comb(len(loops[loop_number]),2*k)
        #         for id_comb in range(combinations):
        #             while model.getConstrByName(f"HP, {loop_number,k,id_comb}")!=None:
        #                 constr_to_remove = model.getConstrByName(f"HP, {loop_number,k,id_comb}")
        #                 model.remove(constr_to_remove) 
        # print("Const bound deleted")
        

        if True:
            gen_var_on_loop(loops)
            ### this creates hyperplanes for c_i and s_i
            for loop_number in range(len(loops)): 
                for k in range(len(loops[loop_number])//2+1):
                    combinations = math.comb(len(loops[loop_number]),2*k)
                    for id_comb in range(combinations):
                        gen_constraints(loop_number,k,id_comb,cs_bounds,cs_bounds_new)

            ### this creates hyperplanes for c_kk
            for loop_number in range(len(loops)): 
                # z_var[(loop_number+len(loops))]
                var_bounds=[]
                for i in range(len(loops[loop_number])):
                    i0=loops[loop_number][i]
                    i1=loops[loop_number][(i+1)%len(loops[loop_number])]
                    if (i0,i1) in fbustbus_ID:
                        id=fbustbus_ID[i0,i1]
                        cLB=0
                        cUB=VMAX[i0]*VMAX[i1]
                        if id<len(case['branch']):
                            var_bounds.append([cLB,cUB])
                        else:
                            id=id-len(case["branch"])
                            # var_bounds.append([-cUB,-cLB])
                            var_bounds.append([cLB,cUB])
                    else:
                        id=fbustbus_ID_new[i0,i1]
                        cLB=0
                        cUB=VMAX[i0]*VMAX[i1]
                        if id>0:
                            var_bounds.append([cLB,cUB])
                        else:
                            # var_bounds.append([-cUB,-cLB])
                            var_bounds.append([cLB,cUB])
                HP=hyperplanes.generate_all_neighbor_hp_local_MAXMIN(var_bounds)
                vars=[]
                for i in range(len(loops[loop_number])):
                    i0=loops[loop_number][i]
                    i1=loops[loop_number][(i+1)%len(loops[loop_number])]
                    if (i0,i1) in fbustbus_ID:
                        id=fbustbus_ID[i0,i1]
                        if id<len(case["branch"]):
                            vars.append(c[id])
                        else:
                            id=id-len(case["branch"])
                            # vars.append(-c[id])
                            vars.append(c[id])
                    else:
                        id=fbustbus_ID_new[i0,i1]
                        if id>0:
                            vars.append(c_new[id])
                        else:
                            id=-id
                            # vars.append(-c_new[id])
                            vars.append(c_new[id])
                for hp in HP : #COMPUTE HYPERPLANES INEQUALITIES
                    hyperplane=hp[0][:-1]
                    constant=hp[0][-1]
                    constraint=np.dot(hyperplane,vars)+constant
                    if hp[1]>0:
                        model.addConstr(constraint >= z_var[(loop_number+len(loops))],name=f"HP_loop_kk, {loop_number}")
                        # print(f"HP_loop_kk, {loop_number}")
                    elif hp[1]<0:
                        model.addConstr(constraint <= z_var[(loop_number+len(loops))],name=f"HP_loop_kk, {loop_number}")
                        # print(f"HP_loop_kk, {loop_number}")
                        


            ### HEART, related loop constraints
            for loop_number in range(len(loops)): 
                constraint=0
                for k in range(len(loops[loop_number])//2):
                    combinations = math.comb(len(loops[loop_number]),2*k)
                    for id_comb in range(combinations):
                        constraint+=((-1)**k)*z_var[(loop_number,k)][id_comb]
                model.addConstr(constraint==z_var[len(loops)+loop_number],name=f"LINKING LOOP {loop_number}")
    
    plt.plot(ecart_tot)
    plt.ylabel("Total bound gap")
    objective = gp.QuadExpr()

    for i in range(len(case['gen'])):
        ca,cb,cc = case['gencost'][i]['costvector']
        objective += ca*Pg[i]*Pg[i]+cb*Pg[i]+cc

    model.setObjective(objective, GRB.MINIMIZE)

    if model.status == GRB.INFEASIBLE:
        print("Model is infeasible.")
        # You can try to find infeasibilities by calling the following:
        model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
        model.write("infeasible_model.ilp")

    print(f'Temps : {time.time()-duration}')

    Cuts=[]

    def JabrCuts(eps,pT,Const_to_add): 
        violation_amount=[]
        for i in range(len(case['branch'])):
            fbus=case['branch'][i]["fbus"]
            tbus=case['branch'][i]["tbus"]
            amount=c[i].X*c[i].X+s[i].X*s[i].X-v2[fbus].X*v2[tbus].X
            if amount>eps:
                violation_amount.append([amount,i,c[i].X,s[i].X,v2[fbus].X,v2[tbus].X])
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
                Const_to_add.append([4*c_*c[branch_i]+4*s_*s[branch_i]+(v2k-v2m-n0)*v2[fbus]+(v2m-v2k-n0)*v2[tbus],f"CutJabr_{fbus},{tbus},{branch_i}"])
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
    nb_gap=0
    gap=0
    
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
        print(f"{T}--------------------CUTS-------------------")
        
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
        nb_const[-1].append(i2Cuts(eps,pi2,Const_to_add))
        nb_const[-1].append(LimitCuts(eps,plim,Const_to_add))
        print(nb_const[-1])
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
            print("Model is infeasible.")
            model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
            model.write("infeasible_model.ilp")
            break
        obj=model.getObjective()
        val.append(obj.getValue())
        if(val[-1]-val[-2]<val[-2]*eps_ftol):
            r+=1
        else:
            r=0
        if val[-1]-initval>10**-5:
            init=False
        if val[-1]-initval>gap:
            gap=val[-1]-initval
            nb_gap=0
        else:
            nb_gap+=1
        if nb_gap>=5:
            T=Tlim
        print(val[-1])
    
        
    print(val)
    print(nb_const_tot)
    print(nb_const)

    fig, ax1 = plt.subplots()

    ax1.plot(temps,val, color='blue', label='Lower Bound')
    ax1.scatter(temps, val, color='blue', marker='|')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Lower Bound', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()
    ax2.plot(temps,nb_const_tot, color='red', label='#Added Constraints')
    ax2.scatter(temps, nb_const_tot, color='red', marker='|')
    ax2.plot(temps,nb_const_deleted, color='green', label='#Deleted Constraints')
    ax2.scatter(temps, nb_const_deleted, color='green', marker='|')
    ax2.set_ylabel('#Constraints', color='black')
    ax2.tick_params(axis='y', labelcolor='black')

    ax2.legend(loc='lower right')

    last_val = val[-1]
    ax1.text(0.95, 0.90, f'{last_val:.5e}', transform=ax1.transAxes, ha='right', va='top', color='blue')

    plt.title(f'Instance : {case["casename"]}')

    print("Jabr = violation")
    for i in range(len(case["branch"])):
        branch=case["branch"][i]
        fbus=branch["fbus"]
        tbus=branch["tbus"]
        aureol=c[i].X**2+s[i].X**2
        new_c=c[i].X/math.sqrt(aureol)
        new_s=s[i].X/math.sqrt(aureol)
        violation=abs(c[i].X**2+s[i].X**2-v2[fbus].X*v2[tbus].X)
        if violation<=1e-04:
            print(f"{i} {violation:.4e}")
        else:
            print(f"{i} {violation:.4e} WARNING")

    ### Computing violation quantity ###

    tot_val_viol=0
    for i in range(len(case['branch'])):
        
        branch_=case["branch"][i]
        fbus=branch_["fbus"]
        tbus=branch_["tbus"]
        G,B=create_admittance_matrix(branch_)
        k=branch_['fbus']
        m=branch_['tbus']
        aureol=c[i].X**2+s[i].X**2
        new_c=c[i].X/math.sqrt(aureol)
        new_s=s[i].X/math.sqrt(aureol)

        val_viol=0
        val_viol+=abs(P[i].X/baseMVA-(G[(k,k)]*v2[k].X+G[(k,m)]*new_c+B[(k,m)]*new_s))
        val_viol+=abs(P[i+N].X/baseMVA-(G[(m,m)]*v2[m].X+G[(m,k)]*new_c-B[(m,k)]*new_s))
        val_viol+=abs(Q[i].X/baseMVA-(-B[(k,k)]*v2[k].X-B[(k,m)]*new_c+G[(k,m)]*new_s))
        val_viol+=abs(Q[i+N].X/baseMVA-(-B[(m,m)]*v2[m].X-B[(m,k)]*new_c-G[(m,k)]*new_s))
        tot_val_viol+=val_viol
        # print(f"Branch {fbus},{tbus} = {val_viol:.3e}")
    print(f"TOT VIOLATION ={tot_val_viol}")

    for bus in range(len(case['bus'])):
        bus_number = int(case['bus'][bus]['bus_i'])
        Pd=case['bus'][bus]['Pd']  # Active load
        Qd=case['bus'][bus]['Qd']  # Reactive load
        GenBus=Bus_Gen.get(bus_number)
        Pg_bus=0
        Qg_bus=0
        if GenBus!=None:
            Pg_bus=sum(Pg[i].X for i in GenBus) # Active power gen on bus k
            Qg_bus=sum(Qg[i].X for i in GenBus) # Reactive power gen on bus k
        P_bus=0
        Q_bus=0
        BranchBus=Branch_Bus.get(bus_number)
        if BranchBus!=None:
            P_bus=sum(P[i].X for i in BranchBus) 
            Q_bus=sum(Q[i].X for i in BranchBus) 
        # print(f"Violation P : {abs(P_bus-(Pg_bus - Pd)):.3e}")
        # print(f"Violation P : {abs(Q_bus-(Qg_bus - Qd)):.3e}")


    # for loop_number in range(len(loops)): 
    #     print(f"z{loop_number}={z_var[loop_number+len(loops)].X}")
    #     constraint=0
    #     for k in range(len(loops[loop_number])//2):
    #         combinations = math.comb(len(loops[loop_number]),2*k)
    #         for id_comb in range(combinations):
    #             print(f"z{loop_number,k,id_comb}={z_var[(loop_number,k)][id_comb].X}")

    loops=init_loops

    sum_tot=0
    for loop in loops:
        sum_=0
        for i in range(len(loop)):
            index=fbustbus_ID[(loop[i],loop[(i+1)%len(loop)])]
            if index>=len(case["branch"]):
                index=index-len(case["branch"])
            aureol=c[index].X**2+s[index].X**2
            new_c=c[index].X/math.sqrt(aureol)
            new_s=s[index].X/math.sqrt(aureol)
            diffTheta=math.asin(new_s)
            sum_+=diffTheta
        sum_tot+=sum_
    print(f"Contraintes de loop : {sum_tot} grad")

    # for j in range(len(loops)):
    #     loop=loops[j]
    #     for k in range(len(loop)//2+1):
    #         z[(j,k)]=[]
    #         z_var[(j,k)]=[]
    #         combinations = list(itertools.combinations(loop, 2*k))
    #         for i in range(len(combinations)):
    #             z[(j,k)].append(combinations[i])
    #             z_var[(j,k)].append(model.addVar(lb=-GRB.INFINITY,ub=GRB.INFINITY,name=f"z{k,i}"))
    loops=loops_3
    bound_gap_tot=0
    for j in range(len(loops)):
        loop=loops[j]
        for k in range(len(loop)//2+1):
            for i in range(len(z[(j,k)])):
                tot_val=1
                s_list=z[(j,k)][i]
                list_str=""
                val_str=""
                for loop_id in loop:
                    index=loop.index(loop_id)
                    if (loop[index],loop[(index+1)%len(loop)]) in fbustbus_ID:
                        id=fbustbus_ID[loop[index],loop[(index+1)%len(loop)]]
                        if id<len(s):
                            if loop_id in s_list:
                                tot_val=tot_val*s[id].X
                                list_str+=str(f"s{id} ")
                                val_str+=f'{s[id].X:.4e} '
                            else:
                                tot_val=tot_val*c[id].X
                                list_str+=str(f"c{id} ")
                                val_str+=f'{c[id].X:.4e} '
                        else:
                            id=id-len(s)
                            if loop_id in s_list:
                                tot_val=tot_val*-s[id].X
                                list_str+=str(f"-s{id} ")
                                val_str+=f'{-s[id].X:.4e} '
                            else:
                                # tot_val=tot_val*-c[id].X
                                # list_str+=str(f"-c{id} ")
                                # val_str+=f'{-c[id].X:.4e} '
                                tot_val=tot_val*c[id].X
                                list_str+=str(f"c{id} ")
                                val_str+=f'{c[id].X:.4e} '
                    elif (loop[index],loop[(index+1)%len(loop)]) in fbustbus_ID_new:
                        id=fbustbus_ID_new[loop[index],loop[(index+1)%len(loop)]]
                        if id>0:
                            if loop_id in s_list:
                                tot_val=tot_val*s_new[id].X
                                list_str+=str(f"s_new{id} ")
                                val_str+=f'{s_new[id].X:.4e} '
                            else:
                                tot_val=tot_val*c_new[id].X
                                list_str+=str(f"c_new{id} ")
                                val_str+=f'{c_new[id].X:.4e} '
                        else:
                            id=-id
                            if loop_id in s_list:
                                tot_val=tot_val*-s_new[id].X
                                list_str+=str(f"-s_new{id} ")
                                val_str+=f'{-s_new[id].X:.4e} '
                            else:
                                # tot_val=tot_val*-c_new[id].X
                                # list_str+=str(f"-c_new{id} ")
                                # val_str+=f'{-c_new[id].X:.4e} '
                                tot_val=tot_val*c_new[id].X
                                list_str+=str(f"c_new{id} ")
                                val_str+=f'{c_new[id].X:.4e} '
                bound_gap_tot+=abs(z_var[(j,k)][i].X-tot_val)
                # if abs(z_var[(j,k)][i].X-tot_val)>0.001:
                #     print(list_str+"   "+val_str)
                #     print(f"Violation : {z_var[(j,k)][i].X-tot_val:.4f}={z_var[(j,k)][i].X:.4f}-{tot_val:.4f}")
    print(f"Bound Gap Violation = {bound_gap_tot}")

    if plot:
        plt.show()
    print(ecart_tot)
    return [val[-1],tot_val_viol]

if __name__=="__main__": 
    
    local_path="cases/"
    name="case9"
    name="case14"
    # name="case30"
    # name="case57"
    # name="case_ACTIVSg500"
    # name="case118"
    path=local_path+name+".m"
    case = conversion.conversion(path)
    print(Solve(case,True,True))