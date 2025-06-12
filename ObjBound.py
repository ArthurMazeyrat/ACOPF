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



def GetBounds(model,s,c,v2):
    s_up_bounds=[]
    s_low_bounds=[]
    c_up_bounds=[]
    c_low_bounds=[]
    
    for type in range(2):
        for direction in range(2):
            for bounds in range(len(s)):
                objective = gp.QuadExpr()
                if type==0:
                    objective=s[bounds]
                else:
                    objective=c[bounds]
                if direction==0:
                    model.setObjective(objective, GRB.MINIMIZE)
                else:
                    model.setObjective(objective, GRB.MAXIMIZE)        
                model.optimize()
                
                if model.status == GRB.INFEASIBLE:
                    print("Model is infeasible.")
                    model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
                    model.write("infeasible_model.ilp")
                    val=["Error"]
                else:
                    obj=model.getObjective()
                    val=[obj.getValue()]

                if type==0:
                    if direction==0:
                        s_low_bounds.append(val[-1])
                    else:
                        s_up_bounds.append(val[-1])
                else:
                    if direction==0:
                        c_low_bounds.append(val[-1])
                    else:
                        c_up_bounds.append(val[-1])
    return s_low_bounds,s_up_bounds,c_low_bounds,c_up_bounds          
            
def GetBoundsnew(model,s_new,c_new,v2):
    s_up_bounds=[]
    s_low_bounds=[]
    c_up_bounds=[]
    c_low_bounds=[]
    
    for type in range(2):
        for direction in range(2):
            for bounds in range(1,len(s_new)+1):
                T=0
                Tlim=0
                pjabr=1
                pi2=0.15
                plim=1
                Tage=5
                Tftol=5
                eps=10E-8
                eps_par=10E-5/2
                eps_ftol=10E-5
                eps_slack=10E-5
                r=0
                objective = gp.QuadExpr()
                if type==0:
                    objective=s_new[bounds]
                else:
                    objective=c_new[bounds]
                if direction==0:
                    model.setObjective(objective, GRB.MINIMIZE)
                else:
                    model.setObjective(objective, GRB.MAXIMIZE)        
                model.optimize()
                
                if model.status == GRB.INFEASIBLE:
                    print("Model is infeasible.")
                    model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
                    model.write("infeasible_model.ilp")
                    val=["Error"]
                else:
                    obj=model.getObjective()
                    val=[obj.getValue()]

                if type==0:
                    if direction==0:
                        s_low_bounds.append(val[-1])
                    else:
                        s_up_bounds.append(val[-1])
                else:
                    if direction==0:
                        c_low_bounds.append(val[-1])
                    else:
                        c_up_bounds.append(val[-1])
    return s_low_bounds,s_up_bounds,c_low_bounds,c_up_bounds          

def GetBoundsTheta(model,case,s,c,v=[]):
    s_up_bounds=[]
    s_low_bounds=[]
    c_up_bounds=[]
    c_low_bounds=[]
    v_up_bounds=[]
    v_low_bounds=[]

    list_var=[]
    if v!=[]:
        list_var=[s,c,v]
    else:
        list_var=[s,c]
    
    for type in range(len(list_var)):
        for direction in range(2):
            for bounds in tqdm(list_var[type],ncols=100,desc=f" {type},{direction}"):
                objective = gp.QuadExpr()
                if type==0:
                    objective=s[bounds]
                elif type==1:
                    objective=c[bounds]
                elif type==2:
                    objective=v[bounds]
                if direction==0:
                    model.setObjective(objective, GRB.MINIMIZE)
                else:
                    model.setObjective(objective, GRB.MAXIMIZE)        
                model.optimize()
                
                if model.status == GRB.INFEASIBLE:
                    print("Model is infeasible.")
                    model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
                    model.write("infeasible_model.ilp")
                    val=["Error"]
                else:
                    obj=model.getObjective()
                    val=[obj.getValue()]

                if type==0:
                    if direction==0:
                        s_low_bounds.append(val[-1])
                    else:
                        s_up_bounds.append(val[-1])
                elif type==1:
                    if direction==0:
                        c_low_bounds.append(val[-1])
                    else:
                        c_up_bounds.append(val[-1])
                elif type==2:
                    if direction==0:
                        v_low_bounds.append(val[-1])
                    else:
                        v_up_bounds.append(val[-1])
    s_b={}
    compteur=0
    for bounds in s:
        s_b[bounds]=[s_low_bounds[compteur],s_up_bounds[compteur]]
        compteur+=1
    c_b={}
    compteur=0
    for bounds in c:
        c_b[bounds]=[c_low_bounds[compteur],c_up_bounds[compteur]]
        compteur+=1
    if v!=[]:
        v_b={}
        compteur=0
        for bounds in v:
            v_b[bounds]=[v_low_bounds[compteur],v_up_bounds[compteur]]
            compteur+=1
        return s_b,c_b,v_b 
    return s_b,c_b

def GetBoundsThetaTheta(model,theta):
    theta_up_bounds=[]
    theta_low_bounds=[]
    model.setParam('MIPGap', 1e-8)  # Tolerance for MIP gap (tighten this for better precision)
    model.setParam('FeasibilityTol', 1e-8)  # Tolerance for feasibility
    model.setParam('OptimalityTol', 1e-8)  # Tolerance for optimality
    model.setParam('IntFeasTol', 1e-8)  # Tolerance for integer feasibility
    for direction in range(2):
        for bounds in tqdm(theta,ncols=100,desc=f" {theta},{direction}"):
            objective = 0
            objective = gp.QuadExpr()
            objective  += theta[bounds]
            if direction==0:
                model.setObjective(objective, GRB.MINIMIZE)
            else:
                model.setObjective(objective, GRB.MAXIMIZE)        
            model.optimize()

            obj=model.getObjective()
            val=obj.getValue()

            if direction==0:
                theta_low_bounds.append(val)
            else:
                theta_up_bounds.append(val)
    theta_b={}
    compteur=0
    for bounds in theta:
        theta_b[bounds]=[theta_low_bounds[compteur],theta_up_bounds[compteur]]
        compteur+=1
    return theta_b

def GetBoundsThetaV(model,v):
    v_up_bounds=[]
    v_low_bounds=[]
    
    for direction in range(2):
        for bounds in tqdm(v,ncols=100,desc=f" {v},{direction}"):
            objective = gp.QuadExpr()
            objective=v[bounds]
            if direction==0:
                model.setObjective(objective, GRB.MINIMIZE)
            else:
                model.setObjective(objective, GRB.MAXIMIZE)        
            model.optimize()

            obj=model.getObjective()
            val=obj.getValue()

            if direction==0:
                v_low_bounds.append(val)
            else:
                v_up_bounds.append(val)
    v_b={}
    compteur=0
    for bounds in v:
        v_b[bounds]=[v_low_bounds[compteur],v_up_bounds[compteur]]
        compteur+=1
    return v_b

def GetBoundsThetai(model,theta,i):
    theta_b=[]

    for direction in range(2):
            objective = 0
            objective = gp.QuadExpr()
            objective  += theta[i]
            if direction==0:
                model.setObjective(objective, GRB.MINIMIZE)
            else:
                model.setObjective(objective, GRB.MAXIMIZE)        
            model.optimize()

            obj=model.getObjective()
            val=obj.getValue()

            theta_b.append(val)
    return theta_b

def GetBoundsVi(model,v,i):
    v_b=[]

    for direction in range(2):
            objective = 0
            objective = gp.QuadExpr()
            objective  += v[i]
            if direction==0:
                model.setObjective(objective, GRB.MINIMIZE)
            else:
                model.setObjective(objective, GRB.MAXIMIZE)        
            model.optimize()

            obj=model.getObjective()
            val=obj.getValue()

            v_b.append(val)
    return v_b

def GetBoundsThetaAll(model,case,theta):
    lambda_={}
    for i in range(len(case["branch"])):
        lambda_[i]=model.addVar(lb=0,ub=1,name=f"lambda {i}")
    Const_added=[]
    Const_added.append(model.addConstr(sum(lambda_[i] for i in lambda_ )== 1))

    objective = gp.QuadExpr()
    for i in range(len(case["branch"])):
        objective+=lambda_[i]*theta[i]
    model.setObjective(objective, GRB.MINIMIZE)
    model.update()
    model.optimize()
    if model.status == GRB.INFEASIBLE:
        print("Model is infeasible.")
        # You can try to find infeasibilities by calling the following:
        model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
        model.write("infeasible_model.ilp")
    obj=model.getObjective()
    Theta_MIN=obj.getValue()

    model.setObjective(objective, GRB.MAXIMIZE)
    model.update()
    model.optimize()
    obj=model.getObjective()
    Theta_MAX=obj.getValue()
    for const in Const_added:
        model.remove(const)
    return [Theta_MIN,Theta_MAX]

def GetBoundsVAll(model,case,v):
    lambda_={}
    for i in v:
        lambda_[i]=model.addVar(lb=0,ub=1,name=f"lambda {i}")
    Const_added=[]
    Const_added.append(model.addConstr(sum(lambda_[i] for i in lambda_ )== 1))

    objective = gp.QuadExpr()
    for i in v:
        objective+=lambda_[i]*v[i]
    model.setObjective(objective, GRB.MINIMIZE)
    model.update()
    model.optimize()
    if model.status == GRB.INFEASIBLE:
        print("Model is infeasible.")
        # You can try to find infeasibilities by calling the following:
        model.computeIIS()  # Compute Irreducible Inconsistent Subsystem
        model.write("infeasible_model.ilp")
    obj=model.getObjective()
    v_MIN=obj.getValue()

    model.setObjective(objective, GRB.MAXIMIZE)
    model.update()
    model.optimize()
    obj=model.getObjective()
    v_MAX=obj.getValue()
    for const in Const_added:
        model.remove(const)
    return [v_MIN,v_MAX]

if __name__=="__main__":
    name="case9"
    
    local_path="cases/"
    case = conversion.conversion(local_path)

    # print(GetBounds(case))