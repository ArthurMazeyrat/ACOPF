### Goal : given a cuboïd, generate all valid neighbors hyperplanes

import random
import itertools
from tqdm import tqdm
from typing import List, Tuple, Optional
import math
import numpy as np
import time

def generate_random_cuboid(dim,LB,UB,symmetric=False):
    L=[]
    for _ in range(dim):
        Lb=random.random()*(UB-LB)+LB
        if not symmetric:
            Ub=random.random()*(UB-LB)+LB
        else:
            Ub=-Lb
        L.append([min(Lb,Ub),max(Lb,Ub)])
    return L


def generate_neighbor_hyperplane(cuboid,point):
    """ An hyperplane is a list (pi_i) and c, reprensenting a fct PI(x)=sum(pi_i*x_i)+c
        The point is a binary = {0,1}^n """
    dim=len(cuboid)
    coeffs=[]
    p=math.prod(cuboid[j][point[j]] for j in range(dim))
    for i in range(dim):
        p_i=math.prod(cuboid[j][point[j]] for j in range(dim) if j!=i)
        coeffs.append(p_i)
    coeffs.append(-p*(dim-1))
    return coeffs

def evaluate(hp,cuboid,point):
    """Evaluate the hyperplane funct on point"""
    c=hp[-1]
    hp=hp[:-1]
    val_point=[]
    for i in range(len(point)):
        val_point.append(cuboid[i][0]+point[i]*(cuboid[i][1]-cuboid[i][0]))
    # val_point = [cuboid[i][point[i]] for i in range(len(cuboid))]
    # print(hp,val_point,np.dot(hp,val_point)+c)
    return np.dot(hp,val_point)+c

def prod(cuboid,point):
    return math.prod(cuboid[i][0]+point[i]*(cuboid[i][1]-cuboid[i][0]) for i in range(len(point)))
    # return math.prod(cuboid[i][point[i]] for i in range(len(point)))

def test_valid_hp(hp,cuboid):
    """Test wether the hyperplane cuts points or not"""
    dim=len(hp)-1
    combinations=[[int(x) for x in bin(i)[2:].zfill(dim)] for i in range(2**dim)]
    maxval=0
    minval=0
    for point in combinations:
       val=evaluate(hp,cuboid,point)-prod(cuboid,point)
       if abs(val)<1e-8:
           val=0
    #    print(hp,point,val)
       maxval=max(maxval,val)
       minval=min(minval,val)
       if maxval*minval<0:
           return False
    return True

def test_valid_random(hp,cuboid,n):
    """Test wether the hyperplane cuts points or not"""
    dim=len(hp)-1
    maxval=0
    minval=0
    for _ in range(n):
       point=[]
       for _ in range(dim):
           point.append(random.random())
       val=evaluate(hp,cuboid,point)-prod(cuboid,point)
       if abs(val)<1e-8:
           val=0
    #    print(hp,point,val)
       maxval=max(maxval,val)
       minval=min(minval,val)
       if maxval*minval<0:
           return False
    return True

def generate_all_neighbor_hp(cuboid):
    """ Generate all valid neighbor hp"""
    HP=[]
    dim=len(cuboid)
    combinations=[[int(x) for x in bin(i)[2:].zfill(dim)] for i in range(2**dim)]
    for point in combinations:
        hp=generate_neighbor_hyperplane(cuboid,point)
        if test_valid_hp(hp,cuboid):
            sign=sign_hp(hp,cuboid)
            HP.append([hp,sign,point])
    return HP

def maxmin(cuboid):
    MAX = math.prod(
    cuboid[i][0] if abs(cuboid[i][0]) >= abs(cuboid[i][1]) else cuboid[i][1]
    for i in range(len(cuboid))
    )
    if MAX<0:
        MIN=MAX
        MAX=0
    else:
        MIN=0
    empty=MAX
    for j in range(len(cuboid)):
        if empty == 0:
            MAX = max(MAX,math.prod(
                (cuboid[i][0] if abs(cuboid[i][0]) >= abs(cuboid[i][1]) else cuboid[i][1])
                if i != j else (cuboid[i][0] if abs(cuboid[i][0]) <= abs(cuboid[i][1]) else cuboid[i][1])
                for i in range(len(cuboid))
            ))
        if empty!=0:
            MIN = min(MIN,math.prod(
                (cuboid[i][0] if abs(cuboid[i][0]) >= abs(cuboid[i][1]) else cuboid[i][1])
                if i != j else (cuboid[i][0] if abs(cuboid[i][0]) <= abs(cuboid[i][1]) else cuboid[i][1])
                for i in range(len(cuboid))
            ))
    return MIN,MAX

def generate_all_neighbor_hp_MAXMIN(cuboid):
    """ Generate all valid neighbor hp using max/min property"""
    HP=[]
    MIN,MAX=maxmin(cuboid)
    dim=len(cuboid)
    combinations=[[int(x) for x in bin(i)[2:].zfill(dim)] for i in range(2**dim)]
    for point in tqdm(combinations,desc="Test...",ncols=100):
        prod_p=prod(cuboid,point)
        if prod_p==MIN or prod_p==MAX:
            hp=generate_neighbor_hyperplane(cuboid,point)
            HP.append([point,hp])
    return HP

def local_maxmin(cuboid,point):
    points=[point]
    for i in range(len(cuboid)):
        copy_point=point.copy()
        copy_point[i]=-point[i]+1
        points.append(copy_point)
    val=[prod(cuboid,points[i]) for i in range(len(points))]
    if max(val)==val[0] or min(val)==val[0]:
        return True
    return False

def sign_hp(hp,cuboid):
    dim=len(hp)-1
    combinations=[[int(x) for x in bin(i)[2:].zfill(dim)] for i in range(2**dim)]
    for point in combinations:
        val=evaluate(hp,cuboid,point)-prod(cuboid,point)
        if val<-1e-8:
           return -1
        elif val>1e-8:
            return 1
    return 0

def generate_all_neighbor_hp_local_MAXMIN(cuboid):
    """ Generate all valid neighbor hp using local max/min property"""
    HP=[]
    dim=len(cuboid)
    combinations=[[int(x) for x in bin(i)[2:].zfill(dim)] for i in range(2**dim)]
    # for point in tqdm(combinations,desc="Test...",ncols=100):
    for point in combinations:
        if local_maxmin(cuboid,point):
            hp=generate_neighbor_hyperplane(cuboid,point)
            if test_valid_hp(hp,cuboid):
                sign=sign_hp(hp,cuboid)
                HP.append([hp,sign])
    return HP


if __name__ == "__main__":
    cuboid=generate_random_cuboid(3,-10,10)
    # cuboid=generate_random_cuboid(6,-10,10,symmetric=True)
    for sublist in cuboid:
        print(" ".join(f"{x:.2f}" for x in sublist), end=" ")
    print("")
    # cuboid=[[1,5],[1,5],[1,-5]]
    # cuboid=[[-2,2],[-2,2]]

    time_=time.time()
    HP=[]

    HP=generate_all_neighbor_hp(cuboid)
    for hp in HP:
        hp_=hp[0]
        if test_valid_random(hp_,cuboid,1000)==False:
            print('ERREUR')
    for i in range(len(HP)):print(HP[i])
    print(f"Nb Hyperplanes = {len(HP)}") 
    print(f"duration = {(time.time()-time_):.2f}")
    time_=time.time()

    HP=generate_all_neighbor_hp_MAXMIN(cuboid)
    for hp in HP:
        hp_=hp[0]
        if test_valid_random(hp_,cuboid,1000)==False:
            print('ERREUR')
    # for i in range(len(HP)):print(HP[i])
    print(f"Nb Hyperplanes = {len(HP)}") 
    print(f"duration = {(time.time()-time_):.2f}")
    time_=time.time()

    HP=generate_all_neighbor_hp_local_MAXMIN(cuboid)
    for hp in HP:
        hp_=hp[0]
        if test_valid_random(hp_,cuboid,1000)==False:
            print('ERREUR')
    # for i in range(len(HP)):print(HP[i])
    print(f"Nb Hyperplanes = {len(HP)}") 
    print(f"duration = {(time.time()-time_):.2f}")
    time_=time.time()