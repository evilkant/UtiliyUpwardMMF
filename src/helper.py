import math
import random
import numpy as np
import utility
import topo
import UIEWF

def scale_utility_functions(utility_functions,times):
    
    new_ufuncs=[]

    for ufunc in utility_functions:
        new_ufunc=[]
        for pt in ufunc:
            new_ufunc.append((pt[0],pt[1]*times))
        new_ufuncs.append(new_ufunc)
    
    return new_ufuncs


def get_bandwidth(u,utility_function):
    for i in range(len(utility_function)-1):
        p_1=utility_function[i]
        p_2=utility_function[i+1]
        if u>=p_1[1] and u<=p_2[1]:
            if p_1[1]==p_2[1]:
                return p_1[0]
            bw=p_1[0]+((u-p_1[1])/(p_2[1]-p_1[1]))*(p_2[0]-p_1[0])
            return bw
    return -1

def raise_utility(current_u,target_u,model):
    cms=model['commodities']
    util_funcs=model['utility_functions']
    splits=model['splits']
    pl_mat=model['pl_matrix']
    c=model['capacities']

    new_c=[c[i] for i in range(len(c))]

    for i in range(len(cms)):
        util_func=util_funcs[i]
        current_bw=get_bandwidth(current_u,util_func)
        target_bw=get_bandwidth(target_u,util_func)
        delta_bw=target_bw-current_bw
        for p in cms[i]:
            delta_bw_p=splits[p]*delta_bw
            for l in range(len(pl_mat[0])):
                if pl_mat[p][l]==1:
                    new_c[l]-=delta_bw_p
    
    return new_c


def is_feasible(current_u,target_u,model):
    new_c=raise_utility(current_u,target_u,model)
    for c in new_c:
        if c<0:
            return False
    return True


def exp_search(u_vec,curr_idx,model):
    vec_len=len(u_vec)
    feas_idx=curr_idx
    if curr_idx==0:
        i=1
    else:
        i=2*curr_idx
        if i>vec_len-1:
            i=vec_len-1

    while is_feasible(u_vec[curr_idx],u_vec[i],model):
        feas_idx=i
        #print(str(i)+"("+str(u_vec[i])+") is feasible")
        if i==vec_len-1:
            return vec_len-1,vec_len
        i=i*2
        if i>vec_len-1:
            i=vec_len-1
    
    return feas_idx,i


def bin_search(u_vec,curr_idx,feas_idx,infeas_idx,model):

    while infeas_idx-feas_idx>1:
        idx=(infeas_idx+feas_idx)//2
        if is_feasible(u_vec[curr_idx],u_vec[idx],model):
            feas_idx=idx
        else:
            infeas_idx=idx

    #print("largest feasible index is "+str(feas_idx))
    return feas_idx


def test():
    ufuncs=utility.read_ufuncs("D:/github/UtiliyUpwardMMF/data/ufuncs/uf_110.txt")
    scaled_ufuncs=scale_utility_functions(ufuncs,100)
    pl_mat,cms,c=topo.read_data("D:/github/UtiliyUpwardMMF/data/topologies/abilene_2_110.txt")
    initial_splits=UIEWF.exp_decay_splits(cms,pl_mat)
    model={'commodities':cms,'utility_functions':scaled_ufuncs,'splits':initial_splits,'pl_matrix':pl_mat,'capacities':c}

    #u_vec=utility.get_util_vec(scaled_ufuncs)
    #print(u_vec)
    #i_feas,i_infeas=exp_search(u_vec,0,model)
    #print("feasible index: "+str(i_feas))
    #print("infeasible index: "+str(i_infeas))
    #print(raise_utility(0,u_vec[4],model))
    #bin_search(u_vec,0,i_feas,i_infeas,model)
    #print(raise_utility(0,u_vec[3],model))


if __name__=='__main__':
    test()

