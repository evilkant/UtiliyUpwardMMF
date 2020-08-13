import math
import random
import numpy as np
import utility
import topo

import time


def Util_IEWF(pl_mat,cms,caps,initial_splits,ufuncs,it=20,th=0.01):

    path_num=pl_mat.shape[0]
    link_num=pl_mat.shape[1]
    cm_num=len(cms)

    link_usage=np.zeros(link_num)

    for p in range(path_num):
        for l in range(link_num):
            if pl_mat[p,l]==1:
                link_usage[l]=1

    path_alloc=np.zeros(path_num)    #flow value allocated to each path
    
    cm_alloc=np.zeros(cm_num)
    
    #th=0.0000001 #threshold for stopping iterating


    splits=initial_splits

    diffs=[]    #difference after each iteration

    current_diff=10000

    u=utility.get_util_vec(ufuncs)
    #u=utility.add_noise_to_uvec(u,9)
    #print(u)


    iteration=0
    while iteration<it and current_diff>th:
        iteration+=1
        print("iteration #"+str(iteration))

        start_time=time.time()
        new_path_alloc=wf(pl_mat,cms,splits,caps,link_usage,ufuncs,u)
        print("time is:")
        print(time.time()-start_time)

        new_cm_alloc=np.zeros(cm_num)
        for i in range(cm_num):
            for p in cms[i]:
                new_cm_alloc[i]+=new_path_alloc[p]

        alloc_diff=new_path_alloc-path_alloc

        if iteration>1:
            current_diff=np.linalg.norm(alloc_diff,2)/np.linalg.norm(path_alloc,2)
            diffs.append(current_diff)

        path_alloc=new_path_alloc
        cm_alloc=new_cm_alloc

        for i in range(cm_num):
            for p in cms[i]:
                splits[p]=path_alloc[p]/cm_alloc[i]

    return cm_alloc,diffs
            




def wf(pl_mat,cms,splits,caps,l_usage,utils,u):

    path_num=pl_mat.shape[0]
    link_num=pl_mat.shape[1]

    j=0

    #start_time=time.time()

    capacity=[caps[i] for i in range(link_num)]
    capacity=np.array(capacity)

    unsat_paths=[i for i in range(path_num)]
    unsat_paths=np.array(unsat_paths)

    link_usage=[l_usage[i] for i in range(len(l_usage))]
    link_usage=np.array(link_usage)

    #print("#1 time:")
    #print(time.time()-start_time)


    alloc=np.zeros(path_num)    #flow value allocated to each path


    unsat_cms=[i for i in range(len(cms))]
    unsat_cms=np.array(unsat_cms)

    cms_state=np.zeros(len(cms))


    weights=np.zeros(path_num)

    path_delta=np.zeros(path_num)

    waterlevel=0

    paths_state=np.zeros(path_num)

    cnt=0

    while len(unsat_paths)>0:

        if j==len(u)-1:
            break

        region=(u[j],u[j+1])

        start_time=time.time()

        for i in range(len(cms)):
            util=utils[i]

            idx=-1
            for k in range(len(util)-1):
                pt_1=util[k]
                pt_2=util[k+1]
                if pt_1[1] <= region[0] and region[1] <= pt_2[1]:
                    idx=k
                    
            weight=(util[idx+1][0]-util[idx][0])/(util[idx+1][1]-util[idx][1])
            for p in cms[i]:
                weights[p]=weight

        #print("#2 time:")
        #print(time.time()-start_time)


        start_time=time.time()

        path_delta=weights*splits
        link_shares=np.zeros(link_num)
        for l in range(link_num):
            l_delta=np.dot(path_delta,pl_mat[:,l])
            if l_delta>0:
                link_shares[l]=capacity[l]/l_delta
            else:
                link_shares[l]=math.inf

        min_share=np.amin(link_shares)

        min_delta=min_share

        
        if waterlevel+min_delta>=region[1]:
            min_delta=region[1]-waterlevel
            j+=1
            

        if math.isinf(min_delta):
            print("something went wrong\n")
            break

        waterlevel+=min_delta

        #print("#3 time:")
        #print(time.time()-start_time)

        start_time=time.time()

        alloc_delta=path_delta*min_delta
        alloc=alloc+alloc_delta

        for l in range(link_num):
            capacity[l]-=np.dot(alloc_delta,pl_mat[:,l])
            if link_usage[l]!=0 and capacity[l]<0.01:
                link_usage[l]=0
                for p in unsat_paths:
                    if pl_mat[p,l]==1:
                        paths_state[p]=1 #sat
                        splits[p]=0

        new_unsat_paths=[]
        for p in unsat_paths:
            if paths_state[p]==0:
                new_unsat_paths.append(p)
        unsat_paths=new_unsat_paths



        for cm_i in unsat_cms:
            unsat_num=0
            unsat_splits=0
            for p in cms[cm_i]:
                if paths_state[p]==0:
                    unsat_num+=1
                    unsat_splits+=splits[p]
            if unsat_num==0:
                cms_state[cm_i]=1 #sat
                continue
            if unsat_splits==0:
                s=1/unsat_num
                for p in cms[cm_i]:
                    if paths_state[p]==0:
                        splits[p]=s
            else:
                for p in cms[cm_i]:
                    if paths_state[p]==0:
                        splits[p]=splits[p]/unsat_splits

        new_unsat_cms=[]
        for cm_i in unsat_cms:
            if cms_state[cm_i]==0:#unsat
                new_unsat_cms.append(cm_i)
        unsat_cms=new_unsat_cms

        #print("#4 time:")
        #print(time.time()-start_time)

        cnt+=1

    print("cnt is "+str(cnt))
    return alloc


def exp_congestion_decay_splits(commodities,pl_mat):

    path_num=pl_mat.shape[0]
    link_num=pl_mat.shape[1]


    link_congestion=np.sum(np.array(pl_mat),axis=0)

    print(link_congestion)

    splits=np.zeros(path_num)

    for cm in commodities:
        
        max_congestion={}

        for p in cm:
            max_congestion[p]=0
            for l in range(link_num):
                if pl_mat[p,l]==1:
                    if link_congestion[l]>max_congestion[p]:
                        max_congestion[p]=link_congestion[l]

        sorted_p=sorted(max_congestion.items(),key=lambda x:x[1],reverse=False)

        total=0
        for i in range(len(sorted_p)):
            total+=1/(10**i)
        for i in range(len(sorted_p)):
            splits[sorted_p[i][0]]=(1/(10**i))/total

    return splits


def congestion_splits(commodities,pl_mat):

    path_num=pl_mat.shape[0]
    link_num=pl_mat.shape[1]


    link_congestion=np.sum(np.array(pl_mat),axis=0)

    print(link_congestion)

    splits=np.zeros(path_num)

    for cm in commodities:
        
        max_congestion={}

        total=0
        for p in cm:
            max_congestion[p]=0
            for l in range(link_num):
                if pl_mat[p,l]==1:
                    if link_congestion[l]>max_congestion[p]:
                        max_congestion[p]=link_congestion[l]
            total+=max_congestion[p]

        for p in cm:
            splits[p]=max_congestion[p]/total

    return splits


def random_splits(commodities):
    paths=[]
    for cm in commodities:
        for p in cm:
            paths.append(p)
    
    splits=np.zeros(len(paths))
    weights=np.zeros(len(paths))
    for cm in commodities:
        total=0
        for p in cm:
            weights[p]=random.randint(1,10)
            total+=weights[p]
        for p in cm:
            splits[p]=weights[p]/total
        
    return splits

def random_decay_splits(commodities):
    paths=[]
    for cm in commodities:
        for p in cm:
            paths.append(p)
    
    splits=np.zeros(len(paths))
    weights={}
    for cm in commodities:
        for p in cm:
            weights[p]=random.randint(1,10)

        sorted_p=sorted(weights.items(),key=lambda x:x[1],reverse=False)

        total=0
        for i in range(len(sorted_p)):
            total+=1/(10**i)
        for i in range(len(sorted_p)):
            splits[sorted_p[i][0]]=(1/(10**i))/total
        
    return splits

def uniform_splits(commodities):
    paths=[]
    for cm in commodities:
        for p in cm:
            paths.append(p)
    
    splits=np.zeros(len(paths))
    for cm in commodities:
        total=len(cm)
        for p in cm:
            splits[p]=1/total
    
    return splits


def exp_decay_splits(commodities,pl_mat):

    path_num=pl_mat.shape[0]
    link_num=pl_mat.shape[1]

    splits=np.zeros(path_num)

    for cm in commodities:
        p_len={}
        for p in cm:
            length=0
            for l in range(link_num):
                if pl_mat[p,l]==1:
                    length+=1
            p_len[p]=length
        sorted_p=sorted(p_len.items(),key=lambda x:x[1])
        total=0
        for i in range(len(sorted_p)):
            total+=1/(10**i)
        for i in range(len(sorted_p)):
            splits[sorted_p[i][0]]=(1/(10**i))/total

    return splits






if __name__=='__main__':
    ufuncs=utility.read_ufuncs("D:/github/UtiliyUpwardMMF/data/ufuncs/uf_110.txt")
    pl_mat,cms,c=topo.read_data("D:/github/UtiliyUpwardMMF/data/topologies/abilene_2_110.txt")

    #initial_splits=exp_decay_splits(cms,pl_mat)
    initial_splits=congestion_splits(cms,pl_mat)
    cm_alloc,diffs=Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,5)

    util_alloc=utility.flow2util(ufuncs,cm_alloc)

    print(sorted(util_alloc))
