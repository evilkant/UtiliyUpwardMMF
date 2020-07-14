import math
import random
import numpy as np
import utility
import topo


def Util_IEWF(pl_mat,cms,caps,initial_splits,ufuncs_,it):

    ufuncs=[]

    for util in ufuncs_:
        new_util=[]
        for pt in util:
            new_util.append((pt[0],pt[1]*100))
        ufuncs.append(new_util)


    # paths of interest
    paths=[]
    for cm in cms:
        for p in cm:
            paths.append(p)

    link_usage=[0 for i in range(len(pl_mat[0]))]

    for p in paths:
        for l in range(len(pl_mat[0])):
            if pl_mat[p][l]==1:
                link_usage[l]=1

    alloc={}    #flow value allocated to each path

    for p in paths:
        alloc[p]=0
    
    cm_alloc=[0 for i in range(len(cms))]


    diff=[10 for i in range(len(paths))]    #difference between two iterations
    
    th=0.05 #threshold for stopping iterating

    splits=initial_splits

    diffs=[]    #difference after each iterationo

    d=1

    u=utility.get_util_vec(ufuncs)


    iteration=0
    while iteration<it and d>th:
        iteration+=1
        print("iteration #"+str(iteration))

        new_alloc=wf(pl_mat,cms,paths,splits,caps,link_usage,ufuncs,u)

        new_cm_alloc=[0 for i in range(len(cms))]
        for i in range(len(cms)):
            for p in cms[i]:
                new_cm_alloc[i]+=new_alloc[p]

        diff=[]
        for p in paths:
            diff.append(new_alloc[p]-alloc[p])

        alloc_ls=[alloc[p] for p in range(len(paths))]  #convert dict to ls

        d=np.linalg.norm(diff,2)/np.linalg.norm(alloc_ls,2)
        diffs.append(d)

        alloc=new_alloc
        cm_alloc=new_cm_alloc


        for i in range(len(cms)):
            for p in cms[i]:
                splits[p]=alloc[p]/cm_alloc[i]

    return cm_alloc,diffs[1:]
            




def wf(pl_mat,cms,paths,splits,caps,l_usage,utils,u):

    j=0


    capacity=[caps[i] for i in range(len(pl_mat[0]))]

    unsat_paths=[paths[i] for i in range(len(paths))]

    link_usage=[l_usage[i] for i in range(len(l_usage))]

    alloc={}    #flow value allocated to each path

    for p in paths:
        alloc[p]=0

    #sat_paths=[]
    sat_cms=[]

    weights={}

    waterlevel=0

    paths_state=[0 for i in range(len(pl_mat))]

    while len(unsat_paths)>0:

        if j==len(u)-1:
            break

        region=(u[j],u[j+1])

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


        min_delta=math.inf
        for l in range(len(pl_mat[0])):
            if link_usage[l]==0:
                continue
            prop=0
            for p in unsat_paths:
                if pl_mat[p][l]==1:
                    prop+=splits[p]*weights[p]
            delta=math.inf
            if prop>0:
                delta=capacity[l]/prop
            if delta<min_delta:
                min_delta=delta
        
        if waterlevel+min_delta>=region[1]:
            min_delta=region[1]-waterlevel
            j+=1
            

        if math.isinf(min_delta):
            print("something went wrong\n")
            break

        waterlevel+=min_delta

        for p in unsat_paths:
            alloc[p]=alloc[p]+min_delta*splits[p]*weights[p]
            for l in range(len(pl_mat[0])):
                if pl_mat[p][l]==1:
                    capacity[l]-=min_delta*splits[p]*weights[p]

        for l in range(len(pl_mat[0])):
            if link_usage[l]!=0 and capacity[l]<0.01:
                link_usage[l]=0
                #print("link "+str(l) +" saturates at "+str(waterlevel))
                for p in unsat_paths:
                    if pl_mat[p][l]==1:
                        paths_state[p]=1 #sat
                        splits[p]=0

                new_unsat_paths=[]
                for p in unsat_paths:
                    if paths_state[p]==0:
                        new_unsat_paths.append(p)
                unsat_paths=new_unsat_paths


        for cm_i in range(len(cms)):
            unsat_num=0
            unsat_splits=0
            for p in cms[cm_i]:
                if paths_state[p]==0:
                    unsat_num+=1
                    unsat_splits+=splits[p]
            if unsat_num==0:
                if cm_i not in sat_cms:
                    sat_cms.append(cm_i)
                    #print(str(cm_i)+' sats in '+str(waterlevel))
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

    return alloc


def congestion_splits(commodities,pl_mat):
    pass

def random_splits(commodities):
    paths=[]
    for cm in commodities:
        for p in cm:
            paths.append(p)
    
    splits={}
    weights={}
    for cm in commodities:
        total=0
        for p in cm:
            weights[p]=random.randint(1,10)
            total+=weights[p]
        for p in cm:
            splits[p]=weights[p]/total
        
    return splits

def uniform_splits(commodities):
    paths=[]
    for cm in commodities:
        for p in cm:
            paths.append(p)
    
    splits={}
    for cm in commodities:
        total=len(cm)
        for p in cm:
            splits[p]=1/total
    
    return splits

        

def exp_decay_splits(commodities,pl_mat):
    paths=[]
    for cm in commodities:
        for p in cm:
            paths.append(p)

    p_len={}
    splits={}
    for p in paths:
        length=0
        for l in range(len(pl_mat[0])):
            if pl_mat[p][l]==1:
                length+=1
        p_len[p]=length
    sorted_p=sorted(p_len.items(),key=lambda x:x[1])
    #print(sorted_p)
    for cm in commodities:
        total=0
        for p in cm:
            index=sorted_p.index((p,p_len[p]))
            #print(1/2**index)
            total+=1/((1.2)**index)
            #print("total is "+str(total))
        for p in cm:
            index=sorted_p.index((p,p_len[p]))
            splits[p]=1/((1.2)**index)/total
    #print(splits)
    return splits




if __name__=='__main__':
    ufuncs=utility.read_ufuncs("D:/github/UtiliyUpwardMMF/data/ufuncs/uf_110.txt")
    pl_mat,cms,c=topo.read_data("D:/github/UtiliyUpwardMMF/data/topologies/abilene_2_110.txt")

    initial_splits=uniform_splits(cms)
    cm_alloc,diffs=Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,5)

    util_alloc=utility.flow2util(ufuncs,cm_alloc)

    print(sorted(util_alloc))
