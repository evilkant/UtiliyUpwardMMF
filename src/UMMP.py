from cvxopt import matrix,solvers
import numpy as np
import utility
import topo

import time

def max_min_program(pl_mat,cms,c,utils,k):

    paths=[]    #paths that are involved
    for cm in cms:
        for p in cm:
            paths.append(p)

    u=utility.get_util_vec(utils)

    u_index=0   #u_index is 0 at first, region [0,...)


    active_cms=[i for i in range(len(cms))] #all commodities are active commodity when starting
    cms_alloc=[0 for i in range(len(cms))]  #allocation is initialized to be 0
    frozen_cms=[]   #no frozen cms yet
    paths_alloc=[0 for i in range(len(paths))]  #path allocation is also 0

    threshold=0.001 #threshold for identifying blocked commodities

    iteration=1    #starts with iteration 1

    #the state of each cm, 0 being active
    cm_state=[0 for i in range(len(cms))]
    
    i=0

    while len(active_cms)>len(cms)-k:
        
        print("iteration #"+str(iteration))

        #start_time=time.time()

        sol=P(pl_mat,cms,c,active_cms,frozen_cms,cms_alloc,utils,u,u_index,paths)

        #print("LP solving takes time:")
        #print(time.time()-start_time)

        if sol['status']!='optimal':
            print("LP solver went wrong\n")
            return cms_alloc,paths_alloc,active_cms

        p_vars=sol['x']
        d_vars=sol['z']
        t=p_vars[-1]

        print('optimal t is:')
        print(t)


        if abs(t-u[u_index+1])<0.00001:  #the right end of region is reached
            u_index+=1
            if u_index!=len(u)-1:   #if utility is not maxium
                iteration+=1
                #print("Iteration takes time:")
                #print(time.time()-start_time)
                print("meet utility end")
                print(str(u[u_index])+" to "+str(u[u_index+1]))
                continue


        cms_alloc=[]
        for cm in cms:
            cm_alloc=0
            for p in cm:
                cm_alloc+=p_vars[paths.index(p)]
            cms_alloc.append(cm_alloc)

        if u_index==len(u)-1 or abs(t-100)<0.01:   #maximum utility is reached
            print('waterlevel reachs 100')
            break

        for i in range(len(active_cms)):
            if d_vars[i]>threshold:
                cm_i=active_cms[i]
                frozen_cms.append(cm_i)
                cm_state[cm_i]=1
                print(str(cm_i)+' is frozen at '+str(t))


        new_active_cms=[]
        for i in range(len(cms)):
            if cm_state[i]==0:
                new_active_cms.append(i)

        active_cms=new_active_cms

        paths_alloc=p_vars[:-1]

        iteration+=1
        
        #print("Iteration takes time:")
        #print(time.time()-start_time)

    return cms_alloc,paths_alloc,[]


def P(pl_mat,cms,caps,active_cms,frozen_cms,cms_alloc,utils,u,u_index,paths):

    # formulate a LP: min cx s.t. Ax<=b  and solve it

    path_num=len(paths)
    link_num=len(pl_mat[0])

    var_num=path_num+1    #number of variables [f1,f2,...,fn,t]

    #objective: max min{} => max t => min -t
    c=[0.0 for i in range(var_num)]
    c[-1]=-1.0

    #constraints
    G=[]
    h=[]


    start_time=time.time()
    for cm_i in active_cms:

        util=utils[cm_i]

        i=-1
        for index in range(len(util)-1):
            pt_1=util[index]
            pt_2=util[index+1]
            if pt_1[1] <= u[u_index] and u[u_index+1] <= pt_2[1]:
                i=index
                break
        k=(util[i+1][1]-util[i][1])/(util[i+1][0]-util[i][0])
        b=util[i][1]-k*util[i][0]

        a=[0.0 for i in range(var_num)]
        for p in cms[cm_i]:
            a[paths.index(p)]=-k
        a[-1]=1.0
        G.append(a)
        h.append(b)
    print("#1 takes time:")
    print(time.time()-start_time)


    start_time=time.time()
    for cm_i in frozen_cms:
        a=[0.0 for i in range(var_num)]
        for p in cms[cm_i]:
            a[paths.index(p)]=-1.0
        G.append(a)
        h.append(-cms_alloc[cm_i]) 

    print("#2 takes time:")
    print(time.time()-start_time)   


    start_time=time.time()
    for l in range(link_num):
        a=[0.0 for i in range(var_num)]
        for i in range(len(paths)):
            if pl_mat[paths[i]][l]==1:
                a[i]=1
        #for p in range(path_num):
        #    if pl_mat[p][l]==1:
        #        a[paths.index(p)]=1
        G.append(a)
        h.append(caps[l])

    print("#3 takes time:")
    print(time.time()-start_time)  

    start_time=time.time()
    for p_i in range(path_num):
        a=[0.0 for i in range(var_num)]
        a[p_i]=-1.0
        G.append(a)
        h.append(0.0)

    print("#4 takes time:")
    print(time.time()-start_time)  


    start_time=time.time()

    a=[0.0 for i in range(var_num)]
    a[-1]=1.0
    G.append(a)
    h.append(u[u_index+1])

    a=[0.0 for i in range(var_num)]
    a[-1]=-1.0
    G.append(a)
    h.append(-u[u_index])

    print("#5 takes time:")
    print(time.time()-start_time) 


    #start_time=time.time()

    #G_np=np.array(G).astype(float)
    #h_np=np.array(h).astype(float)
    #c_np=np.array(c).astype(float)
    #print("#6 takes time:")
    #print(time.time()-start_time)  

    start_time=time.time()
    G=matrix(G).T
    h=matrix(h)
    c=matrix(c)

    print("#7 takes time:")
    print(time.time()-start_time)  

    print(G.size)
    start_time=time.time()
    sol=solvers.lp(c,G,h,solver='glpk')
    print("solving LP takes time:")
    print(time.time()-start_time)

    return sol

def UMMP_alloc(ufuncs,pl_mat,cms,c):
    cms_alloc,paths_alloc,useful_when_error=max_min_program(pl_mat,cms,c,ufuncs,len(cms))
    util_alloc=utility.flow2util(ufuncs,cms_alloc)
    print(sorted(util_alloc))

if __name__=='__main__':
    ufuncs=utility.read_ufuncs("D:/github/UtiliyUpwardMMF/data/ufuncs/uf_110.txt")
    pl_mat,cms,c=topo.read_data("D:/github/UtiliyUpwardMMF/data/topologies/abilene_2_110.txt")
    UMMP_alloc(ufuncs,pl_mat,cms,c)