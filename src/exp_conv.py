import UIEWF
import matplotlib.pyplot as plt
import utility
import topo
import numpy as np


def convergence_abilene():

    plt.rcParams['font.family']='serif'
    plt.rcParams['font.size']=18
    plt.rcParams['font.weight']='bold'
    plt.rcParams['axes.labelweight']='bold'
    plt.rcParams['axes.labelsize']=18

    fig,ax=plt.subplots()
    fig.subplots_adjust(left=0.16,top=0.995,bottom=0.13,right=0.995)

    ufuncs=utility.read_ufuncs("../data/ufuncs/uf_110.txt")
    pl_mat,cms,c=topo.read_data('../data/topologies/abilene_2_110.txt')

    ufuncs = utility.scale_ufuncs(ufuncs)

    new_c = [10000.0 for i in range(len(c))]
    c = new_c

    pl_mat=np.matrix(pl_mat)

    initial_splits=UIEWF.uniform_splits(cms)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,it=20,th=0.00001)
    x=[i for i in range(21)][2:]
    ax.plot(x,diffs,label='uniform',marker='x',linestyle='solid')


    initial_splits=UIEWF.random_splits(cms)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,it=20,th=0.00001)
    ax.plot(x,diffs,label='random',marker='v',linestyle='solid')

    initial_splits=UIEWF.exp_decay_splits(cms,pl_mat)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,it=20,th=0.00001)
    ax.plot(x,diffs,label='len_decay',marker='.',linestyle='solid')

    initial_splits=UIEWF.exp_congestion_decay_splits(cms,pl_mat)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,it=20,th=0.00001)
    ax.plot(x,diffs,label='congestion_decay',marker='^',linestyle='solid')

    ax.set_xticks(np.arange(2,21,step=2))
    ax.set_yticks(np.arange(0.0,0.13,step=0.02))
    ax.set_ylabel('Difference')
    ax.set_xlabel('Iteration number')

    ax.grid(ls='--')
    ax.legend()

    plt.show()



def convergence_waxman_size():

    plt.rcParams['font.family']='serif'
    plt.rcParams['font.size']=18
    plt.rcParams['font.weight']='bold'
    plt.rcParams['axes.labelweight']='bold'
    plt.rcParams['axes.labelsize']=18

    fig,ax=plt.subplots()
    fig.subplots_adjust(left=0.16,top=0.995,bottom=0.13,right=0.995)


    size=[20,30]
    markers=['x','v','.','^']

    x=[i for i in range(2,21)]

    for i in range(len(size)):
        s=size[i]
        ufuncs=utility.read_ufuncs("../data/ufuncs/uf_"+str(s*(s-1))+".txt")
        pl_mat,cms,c=topo.read_data(
            "../data/topologies/waxman_"+str(s)+"_3_"+str(s*(s-1))+".txt")

        ufuncs = utility.scale_ufuncs(ufuncs)

        new_c = [10000.0 for i in range(len(c))]
        c = new_c

        pl_mat=np.matrix(pl_mat)


        initial_splits=UIEWF.exp_decay_splits(cms,pl_mat)
        cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,it=20,th=0.00001)
        ax.plot(x,diffs,label='len_decay',marker=markers[i],linestyle='dotted')

        initial_splits=UIEWF.exp_congestion_decay_splits(cms,pl_mat)
        cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,it=20,th=0.00001)
        ax.plot(x,diffs,label="waxman-"+str(s),marker=markers[i],linestyle='solid')

    ax.set_xticks(np.arange(2,21,step=2))
    ax.set_ylabel('Difference')
    ax.set_xlabel('Iteration number')

    ax.grid(ls='--')
    ax.legend()

    plt.show()

if __name__=='__main__':
    #convergence_abilene()
    convergence_waxman_size()