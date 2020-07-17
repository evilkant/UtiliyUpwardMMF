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

    ufuncs=utility.read_ufuncs("D:/utility-mmf/scripts/dataset/convergence/uf_110.txt")
    pl_mat,cms,c=topo.read_data('D:/utility-mmf/scripts/dataset/convergence/abilene_2_110.txt')


    initial_splits=UIEWF.uniform_splits(cms)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,20)
    x=[i for i in range(21)][2:]
    ax.plot(x,diffs,label='uniform',marker='x',linestyle='solid')


    initial_splits=UIEWF.random_splits(cms)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,20)
    ax.plot(x,diffs,label='random',marker='v',linestyle='solid')

    initial_splits=UIEWF.exp_decay_splits(cms,pl_mat)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,20)
    ax.plot(x,diffs,label='exp_len_decay',marker='.',linestyle='solid')

    initial_splits=UIEWF.congestion_splits(cms,pl_mat)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,20)
    ax.plot(x,diffs,label='congestion',marker='^',linestyle='solid')

    ax.set_xticks(np.arange(2,21,step=2))
    ax.set_ylabel('Difference')
    ax.set_xlabel('Iteration number')

    ax.grid(ls='--')
    ax.legend()

    plt.show()

if __name__=='__main__':
    convergence_abilene()