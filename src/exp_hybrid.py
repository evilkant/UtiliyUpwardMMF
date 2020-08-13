import Hybrid
import UIEWF
import UMMP
import matplotlib.pyplot as plt
import utility
import topo
import numpy as np

def test():

    plt.rcParams['font.family']='serif'
    plt.rcParams['font.size']=18
    plt.rcParams['font.weight']='bold'
    plt.rcParams['axes.labelweight']='bold'
    plt.rcParams['axes.labelsize']=18

    fig,ax=plt.subplots()
    fig.subplots_adjust(left=0.16,top=0.995,bottom=0.13,right=0.995)


    pl_mat, cms, c = topo.read_data(
        '../data/topologies/waxman_30_2_870.txt')
    ufuncs = utility.read_ufuncs(
        "../data/ufuncs/uf_870.txt")

    scaled_ufuncs = utility.scale_ufuncs(ufuncs)


    new_c = [10000.0 for i in range(len(c))]
    c = new_c

    pl_mat=np.matrix(pl_mat)

    x=[i for i in range(0,870)]

    #initial_splits=UIEWF.uniform_splits(cms)
    #cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,scaled_ufuncs,20)
    #util_alloc=utility.flow2util(ufuncs,cm_alloc)
    #ax.plot(x,sorted(util_alloc),label='uniform',linestyle='solid')


    #initial_splits=UIEWF.random_splits(cms)
    #cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,scaled_ufuncs,20)
    #util_alloc=utility.flow2util(ufuncs,cm_alloc)
    #ax.plot(x,sorted(util_alloc),label='random',linestyle='solid')

    #initial_splits=UIEWF.exp_decay_splits(cms,pl_mat)
    #cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,scaled_ufuncs,20)
    #util_alloc=utility.flow2util(ufuncs,cm_alloc)
    #ax.plot(x,sorted(util_alloc),label='len_decay',linestyle='solid')

    cm_alloc, paths_alloc, useful_when_error = UMMP.max_min_program(
        pl_mat, cms, c, scaled_ufuncs, len(cms))

    util_alloc=utility.flow2util(ufuncs,cm_alloc)
    ax.plot(x,sorted(util_alloc),label='Optimal',linestyle='solid')

    initial_splits=UIEWF.exp_congestion_decay_splits(cms,pl_mat)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,scaled_ufuncs,it=20,th=0.01)
    util_alloc=utility.flow2util(ufuncs,cm_alloc)
    ax.plot(x,sorted(util_alloc),label='U-IEWF',linestyle='dotted')

    cm_alloc,diffs=Hybrid.hybrid_alloc(pl_mat,cms,c,scaled_ufuncs,it=20,th=0.01,k=500)
    util_alloc=utility.flow2util(ufuncs,cm_alloc)
    ax.plot(x,sorted(util_alloc),label='Hybrid-500',linestyle='dashed')

    ax.set_xticks(np.arange(0,871,step=100))
    ax.set_ylabel('Utility')
    ax.set_xlabel('Commodity Index')

    ax.grid(ls='--')
    ax.legend()

    plt.show()

if __name__ == '__main__':
    test()