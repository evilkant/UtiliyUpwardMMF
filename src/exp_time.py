import UIEWF
import UIEWF_bin
import matplotlib.pyplot as plt
import utility
import topo
import time
import numpy as np
import UMMP


def running_time():

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = 18
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.labelsize'] = 18

    fig, ax = plt.subplots()
    fig.subplots_adjust(left=0.12, top=0.99, bottom=0.13, right=0.995)

    pl_mat, cms, c = topo.read_data(
        'D:/utility-mmf/scripts/dataset/convergence/waxman_20_4_380_1.txt')

    new_c = [10000.0 for i in range(len(c))]
    c = new_c

    cm_nums = np.arange(100, 381, 40)

    UIEWF_times = [0 for i in range(len(cm_nums))]
    UIEWF_bin_times = [0 for i in range(len(cm_nums))]

    for i in range(3, 4):

        ufuncs = utility.read_ufuncs(
            "D:/utility-mmf/scripts/dataset/convergence/uf_380_time_"+str(i)+".txt")

        for i in range(len(cm_nums)):
            cm_num = cm_nums[i]
            comms = cms[:cm_num]
            util_functions = ufuncs[:cm_num]
            start_time = time.time()
            initial_splits = UIEWF.exp_decay_splits(comms, pl_mat)
            cm_alloc, diffs = UIEWF.Util_IEWF(
                pl_mat, comms, c, initial_splits, util_functions, 5)
            ummf_time = time.time()-start_time
            UIEWF_times[i] += ummf_time

        for i in range(len(cm_nums)):
            cm_num = cm_nums[i]
            comms = cms[:cm_num]
            util_functions = ufuncs[:cm_num]
            start_time = time.time()
            initial_splits = UIEWF.exp_decay_splits(comms, pl_mat)
            cm_alloc, diffs = UIEWF_bin.Util_IEWF(
                pl_mat, comms, c, initial_splits, util_functions, 5)
            ummf_time = time.time()-start_time
            UIEWF_bin_times[i] += ummf_time

    ax.set_ylabel('Time(s)')
    ax.set_xlabel('Number of commodities')

    times = [t/1 for t in UIEWF_times]
    bin_times = [t/1 for t in UIEWF_bin_times]

    ax.plot(cm_nums, times, label='UIEWF', marker='.')
    ax.plot(cm_nums, bin_times, label='UIEWF-bin', marker='v')

    plt.xticks(np.arange(100, 381, 40))

    plt.grid(ls='--')

    plt.legend()

    plt.show()


if __name__ == '__main__':
    # running_time()
    pl_mat, cms, c = topo.read_data(
        '../data/topologies/waxman_50_2_2450.txt')
    ufuncs = utility.read_ufuncs(
        "../data/ufuncs/uf_2450.txt")
    scaled_ufuncs = utility.scale_ufuncs(ufuncs)


    new_c = [10000.0 for i in range(len(c))]
    c = new_c

    pl_mat=np.matrix(pl_mat)

    start_time = time.time()
    cms_alloc, paths_alloc, useful_when_error = UMMP.max_min_program(
        pl_mat, cms, c, scaled_ufuncs, len(cms))
    tot_time = time.time()-start_time
    print("UMMP total time is "+str(tot_time))


    start_time = time.time()
    initial_splits=UIEWF.exp_decay_splits(cms,pl_mat)
    cm_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,initial_splits,ufuncs,50)
    tot_time = time.time()-start_time
    print("UIEWF total time is "+str(tot_time))
    #print(diffs)
    #print(sorted(cms_alloc))
