import UMMP
import UIEWF
import utility
import topo

def hybrid_alloc(pl_mat,cms,c,ufuncs,it,k):

    paths=[]
    for cm in cms:
        for p in cm:
            paths.append(p)
    
    cms_alloc,paths_alloc,active_cms=UMMP.max_min_program(pl_mat,cms,c,ufuncs,k)


    splits={}
    paths_alloc=list(paths_alloc)


    for cm in cms:
        flow_sum=0
        for p in cm:
            flow_sum+=paths_alloc[paths.index(p)]
        if flow_sum==0:
            for p in cm:
                splits[p]=1/len(cm)
        else:
            for p in cm:
                splits[p]=paths_alloc[paths.index(p)]/flow_sum

     
    cms_alloc,diffs=UIEWF.Util_IEWF(pl_mat,cms,c,splits,ufuncs,it)
    return cms_alloc,diffs

if __name__=='__main__':
    ufuncs=utility.read_ufuncs("D:/github/UtiliyUpwardMMF/data/ufuncs/uf_110.txt")
    pl_mat,cms,c=topo.read_data("D:/github/UtiliyUpwardMMF/data/topologies/abilene_2_110.txt")

    cm_alloc,diffs=hybrid_alloc(pl_mat,cms,c,ufuncs,it=5,k=50)

    util_alloc=utility.flow2util(ufuncs,cm_alloc)

    print(sorted(util_alloc))