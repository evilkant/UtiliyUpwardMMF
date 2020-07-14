import networkx as nx
from itertools import islice
from itertools import permutations
import random
import numpy as np


def write_data(fl_mat,cms,c,fp):
    with open(fp,'w',encoding='utf-8') as f:
        for row in fl_mat:
            for i in range(len(row)):
                f.write(str(row[i]))
                if i<len(row)-1:
                    f.write(' ')
                else:
                    f.write('\n')
        f.write('\n')
        for cm in cms:
            for i in range(len(cm)):
                f.write(str(cm[i]))
                if i<len(cm)-1:
                    f.write(' ')
                else:
                    f.write('\n')
        f.write('\n')
        for i in range(len(c)):
            f.write(str(c[i]))
            if i<len(c)-1:
                f.write(' ')


def read_data(fp):

    fl_mat=[]
    cm=[]
    c=[]

    with open(fp,'r',encoding='utf-8') as f:
        lines=f.readlines()
        i=0

        while lines[i]=='\n' or lines[i].startswith('#'):
            i+=1

        # read fl matrix
        while lines[i]!='\n':
            line=lines[i].strip('\n')
            row=line.split(' ')
            fl_mat.append(row)
            i+=1
    
        while lines[i]=='\n' or lines[i].startswith('#'):
            i+=1

        #read commodities
        while lines[i]!='\n':
            line=lines[i].strip('\n')
            comm=line.split(' ')
            cm.append(comm)
            i+=1

        while lines[i]=='\n' or lines[i].startswith('#'):
            i+=1

        #read capacities
        line=lines[i].strip('\n')
        c=line.split(' ')

    return converter(fl_mat,'float'),converter(cm,'int'),converter(c,'float')

def converter(ls,tp):
    # convert list(possibly nested) of strings to list of numbers(float or int)
    res=[]
    for i in ls:
        if isinstance (i,str):
            if tp=='float':
                res.append(float(i))
            else:
                res.append(int(i))
        elif isinstance(i,list):
            res.append(converter(i,tp))
    return res





def k_shortest_paths(G,source,target,k,weight=None):
    return list(islice(nx.shortest_simple_paths(G,source,target,weight=weight),k))


def graph2input(G,path_num,cm_num):

    nodes=list(G.nodes())
    edges=list(G.edges())
    od_pairs=list(permutations(nodes,2))

    random.shuffle(od_pairs)

    demands=od_pairs[:cm_num]

    paths=[]
    comms=[]

    p_no=0
    for i in range(cm_num):
        comm=[]
        ps=k_shortest_paths(G,demands[i][0],demands[i][1],path_num)
        for p in ps:
            comm.append(p_no)
            paths.append(p)
            p_no+=1
        comms.append(comm)

    paths_e=[]  #get paths represented by a set of edges instead of by a sequence of nodes 

    for path in paths:
        path_e=[]
        for i in range(len(path)-1):
            e=(path[i],path[i+1])
            path_e.append(e)
        paths_e.append(path_e)

    pl_mat=np.zeros((len(paths),len(edges)))

    for i in range(len(paths)):
        for j in range(len(edges)):
            if edges[j] in paths_e[i]:
                pl_mat[i][j]=1
    
    return comms,pl_mat



### Abilene network


def write_abilene_data(path_num,cms_num):
    G=nx.read_gml("D:/github/UtiliyUpwardMMF/topologies/Abilene.gml")
    G=G.to_directed()
    cms,pl_mat=graph2input(G,path_num,cms_num)
    caps=[500 for i in range(len(G.edges()))]
    write_data(pl_mat,cms,caps,"D:/github/UtiliyUpwardMMF/data/topologies/abilene_"+str(path_num)+"_"+str(cms_num)+".txt")


### Waxman network

def wm_graph(n):
    G=nx.Graph()
    while True:
        G=nx.waxman_graph(n,0.55,0.55)
        if nx.is_connected(G):
            break
    return G.to_directed()

def write_random_data(node_num,path_num,cms_num):
    G=wm_graph(node_num)
    cms,pl_mat=graph2input(G,path_num,cms_num)
    caps=[500 for i in range(len(G.edges()))]
    write_data(pl_mat,cms,caps,"D:/github/UtiliyUpwardMMF/data/topologies/waxman_"+str(node_num)+'_'+str(path_num)+'_'+str(cms_num)+".txt")

if __name__=="__main__":
    #write_abilene_data(2,110)
    write_random_data(30,2,30*29)