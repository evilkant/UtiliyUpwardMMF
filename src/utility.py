import math
import random

def elastic_ufunc(x,theta=0.1,beta=0):
    return 2/(1+math.exp(-theta*(x-beta)))-1

def realtime_ufunc(x,r=50):
    if x<r:
        return 0
    else:
        return 1


def delayadaptive_ufunc(x,theta=0.2,beta=50):
    return 1/(1+math.exp(-theta*(x-beta)))

def rateadaptive_ufunc(x,theta=0.2,beta=30,r=40):
    if x<r:
        return 0.8*(1/(1+math.exp(-theta*(x-beta))))
    else:
        return 0.8*(1/(1+math.exp(-theta*(x-beta))))+0.5*math.log10(x/r)


def realtime_approx(delta=2,r=50):
    x=[0,r-delta,r+delta,100]
    y=[0,0.0005*(r-delta),1,1]
    return x,y


def get_pl_ufunc(s_type):
    delta=10
    x=[]
    y=[]
    if s_type==1:
        for i in [0,1,2,3,4,6,10]:
            x_i=i*delta
            y_i=elastic_ufunc(i*delta)
            if abs(y_i)<0.01:
                y_i=0
            x.append(x_i)
            y.append(y_i)
        return x,y
    if s_type==2:
        return realtime_approx()
    if s_type==3:
        for i in [0,3,4,6,7,10]:
            x_i=i*delta
            y_i=delayadaptive_ufunc(i*delta)
            if abs(y_i)<0.01:
                y_i=0
            x.append(x_i)
            y.append(y_i)
        return x,y
    if s_type==4:
        for i in [0,1,2,3,4,5,6,10]:
            x_i=i*delta
            y_i=rateadaptive_ufunc(i*delta)
            if abs(y_i)<0.01:
                y_i=0
            x.append(x_i)
            y.append(y_i)
        return x,y


def generate_ufunc():
    s_type=random.randint(1,4)
    x,y=get_pl_ufunc(s_type)
    util=[]
    for i in range(len(x)):
        if abs(1-y[i])<0.01:
            util.append((x[i],1))
            if x[i]!=100:
                util.append((100,1))
            break
        util.append((x[i],y[i]))
    return util

def write_ufuncs(cms_num,fp):
    with open(fp,'w',encoding='utf-8') as f:
        while cms_num>0:
            util=generate_ufunc()
            for pt in util:
                f.write(str(pt[0])+':'+str(pt[1]))
                f.write(' ')
            f.write('\n')
            cms_num-=1
        f.write('\n')

def read_ufuncs(fp):
    with open(fp,'r',encoding='utf-8') as f:
        lines=f.readlines()
        i=0
        utils=[]
        while lines[i]!='\n':
            util=[]
            line=lines[i].strip('\n')
            line=line[:-1]
            pts=line.split(' ')
            for pt in pts:
                p=pt.split(':')
                tup=(float(p[0]),float(p[1]))
                util.append(tup)
            utils.append(util)
            i+=1
    return utils

def flow2util(ufuncs,cms_alloc):
    utils=[]
    for i in range(len(cms_alloc)):
        ufunc=ufuncs[i]
        alloc=cms_alloc[i]

        for idx in range(len(ufunc)-1):
            if ufunc[idx][0]<=alloc and alloc<=ufunc[idx+1][0]:
                k=(ufunc[idx+1][1]-ufunc[idx][1])/(ufunc[idx+1][0]-ufunc[idx][0])
                util=ufunc[idx][1]+(alloc-ufunc[idx][0])*k
                utils.append(util)
                break
            elif alloc>100:
                utils.append(1)
                break
        
    return utils

def get_util_vec(ufuncs):   # sort utilities
    utils=[]
    for ufunc in ufuncs:
        for pt in ufunc:
            utils.append(pt[1])
    sorted_utils=sorted(utils)
    u=[]
    u.append(sorted_utils[0])
    for i in range(1,len(sorted_utils)):
        if sorted_utils[i]!=sorted_utils[i-1]:
            u.append(sorted_utils[i])
    return u

def add_noise_to_uvec(u,times):
    u_len=len(u)
    n=u_len*times
    while n>0:
        a=random.randint(100,200)
        b=random.randint(1,100)
        u.append((a/b)%100)
        n-=1

    sorted_utils=sorted(u)

    new_u=[]
    new_u.append(sorted_utils[0])
    for i in range(1,len(sorted_utils)):
        if sorted_utils[i]!=sorted_utils[i-1]:
            new_u.append(sorted_utils[i])

    return new_u

if __name__ == "__main__":
    write_ufuncs(110,"D:/github/UtiliyUpwardMMF/data/ufuncs/uf_110.txt")
