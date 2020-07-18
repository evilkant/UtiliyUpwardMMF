import matplotlib.pyplot as plt

def draw_demo():
    plt.rcParams['font.family']='serif'
    plt.rcParams['font.size']=18
    plt.rcParams['font.weight']='bold'
    plt.rcParams['axes.labelweight']='bold'
    plt.rcParams['axes.labelsize']=18

    fig,ax=plt.subplots()
    fig.subplots_adjust(left=0.1,top=0.995,bottom=0.13,right=0.995)

    x_1=[0,0,100]
    y_1=[0,1,2]

    x_2=[0,20,100,100]
    y_2=[0,1,1.8,2]

    ax.plot(x_1,y_1,label='low_priority',marker='v',linestyle='solid')
    ax.plot(x_2,y_2,label='high_priority',marker='^',linestyle='solid')

    #ax.set_xticks(np.arange(2,21,step=2))
    ax.set_yticks([0,1,2])
    ax.set_ylabel('Utility')
    ax.set_xlabel('Allocated rate(Mbps)')

    ax.grid(ls='--')
    ax.legend()

    plt.show()

if __name__=='__main__':
    draw_demo()