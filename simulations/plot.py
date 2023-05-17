import matplotlib.pyplot as plt
import numpy as np
import sys

def plot_interdays(sim,lim=(1e-7*0.9,1./0.9),title='',figax=None,alpha=1.0,stylepath=None,d0=1):

    ts_eod=sim[1]
    
    D=len(ts_eod)

    if figax==None:
        # =plt.figure()
        fig,ax=plt.subplots()
    else:   
        fig,ax=figax

    plt.plot(range(1,len(ts_eod)+1),np.sum(ts_eod[:,1:],axis=-1)/np.sum(ts_eod,axis=-1),
    color='black',linestyle='solid', linewidth=.2, marker=None) # hosts
    ax.plot(range(1,len(ts_eod)+1),np.sum(ts_eod[:,1:-1],axis=-1)/np.sum(ts_eod,axis=-1),
    color='green',fillstyle='none',linestyle='solid',linewidth=.2, marker=None,alpha=alpha) # het
    ax.plot(range(1,len(ts_eod)+1),ts_eod[:,-1]/np.sum(ts_eod,axis=-1),
    color='violet',linestyle='solid', linewidth=.2, marker=None,alpha=alpha) # hom mt
    ax.set_yscale('log')
    ax.set_ylim(lim)
    ax.set_xlabel('Days')
    ax.set_ylabel('Frequency')
    ax.set_title(title)
    # ax.grid(b=True, which='major', axis='both')
    #plt.close()
    return figax

def plot_interdays_CI(sims,figax=None,ignorezeros=False):
    if figax==None:
        fig,ax=plt.subplots()
    else:
        fig,ax=figax

    c={ 'het': 'green',
        'hom': 'violet',
    }
    for genotype in ['het','hom']:
        logfreqs=[]
        for i in range(len(sims)):
            ts_eod=sims[i][1]
            with np.errstate(divide='ignore'):
                if genotype=='het':
                    logfreq=np.log10(np.sum(ts_eod[:,1:-1],axis=-1)/np.sum(ts_eod,axis=-1))
                if genotype=='hom':
                    logfreq=np.log10(ts_eod[:,-1]/np.sum(ts_eod,axis=-1))
            logfreqs.append(logfreq)
        logfreqs=np.array(logfreqs)
        if ignorezeros:
            logfreqs[np.isinf(logfreqs)]=np.nan
        else:
            logfreqs[np.isinf(logfreqs)]=-1000

        q=np.nanquantile(logfreqs,q=[(1-.95)/2,(1-.68)/2,1-(1-.68)/2,1-(1-.95)/2],axis=0)

        ax.fill_between(range(1,len(q[0])+1),y1=q[0],y2=q[-1],color=c[genotype],alpha=.2,edgecolor="none")
        ax.fill_between(range(1,len(q[0])+1),y1=q[1],y2=q[-2],color=c[genotype],alpha=.5,edgecolor="none")
    
    ax.grid(linewidth=0.25)
    # ax.set_xlabel('Days')
    # ax.set_ylabel('Frequency')

def plot_exp_medianstd(data,figax=None):
    if figax==None:
        fig,ax=plt.subplots()
    else:
        fig,ax=figax

    c={'hetero_freq':'green','homo_freq':'violet'}

    for genotype in ['hetero_freq','homo_freq']:
        med=data.loc[:,['t',genotype]].groupby(['t']).mean()
        logmed=data.loc[:,['t',genotype]].groupby(['t']).mean().applymap(lambda x: np.log10(x) if x>0 else -1000)
        std=data.loc[:,['t',genotype]].groupby(['t']).std()
        std=0.434*std/med
        std[np.isnan(std)]=0;std[np.isinf(std)]=0;#std[std>1]=0;
        ax.errorbar(range(1,len(med)+1),logmed.to_numpy()[:,0],yerr=std.to_numpy()[:,0],
            color=c[genotype], marker='.',linewidth=.5,markersize=2.,
            capsize=0.75, elinewidth=.5, markeredgewidth=.5)


    return fig

def plot_exp_timetoext(data,figax=None):
    if figax==None:
        fig,ax=plt.subplots()
    else:
        fig,ax=figax

    print(data)

    # hist=np.sum(host_freq==0,axis=0)/6
    # ax.plot(range(1,30+1))
