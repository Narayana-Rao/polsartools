import numpy as np
import matplotlib.pyplot as plt
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = ['Arial']
# plt.rcParams.update({'font.size': 16})


def get_dp_halpha_bounds():
    datal = []
    for m in np.arange(0,1.01,0.01):
    
        data = np.array([[1,0],
                        [0,m]])
        evals, evecs = np.linalg.eig(data)
            
        evals[evals <0] = 0
        evals[evals >1] = 1
        eval_norm1 = (evals[1])/(evals[0] + evals[1])   
        eval_norm2 = (evals[0])/(evals[0] + evals[1])

        
        # # %Alpha 1
        eig_vec_r1 = np.real(evecs[0,1])
        eig_vec_c1 = np.imag(evecs[0,1])
        alpha1 = np.arccos(np.sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1))*180/np.pi
        
        # # %Alpha 2
        eig_vec_r2 = np.real(evecs[0,0])
        eig_vec_c2 = np.imag(evecs[0,0])
        alpha2 = np.arccos(np.sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2))*180/np.pi
        alpha_l = (eval_norm1*alpha1 + eval_norm2*alpha2)
        Hl = - eval_norm1*np.log10(eval_norm1)/np.log10(2) - eval_norm2*np.log10(eval_norm2)/np.log10(2)
        
        
        data = np.array([[m,0],
                        [0,1]])
        evals, evecs = np.linalg.eig(data)
            
            
        evals[evals <0] = 0
        evals[evals >1] = 1
        
        eval_norm1 = (evals[1])/(evals[0] + evals[1])
        eval_norm2 = (evals[0])/(evals[0] + evals[1])
        
        # # %Alpha 1
        eig_vec_r1 = np.real(evecs[0,1])
        eig_vec_c1 = np.imag(evecs[0,1])
        alpha1 = np.arccos(np.sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1))*180/np.pi
        
        # # %Alpha 2
        eig_vec_r2 = np.real(evecs[0,0])
        eig_vec_c2 = np.imag(evecs[0,0])
        alpha2 = np.arccos(np.sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2))*180/np.pi
        alpha_u = (eval_norm1*alpha1 + eval_norm2*alpha2)
        Hu = - eval_norm1*np.log10(eval_norm1)/np.log10(2) - eval_norm2*np.log10(eval_norm2)/np.log10(2)
        
       
        datal.append([alpha_l,Hl,alpha_u,Hu])
        
    datal = np.array(datal)
    datal = np.nan_to_num(datal)
    return datal


def halpha_plot_dp(h,alpha,path=None):
    data = get_dp_halpha_bounds()
    fs = 12
    fig, ax = plt.subplots(figsize=(3.5,3),dpi=300)
    plt.hexbin(h.flatten(), alpha.flatten(), gridsize=300, cmap='viridis',mincnt=1)

    plt.plot(data[:,1],data[:,0],'k-',linewidth=0.3)
    plt.plot(data[:,3],data[:,2],'k-',linewidth=0.3)
    plt.yticks(np.arange(0,100,10),np.arange(0,100,10),fontsize=fs)
    plt.xticks(np.round(np.arange(0,1.2,.2),2),np.round(np.arange(0,1.2,.2),2),fontsize=fs)
    plt.xlim([0,1])
    plt.ylim([0,90])
    plt.xlabel('Entropy, H',fontsize=fs)
    plt.ylabel(r'Alpha(degrees)',fontsize=fs)

    ax.tick_params(axis='both', which='both', direction='in',
                    top=True, bottom=True, left=True, right=True)
    
    ax.minorticks_on()
    ax.tick_params(axis='both', which='both', direction='in',
                   top=True, bottom=True, left=True, right=True,
                   length=4, width=0.5)
    ax.tick_params(axis='both', which='minor', length=2, width=0.2)

    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    if path is not None:
        plt.savefig(path,dpi=300,bbox_inches='tight')

