
import numpy as np
import os
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8})
plt.rcParams["font.family"] = "arial"
# plt.rcParams["mathtext.fontset"] = "dejavuserif"


def plot_sign(data1, data2, title='', pname='',cmap='jet'):
    x, y = np.mgrid[-90: 91 : 1, -45 : 46 : 1]
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=300, subplot_kw={'projection': '3d'})

    for ax, data, subtitle in zip(axes, [data1, data2], ["Co-pol","Cross-pol"]):
        ax.plot_surface(x, y, data, rstride=2, cstride=2, lw=0.2, edgecolor="black", cmap=cmap, alpha=0.3)
        ax.contourf(x, y, data, 100, zdir='z', offset=-0.5, cmap=cmap)
        ax.set_zlim([-0.5, 1])
        
        ax.set_xticks([-90, -45, 0, 45, 90])
        ax.set_xticklabels([f"{tick}°" for tick in [-90, -45, 0, 45, 90]])
        ax.tick_params(axis='y', pad=-3) 
        
        ax.set_yticks([-45, -30, -15, 0, 15, 30, 45])
        ax.set_yticklabels([f"{tick}°" for tick in [-45, -30, -15, 0, 15, 30, 45]])
        ax.tick_params(axis='x', pad=-4) 
        
        ax.set_zticks([0,0.2,0.4,0.6,0.8,1.0])
        ax.set_zticklabels([f"{tick}" for tick in [0,0.2,0.4,0.6,0.8,1.0]])
        ax.tick_params(axis='z', pad=-1) 
        
        ax.set_xlabel(r'Orientation, $\psi$ ', labelpad=-6)
        ax.set_ylabel(r'Ellipticity, $\tau$', labelpad=-6)
        ax.set_zlabel('    Normalized power', labelpad=-4)
        
        for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
            axis._axinfo['grid']['linewidth'] = 0.5
        
        ax.set_title(subtitle)

    plt.suptitle(title)
    plt.tight_layout()
    
    if pname:
        plt.savefig(pname, dpi=300)
    
    plt.show()
def prepare_data(S2):
    
    A0 = (1/4)*np.abs(S2[0,0]+S2[1,1])**2
    B0 = (1/4)*np.abs(S2[0,0]-S2[1,1])**2 + np.abs(S2[0,1])**2
    B_psi = (1/4)*np.abs(S2[0,0]-S2[1,1])**2 - np.abs(S2[0,1])**2
    C_psi = (1/2)*np.abs(S2[0,0]-S2[1,1])**2  
    D_psi = np.imag(S2[0,0]*np.conjugate(S2[1,1]))
    E_psi = np.real(np.conjugate(S2[0,1])*(S2[0,0]-S2[1,1]))
    F_psi = np.imag(np.conjugate(S2[0,1])*(S2[0,0]-S2[1,1]))
    G_psi = np.imag(np.conjugate(S2[0,1])*(S2[0,0]+S2[1,1]))
    H_psi = np.real(np.conjugate(S2[0,1])*(S2[0,0]+S2[1,1]))
    
    K_mat = np.array([
        [A0+B0,C_psi,H_psi,F_psi],
        [C_psi,A0+B_psi,E_psi,G_psi],
        [H_psi,E_psi,A0-B_psi,D_psi],
        [F_psi,G_psi,D_psi,-A0+B0]
        ]) 
    XX_POL = []
    cp_sign = np.zeros((181,91))
    xp_sign = np.zeros((181,91))    
    x, y = np.mgrid[-90: 91 : 1, -45 : 46 : 1]
    
    xx = 0
    for tilt in np.arange(-90,90,1):
        yy=0
        for elep in np.arange(-45,45,1):
            psi = tilt*np.pi/180
            tau = elep*np.pi/180
            
            S_co = np.array([1,np.cos(2*psi)*np.cos(2*tau),
                             np.sin(2*psi)*np.cos(2*tau), np.sin(2*tau)])
            stc = np.conjugate(S_co.T)
            temp_pow_c = np.matmul(stc,np.matmul(K_mat,S_co))
            S_cross = np.array([1,-np.cos(2*psi)*np.cos(2*tau),
                             -np.sin(2*psi)*np.cos(2*tau), -np.sin(2*tau)])
            # S_cr = [1;-cosd(2*tilt)*cosd(2*el);-sind(2*tilt)*cosd(2*el);-sind(2*el)];
            # Cross_pol(ii,jj) = (S_cr'*M_T*S_co); %cross-pol
            stx = np.conjugate(S_cross.T)
            temp_pow_x = np.matmul(stx,np.matmul(K_mat,S_co))
            
            XX_POL.append([tilt,elep,temp_pow_c,temp_pow_x])
            cp_sign[xx,yy] = temp_pow_c
            xp_sign[xx,yy] = temp_pow_x
            
            yy=yy+1  
            # print(xx,yy)
        xx=xx+1
    cp_sign[cp_sign==0]=np.nan
    xp_sign[xp_sign==0]=np.nan
    
    cp_sign = cp_sign/np.nanmax(cp_sign)
    xp_sign = xp_sign/np.nanmax(xp_sign)
    
    return cp_sign, xp_sign

def fp_sign(S2, title='',pname='',cmap='jet'):
    cp_sign, xp_sign = prepare_data(S2)
    plot_sign(cp_sign, xp_sign, title=title,pname=pname,cmap=cmap)
    

# array([[2248.6506   -5218.124j   ,   -4.0451455  -46.991398j],
#     [ -11.6204195  -27.771238j,  800.1831   -5830.6836j  ]],
#     dtype=complex64)