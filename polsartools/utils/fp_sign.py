
import numpy as np
from polsartools.utils.plot_utils import pol_sign, pol_sign2d, poincare_plot,poincare_sphere
def prepare_data(S2):
    
    A0 = (1/4)*np.abs(S2[0,0]+S2[1,1])**2
    B0 = (1/4)*np.abs(S2[0,0]-S2[1,1])**2 + np.abs(S2[0,1])**2
    B_psi = (1/4)*np.abs(S2[0,0]-S2[1,1])**2 - np.abs(S2[0,1])**2
    # C_psi = (1/2)*np.abs(S2[0,0]-S2[1,1])**2 
    # di-hedral fix convert real part of vv to positive
    C_psi = (1/2) * np.abs(S2[0,0] - (np.real(S2[1,1]) * (1 if np.real(S2[1,1]) >= 0 else -1) + 1j * np.imag(S2[1,1])))**2
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
    for tilt in np.arange(-90,91,1):
        yy=0
        for elep in np.arange(-45,46,1):
            psi = tilt*np.pi/180
            tau = elep*np.pi/180
            
            S_co = np.array([1,
                             np.cos(2*psi)*np.cos(2*tau),
                             np.sin(2*psi)*np.cos(2*tau), 
                             np.sin(2*tau)
                             ])
            stc = np.conjugate(S_co.T)
            temp_pow_c = np.matmul(stc,np.matmul(K_mat,S_co))
            S_cross = np.array([1,
                                -np.cos(2*psi)*np.cos(2*tau),
                                -np.sin(2*psi)*np.cos(2*tau), 
                                -np.sin(2*tau)
                                ])
            stx = np.conjugate(S_cross.T)
            temp_pow_x = np.matmul(stx,np.matmul(K_mat,S_co))
            
            
            XX_POL.append([tilt,elep,temp_pow_c,temp_pow_x])
            cp_sign[xx,yy] = temp_pow_c
            xp_sign[xx,yy] = temp_pow_x
            
            yy=yy+1  
            # print(xx,yy)
        xx=xx+1
    # cp_sign[cp_sign==0]=np.nan
    # xp_sign[xp_sign==0]=np.nan
    
    cp_sign = cp_sign/np.nanmax(cp_sign)
    xp_sign = xp_sign/np.nanmax(xp_sign)
    
    return cp_sign, xp_sign



def fp_sign(S2=None, title='',pname='',cmap='jet',plotType = 1):
    if S2 is not None:
        cp_sign, xp_sign = prepare_data(S2)
    if plotType==4:
        poincare_sphere(pname=pname)
    elif S2 is not None and plotType==3:
        poincare_plot(cp_sign, xp_sign, title=title, pname=pname,cmap=cmap)
    elif S2 is not None and plotType==2:
        pol_sign2d(cp_sign, xp_sign, title=title, pname=pname,cmap=cmap)
    else:
        pol_sign(cp_sign, xp_sign, title=title, pname=pname,cmap=cmap)



# S2 = np.array([[1,0],
#                [0,1]])


# # fix this oriented dihedral
# S2 = np.array([[1, 0],
#                [0,-1]])

# S2 = np.array([[0,1],
#                [1,0]])


# S2 = np.array([[0,0],
#                [0,1]])

# S2 = np.array([[1,-1],
#                [-1,1]])

# array([[2248.6506   -5218.124j   ,   -4.0451455  -46.991398j],
#     [ -11.6204195  -27.771238j,  800.1831   -5830.6836j  ]],
#     dtype=complex64)


# fp_sign(S2, plotType = 3)
    

