

import polsartools as pst

#%%
compact_c2 = r'./sample_data/compact_pol/C2_RHV'

cr = pst.cprvi(compact_c2,chi_in=45,window_size=3,  write_flag=0)
dcp = pst.dopcp(compact_c2,chi_in=45,window_size=3,  write_flag=0)
ccp = pst.mf3cc(compact_c2,chi_in=45,window_size=3,  write_flag=0)
socp = pst.mod_is_omega(compact_c2,chi_in=45,psi_in=0,window_size=3,  write_flag=0)
print('compact pol')
#%%
full_T3 = r'./sample_data/full_pol/T3'

cf3 = pst.mf3cf(full_T3,window_size=3,write_flag=0)
cf4 = pst.mf4cf(full_T3,window_size=3,write_flag=0)
dfp = pst.dopfp(full_T3,window_size=3,write_flag=0)
gr = pst.grvi(full_T3,window_size=3,write_flag=0)
rv = pst.rvifp(full_T3,window_size=3,write_flag=0)
pr = pst.prvifp(full_T3,window_size=3,write_flag=0)

print('full pol')

#%%

dxp_C2 = r'./sample_data/dual_pol/C2_VVVH'

dpr = pst.dprvi(dxp_C2,window_size=3,write_flag=0)
rvdp = pst.rvidp(dxp_C2,window_size=3,write_flag=0)
prvdp = pst.prvidp(dxp_C2,window_size=3,write_flag=0)
ddp = pst.dopdp(dxp_C2,window_size=3,write_flag=0)

print('dual cross-pol')

# dpc_T2 = 

