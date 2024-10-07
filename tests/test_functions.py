

import polsartools as pst

#%%
compact_c2 = r'./sample_data/compact_pol/C2_RHV'

cr = pst.polsar.cp.cprvi(compact_c2,chi_in=45,window_size=3,  write_flag=False)
dcp = pst.polsar.cp.dopcp(compact_c2,chi_in=45,window_size=3,  write_flag=False)
ccp = pst.polsar.cp.mf3cc(compact_c2,chi_in=45,window_size=3,  write_flag=False)
socp = pst.polsar.cp.misomega(compact_c2,chi_in=45,psi_in=0,window_size=3,  write_flag=False)
print('compact pol')
#%%
full_T3 = r'./sample_data/full_pol/T3'

cf3 = pst.polsar.fp.mf3cf(full_T3,window_size=3,write_flag=False)
cf4 = pst.polsar.fp.mf4cf(full_T3,window_size=3,write_flag=False)
dfp = pst.polsar.fp.dopfp(full_T3,window_size=3,write_flag=False)
gr = pst.polsar.fp.grvi(full_T3,window_size=3,write_flag=False)
rv = pst.polsar.fp.rvifp(full_T3,window_size=3,write_flag=False)
pr = pst.polsar.fp.prvifp(full_T3,window_size=3,write_flag=False)

print('full pol')

#%%

dxp_C2 = r'./sample_data/dual_pol/C2_VVVH'

dpr = pst.polsar.dxp.dprvi(dxp_C2,window_size=3,write_flag=False)
rvdp = pst.polsar.dxp.rvidp(dxp_C2,window_size=3,write_flag=False)
prvdp = pst.polsar.dxp.prvidp(dxp_C2,window_size=3,write_flag=False)
ddp = pst.polsar.dxp.dopdp(dxp_C2,window_size=3,write_flag=False)

print('dual cross-pol')

