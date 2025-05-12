

import polsartools as pst

def main():
	T3_folder = r'./sample_data/full_pol/T3'
	window_size=3
	
	#%% refined lee polarimetric Speckle filter
	pst.rlee(T3_folder,window_size=5)
	pst.utils.pauliRGB(T3_folder)

	#%%
	compact_c2 = r'./sample_data/compact_pol/C2_RHV'
	
	cprvi = pst.cprvi(compact_c2,chi_in=45,window_size=3)
	dcp = pst.dopcp(compact_c2,chi_in=45,window_size=3)
	ccp = pst.mf3cc(compact_c2,chi_in=45,window_size=3)
	socp = pst.misomega(compact_c2,chi_in=45,psi_in=0,window_size=3)
	print('compact pol')
	#%%
	full_T3 = r'./sample_data/full_pol/T3'
	
	mf3 = pst.mf3cf(full_T3,window_size=3)
	mf4 = pst.mf4cf(full_T3,window_size=3)
	dfp = pst.dopfp(full_T3,window_size=3)
	grvi = pst.grvi(full_T3,window_size=3)
	rvi = pst.rvifp(full_T3,window_size=3)
	prvi = pst.prvifp(full_T3,window_size=3)
	
	print('full pol')
	
	#%%
	
	dxp_C2 = r'./sample_data/dual_pol/C2_VVVH'
	
	dpr = pst.dprvi(dxp_C2,window_size=3)
	rvdp = pst.rvidp(dxp_C2,window_size=3)
	prvdp = pst.prvidp(dxp_C2,window_size=3)
	ddp = pst.dopdp(dxp_C2,window_size=3)
	
	print('dual cross-pol')

if __name__ == "__main__":
    main()