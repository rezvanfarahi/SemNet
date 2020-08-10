"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
=========================================================

"""

# License: BSD (3-clause)


print(__doc__)

import sys

sys.path.append('/imaging/local/software/mne_python/latest_v0.9')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
###

import numpy as np
import mne


###############################################################################
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/UVttest/time_resolved/' # root directory for your MEG data
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
data_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'
inv_fname = 'InverseOperator_EMEG-inv.fif'

# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all = ['meg10_0378/101209/', 
'meg10_0390/101214/',
'meg11_0026/110223/', 
'meg11_0050/110307/', 
'meg11_0052/110307/', 
'meg11_0069/110315/', 
'meg11_0086/110322/', 
'meg11_0091/110328/', 
'meg11_0096/110404/', 
'meg11_0101/110411/', 
'meg11_0102/110411/', 
'meg11_0112/110505/',
'meg11_0104/110412/',  
'meg11_0118/110509/', 
'meg11_0131/110519/', 
'meg11_0144/110602/', 
'meg11_0147/110603/', 


]

# subjects names used for MRI data
subjects=['SS1_Meg10_0378',
'SS2_Meg10_0390',
'SS3_Meg11_0026',
'SS4_Meg11_0050',
'SS5_Meg11_0052',
'SS6_Meg11_0069',
'SS9_Meg11_0086',
'SS10_Meg11_0091',
'SS12_Meg11_0096',
'SS13_Meg11_0101',
'SS14_Meg11_0102',
'SS15_Meg11_0112',
'SS16_Meg11_0104',
'SS18_Meg11_0118',
'SS19_Meg11_0131',
'SS20_Meg11_0144',
'SS21_Meg11_0147'
]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
reject = dict(eeg=120e-6, eog=150e-6, grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3

trange=np.arange(0.005,0.55,0.005)

tmin1=50
tstep1=100
vertices_to = [np.arange(10242), np.arange(10242)]
Matx_c=np.zeros((17,20484,len(trange)))
Matx_a=np.zeros((17,20484,len(trange)))
Matx_s=np.zeros((17,20484,len(trange)))

MTc=np.zeros((20484,len(trange)))
MTa=np.zeros((20484,len(trange)))
MTs=np.zeros((20484,len(trange)))
MPc=np.zeros((20484,len(trange)))
MPa=np.zeros((20484,len(trange)))
MPs=np.zeros((20484,len(trange)))

b='theta'
for ii, meg in enumerate(ll):
	print ii
	fname_c = data_path + meg +  'firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_concrete_meanresampled_5_550_10ms_5overlap_' +b
	stc_c = mne.read_source_estimate(fname_c)

	fname_a = data_path + meg +  'firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_abstract_meanresampled_5_550_10ms_5overlap_' +b
	stc_a = mne.read_source_estimate(fname_a)

	fname_s = data_path + meg +  'firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_subtract_meanresampled_5_550_10ms_5overlap_' +b
	stc_s = mne.read_source_estimate(fname_s)
	
	Matx_c[ii,:,:]=stc_c.copy().data.squeeze()
	Matx_a[ii,:,:]=stc_a.copy().data.squeeze()
	Matx_s[ii,:,:]=stc_s.copy().data.squeeze()

tmin1=5
tstep1=5
vertices_to = [np.arange(10242), np.arange(10242)]
for cnt in range(Matx_c.shape[2]):
    print cnt
    for cntt in range(Matx_c.shape[1]):
        
        MTc[cntt,cnt],MPc[cntt,cnt],h= mne.stats.permutation_t_test(Matx_c[:,cntt:cntt+1,cnt], tail=0, n_permutations=1024);
MTc[np.isnan(MTc)]=0
matx_stc = mne.SourceEstimate(MTc, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')

out_file=out_path + 'Perttest_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_concretebase_meanresampled_5_550_10ms_5overlap_' +b
matx_stc.save(out_file)
MPcl=-10*np.log10(MPc)
matx_stc = mne.SourceEstimate(MPcl, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')

out_file=out_path + 'Perttest_m10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_concretebase_meanresampled_5_550_10ms_5overlap_' +b
matx_stc.save(out_file)

for cnt in range(Matx_a.shape[2]):
    print cnt
    for cntt in range(Matx_a.shape[1]):
        
        MTa[cntt,cnt],MPa[cntt,cnt],h= mne.stats.permutation_t_test(Matx_a[:,cntt:cntt+1,cnt], tail=0, n_permutations=1024);
MTa[np.isnan(MTa)]=0
matx_stc = mne.SourceEstimate(MTa, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')

out_file=out_path + 'Perttest_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_abstractbase_meanresampled_5_550_10ms_5overlap_' +b
matx_stc.save(out_file)

MPal=-10*np.log10(MPa)
matx_stc = mne.SourceEstimate(MPal, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')

out_file=out_path + 'Perttest_m10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_abstractbase_meanresampled_5_550_10ms_5overlap_' +b
matx_stc.save(out_file)

for cnt in range(Matx_s.shape[2]):
    print cnt
    for cntt in range(Matx_s.shape[1]):
        
        MTs[cntt,cnt],MPs[cntt,cnt],h= mne.stats.permutation_t_test(Matx_s[:,cntt:cntt+1,cnt], tail=0, n_permutations=1024);
MTs[np.isnan(MTs)]=0
matx_stc = mne.SourceEstimate(MTs, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')

out_file=out_path + 'Perttest_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_subtract_meanresampled_5_550_10ms_5overlap_' +b
matx_stc.save(out_file)

MPsl=-10*np.log10(MPs)
matx_stc = mne.SourceEstimate(MPsl, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')

out_file=out_path + 'Perttest_m10logPval_firstMorphed_SemLoc_icaclean_degree_proportional_25_pca_ppc_subtract_meanresampled_5_550_10ms_5overlap_' +b
matx_stc.save(out_file)




"""
import matplotlib.pyplot as plt
in_file1=out_path + 'GrandAverage_SemLoc_icaclean_Subtract_SourcePower_ratio_m500_700_gamma'  
stc1=mne.read_source_estimate(in_file1)
in_file2=out_path + 'GrandAverage_SemLoc_icaclean_Subtract_SourcePower_ratio_m500_700_gamma' 
stc2=mne.read_source_estimate(in_file2) 
in_file3=out_path + 'GrandAverage_SemLoc_icaclean_Subtract_SourcePower_ratio_m500_700_gamma'  
stc3=mne.read_source_estimate(in_file3)
in_file4=out_path + 'GrandAverage_SemLoc_icaclean_Subtract_SourcePower_ratio_m500_700_gamma'  
stc4=mne.read_source_estimate(in_file4)
plt.plot(stc1.times[100:1100], stc1.data[:,100:1100].mean(axis=0), label='gamma')
plt.plot(stc2.times[100:1100], stc2.data[:,100:1100].mean(axis=0), label='gamma')
plt.plot(stc3.times[100:1100], stc3.data[:,100:1100].mean(axis=0), label='gamma')
plt.plot(stc4.times[100:1100], stc4.data[:,100:1100].mean(axis=0), label='gamma')
plt.xlabel('Time (ms)')
plt.ylabel('Power')
plt.legend()
plt.title('Mean source induced power Subtract')
plt.show()
"""

"""
out_file1=out_path + 'GrandAverage_SemLoc_icaclean_Concrete_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs1 = mne.read_source_estimate(out_file1)
out_file2=out_path + 'GrandAverage_SemLoc_icaclean_Concrete_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs2 = mne.read_source_estimate(out_file2)
stc_subtract=np.subtract(stcs1,stcs2)
out_file=out_path + 'GrandAverage_SemLoc_icaclean_Subtract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stc_subtract.save(out_file)
"""

"""
out_file1=out_path + 'GrandZscore_SemLoc_icaclean_Concrete_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs1 = mne.read_source_estimate(out_file1)
out_file2=out_path + 'GrandZscore_SemLoc_icaclean_Abstract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs2 = mne.read_source_estimate(out_file2)
stc_subtract=np.subtract(stcs1,stcs2)
out_file=out_path + 'GrandZscore_SemLoc_icaclean_Subtract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stc_subtract.save(out_file)
"""
