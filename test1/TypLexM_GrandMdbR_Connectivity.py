"""
=========================================================
TF representation of SemLoc for frequency bands in source labels
=========================================================

"""
# Authors: Alexandre Gramfort <gramfort@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)


print(__doc__)

import sys
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import numpy as np
import mne


###############################################################################
data_path = '/imaging/rf02/TypLexMEG/' # root directory for your MEG data
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/GrandAverage/connectivity/'
subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files are
event_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/fsaverage/label'


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
reject = dict(eeg=120e-6,grad=200e-12, mag=4e-12)
event_ids = {'cncrt_wrd': 1, 'abs_wrd': 2}
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s
frequencies = np.arange(8, 45, 1)
n_cycles = frequencies / float(7)
n_cycles[frequencies<=20] = 2
n_cycles[frequencies<=20] += np.arange(0,1,0.077)
    # n_cycles[np.logical_and((frequencies>10), (frequencies<=20))] = 3
labellist = ['leftATL-lh','rightATL-rh']#, 'rightATL-rh']#, 'medialtempright-rh','medialtempleft-lh']
Matx_c=np.zeros((17,20484,4))
Matx_a=np.zeros((17,20484,4))
Matx_cb=np.zeros((17,20484,4))
Matx_ab=np.zeros((17,20484,4))
trange='_250_450_'
conn_method='PPC'
for label_name in labellist:
	
	for ii, meg in enumerate(ll):
		print ii
		fname_c = data_path + meg +  'firstMorphed_SemLoc_icaclean_Concrete_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3]
		stc_c = mne.read_source_estimate(fname_c)
		fname_cb = data_path + meg +  'firstMorphed_SemLoc_icaclean_ConcreteBaseline_' + conn_method + '_ThetaAlphaBetaGamma_m300_m100_' + label_name[0:-3]
		stc_cb = mne.read_source_estimate(fname_cb)
  
		fname_a = data_path + meg +  'firstMorphed_SemLoc_icaclean_Abstract_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3]
		stc_a = mne.read_source_estimate(fname_a)
		fname_ab = data_path + meg +  'firstMorphed_SemLoc_icaclean_AbstractBaseline_' + conn_method + '_ThetaAlphaBetaGamma_m300_m100_' + label_name[0:-3]
		stc_ab = mne.read_source_estimate(fname_ab)
  
		stc_s=np.subtract(stc_c,stc_a)   
		stccb_s=np.subtract(stc_c,stc_cb)           
		stcab_s=np.subtract(stc_a,stc_ab)           
        
		fname_s = data_path + meg +  'firstMorphed_SemLoc_icaclean_Subtract_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3]
		stc_s.save(fname_s)           	
		Matx_c[ii,:,:]=stc_c.copy().data
		Matx_a[ii,:,:]=stc_a.copy().data
		Matx_cb[ii,:,:]=stc_cb.copy().data
		Matx_ab[ii,:,:]=stc_ab.copy().data
	Matx=Matx_c-Matx_a
	M_mean=Matx.mean(axis=0)
	#M_range=Matx.max(axis=0)-Matx.min(axis=0)
	#MT=M_mean/M_range
	tmin1=1
	tstep1=1
	vertices_to = [np.arange(10242), np.arange(10242)]
	matx_stc = mne.SourceEstimate(M_mean, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
	out_file=out_path + 'Grandaverage_firstMorphed_SemLoc_icaclean_Subtract_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3]
	matx_stc.save(out_file)

	vertices_to = [np.arange(10242), np.arange(10242)]
	Matx=Matx_c-Matx_cb
	M_mean=Matx.mean(axis=0)
	matx_stc = mne.SourceEstimate(M_mean, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
	out_file=out_path + 'Grandaverage_firstMorphed_SemLoc_icaclean_ConcreteBase_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3]
	matx_stc.save(out_file)

	Matx=Matx_a-Matx_ab
	M_mean=Matx_a.mean(axis=0)
	vertices_to = [np.arange(10242), np.arange(10242)]
	matx_stc = mne.SourceEstimate(M_mean, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
	out_file=out_path + 'Grandaverage_firstMorphed_SemLoc_icaclean_AbstractBase_' + conn_method + '_ThetaAlphaBetaGamma' + trange + label_name[0:-3]
	matx_stc.save(out_file)    	
	
