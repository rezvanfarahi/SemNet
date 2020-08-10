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
sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed) 
sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
###

import matplotlib.pyplot as plt
import numpy as np
import mne
from mne import io
from mne.io import Raw
import pylab as pl
from mne import (read_evokeds, equalize_channels,read_proj, read_selection)
from mne.preprocessing.ica import ICA
import scipy.io as sio
import operator

from mne.minimum_norm import (read_inverse_operator, make_inverse_operator, write_inverse_operator, apply_inverse, read_inverse_operator,apply_inverse_epochs, apply_inverse, source_induced_power,source_band_induced_power)
from mne.minimum_norm.inverse import (prepare_inverse_operator, _assemble_kernel)

from surfer import Brain
from surfer.io import read_stc
import logging
import os
import sklearn
import scipy.io
from mne import filter
from mne import find_events
from mne.epochs import combine_event_ids


from mne.layouts import read_layout
from mne.time_frequency import compute_raw_psd
from mne.connectivity import seed_target_indices, spectral_connectivity

###############################################################################
out_path = '/imaging/rf02/TypLexMEG/icaanalysis_results/stc/GrandAverage/power/' # root directory for your MEG data
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


tmin1=50
tstep1=100
vertices_to = [np.arange(10242), np.arange(10242)]
Matx=np.zeros((17,20484,5))


fname0 = data_path + ll[0] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc0 = mne.read_source_estimate(fname0)
Matx0=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx0[:,cc]=stc0.copy().crop(b1,b2).mean().data.squeeze()
matx_stc0 = mne.SourceEstimate(Matx0, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file0=data_path + ll[0] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc0.save(out_file0)
print "0 finished!"

fname1 = data_path + ll[1] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc1 = mne.read_source_estimate(fname1)
Matx1=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx1[:,cc]=stc1.copy().crop(b1,b2).mean().data.squeeze()
matx_stc1 = mne.SourceEstimate(Matx1, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file1=data_path + ll[1] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc1.save(out_file1)
print "1 finished!"


fname2 = data_path + ll[2] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc2 = mne.read_source_estimate(fname2)
Matx2=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx2[:,cc]=stc2.copy().crop(b1,b2).mean().data.squeeze()
matx_stc2 = mne.SourceEstimate(Matx2, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file2=data_path + ll[2] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc2.save(out_file2)
print "2 finished!"


fname3 = data_path + ll[3] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc3 = mne.read_source_estimate(fname3)
Matx3=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx3[:,cc]=stc3.copy().crop(b1,b2).mean().data.squeeze()
matx_stc3 = mne.SourceEstimate(Matx3, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file3=data_path + ll[3] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc3.save(out_file3)
print "3 finished!"


fname4 = data_path + ll[4] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc4 = mne.read_source_estimate(fname4)
Matx4=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx4[:,cc]=stc4.copy().crop(b1,b2).mean().data.squeeze()
matx_stc4 = mne.SourceEstimate(Matx4, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file4=data_path + ll[4] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc4.save(out_file4)
print "4 finished!"


fname5 = data_path + ll[5] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc5 = mne.read_source_estimate(fname5)
Matx5=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx5[:,cc]=stc5.copy().crop(b1,b2).mean().data.squeeze()
matx_stc5 = mne.SourceEstimate(Matx5, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file5=data_path + ll[5] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc5.save(out_file5)
print "5 finished!"


fname6 = data_path + ll[6] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc6 = mne.read_source_estimate(fname6)
Matx6=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx6[:,cc]=stc6.copy().crop(b1,b2).mean().data.squeeze()
matx_stc6 = mne.SourceEstimate(Matx6, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file6=data_path + ll[6] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc6.save(out_file6)
print "6 finished!"

fname7 = data_path + ll[7] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc7 = mne.read_source_estimate(fname7)
Matx7=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx7[:,cc]=stc7.copy().crop(b1,b2).mean().data.squeeze()
matx_stc7 = mne.SourceEstimate(Matx7, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file7=data_path + ll[7] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc7.save(out_file7)
print "7 finished!"


fname8 = data_path + ll[8] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc8 = mne.read_source_estimate(fname8)
Matx8=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx8[:,cc]=stc8.copy().crop(b1,b2).mean().data.squeeze()
matx_stc8 = mne.SourceEstimate(Matx8, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file8=data_path + ll[8] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc8.save(out_file8)
print "8 finished!"


fname9 = data_path + ll[9] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc9 = mne.read_source_estimate(fname9)
Matx9=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx9[:,cc]=stc9.copy().crop(b1,b2).mean().data.squeeze()
matx_stc9 = mne.SourceEstimate(Matx9, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file9=data_path + ll[9] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc9.save(out_file9)
print "9 finished!"


fname10 = data_path + ll[10] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc10 = mne.read_source_estimate(fname10)
Matx10=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx10[:,cc]=stc10.copy().crop(b1,b2).mean().data.squeeze()
matx_stc10 = mne.SourceEstimate(Matx10, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file10=data_path + ll[10] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc10.save(out_file10)
print "10 finished!"


fname11 = data_path + ll[11] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc11 = mne.read_source_estimate(fname11)
Matx11=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx11[:,cc]=stc11.copy().crop(b1,b2).mean().data.squeeze()
matx_stc11 = mne.SourceEstimate(Matx11, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file11=data_path + ll[11] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc11.save(out_file11)
print "11 finished!"


fname12 = data_path + ll[12] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc12 = mne.read_source_estimate(fname12)
Matx12=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx12[:,cc]=stc12.copy().crop(b1,b2).mean().data.squeeze()
matx_stc12 = mne.SourceEstimate(Matx12, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file12=data_path + ll[12] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc12.save(out_file12)
print "12 finished!"

fname13 = data_path + ll[13] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc13 = mne.read_source_estimate(fname13)
Matx13=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx13[:,cc]=stc13.copy().crop(b1,b2).mean().data.squeeze()
matx_stc13 = mne.SourceEstimate(Matx13, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file13=data_path + ll[13] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc13.save(out_file13)
print "13 finished!"

fname14 = data_path + ll[14] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc14 = mne.read_source_estimate(fname14)
Matx14=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx14[:,cc]=stc14.copy().crop(b1,b2).mean().data.squeeze()
matx_stc14 = mne.SourceEstimate(Matx14, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file14=data_path + ll[14] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc14.save(out_file14)
print "14 finished!"

fname15 = data_path + ll[15] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc15 = mne.read_source_estimate(fname15)
Matx15=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx15[:,cc]=stc15.copy().crop(b1,b2).mean().data.squeeze()
matx_stc15 = mne.SourceEstimate(Matx15, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file15=data_path + ll[15] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc15.save(out_file15)
print "15 finished!"


fname16 = data_path + ll[16] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_m500_700_gamma'  
stc16 = mne.read_source_estimate(fname16)
Matx16=np.zeros((20484,5))
b2=0.05
for cc in range(5):
		b1=b2+0.001; b2=b1+0.1-0.001
	 	Matx16[:,cc]=stc16.copy().crop(b1,b2).mean().data.squeeze()
matx_stc16 = mne.SourceEstimate(Matx16, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='avgsubj')

out_file16=data_path + ll[16] + 'Morphed_SemLoc_icaclean_Abstract_Power_ratio_gamma_meanresampled_50_550_100ms'
matx_stc16.save(out_file16)
print "16 finished!"


stc_mean=np.mean([matx_stc0,matx_stc1,matx_stc2,matx_stc3,matx_stc4,matx_stc5,matx_stc6,matx_stc7,matx_stc8,matx_stc9,matx_stc10,matx_stc11,matx_stc12,matx_stc13,matx_stc14,matx_stc15,matx_stc16])
import scipy
stc_min=np.amin([matx_stc0,matx_stc1,matx_stc2,matx_stc3,matx_stc4,matx_stc5,matx_stc6,matx_stc7,matx_stc8,matx_stc9,matx_stc10,matx_stc11,matx_stc12,matx_stc13,matx_stc14,matx_stc15,matx_stc16])
stc_max=np.amax([matx_stc0,matx_stc1,matx_stc2,matx_stc3,matx_stc4,matx_stc5,matx_stc6,matx_stc7,matx_stc8,matx_stc9,matx_stc10,matx_stc11,matx_stc12,matx_stc13,matx_stc14,matx_stc15,matx_stc16])
stc_sub=np.subtract(stc_max,stc_min)
stc_grand=np.divide(stc_mean,stc_sub)

#datag=np.log(stc_grand.data)
#datag_stc = mne.SourceEstimate(datag, vertices=stc_grand.vertno, tmin=stc_grand.tmin, tstep=stc_grand.tstep, subject='avgsubject')
#out_file=out_path + 'GrandAverage_SemLoc_icaclean_Abstract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
#stc_grand.save(out_file)

out_file=out_path + 'GrandZscore_SemLoc_icaclean_Abstract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stc_zscore.save(out_file)




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
out_file1=out_path + 'GrandAverage_SemLoc_icaclean_Abstract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs1 = mne.read_source_estimate(out_file1)
out_file2=out_path + 'GrandAverage_SemLoc_icaclean_Abstract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs2 = mne.read_source_estimate(out_file2)
stc_subtract=np.subtract(stcs1,stcs2)
out_file=out_path + 'GrandAverage_SemLoc_icaclean_Subtract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stc_subtract.save(out_file)
"""

"""
out_file1=out_path + 'GrandZscore_SemLoc_icaclean_Abstract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs1 = mne.read_source_estimate(out_file1)
out_file2=out_path + 'GrandZscore_SemLoc_icaclean_Abstract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stcs2 = mne.read_source_estimate(out_file2)
stc_subtract=np.subtract(stcs1,stcs2)
out_file=out_path + 'GrandZscore_SemLoc_icaclean_Subtract_SourcePower_ratio_gamma_meanresampled_50_550_100ms'
stc_subtract.save(out_file)
"""
