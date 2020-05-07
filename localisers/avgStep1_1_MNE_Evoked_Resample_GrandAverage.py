"""
=========================================================
GrandAvg of SemDec for different conditions 
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

import numpy as np
import mne
import os

###############################################################################
out_path = '/imaging/rf02/Semnet/stc/GrandAverage/evoked/localisers/' # root directory for your MEG data
if not os.path.exists(out_path):
    os.makedirs(out_path)
subjects_dir = '/imaging/rf02/Semnet/'    # where your MRI subdirectories are
# where event-files are
data_path = '/imaging/rf02/Semnet/'    # where event files are


# get indices for subjects to be processed from command line input
# 
print sys.argv
subject_inds = []
for ss in sys.argv[1:]:
   subject_inds.append( int( ss ) )

subject_inds=[0,1,2,3,4,5,6,7,8,9,10,11,12,13]
print "subject_inds:"
print subject_inds
print "No rejection"

list_all =  ['/meg16_0045/160303/', #8
            '/meg16_0052/160310/', #10
            '/meg16_0056/160314/',#11
            '/meg16_0069/160405/',#12 
            '/meg16_0070/160407/', #13
            '/meg16_0072/160408/', #15
            '/meg16_0073/160411/', #16
            '/meg16_0075/160411/', #17
            '/meg16_0078/160414/', #18
            '/meg16_0082/160418/', #19
            '/meg16_0086/160422/', #20
            '/meg16_0097/160512/', #21 
            '/meg16_0122/160707/', #22 
            '/meg16_0125/160712/', #24 
            ]


# subjects names used for MRI data
subjects=  [
            'MRI_meg16_0045' ,#8
            'MRI_meg16_0052' ,#10
            'MRI_meg16_0056' ,#11
    	       'MRI_meg16_0069' ,#12
 	       'MRI_meg16_0070' ,#13
            'MRI_meg16_0072' ,#15
            'MRI_meg16_0073' ,#16
            'MRI_meg16_0075' ,#17
            'MRI_meg16_0078' ,#18
            'MRI_meg16_0082' ,#19
            'MRI_meg16_0086' ,#20
            'MRI_meg16_0097' ,#21
            'MRI_meg16_0122' ,#22
            'MRI_meg16_0125' ,#24
            ]

ll = []
for ss in subject_inds:
 ll.append(list_all[ss])

print "ll:"
print ll



# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.3, 0.6


tmin1=50
tstep1=100
stc_allc=range(len(list_all))
stc_allcr=range(len(list_all))
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(len(list_all))
stc_allar=range(len(list_all))
vertices_to = [np.arange(10242), np.arange(10242)]
event_names = ['button','audio',  'colour', 'grey', 'shapes', 'scrambled']
method='MNE'
stc_grandc=range(len(event_names))
stc_grandcr=range(len(event_names))
for event_no in range(len(event_names)):
    print event_no
    for ii, meg in enumerate(ll):
        print ii
        fname = data_path + meg +method+ '_firstMorphed_ico_newreg_localisers_ica_'+event_names[event_no]+'_Source_Evoked_m100_500'  
        stc_allc[ii] = mne.read_source_estimate(fname)
#        Matx0=np.zeros((20484,5))
#        b2=0.05
#        for cc in range(5):
#            b1=b2+0.0001; b2=b1+0.1-0.0001
#            Matx0[:,cc]=stc_allc[ii].copy().crop(b1,b2).mean().data.squeeze()
#        matx_stcc = mne.SourceEstimate(Matx0, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#        stc_allcr[ii]=matx_stcc
#        out_file0=data_path + meg + 'firstMorphed_ico_localisers_ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
#        matx_stcc.save(out_file0)
        
    stc_grandc[event_no]=np.mean(stc_allc)
#    stc_grandcr[event_no]=np.mean(stc_allcr)
    
    #datag=np.log(stc_grand.data)
    #datag_stc = mne.SourceEstimate(datag, vertices=stc_grand.vertno, tmin=stc_grand.tmin, tstep=stc_grand.tstep, subject='fsaverageect')
    out_file=out_path + method+'_GrandAverage_firstMorphed_ico_newreg_localisers_ica_'+event_names[event_no]+'_Source_Evoked_m100_500'
    out_filer=out_path + method+'_GrandAverage_firstMorphed_ico_newreg_localisers_ica_'+event_names[event_no]+'_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
    stc_grandc[event_no].save(out_file)
#    stc_grandcr[event_no].save(out_filer)

out_file=out_path + method+'_GrandAverage_firstMorphed_ico_newreg_localisers_ica_'+event_names[2]+'_'+event_names[3]+'_Source_Evoked_m100_500'
stc_col=np.subtract(stc_grandc[2],stc_grandc[3])
stc_col.save(out_file)
out_file=out_path + method+'_GrandAverage_firstMorphed_ico_newreg_localisers_ica_'+event_names[4]+'_'+event_names[5]+'_Source_Evoked_m100_500'
stc_shp=np.subtract(stc_grandc[4],stc_grandc[5])
stc_shp.save(out_file)

#out_filef=out_path + 'GrandAverage_firstMorphed_ico_localisers_ica_all_Source_Evoked_m300_600'
#stc_grand=np.mean(stc_grandc)
#stc_grand.save(out_filef)
#out_filefr=out_path + 'GrandAverage_firstMorphed_ico_localisers_ica_all_Source_Evoked_mnrsmpl_50_550_100ms_0ov'
#stc_grandr=np.mean(stc_grandcr)
#stc_grandr.save(out_filefr)
"""
import matplotlib.pyplot as plt
in_file1=out_path + 'GrandAverage_localisers_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_theta'  
stc1=mne.read_source_estimate(in_file1)
in_file2=out_path + 'GrandAverage_SemDec_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_alpha' 
stc2=mne.read_source_estimate(in_file2) 
in_file3=out_path + 'GrandAverage_SemDec_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_beta'  
stc3=mne.read_source_estimate(in_file3)
in_file4=out_path + 'GrandAverage_SemDec_icaclean_oldreg_Subtract_SourcePower_ratio_m500_700_gamma'  
stc4=mne.read_source_estimate(in_file4)
plt.plot(stc1.times[100:1100], stc1.data[:,100:1100].mean(axis=0), label='theta')
plt.plot(stc2.times[100:1100], stc2.data[:,100:1100].mean(axis=0), label='alpha')
plt.plot(stc3.times[100:1100], stc3.data[:,100:1100].mean(axis=0), label='beta')
plt.plot(stc4.times[100:1100], stc4.data[:,100:1100].mean(axis=0), label='gamma')
plt.xlabel('Time (ms)')
plt.ylabel('Power')
plt.legend()
plt.title('Mean source induced power Subtract')
plt.show()
"""


