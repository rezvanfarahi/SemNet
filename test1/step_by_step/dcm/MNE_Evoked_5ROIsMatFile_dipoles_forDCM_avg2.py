# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 06:26:19 2017

@author: rf02
"""

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
#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.3.1')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
sys.path.insert(1,'/imaging/rf02/TypLexMEG/mne_python0.9')

# for qsub
# (add ! here if needed) 
#sys.path.append('/imaging/local/software/anaconda/latest/x86_64/bin/python')
#sys.path.append('/imaging/local/software/mne_python/git-master-0.8/mne-python/')
#sys.path.insert(1,'/imaging/local/software/anaconda/2.4.1/3/lib/python3.5/site-packages/scipy')
###

#sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
#sys.path.insert(1,'/imaging/rf02/scikit-learn-0.15.0')
sys.path.insert(1,'/imaging/rf02/scikit-learn-0.16.1')
sys.path.insert(1,'/home/rf02/.local/lib/python2.7/site-packages')

#import joblib
###
import sklearn

import os
import numpy as np
import mne
import scipy
from scipy.stats.mstats import zscore
import scipy.io as scio
#import sklearn
from mne.io import Raw
from mne.epochs import equalize_epoch_counts
from mne.minimum_norm import (read_inverse_operator, psf_ctf_28Oct15, apply_inverse)
from mne.stats import f_threshold_mway_rm, f_mway_rm
from mne import (io, spatial_tris_connectivity, grade_to_tris)
from mne.epochs import equalize_epoch_counts
from mne.stats import (permutation_cluster_test,summarize_clusters_stc)
from mne.minimum_norm import apply_inverse, read_inverse_operator
from mne.stats import f_threshold_mway_rm, f_mway_rm, fdr_correction, spatio_temporal_cluster_test
from mne.datasets import sample
from mne.viz import mne_analyze_colormap
import copy
from copy import deepcopy
from math import sqrt
import numpy as np
from scipy import linalg

from mne.io.constants import FIFF
from mne.io.open import fiff_open
from mne.io.tag import find_tag
from mne.io.matrix import (_read_named_matrix, _transpose_named_matrix,
                         write_named_matrix)
from mne.io.proj import _read_proj, make_projector, _write_proj
from mne.io.proj import _needs_eeg_average_ref_proj
from mne.io.tree import dir_tree_find
from mne.io.write import (write_int, write_float_matrix, start_file,
                        start_block, end_block, end_file, write_float,
                        write_coord_trans, write_string)

from mne.io.pick import channel_type, pick_info, pick_types
from mne.cov import prepare_noise_cov, _read_cov, _write_cov, Covariance
from mne.forward import (compute_depth_prior, _read_forward_meas_info,
                       write_forward_meas_info, is_fixed_orient,
                       compute_orient_prior, convert_forward_solution)
from mne.source_space import (_read_source_spaces_from_tree,
                            find_source_space_hemi, _get_vertno,
                            _write_source_spaces_to_fid, label_src_vertno_sel)
#from mne.transforms import _ensure_trans, transform_surface_to
#from mne.source_estimate import _make_stc
from mne.utils import check_fname, logger, verbose
#from functools import reduce

class InverseOperatorMine(dict):
    """InverseOperator class to represent info from inverse operator."""

    def copy(self):
        """Return a copy of the InverseOperator."""
        return InverseOperator(deepcopy(self))

    def __repr__(self):  # noqa: D105
        """Summarize inverse info instead of printing all."""
        entr = '<InverseOperator'

        nchan = len(pick_types(self['info'], meg=True, eeg=False))
        entr += ' | ' + 'MEG channels: %d' % nchan
        nchan = len(pick_types(self['info'], meg=False, eeg=True))
        entr += ' | ' + 'EEG channels: %d' % nchan

        # XXX TODO: This and the __repr__ in SourceSpaces should call a
        # function _get_name_str() in source_space.py
        if self['src'][0]['type'] == 'surf':
            entr += (' | Source space: Surface with %d vertices'
                     % self['nsource'])
        elif self['src'][0]['type'] == 'vol':
            entr += (' | Source space: Volume with %d grid points'
                     % self['nsource'])
        elif self['src'][0]['type'] == 'discrete':
            entr += (' | Source space: Discrete with %d dipoles'
                     % self['nsource'])

        source_ori = {FIFF.FIFFV_MNE_UNKNOWN_ORI: 'Unknown',
                      FIFF.FIFFV_MNE_FIXED_ORI: 'Fixed',
                      FIFF.FIFFV_MNE_FREE_ORI: 'Free'}
        entr += ' | Source orientation: %s' % source_ori[self['source_ori']]
        entr += '>'

        return entr
def _prepare_forward_mine(forward, info, vertsub_vect,noise_cov, pca=False, rank=None,
                     verbose=None):
    """Prepare forward solution for inverse solvers."""
    # fwd['sol']['row_names'] may be different order from fwd['info']['chs']
    fwd_sol_ch_names = forward['sol']['row_names']
    ch_names = [c['ch_name'] for c in info['chs']
                if ((c['ch_name'] not in info['bads'] and
                     c['ch_name'] not in noise_cov['bads']) and
                    (c['ch_name'] in fwd_sol_ch_names and
                     c['ch_name'] in noise_cov.ch_names))]

    if not len(info['bads']) == len(noise_cov['bads']) or \
            not all(b in noise_cov['bads'] for b in info['bads']):
        logger.info('info["bads"] and noise_cov["bads"] do not match, '
                    'excluding bad channels from both')

    n_chan = len(ch_names)
    logger.info("Computing inverse operator with %d channels." % n_chan)

    #
    #   Handle noise cov
    #
    noise_cov = prepare_noise_cov(noise_cov, info, ch_names, rank)

    #   Omit the zeroes due to projection
    eig = noise_cov['eig']
    nzero = (eig > 0)
    n_nzero = sum(nzero)

    if pca:
        #   Rows of eigvec are the eigenvectors
        whitener = noise_cov['eigvec'][nzero] / np.sqrt(eig[nzero])[:, None]
        logger.info('Reducing data rank to %d' % n_nzero)
    else:
        whitener = np.zeros((n_chan, n_chan), dtype=np.float)
        whitener[nzero, nzero] = 1.0 / np.sqrt(eig[nzero])
        #   Rows of eigvec are the eigenvectors
        whitener = np.dot(whitener, noise_cov['eigvec'])

    gain = forward['sol']['data']

    # This actually reorders the gain matrix to conform to the info ch order
    fwd_idx = [fwd_sol_ch_names.index(name) for name in ch_names]
    gain = gain[fwd_idx]
    # Any function calling this helper will be using the returned fwd_info
    # dict, so fwd['sol']['row_names'] becomes obsolete and is NOT re-ordered

    info_idx = [info['ch_names'].index(name) for name in ch_names]
    fwd_info = pick_info(info, info_idx)

    logger.info('Total rank is %d' % n_nzero)
    #Rez pick sub leadfield
#    gainpick=np.hstack((vertsub_vect,vertsub_vect+int(gain.shape[1]/3),vertsub_vect+2*int(gain.shape[1]/3)))
    gainpick=np.sort(np.hstack((3*vertsub_vect,3*vertsub_vect+1,3*vertsub_vect+2)))

    gain=gain[:,gainpick].copy()
    return fwd_info, gain, noise_cov, whitener, n_nzero


@verbose

def make_inverse_operator_mine(info, forward, vertsub_vect, noise_cov, loose=0.2, depth=0.8,
                          fixed=False, limit_depth_chs=True, rank=None,
                          verbose=None):
    """Assemble inverse operator.

    Parameters
    ----------
    info : dict
        The measurement info to specify the channels to include.
        Bad channels in info['bads'] are not used.
    forward : dict
        Forward operator.
    noise_cov : instance of Covariance
        The noise covariance matrix.
    loose : None | float in [0, 1]
        Value that weights the source variances of the dipole components
        defining the tangent space of the cortical surfaces. Requires surface-
        based, free orientation forward solutions.
    depth : None | float in [0, 1]
        Depth weighting coefficients. If None, no depth weighting is performed.
    fixed : bool
        Use fixed source orientations normal to the cortical mantle. If True,
        the loose parameter is ignored.
    limit_depth_chs : bool
        If True, use only grad channels in depth weighting (equivalent to MNE
        C code). If grad chanels aren't present, only mag channels will be
        used (if no mag, then eeg). If False, use all channels.
    rank : None | int | dict
        Specified rank of the noise covariance matrix. If None, the rank is
        detected automatically. If int, the rank is specified for the MEG
        channels. A dictionary with entries 'eeg' and/or 'meg' can be used
        to specify the rank for each modality.
    verbose : bool, str, int, or None
        If not None, override default verbose level (see :func:`mne.verbose`).

    Returns
    -------
    inv : instance of InverseOperator
        Inverse operator.

    Notes
    -----
    For different sets of options (**loose**, **depth**, **fixed**) to work,
    the forward operator must have been loaded using a certain configuration
    (i.e., with **force_fixed** and **surf_ori** set appropriately). For
    example, given the desired inverse type (with representative choices
    of **loose** = 0.2 and **depth** = 0.8 shown in the table in various
    places, as these are the defaults for those parameters):

        +---------------------+-----------+-----------+-----------+-----------------+--------------+
        | Inverse desired                             | Forward parameters allowed                 |
        +=====================+===========+===========+===========+=================+==============+
        |                     | **loose** | **depth** | **fixed** | **force_fixed** | **surf_ori** |
        +---------------------+-----------+-----------+-----------+-----------------+--------------+
        | | Loose constraint, | 0.2       | 0.8       | False     | False           | True         |
        | | Depth weighted    |           |           |           |                 |              |
        +---------------------+-----------+-----------+-----------+-----------------+--------------+
        | | Loose constraint  | 0.2       | None      | False     | False           | True         |
        +---------------------+-----------+-----------+-----------+-----------------+--------------+
        | | Free orientation, | None      | 0.8       | False     | False           | True         |
        | | Depth weighted    |           |           |           |                 |              |
        +---------------------+-----------+-----------+-----------+-----------------+--------------+
        | | Free orientation  | None      | None      | False     | False           | True | False |
        +---------------------+-----------+-----------+-----------+-----------------+--------------+
        | | Fixed constraint, | None      | 0.8       | True      | False           | True         |
        | | Depth weighted    |           |           |           |                 |              |
        +---------------------+-----------+-----------+-----------+-----------------+--------------+
        | | Fixed constraint  | None      | None      | True      | True            | True         |
        +---------------------+-----------+-----------+-----------+-----------------+--------------+

    Also note that, if the source space (as stored in the forward solution)
    has patch statistics computed, these are used to improve the depth
    weighting. Thus slightly different results are to be expected with
    and without this information.
    """  # noqa: E501
    is_fixed_ori = is_fixed_orient(forward)

#    if fixed and loose is not None:
#        #        mne.utils.warn('When invoking make_inverse_operator with fixed=True, the loose '
#        #             'parameter is ignored.')
#        loose = None

    if is_fixed_ori and not fixed:
        raise ValueError('Forward operator has fixed orientation and can only '
                         'be used to make a fixed-orientation inverse '
                         'operator.')
    if fixed:
        if depth is not None:
            if is_fixed_ori or not forward['surf_ori']:
                raise ValueError('For a fixed orientation inverse solution '
                                 'with depth weighting, the forward solution '
                                 'must be free-orientation and in surface '
                                 'orientation')
        elif forward['surf_ori'] is False:
            raise ValueError('For a fixed orientation inverse solution '
                             'without depth weighting, the forward solution '
                             'must be in surface orientation')

    # depth=None can use fixed fwd, depth=0<x<1 must use free ori
    if depth is not None:
        if not (0 < depth <= 1):
            raise ValueError('depth should be a scalar between 0 and 1')
        if is_fixed_ori:
            raise ValueError('You need a free-orientation, surface-oriented '
                             'forward solution to do depth weighting even '
                             'when calculating a fixed-orientation inverse.')
        if not forward['surf_ori']:
            forward = convert_forward_solution(forward, surf_ori=True)
        assert forward['surf_ori']
    if loose is not None:
        if not (0 <= loose <= 1):
            raise ValueError('loose value should be smaller than 1 and bigger '
                             'than 0, or None for not loose orientations.')
        if loose < 1 and not forward['surf_ori']:
            raise ValueError('Forward operator is not oriented in surface '
                             'coordinates. A loose inverse operator requires '
                             'a surface-based, free orientation forward '
                             'operator.')

    #
    # 1. Read the bad channels
    # 2. Read the necessary data from the forward solution matrix file
    # 3. Load the projection data
    # 4. Load the sensor noise covariance matrix and attach it to the forward
    #
    n_src=len(forward['src'][0]['vertno'])+len(forward['src'][1]['vertno'])
    gainpick=np.sort(np.hstack((3*vertsub_vect,3*vertsub_vect+1,3*vertsub_vect+2)))

    gain_info, gain, noise_cov, whitener, n_nzero = _prepare_forward_mine(forward, info, vertsub_vect,noise_cov, rank=rank)
#    forward['info']._check_consistency()

    #
    # 5. Compose the depth-weighting matrix
    #

    if depth is not None:
        patch_areas = forward.get('patch_areas', None)
        depth_prior = compute_depth_prior(gain, gain_info, is_fixed_ori,
                                          exp=depth, patch_areas=patch_areas,
                                          limit_depth_chs=limit_depth_chs)
    else:
        depth_prior = np.ones(gain.shape[1], dtype=gain.dtype)

    # Deal with fixed orientation forward / inverse
    if fixed:
        if depth is not None:
            # Convert the depth prior into a fixed-orientation one
            logger.info('    Picked elements from a free-orientation '
                        'depth-weighting prior into the fixed-orientation one')
        if not is_fixed_ori:
            # Convert to the fixed orientation forward solution now
            depth_prior = depth_prior[2::3]
            forward = convert_forward_solution(
                forward, surf_ori=forward['surf_ori'], force_fixed=True)
            is_fixed_ori = is_fixed_orient(forward)
            gain_info, gain, noise_cov, whitener, n_nzero = \
                _prepare_forward_mine(forward, info, vertsub_vect,noise_cov, verbose=False)

    logger.info("Computing inverse operator with %d channels."
                % len(gain_info['ch_names']))

    #
    # 6. Compose the source covariance matrix
    #

    logger.info('Creating the source covariance matrix')
    source_cov = depth_prior.copy()
    depth_prior = dict(data=depth_prior, kind=FIFF.FIFFV_MNE_DEPTH_PRIOR_COV,
                       bads=[], diag=True, names=[], eig=None,
                       eigvec=None, dim=depth_prior.size, nfree=1,
                       projs=[])

    # apply loose orientations
    if not is_fixed_ori:
        orient_prior = compute_orient_prior(forward, loose=loose)
        orient_prior_all = compute_orient_prior(forward, loose=loose)
#        gainpick=np.hstack((vertsub_vect,vertsub_vect+int(gain.shape[1]/3),vertsub_vect+2*int(gain.shape[1]/3)))
        gainpick=np.sort(np.hstack((3*vertsub_vect,3*vertsub_vect+1,3*vertsub_vect+2)))

        orient_prior=orient_prior[gainpick].copy()
        source_cov *= orient_prior
#        iop=np.zeros((n_src*3,))
#        iop[gainpick]=orient_prior.copy()
        orient_prior=orient_prior_all.copy()
        orient_prior = dict(data=orient_prior,
                            kind=FIFF.FIFFV_MNE_ORIENT_PRIOR_COV,
                            bads=[], diag=True, names=[], eig=None,
                            eigvec=None, dim=orient_prior.size, nfree=1,
                            projs=[])
    else:
        orient_prior = None

    # 7. Apply fMRI weighting (not done)

    #
    # 8. Apply the linear projection to the forward solution
    # 9. Apply whitening to the forward computation matrix
    #
    
    logger.info('Whitening the forward solution.')
    gain = np.dot(whitener, gain)

    # 10. Exclude the source space points within the labels (not done)

    #
    # 11. Do appropriate source weighting to the forward computation matrix
    #

    # Adjusting Source Covariance matrix to make trace of G*R*G' equal
    # to number of sensors.
    logger.info('Adjusting source covariance matrix.')
    source_std = np.sqrt(source_cov)
    gain *= source_std[None, :]
    trace_GRGT = linalg.norm(gain, ord='fro') ** 2
    scaling_source_cov = n_nzero / trace_GRGT
    source_cov *= scaling_source_cov
    gain *= sqrt(scaling_source_cov)
    depth_prior = np.ones(gain.shape[1], dtype=gain.dtype)
    scv=np.zeros((n_src*3,))
    scv[gainpick]=source_cov.copy()
    source_cov=scv.copy()
    source_cov = dict(data=source_cov, dim=source_cov.size,
                      kind=FIFF.FIFFV_MNE_SOURCE_COV, diag=True,
                      names=[], projs=[], eig=None, eigvec=None,
                      nfree=1, bads=[])

    # now np.trace(np.dot(gain, gain.T)) == n_nzero
    # logger.info(np.trace(np.dot(gain, gain.T)), n_nzero)

    #
    # 12. Decompose the combined matrix
    #

    logger.info('Computing SVD of whitened and weighted lead field '
                'matrix.')
    eigen_fields, sing, eigen_leads_pre = linalg.svd(gain, full_matrices=False)#el: nchanxnsource
    eigen_leads=np.zeros((gain.shape[0],n_src*3))
    eigen_leads[:,gainpick]=eigen_leads_pre.copy()
    logger.info('    largest singular value = %g' % np.max(sing))
    logger.info('    scaling factor to adjust the trace = %g' % trace_GRGT)

    eigen_fields = dict(data=eigen_fields.T, col_names=gain_info['ch_names'],
                        row_names=[], nrow=eigen_fields.shape[1],
                        ncol=eigen_fields.shape[0])#size nchanxnchan
    eigen_leads = dict(data=eigen_leads.T, nrow=eigen_leads.shape[1],
                       ncol=eigen_leads.shape[0], row_names=[],
                       col_names=[])#size nsourcexnchan
    nave = 1.0

    # Handle methods
    has_meg = False
    has_eeg = False
    ch_idx = [k for k, c in enumerate(info['chs'])
              if c['ch_name'] in gain_info['ch_names']]
    for idx in ch_idx:
        ch_type = channel_type(info, idx)
        if ch_type == 'eeg':
            has_eeg = True
        if (ch_type == 'mag') or (ch_type == 'grad'):
            has_meg = True
    if has_eeg and has_meg:
        methods = FIFF.FIFFV_MNE_MEG_EEG
    elif has_meg:
        methods = FIFF.FIFFV_MNE_MEG
    else:
        methods = FIFF.FIFFV_MNE_EEG

    # We set this for consistency with mne C code written inverses
    if depth is None:
        depth_prior = None
    inv_op = dict(eigen_fields=eigen_fields, eigen_leads=eigen_leads,
                  sing=sing, nave=nave, depth_prior=depth_prior,
                  source_cov=source_cov, noise_cov=noise_cov,
                  orient_prior=orient_prior, projs=deepcopy(info['projs']),
                  eigen_leads_weighted=False, source_ori=forward['source_ori'],
                  mri_head_t=deepcopy(forward['mri_head_t']),
                  methods=methods, nsource=forward['nsource'],
                  coord_frame=forward['coord_frame'],
                  source_nn=forward['source_nn'].copy(),
                  src=deepcopy(forward['src']), fmri_prior=None)
    inv_info = deepcopy(forward['info'])
    inv_info['bads'] = [bad for bad in info['bads']
                        if bad in inv_info['ch_names']]
#    inv_info._check_consistency()
    inv_op['units'] = 'Am'
    inv_op['info'] = inv_info

    return inv_op, InverseOperatorMine(inv_op)

###############################################################################
#out_path = '/home/rf02/rezvan/test1/step_by_step/dcm/' # root directory for your MEG data
out_path = '/imaging/rf02/TypLexMEG/dcm/dipole/'
orig_path='/group/erp/data/olaf.hauk/MEG/TypLexM/data/MF22_new/'

subjects_dir = '/imaging/rf02/TypLexMEG/'    # where your MRI subdirectories are
# where event-files aremean
data_path = '/imaging/rf02/TypLexMEG/'    # where event files are

label_path = '/home/rf02/rezvan/TypLexMEG/createdlabels_SL/'

inv_path = '/imaging/rf02/TypLexMEG/'

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
labelss = mne.read_labels_from_annot(subject='fsaverage', parc='rez_aparc_25oct16', subjects_dir=data_path)
#labelss.remove(labelss[-3])
#labelss.remove(labelss[-2])

label_names=[label.name for label in labelss]
label_index=np.array([label_names.index('ATL_latmed-lh'),label_names.index('ATL_latmed-rh'),label_names.index('posterior_inferiorfrontal-lh'),label_names.index('posterior_inferiorfrontal-rh'),label_names.index('middle_temporal-lh'),label_names.index('middle_temporal-rh'),label_names.index('AG_new-lh'),label_names.index('AG_new-rh'),label_names.index('vWFA2-lh'),label_names.index('vWFA2-rh'),label_names.index('lateraloccipital-lh'),label_names.index('lateraloccipital-rh')])#label_names.index('ATL_latmed-lh'),label_names.index('supramarginal-lh'),
lfname=data_path+'fsaverage/label/AG_manual-lh.label'
labelaglh=mne.read_label(lfname,subject='fsaverage')
lfname=data_path+'fsaverage/label/AG_manual-rh.label'
labelagrh=mne.read_label(lfname,subject='fsaverage')
#lfname=data_path+'fsaverage/label/ATL_manual-lh.label'
#labelatllh=mne.read_label(lfname,subject='fsaverage')
#lfname=data_path+'fsaverage/label/ATL_manual-rh.label'
#labelatlrh=mne.read_label(lfname,subject='fsaverage')

lfname=data_path+'fsaverage/label/lh.ATL_DCM-lh.label'
labelatllh=mne.read_label(lfname,subject='fsaverage')
lfname=data_path+'fsaverage/label/lh.parsorbitalis.label'
labelifglh=mne.read_label(lfname,subject='fsaverage')
lfname=data_path+'fsaverage/label/rh.parsorbitalis.label'
labelifgrh=mne.read_label(lfname,subject='fsaverage')
#lfname=data_path+'fsaverage/label/AG_hsmetafinal2-lh.label'
#labelag=mne.read_label(lfname,subject='fsaverage')
#lfname=data_path+'fsaverage/label/lh.ATL_DCM-lh.label'
#labelatl=mne.read_label(lfname,subject='fsaverage')
#lfname=data_path+'fsaverage/label/lh.parsorbitalis.label'
#labelifg_por=mne.read_label(lfname,subject='fsaverage')
#lfname=data_path+'fsaverage/label/lh.parstriangularis.label'
#labelifg_ptr=mne.read_label(lfname,subject='fsaverage')
#labelifg=labelifg_por+labelifg_ptr
labellist=[labelss[ii] for ii in label_index]
labellist[0]=labelatllh.copy()
#labellist[1]=labelatlrh.copy()
labellist[2]=labelifglh.copy()
labellist[3]=labelifgrh.copy()
labellist[6]=labelaglh.copy()
labellist[7]=labelagrh.copy()
#labellist=[labelatllh,labelatlrh,labelss[label_index[2]],labelss[label_index[3]],labelss[label_index[2]],labelss[label_index[3]],labelag,labelss[label_index[4]]]
for lbl in labellist:
    lbl.values.fill(1.0)


#labellist = ['lh.ATL_latmed.label']#['atlleft-lh.label']#['lh.ATL_latmed.label']
#labellist2 = ['AG_new-lh.label','vWFA2-lh.label']#'AGSMG_meta2-lh.label'['lh.ATL_DCM-lh.label','lh.WFA_DCM-lh.label']
#lfname=data_path+'fsaverage/label/'+labellist2[0]#label_path+labellist[0]#data_path+'fsaverage/label/'+labellist[0]
#label0=mne.read_label(lfname,subject='fsaverage')
#labellist.append(label0)
#lfname=data_path+'fsaverage/label/'+labellist2[1]#label_path+labellist[0]#data_path+'fsaverage/label/'+labellist[0]
#label1=mne.read_label(lfname,subject='fsaverage')
#labellist.append(label1)
#print labellist
#labellist = ['lh.WFA_DCM-lh.label','lh.ATL_DCM-lh.label']

# Labels/ROIs
#labellist = ['atlleft-lh']#, 'atlleftt-rh'
tmin, tmax = -0.5, 0.7
# frequency bands with variable number of cycles for wavelets
stim_delay = 0.034 # delay in s

tmin1=-500
tstep1=700
stc_allc=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
stc_alla=range(17)
vertices_to = [np.arange(10242), np.arange(10242)]
n_ROIs=len(labellist)
n_ROIs_final=5

conabs_mat=np.zeros((n_ROIs_final,1201,2))
sfrq=1.
n_subjects=17
n_times=int(1201/sfrq)
#n_levels_facts=4
#X=np.zeros((n_subjects,5,n_times,2))
sonset=int(abs(tmin*1000)/sfrq)
Matx_all=np.zeros((20484,n_ROIs,n_subjects,2))
stc_con=list()
stc_abs=list()
roivert=-50
roivert_upperbad=-50
for ii, meg in enumerate(ll):
    conabs_mat=np.zeros((n_ROIs_final,1201,2))
    subject_no=subject_inds[ii]
#    conind=0
    for event_no,event_name in enumerate(['Concrete','Abstract']):
#        conind=conind+1
#        print conind
        snr = 3.0
        print snr
        lambda2 = 1.0 / snr ** 2
        fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
        forward = mne.read_forward_solution(fname_fwd,surf_ori=True)
        print "load source estimates"
        method = "MNE"  # use dSPM method (could also be MNE or sLORETA)
        fname_inv = data_path + meg + 'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'#'InvOp_ico5diagreg_fft_1_48_clean_ica_EMEG-inv.fif'#
        # Load data
        inverse_operator = read_inverse_operator(fname_inv)
#        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
#        src = mne.read_source_spaces(srcin)
#        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_'+event_name+'_Source_Evoked_m500_700' 
#        this_stc = mne.read_source_estimate(fname)  
#        tmin1=-500
#        tstep1=1
#        thisstc_data=np.ones((this_stc.data.shape[0],1))
#        vertices_to = [np.arange(10242), np.arange(10242)]
#        thisstc_tolabel = mne.SourceEstimate(thisstc_data, vertices=vertices_to,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#        thisstc_label=mne.stc_to_label(thisstc_tolabel, src)

#        ntimes=this_stc.data.shape[1]
        print "select 100 most sensitive labels"
        
        method = 'MNE'  # can be 'MNE' or 'sLORETA'
        mode = 'svd'
        n_svd_comp = 1
        labellist_temp=copy.deepcopy(labellist)
        labels = [labellist_temp[jj].morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path) for jj in range(len(labellist))]
        stc_ctf_mne=psf_ctf_28Oct15.cross_talk_function(inverse_operator, forward, labels,method=method, lambda2=lambda2, signed=False, mode=mode, n_svd_comp=1, verbose=None)
        
        subject_from = subjects[subject_no]
        subject_to = 'fsaverage'
        vertices_to = [np.arange(10242), np.arange(10242)]
        morphmat1=mne.compute_morph_matrix(subject_from, subject_to, stc_ctf_mne.vertices, vertices_to, subjects_dir=data_path)
        stc_to=stc_ctf_mne.morph_precomputed(subject_to,vertices_to,morphmat1, subject_from)
        if event_no==0:
            stc_con.append(stc_to)
        else:
            stc_abs.append(stc_to)
        Matx_all[:,:,ii,event_no]=stc_to.data[:,:-1]# Nv x Nlabels+1
        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
        src_avg = mne.read_source_spaces(srcin)
stc_conm=np.mean(stc_con)
stc_absm=np.mean(stc_abs)
stc_bothm=np.mean([stc_conm,stc_absm])
#ctfstc_path=data_path+'dcm/stc_ctf_bothm'
#stc_bothm.save(ctfstc_path)
#stc_bothm = mne.read_source_estimate(ctfstc_path) 
#Matx=np.mean(Matx_all,3)
#Matx=np.mean(Matx,2)  
allwin_verts=np.zeros((n_subjects,n_ROIs,np.abs(roivert)))
for ii, meg in enumerate(ll):
    conabs_mat=np.zeros((n_ROIs_final,1201,2))
    subject_no=subject_inds[ii]
    print "reading raw file"
    raw_fname = data_path + meg + 'semloc_ssstf_fft_1_48_clean_ica_raw.fif'
    raw = mne.io.Raw(raw_fname, preload=True)

    print "Reading events"   
    events_fname = orig_path + meg + 'semloc_raw_ssstnew.txt'
    events = mne.read_events(events_fname)

    print "Reading evevokeds"   
    #fname_evoked = directory + 'semloc_ssstf_fft_1_48_clean_ica_raw-ave.fif'
    #evoked = mne.read_evokeds(fname_evoked, condition=0, baseline=(None, 0))
    stim_delay = 0.034 # delay in s
    events[:,0] = events[:,0]+np.round( raw.info['sfreq']*stim_delay )

##################################################################################
    #extracting epochs
    tmin, tmax = -0.5, 0.7
    event_id={'cncrt_wrd': 1, 'abs_wrd': 2}
    exclude = raw.info['bads']
    include = []
    picks = mne.pick_types(raw.info, meg=True, eeg=True, eog=True, stim=False, include=include, exclude=exclude)
    print "Epoching"
    if ii==3:
        reject = dict(eeg=150e-6, grad=200e-12, mag=4e-12)# eog=150e-6,    
    else:
        reject = dict(eeg=120e-6, grad=200e-12, mag=4e-12)# eog=150e-6,
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, proj=True, baseline=(None, 0), reject=reject)
    evoked=epochs.average()

   ##################################################################################


    print "making noise covariance matrix"   

   #Compute the noise covariance on baseline
    
    method_params=dict(diagonal_fixed=dict(grad=0.1, mag=0.1, eeg=0.1))
#   noise_cov = mne.compute_covariance(epochs, method=['shrunk','empirical', 'diagonal_fixed'],  tmin=None, tmax=0., return_estimators=True, method_params=method_params)
    noise_cov = mne.compute_covariance(epochs, method=['diagonal_fixed'],  tmin=None, tmax=0., method_params=method_params, projs=raw.info['projs'])
    fname_fwd = data_path + meg + 'ico5_forward_5-3L-EMEG-fwd.fif'
    forward = mne.read_forward_solution(fname_fwd,surf_ori=True)
    vertices_avg = [np.arange(10242), np.arange(10242)]
#    vertices_sub=[np.arange(forward['src'][0]['nuse']),np.arange(forward['src'][1]['nuse'])]
    vertices_sub=[forward['src'][0]['vertno'],forward['src'][1]['vertno']]

#    morphmat2=mne.compute_morph_matrix(subject_from='fsaverage', subject_to=subjects[subject_no], vertices_from=vertices_avg, vertices_to=vertices_sub, subjects_dir=data_path)
#    stc_bothm_sub=stc_bothm.copy().morph_precomputed(subject_to=subjects[subject_no],vertices_to=vertices_sub,morph_mat=morphmat2, subject_from='fsaverage')
    Matx=stc_bothm.data[:,:-1].copy()
#    Matxs=np.sum(Matx,1)
#    Matx=Matx/np.vstack((Matxs,Matxs,Matxs,Matxs,Matxs,Matxs,Matxs,Matxs,Matxs,Matxs,Matxs,Matxs)).T
    
    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage_dist-ico-5-src.fif'
    src_avg = mne.read_source_spaces(srcin)
    src_sub=forward['src']
    print "select 100 most sensitive labels"
    sub_labellist=list()#range(len(labellist))  
    winvert_sub=list()
    for this_labelc, this_label in enumerate(labellist):
        this_label_temp=copy.deepcopy(this_label)
        if this_label.hemi=='lh':
            label_verts=stc_bothm.in_label(this_label_temp).lh_vertno
        else:
            label_verts=stc_bothm.in_label(this_label_temp).rh_vertno+len(vertices_avg[0])
#        label_verts=np.asarray([np.where(forward['src'][0]['vertno']==lv0)[0][0] for lv0 in label_verts_pre])
        #            np.intersect1d(label_verts,vertices_to[0])
        avgvert_select=-np.min([len(label_verts),np.abs(roivert_upperbad)])
        vert50=np.sort(Matx[label_verts,this_labelc])[avgvert_select]
        winner_verts=label_verts[np.where(Matx[label_verts,this_labelc]>=vert50)[0]][:np.abs(avgvert_select)]#label_verts[np.where(Matx[label_verts,this_labelc]>=np.max(Matx[label_verts,this_labelc])/2.)[0]]#np.where(Matx[:,this_labelc]>=np.max(Matx[:,this_labelc])/2.)[0]#
        
        thisstc_testdata=np.zeros((len(vertices_avg[0])+len(vertices_avg[1]),1))
        thisstc_testdata[winner_verts,0]=1
        thisstc_tolabel_test = mne.SourceEstimate(thisstc_testdata, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
        morphmat2=mne.compute_morph_matrix(subject_from='fsaverage', subject_to=subjects[subject_no], vertices_from=vertices_avg, vertices_to=vertices_sub, subjects_dir=data_path)
        thisstc_tolabel_test_sub=thisstc_tolabel_test.copy().morph_precomputed(subject_to=subjects[subject_no],vertices_to=vertices_sub,morph_mat=morphmat2, subject_from='fsaverage')
        label_verts_subm=np.where(thisstc_tolabel_test_sub.data[:,0]>0)[0]
#        subvert_select=np.random.permutation(len(label_verts_subm))[:np.abs(roivert)]
#        label_verts_subm=np.sort(label_verts_subm[subvert_select])
        thisstc_testdata2=np.zeros((len(vertices_sub[0])+len(vertices_sub[1]),1))
        thisstc_testdata2[label_verts_subm,0]=1
        thisstc_tolabel_test2 = mne.SourceEstimate(thisstc_testdata2, vertices=vertices_sub,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject=subjects[subject_no])
        if this_label.hemi=='lh':
            label_temp_sub=mne.stc_to_label(thisstc_tolabel_test2, forward['src'],smooth=False)[0]
        else:
            label_temp_sub=mne.stc_to_label(thisstc_tolabel_test2, forward['src'],smooth=False)[1]

        ##all ones label
#        thisstc_allone=np.ones((len(vertices_sub[0])+len(vertices_sub[1]),1))
#        thisstc_allone_tolabel = mne.SourceEstimate(thisstc_allone, vertices=vertices_sub,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject=subjects[subject_no])
#        allone_rhlabel_temp_sub=mne.stc_to_label(thisstc_allone_tolabel, forward['src'],smooth=False)[1]
#        label_temp_sub_rhpos=label_temp_sub.pos.copy()
#        label_temp_sub_rhpos[:,0]=-label_temp_sub_rhpos[:,0]
#        allone_rhlabel_temp_sub_rhpos=allone_rhlabel_temp_sub.pos.copy()
#        thisstc_rhdata=np.zeros((len(vertices_sub[0])+len(vertices_sub[1]),1))
#        
#        for rhcnt in range(label_temp_sub_rhpos.shape[0]):#for each vert in lh label find the closest in the rh
#            taken_verts=np.where(thisstc_rhdata==1)[0]
#            aa=label_temp_sub_rhpos[rhcnt,:]
#            bb=np.tile(aa,(allone_rhlabel_temp_sub_rhpos.shape[0],1))
#            cc=allone_rhlabel_temp_sub_rhpos-bb
#            dd=np.argsort(np.sqrt(np.sum(cc**2,1)))[0]+len(vertices_sub[0])
#            sortcnt=0
#            while dd in taken_verts:
#                sortcnt=sortcnt+1
#                dd=np.argsort(np.sqrt(np.sum(cc**2,1)))[sortcnt]+len(vertices_sub[0])
#            thisstc_rhdata[dd,0]=1
#            
#        thisstc_tolabel_rh = mne.SourceEstimate(thisstc_rhdata, vertices=vertices_sub,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject=subjects[subject_no])
#        labelrh_temp_sub=mne.stc_to_label(thisstc_tolabel_rh, forward['src'],smooth=False)[1]
#        rhlabel_verts_subm=np.where(thisstc_rhdata>0)[0]
#        label_temp_sub=thisstc_labeltest.morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path) 
#        allwin_verts[ii,this_labelc,:]=label_verts_subm.copy()
        if this_labelc==0:
            winvert_sub=label_verts_subm.copy()
        else:
            winvert_sub=np.hstack((winvert_sub,label_verts_subm))
#        allwin_verts[ii,2*this_labelc+1,:]=rhlabel_verts_subm.copy()
        testlabel_path=data_path+'dcm/dipoleDCMsublabs/'+'subject'+str(ii)+'_'+this_label.name
        label_temp_sub.save(testlabel_path)
#        testlabel_path=data_path+'dcm/dipoleDCMsublabs/'+'subject'+str(ii)+'_'+this_label.name[:-2]+'rh'
#        labelrh_temp_sub.save(testlabel_path)
        sub_labellist.append(label_temp_sub)
#        sub_labellist.append(labelrh_temp_sub)
#    winvert_sub=np.squeeze(allwin_verts[ii,:,:]).copy()
    vertsub_vect=winvert_sub.astype(int)#np.reshape(winvert_sub,(n_ROIs*np.abs(roivert),)).astype(int)
    print len(vertsub_vect)
    info = evoked.info
    print "a"
    invop, inverse_operator_meeg = make_inverse_operator_mine(info, forward, vertsub_vect, noise_cov,loose=0.2, depth=None)
    event_id=1
    epochs1 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,
    		             baseline=(None, 0), proj=True, reject=reject, preload=True)
    
    event_id = 2  # abstract
    epochs2 = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks, baseline=(None, 0), proj=True, reject=reject, preload=True)
    event_id={'cncrt_wrd': 1, 'abs_wrd': 2}
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, picks=picks,baseline=(None, 0), proj=True, reject=reject, preload=True)
    equalize_epoch_counts([epochs1, epochs2],method='mintime')
    snr = 3.0#np.mean(ediv)
    print snr
    lambda2 = 1.0 / snr ** 2
    method = 'MNE'
    evoked1 = epochs1.average()
    
    #evoked1 = whiten_evoked(evoked11, noise_cov)
    #evoked1.resample(50)
    condition1 = apply_inverse(evoked1, inverse_operator_meeg, lambda2, method,pick_ori="normal")
    
    evoked2 = epochs2.average()
    #evoked2 = whiten_evoked(evoked21, noise_cov)
    #evoked2.resample(50)
    condition2 = apply_inverse(evoked2, inverse_operator_meeg, lambda2, method,pick_ori="normal")
#    data_out = data_path + meg + 'notMorphed_ico_signed_dipDCM_SemLoc_ica_Concrete_Source_Evoked_m500_700'
#    condition1.save(data_out)
#    data_out = data_path + meg + 'notMorphed_ico_signed_dipDCM_SemLoc_ica_Abstract_Source_Evoked_m500_700'
#    condition2.save(data_out)
    for tlc, tlab in enumerate(sub_labellist[:-2:2]):
        thisscale=mne.label.label_sign_flip(tlab,src_sub)#np.hstack((mne.label.label_sign_flip(tlab,src_sub),mne.label.label_sign_flip(thisstc_label[1],src_sub)))[winner_verts]#mne.label.label_sign_flip(this_label,src_avg)[winner_verts]#1#Matx[winner_verts,this_labelc]/np.max(Matx[winner_verts,this_labelc])
#        print thisscale
        label_predata1=condition1.in_label(tlab).data
        #            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
        label_data1=np.tile(thisscale[:,np.newaxis],(1,label_predata1.shape[1]))*label_predata1
        ini_tc1=np.mean(label_data1, axis=0)
        label_predata2=condition2.in_label(tlab).data
        #            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
        label_data2=np.tile(thisscale[:,np.newaxis],(1,label_predata2.shape[1]))*label_predata2
        ini_tc2=np.mean(label_data2, axis=0)
        conabs_mat[tlc,:,0]=mne.extract_label_time_course(condition1,tlab,src_sub,mode='mean_flip')
        conabs_mat[tlc,:,1]=mne.extract_label_time_course(condition2,tlab,src_sub,mode='mean_flip')
#    X[ii,:,:,:]=conabs_mat
    out_file=out_path + meg[:10] + '_SemLoc_Evoked_5ROIs_dipole_exttc_'+str(np.abs(roivert))+'verts_avg.mat'
    scio.savemat(out_file,{'conabs_mat':conabs_mat})

#    conind=0
#    for event_no,event_name in enumerate(['Concrete','Abstract']):
##        conind=conind+1
##        print conind
#                # Load data
#        
#        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_'+event_name+'_Source_Evoked_m500_700' 
#        this_stc = mne.read_source_estimate(fname)  
#        tmin1=-500
#        tstep1=1
#        thisstc_data=np.ones((this_stc.data.shape[0],1))
#        thisstc_tolabel = mne.SourceEstimate(thisstc_data, vertices=vertices_avg,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#        thisstc_label=mne.stc_to_label(thisstc_tolabel, src_avg)
#
#        ntimes=this_stc.data.shape[1]
#        print "select 100 most sensitive labels"      
#        for this_labelc, this_label in enumerate(labellist): 
#            label_origdata=this_stc.in_label(this_label).data
#            label_verts=this_stc.in_label(this_label).lh_vertno
#            #            np.intersect1d(label_verts,vertices_to[0])
#            vert50=np.sort(Matx[label_verts,this_labelc])[roivert]
#            winner_verts=label_verts[np.where(Matx[label_verts,this_labelc]>=vert50)[0]]#label_verts[np.where(Matx[label_verts,this_labelc]>=np.max(Matx[label_verts,this_labelc])/2.)[0]]#np.where(Matx[:,this_labelc]>=np.max(Matx[:,this_labelc])/2.)[0]#
#            allwin_verts[ii,this_labelc,:]=winner_verts.copy()
#        winvert_sub=np.squeeze(allwin_verts[ii,:,:]).copy()
#        vertsub_vect=np.reshape(winvert_sub,(n_ROIs*np.abs(roivert),))
#    info = evoked.info
#    print "a"
#    inverse_operator_meeg = make_inverse_operator_mine(info, forward, noise_cov,loose=0.2, depth=None)












            
            #            print winner_verts
#            thisscale=np.hstack((mne.label.label_sign_flip(thisstc_label[0],src_avg),mne.label.label_sign_flip(thisstc_label[1],src_avg)))[winner_verts]#mne.label.label_sign_flip(this_label,src_avg)[winner_verts]#1#Matx[winner_verts,this_labelc]/np.max(Matx[winner_verts,this_labelc])
#            label_predata=this_stc.data[winner_verts,:]
#            #            label_data1=np.tile(thisscale,[ntimes,1]).transpose()*label_predata
#            label_data1=thisscale[:,np.newaxis]*label_predata
#            active_ones=np.where(np.mean(np.abs(label_data1[:,sonset:]),1)>=np.min(np.min(np.abs(label_data1[:,sonset:]),1)))[0]
#            print len(active_ones)
#            label_data=label_data1[active_ones,:]
#            ini_tc=np.mean(label_data[:,:], axis=0)
#            thisstc_testdata=np.zeros((this_stc.data.shape[0],1))
#            thisstc_testdata[winner_verts[active_ones],0]=1
#            vertices_test = [np.arange(10242), np.arange(10242)]
#            thisstc_tolabel_test = mne.SourceEstimate(thisstc_testdata, vertices=vertices_test,tmin=1e-3 * tmin1, tstep=1e-3 * tstep1, subject='fsaverage')
#            thisstc_labeltest=mne.stc_to_label(thisstc_tolabel_test, src_avg,smooth=False)[0]
#            thisstc_labelsubj = thisstc_labeltest.morph(subject_from='fsaverage', subject_to=subjects[subject_no], subjects_dir=data_path)
#
#            testlabel_path=data_path+'dcm/testDCMoldlabs/'+'subject'+str(ii)+'_'+event_name+'_'+this_label.name
#            thisstc_labeltest.save(testlabel_path)
##            conabs_mat[this_labelc,:,event_no]=ini_tc.copy()#(ini_tc-np.min(ini_tc))/(np.max(ini_tc)-np.min(ini_tc))zscore(ini_tc)#
#            conabs_mat[this_labelc,:,event_no]=mne.extract_label_time_course(this_stc,thisstc_labeltest,src_avg,mode='mean_flip')
##        wnw_mat=np.concatenate((np.expand_dims(cncrt_mat,2),np.expand_dims(abs_mat,2)),2)
##        con=mne.connectivity.spectral_connectivity(data=cncrt_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
##        conc=np.corrcoef(cncrt_mat[:,550:750])#con[0]+con[0].transpose(1,0,2)
##        con2=mne.connectivity.spectral_connectivity(data=abs_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
##        cona=np.corrcoef(abs_mat[:,550:750])#con2[0]+con2[0].transpose(1,0,2)
##        cons[ii,:,:]=np.squeeze(conc-cona)
#    X[ii,:,:,:]=conabs_mat
#    out_file=out_path + meg[:10] + '_SemLoc_Evoked_5ROIs_meanCTF_50verts_forDCM_avg.mat'
#    scio.savemat(out_file,{'conabs_mat':conabs_mat})
#xms=np.mean(X,0)
#xms=np.mean(xms,2)
#import matplotlib.pyplot as plt
#plt.figure(1),plt.plot(np.linspace(-500,700,1201),xms[0,:])
#plt.figure(2),plt.plot(np.linspace(-500,700,1201),xms[1,:])
#plt.figure(3),plt.plot(np.linspace(-500,700,1201),xms[2,:])
#plt.figure(4),plt.plot(np.linspace(-500,700,1201),xms[3,:])
#plt.figure(5),plt.plot(np.linspace(-500,700,1201),xms[4,:])
#            conind=conind+1
#cl_ideal[cntcl]=simstc_morphed_ideal[cntlts].in_label(labels[cntcl]).data
#clw_ideal[cntcl]=np.argmax(np.mean(np.abs(cl_ideal[cntcl][:,125:]),axis=1))
#thislabelsign2=mne.label.label_sign_flip(labels[cntcl],src_avg)[clw_ideal[cntcl]]
#label_ts_ideal[cntlts][cntcl,:]=thislabelsign2*np.expand_dims(cl_ideal[cntcl][clw_ideal[cntcl],:].transpose(),0)

#    print ii
#    srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
#    src_avg = mne.read_source_spaces(srcin)
#    fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Concrete_Source_Evoked_m500_700' 
#    stc_allc = mne.read_source_estimate(fname)
#    thislabel_tcc=mne.extract_label_time_course(stc_allc, thislabel,src, mode='mean_flip') 
#    labels_tcc[ii,:,thislabelc]=thislabel_tcc.copy()
#    
##        cl1=stc_allc.in_label(thislabel).data
##    cl1w=np.argmax(np.mean(np.abs(cl1[:,550:950]),axis=1))
##    print cl1w
#    
#    
#    fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Abstract_Source_Evoked_m500_700'  
#    stc_alla = mne.read_source_estimate(fname)
#    thislabel_tca=mne.extract_label_time_course(stc_alla, thislabel,src, mode='mean_flip')
#    labels_tca[ii,:,thislabelc]=thislabel_tca.copy()
#
#labelsa_forsign=np.squeeze(np.mean(labels_tca,axis=0)) #size time x ROI
#labelsc_forsign=np.squeeze(np.mean(labels_tcc,axis=0)) #size time x ROI
#    
#for thislabelc,thislabel in enumerate(labellist):
#    print thislabel
#    for ii, meg in enumerate(ll):
#        print ii
#        srcin='/imaging/rf02/TypLexMEG/fsaverage/bem/fsaverage-ico-5-src.fif'
#        src = mne.read_source_spaces(srcin)
#        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Concrete_Source_Evoked_m500_700' 
#        stc_allc = mne.read_source_estimate(fname)
#        thislabel_tcc=mne.extract_label_time_course(stc_allc, thislabel,src, mode='pca_flip') 
#        #        thislabel_signtc=mne.extract_label_time_course(stc_allc, thislabel,src, mode='mean_flip')
#        Pc=np.polyfit(np.linspace(50,200,150),thislabel_tcc[0,550:700],1)[0]
#        #        Pct=np.polyfit(np.linspace(50,200,150),thislabel_signtc[0,550:700],1)[0]
#        if np.corrcoef(np.squeeze(thislabel_tcc),labelsc_forsign[:,thislabelc])[0,1]>=0:#np.corrcoef(thislabel_tcc,thislabel_signtc)[0,1]>=0:#np.sign(Pc)==np.sign(Pct):
#            cncrt_mat[thislabelc,:]=thislabel_tcc.copy()
#            labels_tcc[ii,:,thislabelc]=thislabel_tcc.copy()
#        else:
#            cncrt_mat[thislabelc,:]=-thislabel_tcc.copy()
#            labels_tcc[ii,:,thislabelc]=-thislabel_tcc.copy()
#
#        #        cl1=stc_allc.in_label(thislabel).data
#        #    cl1w=np.argmax(np.mean(np.abs(cl1[:,550:950]),axis=1))
#        #    print cl1w
#        
#        
#        fname = data_path + meg + 'firstMorphed_ico_signed_SemLoc_ica_Abstract_Source_Evoked_m500_700'  
#        stc_alla = mne.read_source_estimate(fname)
#        thislabel_tca=mne.extract_label_time_course(stc_alla, thislabel,src, mode='pca_flip')
#        #        thislabel_signta=mne.extract_label_time_course(stc_alla, thislabel,src, mode='mean_flip')
#        Pa=np.polyfit(np.linspace(50,200,150),thislabel_tca[0,550:700],2)[0]
#        #        Pat=np.polyfit(np.linspace(50,200,150),thislabel_signta[0,550:700],2)[0]
#        if np.corrcoef(np.squeeze(thislabel_tca),labelsa_forsign[:,thislabelc])[0,1]>=0:#np.corrcoef(thislabel_tca,thislabel_signta)[0,1]>=0:#np.sign(Pa)==np.sign(Pat):
#            labels_tca[ii,:,thislabelc]=thislabel_tca.copy()
#            abs_mat[thislabelc,:]=thislabel_tca.copy()
#        else:
#            labels_tca[ii,:,thislabelc]=-thislabel_tca.copy()
#            abs_mat[thislabelc,:]=-thislabel_tca.copy()
        
#        wnw_mat=np.concatenate((np.expand_dims(cncrt_mat,2),np.expand_dims(abs_mat,2)),2)
#        con=mne.connectivity.spectral_connectivity(data=cncrt_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
#        conc=np.corrcoef(cncrt_mat[:,550:750])#con[0]+con[0].transpose(1,0,2)
#        con2=mne.connectivity.spectral_connectivity(data=abs_mat[np.newaxis,:,550:750],method='coh',sfreq=1000.,mode='multitaper',fmin=1.,fmax=45.,mt_adaptive=True,faverage=True)
#        cona=np.corrcoef(abs_mat[:,550:750])#con2[0]+con2[0].transpose(1,0,2)
#        cons[ii,:,:]=np.squeeze(conc-cona)
#        
#        out_file=out_path + meg[:10] + '_SemLoc_Evoked_5ROIs_PCA_forDCM.mat'
#        scio.savemat(out_file,{'wnw_mat':wnw_mat})

#consr=np.reshape(cons,(17,25))
#import scipy.stats as scist
#MTc= mne.stats.ttest_1samp_no_p(consr, sigma=1e-3, method='relative')#Tc,pTc,tTc= mne.stats.permutation_t_test(consr)
#pTc=scist.t.sf(np.abs(MTc),16)*2
#pfinal=np.reshape(pTc,(5,5))
#print pfinal
#import matplotlib.pyplot as plt
#plt.figure(3)
#for ii in range(17):
#    plt.subplot(4,5,ii+1) 
#    plt.plot(np.linspace(-500,700,1201),labels_tcc[ii,:,4].T)
#    
#plt.figure(4)
#for ii in range(17):
#    plt.subplot(4,5,ii+1) 
#    plt.plot(np.linspace(-500,700,1201),labels_tca[ii,:,4].T)