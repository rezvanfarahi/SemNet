# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:55:11 2015

@author: rf02
"""

# Author: Martin Luessi <mluessi@nmr.mgh.harvard.edu>
#
# License: BSD (3-clause)

print(__doc__)
print "hi! this script is for spectral functional connectivity. Very nice maps are fruits which worth waiting for I bet!"
import sys

#sys.path.append('/imaging/local/software/python_packages/nibabel/1.3.0')
#sys.path.append('/imaging/local/software/python_packages/pysurfer/v0.4')
#sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.append('/imaging/local/software/mne_python/latest')

# for qsub
# (add ! here if needed)
sys.path.append('/imaging/local/software/EPD/latest/x86_64/bin/python')
sys.path.append('/imaging/local/software/mne_python/latest_v0.9full')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.16.1')

import sklearn 