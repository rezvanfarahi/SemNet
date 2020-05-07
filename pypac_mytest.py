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
sys.path.append('/imaging/local/software/python_packages/scikit-learn/v0.14.1')
#sys.path.insert(1,'/imaging/local/software/mne_python/latest_v0.9')#/imaging/rf02/TypLexMEG/mne_python0.9
sys.path.insert(1,'/imaging/local/software/python')
sys.path.insert(1,'/imaging/local/software/python_packages/scikit-learn/v0.13.1/lib.linux-x86_64-2.7/sklearn/externals')
sys.path.insert(1,'/imaging/rf02/Semnet/pacpy-1.0.3.1')
import numpy as np
from scipy.signal import hilbert
from pacpy.pac import plv

t = np.arange(0, 10, .001) # Define time array
lo = np.sin(t * 2 * np.pi * 6) # Create low frequency carrier
hi = np.sin(t * 2 * np.pi * 100) # Create modulated oscillation
hi[np.angle(hilbert(lo)) > -np.pi*.5] = 0 # Clip to 1/4 of cycle

pp=plv(lo, hi, (4,8), (80,150)) # Calculate PAC