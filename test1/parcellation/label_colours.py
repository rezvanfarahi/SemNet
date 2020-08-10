# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 14:47:03 2017

@author: rf02
"""

from __future__ import division

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors


#colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

# Sort colors by hue, saturation, value and name.
#by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
#                for name, color in colors.items())
#sorted_names = [name for hsv, name in by_hsv]
sorted_names=copy.deepcopy(label_names)
colors=copy.deepcopy(label_colors)
n = len(sorted_names)
ncols = 1
nrows =  1

fig, ax = plt.subplots(figsize=(8, 9.5))
fig.set_facecolor('black')

# Get height and width
X, Y = fig.get_dpi() * fig.get_size_inches()
h = 10#Y / (nrows + 1)
w = 400#X / ncols

for i, name in enumerate(sorted_names):
    col = i % ncols
    row = i // ncols
    y = Y - (row * h) - h

    xi_line = w * (col )
    xf_line = w * (col + 0.1)
    xi_text = w * (col + 0.3)

#    ax.text(xi_text, y, name, fontsize=11,
#            horizontalalignment='left',
#            verticalalignment='center')

    ax.hlines(y , xi_line, xf_line,
              color=colors[i], linewidth=(h * 0.8))

ax.set_xlim(0, X)
ax.set_ylim(0, Y)
ax.set_axis_off()

fig.subplots_adjust(left=0, right=1,
                    top=1, bottom=0,
                    hspace=0, wspace=0)
plt.show()