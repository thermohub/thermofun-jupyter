# Copyright (C) 2022 dmiron
# 
# This file is part of pitzer-plots.
# 
# pitzer-plots is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# pitzer-plots is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with pitzer-plots.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib as mpl
from collections import OrderedDict

def set_plot_dimensions():
    # plot dimensions and settings
    mpl.rcParams['lines.linewidth']=2
    mpl.rcParams['axes.labelsize']=17
    mpl.rcParams['axes.linewidth']=2
    mpl.rcParams['lines.linewidth']=2
    mpl.rcParams['lines.markersize']=10
    mpl.rcParams['font.size']=17
    mpl.rcParams['legend.labelspacing']=0.2
    # linewidth=2.5
    # markersize=15
    
    mpl.rcParams['figure.figsize']=[7,7]
    mpl.rcParams['figure.subplot.left'] = 0.13
    mpl.rcParams['figure.subplot.bottom'] = 0.12
    mpl.rcParams['figure.subplot.right'] = 0.96
    mpl.rcParams['figure.subplot.top'] = 0.95

def legend_(plt_, location, tit=''):
    handles, labels = plt_.gca().get_legend_handles_labels()
    #handles = reversed(handles)
    #labels = reversed(labels)
    by_label = OrderedDict(zip(labels, handles))
    legend1 = plt_.legend(by_label.values(), by_label.keys(), ncol=1, fontsize=16, loc=location, title=tit)