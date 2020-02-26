#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# This file is part of PVschedule.
#
# Copyright 2019-2020 Patrick Lilienthal, Manuel Tetschke and Sebastian Sager
#
# PVschedule is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PVschedule is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with PVschedule. If not, see <http://www.gnu.org/licenses/>.


This file contains the standard plotting routine for PVschedule and a separation 
algorithm for Casadi NLP output

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm


"""
Routine for plotting single or multiple outputs returned by PVschedule methods.
All inputs with (*):
for num_inputs == 1 its one single numpy array,
for num_inputs > 1 its a list of numpy arrays
 
Inputs:
x1_opt:             Optimal trajectory for x1 (*)
x2_opt:             Optimal trajectory for x2 (*)
x3_opt:             Optimal trajectory for x3 (*)
q_opt:              Objective trajectory (*)
u_opt:              Optimal control trajectory (*)
tgrids:             Time grid of trajectories (*)
allowed_arr:        List of allowed treatment times (*)
num_inputs:         Number of different plots
input_names:        List with names of different plots (optional)
title:              Super title of plot (optional)
color_inputs:       List of colors of trajectories in matplotlib.pyplot syntax (optional)
colorbar_type:      List of types of colorbar (optional):
    None:               Do not create any colorbar (default)
    'linear':           Create color bar directly using input values
    'log':              Create color bar using logarithm with base 10 on input values
show_forbidden:     Show forbidden treatment times in control plot (Optional)
show_plot:          Show plot at end of execution (optional)
limit:              Number for upper bound which is displayed in x3 plot (optional)
marker_styles:      List of marker styles of trajectories (optional)
line_width:         List of line widths of trajectories (optional)

Outputs:
Function has no return values. 
"""


def plot_sol(x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrids, allowed_arr, num_inputs=1,
             input_names=None, title=None, color_inputs=None, colorbar_type=None, show_forbidden=True,
             show_plot=False, limit=None, marker_styles=None, line_width=None):
    
    # for consistency: single array input will be put into list 
    if num_inputs == 1:
        x1_opt = [x1_opt]
        x2_opt = [x2_opt]
        x3_opt = [x3_opt]
        q_opt = [q_opt]
        u_opt = [u_opt]
        tgrids= [tgrids]
        allowed_arr= [allowed_arr]
        
        
        
    # Get marker styles and line width if specified
    if marker_styles is None:
        marker_styles = ['None'] * num_inputs
    if line_width is None:
        line_width = [1.0] * num_inputs

    # Set legend by inputs or by default
    if input_names is not None:
        legend_arr = input_names
    else:
        legend_arr = []
        for i in range(num_inputs):
            legend_arr.append(str(i+1))


    # create multi plot 
    f, [p1, p2, p3, p4, p5] = plt.subplots(5, 1, sharex=True)
    p1.grid(True)
    p1.set_title('State x1')    
    p2.grid(True)
    p2.set_title('State x2')    
    p3.grid(True)
    p3.set_title('State x3')    
    p4.grid(True)
    p4.set_title('Control')
    p5.grid(True)
    p5.set_title('Objective')
    
    # set super title if specified
    if title is not None:
        f.suptitle(title)
    
    # plot upper bound of x3 if existing
    if limit != None:
        p3.axhline(limit, color='gray', linestyle='dashed')
    
    # Set color options and color bar if specified
    if color_inputs is not None:
        colorparams = color_inputs
        colormap = cm.jet
        if colorbar_type == 'linear':
            normalize = mcolors.Normalize(vmin=np.min(colorparams), vmax=np.max(colorparams))    
        else: # colorbar_type = log
            normalize = mcolors.LogNorm(vmin=np.min(colorparams), vmax=np.max(colorparams))
        s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
        s_map.set_array(colorparams)

    # Convert and plot each output
    for i in range(num_inputs):
        if color_inputs is not None:
            color = colormap(normalize(color_inputs[i]))

        x1_i = x1_opt[i]
        x2_i = x2_opt[i]
        x3_i = x3_opt[i]
        u_i = u_opt[i]
        tgrid = tgrids[i]
        q_i = q_opt[i]
        
        # Get suitable time grids for allowed and forbidden times
        u_time = [tgrid[j] for j in range(0, len(allowed_arr[i])) if allowed_arr[i][j] >= 1e-8]
        forbid_time = [tgrid[j] for j in range(0, len(allowed_arr[i])) if allowed_arr[i][j] == 0]

        # Plot trajectories
        if color_inputs is not None:
            p1.plot(tgrid, x1_i, color=color, marker = marker_styles[i], linewidth=line_width[i])
            p2.plot(tgrid, x2_i, color=color, marker = marker_styles[i], linewidth=line_width[i])
            p3.plot(tgrid, x3_i, color=color, marker = marker_styles[i], linewidth=line_width[i])
            p5.plot(tgrid, q_i, color=color, marker = marker_styles[i], linewidth=line_width[i])
        else:
            p1.plot(tgrid, x1_i, marker = marker_styles[i], linewidth=line_width[i])
            p2.plot(tgrid, x2_i, marker = marker_styles[i], linewidth=line_width[i])
            p3.plot(tgrid, x3_i, marker = marker_styles[i], linewidth=line_width[i])
            p5.plot(tgrid, q_i, marker = marker_styles[i], linewidth=line_width[i])

        
        # Check: due to rounding last entry of u_time might not be needed
        if len(u_time) > len(u_i):
            u_time = u_time[:-1]
            
        # Plot control
        if color_inputs is not None:
            p4.plot(u_time, u_i, color=color, marker='x', linestyle='None')
        else:
            p4.plot(u_time, u_i, marker='x', linestyle='None')
        
        # Plot forbidden times
        if i == 0 and show_forbidden:
            for f in forbid_time:
                p4.axvline(f, color='red', alpha=0.3 )
    # Plot color bar
    if color_inputs is not None and colorbar_type is not None:
        plt.colorbar(s_map, orientation='horizontal')
        
    # Plot legend in plot 1
    p1.legend(legend_arr)
    
    # Show plot 
    if show_plot:
        plt.show()
    
"""
Separate state x obtained from casadi solution into plotable trajectories

Input:
x_in:   List of Casadi state vectors

Outputs:
x1_out: Numpy array with trajectory for x1
x2_out: Numpy array with trajectory for x2
x3_out: Numpy array with trajectory for x3
"""   

def state_separator(x_in):
    x1_out = []
    x2_out = []
    x3_out = []

    for x in x_in:
        x1_out.append(x[0])
        x2_out.append(x[1])
        x3_out.append(x[2])

    # conversion to numpy array
    x1_out = np.array(x1_out)
    x2_out = np.array(x2_out)
    x3_out = np.array(x3_out)    

    return x1_out, x2_out, x3_out








