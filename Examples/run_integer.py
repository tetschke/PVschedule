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


This routine uses the integer end point method to compute integer pv schedules
based on an iterative scheme using BONMIN with Casadi
"""



import sys
# path to casadi if not in PYTHONPATH
sys.path.append(r'/home/lilienthal/Programmieren/casadi-linux-py36-v3.4.5-64bit/')
# path to pv_schedule
sys.path.append('../')

from pv_schedule import pv_schedule
from Modules.Tools.patient_parameters import return_parameters
from Modules.Tools.plot_tools import plot_sol

# global settings
Tf = 103.0 
Nperday = 6 

# get patient parameters; for more information see Modules/Tools/patient_parameters.py
patient_index = 'F20'
lambda_version = 1   
gamma, beta, base, patient_volume, x0, pv_lambda = return_parameters(patient_index, lambda_version = lambda_version)

# lists for storing trajectories and other information for plotting
x1 = []
x2 = []
x3 = []
q = []
u = []
tgrids = []
allowed_arrs = []
input_names = []
marker_styles = []

# option snipplets for later use
default_options = {'obj_factor':1}
# experiment configurations
allowed_days = [1, 0, 0, 0, 0, 0, 0]
allowed_hours = [1, 0, 0, 0, 0, 0]   
forbidden_days = []
allowed_options = {'allowed_hours': allowed_hours, 'allowed_days': allowed_days, 'forbidden_days': forbidden_days}

# parameter array
p_in = [beta, gamma]    

# length of time horizon for heuristic
Tf_heuristic = int(Tf*2.5) 

########################### HEURISTIC with constraints ####################
options = {'objective':'heuristic'}
options.update(default_options)
options.update(allowed_options)
_,_,_,_,_,_,_,allowed_old, _ = pv_schedule(Tf, Nperday, x0, base, pv_lambda, p_in, patient_volume, options)
x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid, sol, allowed_arr, error_flag = pv_schedule(Tf_heuristic, Nperday, x0, base, pv_lambda, p_in, patient_volume, options)


if error_flag == 0:
    x1.append(x1_opt)
    x2.append(x2_opt)
    x3.append(x3_opt)
    q.append(q_opt)
    u.append(u_opt)
    tgrids.append(tgrid)
    allowed_arrs.append(allowed_arr)
    input_names.append('Heuristic grid')
    marker_styles.append('+')
    
    # get number of treatments until Tf
    u_till_tf = u_opt[:sum(allowed_old)]
    u_up = int(sum(u_till_tf) + 1e-4)
    valid_k = True

else:
    print('Heuristic on grid failed: treatment density too high')
    Tf2_heur_grid = 9e5
    valid_k = False

####################### HEURISTIC w/o restrictions ########################

options = {'objective':'heuristic'}
options.update(default_options)
_,_,_,_,_,_,_,allowed_old, error_flag = pv_schedule(Tf, Nperday, x0, base, pv_lambda, p_in, patient_volume, options)
x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid, sol, allowed_arr, error_flag = pv_schedule(Tf_heuristic, Nperday, x0, base, pv_lambda, p_in, patient_volume, options)

x1.append(x1_opt)
x2.append(x2_opt)
x3.append(x3_opt)
q.append(q_opt)
u.append(u_opt)
tgrids.append(tgrid)
allowed_arrs.append(allowed_arr)
input_names.append('Heuristic free')
marker_styles.append('None')

# get number of treatments until Tf
u_till_tf = u_opt[:sum(allowed_old)]
u_lo = int(sum(u_till_tf) + 1e-4)

# Start integer end point method at k = u_up
Tf2_bonmin = -9001
u_min = 9001
k = u_up 
print('Using upper bound of ', k)

while valid_k:
    print('Using k=', k)

    options = {'objective': 'integer_end_point', 'u_max':k}
    options.update(default_options)
    options.update(allowed_options)

    x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid, sol, allowed_arr, error_flag = pv_schedule(Tf, Nperday, x0, base, pv_lambda, p_in, patient_volume, options)
    
    x1.append(x1_opt)
    x2.append(x2_opt)
    x3.append(x3_opt)
    q.append(q_opt)
    u.append(u_opt)     
    tgrids.append(tgrid)
    allowed_arrs.append(allowed_arr)
    

    # update minimal number of donations; u_opt > 0 should catch infeasible bonmin solutions as no control is applied there
    print(tgrid[-1], sum(u_opt))
    if tgrid[-1] >= Tf2_bonmin and sum(u_opt) > 0:
        u_min = sum(u_opt)
        Tf2_bonmin = tgrid[-1]
    
    # if sum(u_opt) != k, then the try with lower k failed and the loop should break
    if sum(u_opt) != k:
        valid_k = False
        print("Break process, as control is not valid anymore")
    else:
        k = k - 1
    input_names.append('Bonmin_' + str(patient_index) + str(lambda_version) + '_' + str(k))
    marker_styles.append('None')


# plot results 
plot_sol(x1, x2, x3, q, u, tgrids, allowed_arrs, len(input_names), input_names = input_names,
         title='Bonmin', show_forbidden=True, show_plot=True, limit=1.1*base, marker_styles=marker_styles)






