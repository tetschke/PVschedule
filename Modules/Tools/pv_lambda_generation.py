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


Routine for generation of plausible pv_lambda parameters for given healthy subjects data.
"""

import sys
# path to pv_schedule
sys.path.append('../../')
from pv_schedule import pv_schedule
from Modules.Tools.patient_parameters import return_parameters
import numpy as np


Tf = 365    
Nperday = 6

# amount of pv_lambdas per subject generated
num_pv_lambdas = 5

# Patient indices
freiburg_indices = ['F01', 'F02', 'F03', 'F04', 'F05',
            'F06', 'F07', 'F08', 'F09', 'F10',
            'F11', 'F12', 'F13', 'F14', 'F15',
            'F16', 'F17', 'F18', 'F19', 'F20',
            'F21', 'F23', 'F24', 'F25',
            'F26', 'F27', 'F28', 'F29']
freiburg_indices = freiburg_indices * num_pv_lambdas

# lower and upper bound of treatments
u_lim_lo = 1
u_lim_up = int(Tf/14)

# dictionary options are those from the pv_schedule function
dict_opts = {'objective': 'heuristic'}

# Lists for storing number of Treatments with according pv_lambdas
u_pat = []
pv_lambdas = []

for idx in freiburg_indices:

    num_iter = 0
    is_successful = False
    
    while not is_successful and num_iter < 50:
        
        # choose random lambda between 0 and 1
        pv_lambda = np.random.rand()
        
        # Get patient parameters
        gamma, beta, Base, patient_volume, x0, _ = return_parameters(idx)
        
        p_in = [beta, gamma]
        
        # upper bound for x3
        x3_up = 1.1 * Base
        x3_lo = 0.8 * Base
       
        # call code for given setting and subject
        x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid, sol, allowed_arr, error_flag = \
        pv_schedule(Tf, Nperday, x0, Base, pv_lambda, p_in, patient_volume, dict_opts)
        
        # Count and check calculated number of treatments
        num_treats = sum(u_opt)
        if num_treats < u_lim_lo:
            print('pv_lambda failed at lower bound: ', pv_lambda, num_treats)
            num_iter += 1
        elif num_treats > u_lim_up:
            print('pv_lambda failed at upper bound: ', pv_lambda, num_treats)
            num_iter += 1
        else:
            is_successful = True
    
    print('pv_lambda successful', pv_lambda, num_treats)
    pv_lambdas.append(pv_lambda)
    u_pat.append(num_treats)


# =============================================================================
# # count entries in num_treats and plot distribution
# count_array = np.zeros(u_lim_up - u_lim_lo + 1)
# value_arr = np.arange(u_lim_lo, u_lim_up+1e-8, 1)
# 
# for n in u_pat:
#     count_array[n - u_lim_lo] += 1
# 
# plt.plot(value_arr, count_array)
# plt.xlabel('Number of Treatments')
# plt.ylabel('Number of Konfigurations')
# plt.grid()
# =============================================================================

