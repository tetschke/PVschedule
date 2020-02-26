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


This routine uses the dynamic programming method for generation of treatment schedules for pv patients
"""



import sys
# path to casadi if not in PYTHONPATH
sys.path.append(r'/home/lilienthal/Programmieren/casadi-linux-py36-v3.4.5-64bit/')
# path to pv_schedule
sys.path.append('../')

from Modules.Tools.patient_parameters import return_parameters
from Modules.Model.model_integrator import model_integrator
from Modules.Tools.plot_tools import plot_sol
from Modules.DynamicProg.dynamic_programming import dynamic_prog_pv


# length of time horizon (days)
Tf = 365   
# integration points per day
Nperday = 6


# Index of Patient or Subject; see Modules/Tools/patient_parameters.py for more details
index = 'F02' 
pv_lambda_index = 1

# allowed array options similar to pv_schedule
allowed_opts = {'allowed_hours': [1, 0, 0, 0, 0, 0], 'allowed_days': [1, 1, 1, 1, 1, 0, 0],
                'forbidden_days': [range(81, 96), range(280, 302)]
                }

# dynamic programming options
dp_options = {'NX':101, 'NU':2, 'NK':6}
#dp_options = {'NX':101, 'NU':2, 'NK':6, 'p_shift': 0.4}  # shift of rounding rule

# Get patient parateters
gamma, beta, Base, patient_volume, x0, pv_lambda = return_parameters(index, lambda_version=pv_lambda_index)
p_in = [beta, gamma]

# upper bound for x3
x3_up = 1.1 * Base

# number of integration points for integrator function
N = int(Tf * Nperday)    

# stepsize of equidistant integrator
dt = 1. / Nperday    

# maximal fractional blood loss
max_donation_volume = 500.
max_fraction = max_donation_volume / patient_volume

# get integrator function
integrator_function = model_integrator(N, dt, Tf, Nperday, Base, max_fraction, pv_lambda, False)

# execute dynamic programming
x1_dp, x2_dp, x3_dp, q_dp, dp_u_plot, tgrid = dynamic_prog_pv(Tf, Nperday, x0, Base, pv_lambda, p_in, patient_volume, allowed_opts,
                     max_fraction, x3_up,  dp_options)

# plot solution
plot_sol(x1_dp, x2_dp, x3_dp, q_dp, dp_u_plot, tgrid, [1]*N, 1, input_names=['dp'], title='DynamicProgramming',
         show_forbidden=True, show_plot=True, limit=x3_up)


