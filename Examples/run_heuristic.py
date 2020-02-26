#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 10:06:53 2019

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


This routine generates pv schedules of a chosen patient using
- Heuristic algorithm with restrictions
- Relaxed NLP algorithm 
- Sumup method applied on solution of the relaxed NLP algorithm

"""

import sys
# path to casadi if not in PYTHONPATH
sys.path.append(r'/home/lilienthal/Programmieren/casadi-linux-py36-v3.4.5-64bit/')
# path to pv_schedule
sys.path.append('../')

from pv_schedule import pv_schedule
from Modules.Tools.patient_parameters import return_parameters
from Modules.Model.model_integrator import model_integrator
from Modules.NLP.integrate_nlp_sol import integrate_nlp_sol
from Modules.Sumup.sumup_alg import sumup
from Modules.Tools.plot_tools import plot_sol

# length of time horizon (days)
Tf = 365
Nperday = 6

# Patient parameters; for more information see Modules/Tools/patient_parameters.py
index = 'F02'
pv_lambda_index = 1
gamma, beta, Base, patient_volume, x0, pv_lambda = return_parameters(index, lambda_version=pv_lambda_index)
p_in = [beta, gamma]

# dictionary options 
dict_opts = {}

dict_allowed = {'allowed_hours': [0, 0, 1, 1, 1, 0], 'allowed_days': [1, 1, 1, 1, 1, 0, 0],
                'forbidden_days': [range(81, 96), range(280, 302)]
                }
dict_opts.update(dict_allowed)

# upper bound for x3
x3_up = 1.1 * Base
x3_lo = 0.8 * Base

# number of integration points;
N = int(Tf * Nperday)    

# stepsize of equidistant integrator
dt = 1. / Nperday    

# max fractional blood loss
max_donation_volume = 500.
max_fraction = max_donation_volume / patient_volume

# get integrator function
integrator_function = model_integrator(N, dt, Tf, Nperday, Base, max_fraction, pv_lambda)



# Lists for storing trajectories and other plot information
x1 = []
x2 = []
x3 = []
q = []
u = []
tgrids = []
allowed_arrs = []
input_names = []
marker_styles = []
line_width = []

 # -------------------------------------------------------------------
 ################## HEURISTIC #####################
heur_options = {'objective': 'heuristic'}
heur_options.update(dict_allowed)


x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid, sol, allowed_arr, err_flag = \
    pv_schedule(Tf, Nperday, x0, Base, pv_lambda, p_in, patient_volume, heur_options)
    
if not err_flag:
    x1.append(x1_opt)
    x2.append(x2_opt)
    x3.append(x3_opt)
    q.append(q_opt)
    u.append(u_opt)
    tgrids.append(tgrid)
    allowed_arrs.append(allowed_arr)
    input_names.append('Heuristic method')
    marker_styles.append(None)
    line_width.append(3.0)

################## RELAXED #####################
# Relaxed solution without restrictions with sumup rule

heur_options = {'objective': 'relaxed_int_u'}
heur_options.update(dict_allowed)
x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid, sol, allowed_arr, err_flag = \
    pv_schedule(Tf, Nperday, x0, Base, pv_lambda, p_in, patient_volume, dict_opts)

# Append entries
x1.append(x1_opt)
x2.append(x2_opt)
x3.append(x3_opt)
q.append(q_opt)
u.append(u_opt)
tgrids.append(tgrid)
allowed_arrs.append(allowed_arr)
input_names.append('Relaxed method')
marker_styles.append(None)
line_width.append(1.0)


################## SUMUP #####################

# Generation of sumup control from relaxed control
sum_u_opt = sumup(u_opt)

# Integration of system using sumup control
x1_sumup, x2_sumup, x3_sumup, q_sumup, u_sumup, tgrid = \
    integrate_nlp_sol(sum_u_opt, x0, allowed_arr, N, dt, Nperday, Tf, integrator_function, None,
                      max_fraction, p_in, sol_is_u=True)

# Append entries
x1.append(x1_sumup)
x2.append(x2_sumup)
x3.append(x3_sumup)
q.append(q_sumup)
u.append(u_sumup)
tgrids.append(tgrid)
allowed_arrs.append(allowed_arr)
input_names.append('Sumup method')
marker_styles.append(None)
line_width.append(1.0)

# Plot solution
plot_sol(x1, x2, x3, q, u, tgrids, allowed_arrs, len(input_names), input_names=input_names, title='pyCombina',
         show_forbidden=True, show_plot=True, limit=x3_up, marker_styles=marker_styles, line_width=line_width)

