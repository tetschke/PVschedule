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


"""

from Modules.NLP.nlp_builder import nlp_builder
from Modules.NLP.integrate_nlp_sol import integrate_nlp_sol
from Modules.Heuristic.heuristic_alg import pv_heuristic_alg
from Modules.Model.allowed_generator import allowed_generator
from Modules.Model.model_integrator import model_integrator

import numpy as np
import casadi as ca


"""
Main routine for generation of pv treatment schedules using the heuristic approach,
a relaxed NLP formulation minimizing the control integral and an integer end point method maximizing Tf.

Inputs:
Tf:             End point of observed interval [0, Tf]
Nperday:        Number of integration points per day
x0:             Initial value of x
B:              Steady state value of x3
pv_lambda:      Patient parameter pv_lambda
p_in:           Patient parameters [beta, gamma]
patient_volume: Total blood volume of patients
dict_opts:      Dictionary specifiying more options:
    objective:              Determines the approach for generation of pv schedules:
        'heuristic' :           Using heuristic algorithm
        'relaxed_int_u' :       Using the integral of the control u as objective of an NLP formulation of the problem.
                                This problem is solved as a relaxed problem using IPOPT
        'integer_end_point':    Using end point optimzation on an NLP maximizing Tf. This problem is solved using BONMIN
    max_treatment_volume:   Maximal treatment volume in ml                      -> default: 500             ,other: any float > 0
    obj_factor:             Scaling factor of objective                         -> default: 10              ,other: any float > 0
    allowed_hours:          Daily allowed treatment time sections               -> default: [1]*Nperday     ,other: 0/1 list with length = Nperday
    allowed_days:           Weekly allowed treatment days                       -> default: [1]*7           ,other: 0/1 list with length = 7
    forbidden_days:         Absolute days in which a treatment is not allowed   -> default: None            ,other: 0/1 list with integers (single days) or range(i, j+1) (forbidden from day i to day j)
    u_start:                Initial control trajectory for optimization         -> default: zero control    ,other: list of numbers valid for control
                                    
    'u_max':                number of treatments used for integer end point method

Outputs:
x1_opt:         Optimal trajectory for x1
x2_opt:         Optimal trajectory for x2
x3_opt:         Optimal trajectory for x3
q_opt:          Objective value
u:              Optimal control
tgrid:          Time grid of trajectories for plotting
sol:            Casadi NLP solution object, {} for 'heuristic'
allowed_arr:    Returns list of allowed times formatted for plotting
error_flag:     Error flag of heuristic approach, 0 for other methods
 
"""    

def pv_schedule(Tf, Nperday, x0, B, pv_lambda, p_in, patient_volume, dict_opts):
    
    N = int(Tf * Nperday)
    dt = Tf/N
    
    # set maximal treatment volume, default is 500 ml
    if 'max_treatment_volume' in dict_opts.keys():
        max_treatment_volume = dict_opts['max_treatment_volume']
    else:
        max_treatment_volume = 500
    # maximal fractional blood removal
    max_fraction = max_treatment_volume / patient_volume
    

    # Allowed treatment times
    if 'allowed_days' in dict_opts.keys() or 'allowed_hours' in dict_opts.keys() or 'forbidden_days' in dict_opts.keys():
        allowed_arr = allowed_generator(Nperday, dict_opts, N, dt)
    else:
        allowed_arr = [1]*N

    # initial value for control, default is zero control
    if 'u_start' in dict_opts.keys():
        u_start = dict_opts['u_start']
    else:
        u_start = [ca.DM(0.)] * len(allowed_arr)
    
    # obtaining objective
    if 'objective' in dict_opts.keys() and dict_opts['objective'] == 'integer_end_point':
        objective = 'integer_end_point'
    elif 'objective' in dict_opts.keys() and dict_opts['objective'] == 'heuristic':
        objective = 'heuristic'
    else: # default
        objective = 'relaxed_int_u'
    
    if objective == 'integer_end_point':
        # model formulation and integrator
        integrator_function, integrator_function_2 = model_integrator(N, dt, Tf, Nperday, B, max_fraction, pv_lambda,  two_stage=True)    
    else:   
        # model formulation and integrator
        integrator_function = model_integrator(N, dt, Tf, Nperday, B, max_fraction, pv_lambda)    
        integrator_function_2 = None
        
    if objective == 'heuristic':
        x1_opt, x2_opt, x3_opt, q_opt, u, tgrid, error_flag = pv_heuristic_alg(Tf, N, np.array(x0), integrator_function, p_in, max_fraction, B, allowed_arr)
        return x1_opt, x2_opt, x3_opt, q_opt, u, tgrid, {}, allowed_arr, error_flag        
    else:
        # NLP based problem
        # maximal number of donations, default is none
        if 'u_max' in dict_opts.keys():
            u_max = dict_opts['u_max']
        else:
            u_max = None    
        Q, w, w0, g, lbw, ubw, lbg, ubg, discrete = nlp_builder(N, Nperday, 0.05, B, x0, max_fraction, allowed_arr, integrator_function, integrator_function_2, p_in, u_start, u_max)
            
        # build NLP
        # Solve problem using NLP solver
        if 'obj_factor' in dict_opts.keys():
            obj_factor = dict_opts['obj_factor']
            print('Using objective factor')
        else:
            obj_factor = 10
            
        # problem formulation
        nlp_prob = {'f': obj_factor*Q, 'x': w, 'g': g}

        # define and run NLP solver
        if objective == 'integer_end_point':
            bonmin_options = {'variable_selection': 'most-fractional', 'tree_search_strategy':'dive'} # options used in paper
            # bonmin_options = {'variable_selection': 'nlp-strong-branching', 'tree_search_strategy':'dive'}  
            # bonmin_options = {}   
            nlp_solver = ca.nlpsol('nlp_solver', 'bonmin', nlp_prob, {"discrete": discrete, "bonmin": bonmin_options});            

        else: # objective == 'relaxed_int_u'
            ipopt_opts = {}  
            nlp_solver = ca.nlpsol('nlp_solver', 'ipopt', nlp_prob, {"ipopt": ipopt_opts});
        sol = nlp_solver(x0=ca.vertcat(*w0), lbx=lbw, ubx=ubw, lbg=lbg, ubg=ubg)
        
        # Integrate solution object to obtain trajectories
        x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid = integrate_nlp_sol(sol, x0, allowed_arr,  N, dt, Nperday, Tf, integrator_function, integrator_function_2, max_fraction, p_in)
        return x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid, sol, allowed_arr, 0 


