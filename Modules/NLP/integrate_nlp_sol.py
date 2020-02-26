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

import numpy as np
from Modules.Tools.plot_tools import state_separator

"""
Integration of the casadi NLP solution using the casadi 'sol' object. This function also can be used for a 
control u directly, if the according option is used.

Inputs: 
sol:                    Casadi 'sol' object, or list with control values u, if 'sol_is_u' is True
x0:                     Initial values of the dynamical system
allowed_arr:            Array of length N indication integration points in which a treatment is allowed
N:                      Number of grid points
dt:                     Integration stepsize
Nperday:                Number of integration points per day
Tf:                     End time of observed time horizon [0, Tf]
integrator_function:    Casadi Integrator function for NLP
integrator_function_2:  Casadi Integrator function used for the second part
                        of the end time optimization, 'None' if other objective is used
max_fraction:           Maximal fractional blood loss
p_in:                   Patient parameters (beta, gamma)
sol_is_u:               Using this option a known control u can be used instead of 'sol'.
                        Useful for integration outside of pv_schedule

Outputs: 
x1_opt:                 Optimal trajectory of x1
x2_opt:                 Optimal trajectory of x2
x3_opt:                 Optimal trajectory of x3
q_opt:                  Objective value of optimal solution
u_opt:                  Optimal control function for allowed time points
tgrid:                  Time grid of trajectories for ploting

"""

def integrate_nlp_sol(sol, x0, allowed_arr, N, dt, Nperday, Tf, integrator_function, integrator_function_2, max_fraction, p_in, sol_is_u=False):
    # check whether sol is the casadi solution object or the control directly
    if sol_is_u:
        u_opt = sol
    else:
        w1_opt = sol['x']
        num_controls = np.sum(allowed_arr)
        w1_opt = w1_opt.full().flatten()
        u_opt = w1_opt[-num_controls:]
    
    tgrid = [Tf/N*k for k in range(N+1)]

    # integrate solution to obtain trajectories
    x_opt = [x0]
    q_opt = [0]
    u_idx = 0

    hour_idx = 0
    day_idx = 0
    current_time = 0

    for k in range(N):
        allowed = allowed_arr[k]

        if allowed>=1e-8:
            i_out = integrator_function(x0 = x_opt[-1], q0 = q_opt[-1], u = u_opt[u_idx], p=p_in)
            # include jump if control > 0
            i_out['xf'][5] = i_out['xf'][5]*(1 - u_opt[u_idx]*max_fraction)
            if u_idx < len(u_opt) - 1:
                u_idx+=1
        else:
            i_out = integrator_function(x0 = x_opt[-1], q0 = q_opt[-1], u = 0, p=p_in)

        # update times
        hour_idx +=1
        if hour_idx == Nperday:
            hour_idx = 0
        day_idx +=1
        if day_idx == 7:
            day_idx = 0
        current_time += dt

        x_opt.append(i_out['xf'][3:6])
        q_opt.append(i_out['li'][1])

    
    # integer end point extension if applicable
    if integrator_function_2 is not None:
        scale = w1_opt[-num_controls-1]
        dt_two_stage = 1./20
        retransformed_stepsize = dt_two_stage / scale
        for k in range(20):
            i_out = integrator_function_2(x0 = x_opt[-1], q0 = q_opt[-1], u = 0, p=p_in, scale = 1./scale)
            x_opt.append(i_out['xf'][3:6])
            q_opt.append(i_out['li'][1])    
            tgrid.append(tgrid[-1] + retransformed_stepsize)   
            
    # separation of states in suitable format for plotting routine
    x1_opt, x2_opt, x3_opt = state_separator(x_opt)   
    return x1_opt, x2_opt, x3_opt, q_opt, u_opt, tgrid

