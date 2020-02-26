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
import sys
sys.path.append('./Modules/Integrator')

import numpy as np
from Modules.Tools.plot_tools import state_separator

"""
Heuristic algorithm for generation of treatment schedules according to Algorithm 1
This algorithm uses a forward mode for integration of the dynamic system until the upper contraint is violated.
After violation the closest time point for a treatment is searched via backtracking (backward mode)


Inputs:     
Tf:                     End of time interval [0, Tf]
N:                      Absolute number of integrator steps
x0:                     Initial value of dynamic system
integrator_function:    Casadi integrator function defined in pv_schedule routine by model_integrator
p_in:                   Parameter vector (beta, gamma) of subject / patient
max_fraction:           Maximal allowed fractional blood removal 
Base:                   steady state value of x3 
allowed_arr:            list indicating whether treatment is allowed at time point N

Outputs: 
x1_opt:                 Optimal trajectory (based on algorithm) of x1
x2_opt:                 Optimal trajectory (based on algorithm) of x2
x3_opt:                 Optimal trajectory (based on algorithm) of x3
q_opt:                  Trajectory of objective value
u:                      Optimal control with control value for each valid grid point
tgrid:                  Time grid of trajectories for plotting
error_flag:             Indication whether the heuristic algorithm found a valid solution

"""


def pv_heuristic_alg(Tf, N, x0, integrator_function, p_in, max_fraction, Base, allowed_arr):
    
    # time grid
    tgrid = [Tf / N * k for k in range(N + 1)]

    # Start values
    x_opt = [x0]
    q_opt = [0]
    u = [0]*N
    
    # Integration starts in forward mode
    back_flag = False
    k = 0

    count = 0
    
    # Start integration
    while k < N:
        count += 1
        # FORWARD MODE
        if not back_flag:        
            i_out = integrator_function(x0=x_opt[-1], q0=q_opt[-1], u=u[k], p=p_in)   # forward integration
            i_out['xf'][5] = i_out['xf'][5] * (1 - u[k]*max_fraction)                 # apply treatment if scheduled
            
            
            # check if constraint in x3 is violated: start back search if yes
            if i_out['xf'][5] > 1.1*Base:
                back_flag = True
                first_step_back = True
                q_opt_tmp = np.copy(q_opt)
                x_opt_tmp = np.copy(x_opt)
                # print('Backtracking started at time ', tgrid[k])
   
            else: # Append solution of integration to solution vector
                
                x_opt.append(i_out['xf'][3:6])
                q_opt.append(i_out['li'][1])
                k = k + 1
        else:
            # BACKWARD MODE
            # Exceptional case: No valid treatment time is found (e.g. because of too sparse grid)
            # Count >= number ensures that this happens at the beginning of the integration
            if count >=10 and k<=0:
                print('Allowed grid is too sparse. Backward mode does not find a valid solution anymore')
                # Store and return solution up to this point
                x1_opt, x2_opt, x3_opt = state_separator(x_opt_tmp)
                
                # extend solutions to correct length by zero entries
                x1_opt = np.concatenate([x1_opt, np.zeros(len(tgrid) - len(x1_opt))])
                x2_opt = np.concatenate([x2_opt, np.zeros(len(tgrid) - len(x2_opt))])
                x3_opt = np.concatenate([x3_opt, np.zeros(len(tgrid) - len(x3_opt))])
                q_opt = np.concatenate([q_opt_tmp, np.zeros(len(tgrid) - len(q_opt_tmp))])
                
                u = [u[i] for i in range(0, len(allowed_arr)) if allowed_arr[i]>=1e-8]
                error_flag = 1  # error occured
                return x1_opt, x2_opt, x3_opt, q_opt, u, tgrid, error_flag
            
            # first step back is actually no real step back, therefore it should not be deleted
            if (k != 0 and not first_step_back):
                del x_opt[-1]
                del q_opt[-1]
            elif first_step_back:
                first_step_back = False

            if allowed_arr[k]>=1e-8 and u[k] == 0:     
                # if allowed time point and no treatment was applied, apply treatment
                i_out = integrator_function(x0=x_opt[-1], q0=q_opt[-1], u=1, p=p_in)
                i_out['xf'][5] = i_out['xf'][5] * (1 - max_fraction)
                #print('Applied treatment at time', tgrid[k])
                if i_out['xf'][5] > 0.8*Base:    # if lower constraint is not violated, donate here
                    back_flag = False
                    u[k] = 1
                    x_opt.append(i_out['xf'][3:6])
                    q_opt.append(i_out['li'][1])
                    k = k + 1
                else:
                    # print('Lower constraint is violated. Going further back.')
                    k = k - 1
            else:
                k = k - 1    # if time point not allowed, go further back to find last valid time point
    
    # split solution returned by integrator into a format suitable for evaluation
    x1_opt, x2_opt, x3_opt = state_separator(x_opt)
        
    # Plot routine needs control only on valid time points
    u = [u[i] for i in range(0, len(allowed_arr)) if allowed_arr[i]>=1e-8]
    
    error_flag = 0
    return x1_opt, x2_opt, x3_opt, q_opt, u, tgrid, error_flag




