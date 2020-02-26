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

from Modules.Integrator.integrator_rk4 import rk4_integrator
from Modules.Integrator.integrator_rk4_scaled import rk4_scaled_integrator

import casadi as ca

"""
Generation of casadi integrator function for use with the dynamic pv model. 
Also includes a second integrator function on [0, 1] for the two_stage extention.

Inputs:
N:              Number of integration points
dt:             Step size on the integration grid
Tf:             End time of integration
Nperday:        Number of integration points per day
B:              Steady state value of x3
max_fraction:   Maximal fractional blood removal by treatment
pv_lambda:      Patient parameter pv_lambda
two_stage:      Option for generation of a second integrator function for the two_stage extension


Output:
integrator_function:    Integrator function for the dynamic model
Additional Output for two_stage == True:
integrator_function_2:  Additional integrator function for the two_stage process

"""
def model_integrator(N, dt, Tf, Nperday, B, max_fraction, pv_lambda, two_stage=False):
    k1 = 1./8
    k2 = 1./6  
    alpha = 1./120    
    p = ca.SX.sym('p', 2)  
    x = ca.SX.sym('x', 3)
    u = ca.SX.sym('u', 1)

    # Dynamic model
    X0_const = alpha * B
    gamma_pv = p[0] * 0.1
    
    # New stepsize for two_stage process
    dt_two_stage = 1./20
    
    # Model equations
    ode_rhs = ca.vertcat(p[0] * (X0_const - k1 * x[0]) +
               p[1] * (1 - pv_lambda) * (1 - x[2]/B) * x[0] +
               pv_lambda * gamma_pv * x[0],
               p[0] * (k1 * x[0] - k2 * x[1]),
               p[0] * (k2 * x[1] - alpha * x[2]))
    
    # Objective observed by integrator
    objective = u
    
    # Casadi function for integration
    f = ca.Function('f', [x, u, p], [ode_rhs, objective])
    integrator_function = rk4_integrator(f, dt, 1, 2, 3)   
    
    if two_stage:
        # second function for two_stage extension
        integrator_function_2 = rk4_scaled_integrator(f, dt_two_stage, 1, 2, 3)
        
    if two_stage:
        return integrator_function, integrator_function_2
    else:
        return integrator_function






