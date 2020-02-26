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


import casadi as ca

"""
Generates a Casadi Function for integration of the dynamic model using the Runge-Kutta method of order 4


Inputs: 
    f:              Casadi function for the ode rhs; function inputs are (X, U, P ...),
                    where X is state vector, U is control variable and P is parameter vector
    dt:             Length of the control interval
    M:              Number of integrator steps per control interval
    n_params:       Number of input parameters
    n_states:       Number of state variables
    
Outputs: 
    RK_func:        Integrator function for stepwise integration or stepwise nlp formulation
                    This function returns all M steps as casadi variables            
"""



# RK4 integrator
def rk4_integrator(f, dt, M, n_params, n_states):

    dt_m = dt/M     # integration step size
    
    # define symbolic variables
    U = ca.MX.sym('U')             # 1D symb variable for control
    P = ca.MX.sym('P', n_params)   # symb variable for input parameters
    
    X0 = ca.MX.sym('X0', n_states) # symb variable for initial value of integration  
    X = X0                         # set temporary X to initial value

    Q0 = ca.MX.sym('Q0', 1)        # symb variable for the integral of the control function
    Q = Q0                         # set temporary Q to initial value
    
    X_sol = X0                     # store solutions as concatinated variable, start with initial values
    Q_sol = Q
    
    # RK4 integration scheme
    for k in range(M):
        k1, l1 = f(X, U, P)
        k2, l2 = f(X + dt_m/2 * k1, U, P)
        k3, l3 = f(X + dt_m/2 * k2, U, P)
        k4, l4 = f(X + dt_m * k3, U, P)
        
        # update solution
        X = X +  dt_m / 6 * (k1 + 2*k2 + 2*k3 + k4)
        Q = Q +  dt_m / 6 * (l1 + 2*l2 + 2*l3 + l4)
        
        # store solution step 
        X_sol= ca.vertcat(X_sol, X)
        Q_sol = ca.vertcat(Q_sol, Q)
    
    # return scheme as casadi function
    RK_func = ca.Function('RK4',  [X0, Q0, U, P], [X_sol, Q_sol], ['x0', 'q0', 'u', 'p'], ['xf', 'li'])
    
    return RK_func

