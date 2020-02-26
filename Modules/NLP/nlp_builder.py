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
import casadi as ca


"""
Generation of an NLP for the dynamic pv system in casadi with multiple shooting structure

Inputs:
N:                      Total number of integration points
Nperday:                Number of integration points per day
dt:                     Integration step size
B:                      Steady state value of x3
x0:                     Steady state value of x
max_fraction:           Maximal fractional blood removal per treatment
allowed_arr:            List of allowed treatment integration points
integrator_function:    Casadi integrator function for NLP
integrator_function_2:  Second integrator function for integer end point method, use 'None' if other method is used
p_in:                   Patient parameter [beta, gamma]
u_start:                Initialization of control values
u_max:                  Maximal number of allowed donations, only used for integer end point method

Outputs:
Q:                      Cumulated objective function
w:                      Optimization variables including x and u at certain time points
w0:                     Initialization for w
g:                      Constraints g
lbw:                    Lower bound on w
ubw:                    Upper bound on w
lbg:                    Lower bound of g
ubg:                    Upper bound of g
discrete:               Specifies, whether parts of w are discrete by definition. Here always 'False', needed for NLP formulation in casadi
"""



def nlp_builder(N, Nperday, dt, B, x0, max_fraction, allowed_arr,
                integrator_function, integrator_function_2, p_in, u_start, u_max):
    
    # Bounds in u and x
    lbu = 0
    ubu = 1
    lbx = [0, 0, 0]
    ubx = [np.Infinity]*3
    # ubx = [1000, 1000, 2000 ]    
   
    # get a feasible trajectory as initial guess
    # Incorporate initial value x0
    xk = ca.DM(x0)
    x_start = [xk]
    
    # Integration for rest of x trajectory from given u_start
    u_idx = 0
    for k in range(N):
        if allowed_arr[k]>=1e-8:
            xk = integrator_function(x0=x_start[-1], q0 = 0, u = u_start[u_idx], p = p_in)
            u_idx+=1
        else:
            xk = integrator_function(x0=x_start[-1], q0 = 0, u = 0, p = p_in)
        x_start +=[xk['xf'][3:6]]
        
    # Empty NLP
    w = []
    w0 = []
    lbw = []
    ubw = []
    discrete = []
    
    # Decouple parts given by Uk first from w
    # Needed for evaluation routine
    w_u = []
    w0_u = []
    lbw_u = []
    ubw_u = []
    discrete_u = []
    
    # Constraints
    g = []
    lbg = []
    ubg = []
    
    # "Lift" initial conditions
    X0 = ca.MX.sym('X0', 3)
    w += [X0]
    lbw += x0
    ubw += x0
    w0 += [x0]
    discrete += [False]*3
    
    # variable for objective
    Q = 0
    Xk = X0
    
    # allowed control times
    u_idx = 0

    # variable for number of controls
    u_sum = 0
    
    # Integration
    for k in range(N):
        if allowed_arr[k]>1e-8:
            # New NLP variable for the control
            Uk = ca.MX.sym('U_' + str(k))
            w_u   += [Uk]
            lbw_u += [lbu]
            ubw_u += [ubu]
            w0_u  += [u_start[u_idx]]
            discrete_u += [True]
            u_sum+=Uk
            
            u_idx +=1
                        
            F_output = integrator_function(x0=Xk, q0=Q , u=Uk , p=p_in)
        else:
            F_output = integrator_function(x0=Xk, q0=Q , u=0 , p=p_in)

        Xk_end = F_output['xf'][3:6]
        Q = F_output['li'][1]
        
        # Multiple shooting
        Xk = ca.MX.sym('X_' + str(k+1), 3)
        w   += [Xk]
        lbw += lbx
        ubw += ubx
        w0  += [x_start[k+1]]
        discrete += [False, False, False]
        
        # include jump if control > 0
        if allowed_arr[k]>1e-8:  
            g   += [Xk_end[0]-Xk[0], Xk_end[1]-Xk[1], Xk_end[2]*(1 - Uk * max_fraction)-Xk[2]]
            lbg += [0, 0, 0]
            ubg += [0, 0, 0]        
        else:
            g   += [Xk_end-Xk]
            lbg += [0, 0, 0]
            ubg += [0, 0, 0]  
        
        # Strict constraints on x3, enforce only after a few integration steps, as it is possible that a patient starts in region with too high thb
        if k > 3:        
            g += [Xk[2]]
            lbg += [0.8*B]
            ubg += [1.1*B]
    
    
    # NLP extension for the second stage of the integer end point process if applied
    if integrator_function_2 != None:
        # optimization variable for variable end time
        Tf_inv = ca.MX.sym('Tf_inv', 1)      
    
        # Integration on [0, Tf2] with 20 steps
        for k in range(20): 
            
            # Integration
            F_output = integrator_function_2(x0=Xk, q0=Q , u=0 , p=p_in, scale=1./Tf_inv)
            Xk_end = F_output['xf'][3:6]
            Q = F_output['li'][1]
            
            
            # Multiple shooting
            Xk = ca.MX.sym('X_' + str(N+k+1), 3)
            w   += [Xk]
            lbw += lbx
            ubw += ubx
            w0  += [x_start[k+1]]
            discrete += [False, False, False]
            
            # Contraint for connection of multiple shooting intervals
            g   += [Xk_end-Xk]
            lbg += [0, 0, 0]
            ubg += [0, 0, 0]  
            
            # Constraints on x3
            g += [Xk[2]]
            ubg += [1.1*B]
            if k == 20-1:
                lbg += [1.1*B] # endpoint condition x_3(t=1) = x_up
            else:
                lbg += [0.8*B]

        # include Tf_inv into optimization problem
        w += [Tf_inv]
        lbw +=[0.01]      # heuristic: treatment should not take more than 100 days
        ubw +=[1e8]       
        w0 += [0.2]
        Q = -1./Tf_inv
        discrete += [False]
    
    # Include number of treatments into objective for end point approach
    if u_max != None: 
        g += [u_sum]
        ubg += [u_max]
        lbg += [0]

    # concatenate w_u with w
    w += w_u
    w0 += w0_u
    lbw += lbw_u
    ubw += ubw_u
    discrete += discrete_u

    # Concatenate decision variables and constraint terms
    w = ca.vertcat(*w)
    g = ca.vertcat(*g)
    
    return Q, w, w0, g, lbw, ubw, lbg, ubg, discrete

