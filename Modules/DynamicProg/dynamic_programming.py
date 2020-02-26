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
# path to pv_schedule
sys.path.append('../../')

from Modules.NLP.integrate_nlp_sol import integrate_nlp_sol
from Modules.Model.model_integrator import model_integrator
from Modules.Tools.patient_parameters import *

from pylab import * 

import sparse
import time
import resource
"""
Dynamic programming approach for generation of schedules for PV patients. 
!!! Caution: very high storage demand, especially with highly increasing NX !!!

Inputs:
Tf:             End Time of optimization
Nperday:        Number of integration points per day, only used for integration 
x0:             Initial value of dynamical variables
Base:           Steady state value of x3
pv_lambda:      pv_lambda value of patient 
patient_volume: patient blood volume
allowed_opts:   dictionary options for allowed days in format similar to pv_schedule
max_fraction:   maximal fractional blood loss
x3_up:          upper bound for x3
dp_options:     dictionary with options for dynamic programming including:
        NK:         number of RK4 steps per interval
        NU:         number of discrete control values
        NX:         number of state values (highly increases precision, but also storage demand!)
        p_shift:    shifting factor in rounding rule

Outputs:
    x1_opt:     Trajectory for computed x1
    x2_opt:     Trajectory for computed x2
    x3_opt:     Trajectory for computed x3
    q_opt:      Trajectory for objective function
    u_dp:       Optimal control
    tgrid:      Time grid for plotting of trajectories

"""


def dynamic_prog_pv(Tf, Nperday, x0, Base, pv_lambda, p_in, patient_volume,
                    allowed_opts, max_fraction, x3_up,  dp_options):

    # start of time measurement
    t = time.clock()

    # Number of control intervals = number of days here
    N = int(Tf)

    # necessary transformation of allowed configuration for dp algorithm
    allowed_days = allowed_opts['allowed_days']
    allowed_idx = [i for i in range(7) if allowed_days[i]==1]
    t_blocked = [j for j in range(N) if int(j/1) % 7 not in allowed_idx]
    
    if 'forbidden_days' in allowed_opts.keys():
        forbidden_days = allowed_opts['forbidden_days']
        for f in forbidden_days:
            t_blocked += list(f)

    # reading options from dp_options
    # Number of Runge-Kutta 4 steps per interval and step size
    NK = dp_options['NK']  
    DT = Tf/(N*NK)

    # Number of discrete control values
    NU = dp_options['NU'] 

    # Number of discrete state values
    NX = dp_options['NX'] 

    # System dynamics, can be called with matricex
    def f(x1, x2, x3, u):
        k1 = 1./8
        k2 = 1./6
        alpha = 1./120
        gamma_pv = 0.1*p_in[0]
        X0_const = alpha * Base
        x1_dot = p_in[0] * (X0_const - k1*x1) + p_in[1] * (1 - pv_lambda) * (1 - x3/Base) * x1 + pv_lambda * gamma_pv * x1
        x2_dot = p_in[0] * (k1 * x1 - k2 * x2)
        x3_dot = p_in[0] * (k2 * x2 - alpha * x3)
        q_dot = u
        return (x1_dot, x2_dot, x3_dot, q_dot)

    # Control enumeration
    U  = linspace(0,1,NU)

    # State space enumeration
    # This is a hard coded reasonable domain for x
    # Can be reduced for higher precision without increased storage demand if 
    # tighter bounds are known
    x1 = linspace(50, 170, NX)
    x2 = linspace(35, 150, NX)
    x3 = linspace(0.8*Base, 1.1*Base, NX)
    X1, X2, X3 = meshgrid(x1, x2, x3, indexing='ij')

    # For each control action and state, precalculate next state and stage cost
    stage_J = []
    next_x1 = []
    next_x2 = []
    next_x3 = []
    for u in U:
      # Take number of integration steps
      X1_k = copy(X1)
      X2_k = copy(X2)
      X3_k = copy(X3)
      Q_k = zeros(X1.shape)
      for k in range(NK):
        # RK4 integration for x1, x2 and q
        k1_x1, k1_x2, k1_x3, k1_q = f(X1_k, X2_k, X3_k, u)
        k2_x1, k2_x2, k2_x3, k2_q = f(X1_k + DT/2 * k1_x1, X2_k + DT/2 * k1_x2, X3_k + DT/2 * k1_x3, u)
        k3_x1, k3_x2, k3_x3, k3_q = f(X1_k + DT/2 * k2_x1, X2_k + DT/2 * k2_x2, X3_k + DT/2 * k2_x3, u)
        k4_x1, k4_x2, k4_x3, k4_q = f(X1_k + DT * k3_x1, X2_k + DT * k3_x2, X3_k + DT * k3_x3, u)
        X1_k += DT/6*(k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1)
        X2_k += DT/6*(k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2)
        X3_k += DT/6*(k1_x3 + 2*k2_x3 + 2*k3_x3 + k4_x3)
        Q_k += DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q)

      if u == 1:
        X3_k *= (1-500./patient_volume)

      # Save the stage cost and next state
      next_x1.append(X1_k)
      next_x2.append(X2_k)
      next_x3.append(X3_k)
      stage_J.append(Q_k)

    X1_before_shift = copy(next_x1)
    X2_before_shift = copy(next_x2)
    X3_before_shift = copy(next_x3)
    Q_before_shift = copy(stage_J)
    
    # Use shift in rounding rule if stated
    if 'p_shift' in dp_options.keys():
        p_shift = dp_options['p_shift']
    else:
        # default shift is 0
        p_shift = 0.0

#    for p_shift in P:
    stage_J = []
    next_x1 = []
    next_x2 = []
    next_x3 = []
    X1_k = copy(X1)
    X2_k = copy(X2)
    X3_k = copy(X3)
    Q_k = zeros(X1.shape)
    for u in {0, 1}:
        # Find out which state comes next (index)
        X1_k = matrix.round(p_shift + (X1_before_shift[u] - x1[0]) / (x1[-1] - x1[0]) * (NX - 1)).astype(int)  
        X2_k = matrix.round(p_shift + (X2_before_shift[u] - x2[0]) / (x2[-1] - x2[0]) * (NX - 1)).astype(int)
        X3_k = matrix.round(p_shift + (X3_before_shift[u] - x3[0]) / (x3[-1] - x3[0]) * (NX - 1)).astype(int)

        Q_k = Q_before_shift[u]
        # Infinite cost if state gets out-of-bounds
        I = X1_k < 0;
        Q_k[I] = inf;
        X1_k[I] = 0
        I = X2_k < 0;
        Q_k[I] = inf;
        X2_k[I] = 0
        I = X3_k < 0;
        Q_k[I] = inf;
        X3_k[I] = 0
        I = X1_k >= NX;
        Q_k[I] = inf;
        X1_k[I] = 0
        I = X2_k >= NX;
        Q_k[I] = inf;
        X2_k[I] = 0
        I = X3_k >= NX;
        Q_k[I] = inf;
        X3_k[I] = 0

        # Save the stage cost and next state
        next_x1.append(X1_k)
        next_x2.append(X2_k)
        next_x3.append(X3_k)
        stage_J.append(Q_k)

    print("Initial table ready. This took ", time.clock() - t, " seconds.")
    print("Memory: ", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000000, " GB")
    t = time.clock()
    for k in range(N):
      print('-',end="")
    print(' ')
    
    # Calculate cost-to-go (no end cost) and optimal control
    J = zeros(X1.shape)
    U_opt = []
    for k in reversed(list(range(N))):
        # Cost to go for the previous step, optimal control action
        print('+',end="", flush=True)
        J_prev = inf*ones(X1.shape)
        u_prev = zeros(X1.shape, dtype=bool)
        # if time in t_blocked fix control to 0, else allow 1
        if k in t_blocked:
            U_feas = [0]
        else:
            U_feas = [0, 1]
            
        # Test all control actions
        for uind in U_feas:
            J_prev_test = J[next_x1[uind], next_x2[uind], next_x3[uind]]+stage_J[uind]*1 # u_weights[k]
            better = J_prev_test<J_prev
            u_prev[better] = uind
            J_prev[better] = J_prev_test[better]
             
        # Update cost-to-go and save optimal control
        J = J_prev
        u_sparse = sparse.COO(u_prev)
        U_opt.append(u_sparse)
        
    # Reorder U_opt by stage
    U_opt.reverse()
    print(" ")
    print("Computation of optimal control took ", time.clock() - t, " seconds.")
    # Find optimal control
    u_opt = []
    x1_opt = [x0[0]] 
    x2_opt = [x0[1]] 
    x3_opt = [x0[2]] 
    i1 = int(round((x1_opt[0] - x1[0]) / (x1[-1] - x1[0]) * (NX - 1))) 
    i2 = int(round((x2_opt[0] - x2[0]) / (x2[-1] - x2[0]) * (NX - 1)))
    i3 = int(round((x3_opt[0] - x3[0]) / (x3[-1] - x3[0]) * (NX - 1)))
    cost = 0
    for k in range(N):
        # Get the optimal control and go to next step
        u_ind = U_opt[k][i1, i2, i3]
        cost += stage_J[u_ind][ i1, i2, i3]
        i1, i2, i3 = next_x1[u_ind][i1, i2, i3], next_x2[u_ind][i1, i2, i3], next_x3[u_ind][i1, i2, i3]

        # Save the trajectories
        if u_ind:
            u_opt.append(1)
        else:
            u_opt.append(0)
        x1_opt.append(x1[i1])
        x2_opt.append(x2[i2])
        x3_opt.append(x3[i3])

    # Optimal cost
    print("\n Minimal cost: ", cost)

    
    memory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1000000
    print("Memory: ", memory, " GB")


    # integration of solution for external plotting
    u_dp = zeros(Tf*Nperday)
    for i in range(len(u_opt)):
        u_dp[int(i * Nperday)] = u_opt[i]
    dt = 1./Nperday
    N_int = int(Tf * Nperday)

    integrator_function = model_integrator(N_int, dt, Tf, Nperday, Base, max_fraction, pv_lambda)

    x1_opt, x2_opt, x3_opt, q_opt, u_dp, tgrid = integrate_nlp_sol(u_dp, x0, [1] * N_int, N_int, dt, Nperday,
                                                                    Tf, integrator_function, None, max_fraction, p_in,
                                                                    sol_is_u=True)

    return x1_opt, x2_opt, x3_opt, q_opt, u_dp, tgrid