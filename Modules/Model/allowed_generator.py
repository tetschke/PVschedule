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



"""
Generation of an array of allowed times for each integration point using options
'allowed_hours', 'allowed_days', 'forbidden_days' in 'dict_opts'

Inputs: 
Nperday:        Number of steps per day
dict_opts:      Dictionary options given to pv_schedule
N:              Length of integration grid
dt:             Integration stepsize

Output:
allowed_arr:    List with 0-1 entries for each grid point indicating whether a treatment is allowed or not
"""
def allowed_generator(Nperday, dict_opts, N, dt):
    
    if 'allowed_hours' in dict_opts.keys():
        allowed_hours = dict_opts['allowed_hours']
    else:
        allowed_hours = [1]*Nperday
        
    if 'allowed_days' in dict_opts.keys():
        allowed_days = dict_opts['allowed_days']
    else:
        allowed_days = [1]*7
    
    if 'forbidden_days' in dict_opts.keys():
        tmp = dict_opts['forbidden_days']
        forbidden_days = []
        
        for t in tmp:
            if type(t) == int:
                forbidden_days.append(t)
            elif type(t) == range:
                for i in t:
                    forbidden_days.append(i)
    else: 
        forbidden_days = []   
        
    # Weave all dictionary options together 
    hour_idx = 0
    day_idx = 0
    current_time = 0  

    allowed_arr = []
    for k in range(N):
        # evaluate current time step
        allowed = allowed_hours[hour_idx] > 0 and allowed_days[day_idx] > 0 and int(current_time) not in forbidden_days
        if allowed:
            allowed_arr.append(1)
        else:
            allowed_arr.append(0)    
    
        # update indices
        hour_idx +=1
        if hour_idx == Nperday:
            hour_idx = 0
            day_idx +=1
            if day_idx == 7:
                day_idx = 0
        current_time += dt
        
    return allowed_arr

