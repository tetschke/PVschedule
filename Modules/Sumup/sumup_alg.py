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

"""
Conversion of relaxed control to integer control using the sum up method
Inputs:
u_opt:      Optimal control of the relaxed problem

Output:
u_sumup:    Optimal control using the sum up method

"""
def sumup(u_opt):
    u_len = len(u_opt)
    u_sumup = np.zeros(u_len)
    # minimal value of u_opt which is recognized as non-zero, prevents incorporation of numerical artifacts
    epsilon = 1e-5  
    w = 0
    for k in range(u_len):
        if u_opt[k] > epsilon:
            w += 1.5 * u_opt[k]
        if w > 0:
            u_sumup[k] = 1
            w -= 1
    return u_sumup



