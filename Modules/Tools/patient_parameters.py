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



This file contains all 28 patients from the study by Pottgiesser et al. 2008, 
which was cited and used for generation of the results in Tetschke et al. 2018.
Also the resulting parameters of the three subsets of clinical data are given.

Subjects from the Pottgiesser study have indices 'F01', 'F02', ... , 'F29', where number 22 was excluded.
Patients datasets are referenced by indices 'P01', 'P02' and 'P03'

"""

# =============================================================================
"""
Calculation of patients total blood volume based on tHb difference before and
after removal of 500 ml blood
Inputs: 
thb_before:     tHb before blood removal
thb_after:      tHb after blood removal

Output:
blood_volume: total blood volume of patients

"""
def calc_blood_volume(thb_before, thb_after):
    thb_loss = thb_before - thb_after
    blood_volume = thb_before / thb_loss * 500
    return blood_volume

# =============================================================================
"""
Calculation of steady state value
Input:
B: Steady state value of x3
Output:
Vector of steady state values of x1, x2 and x3
"""

# =============================================================================
def calc_x0(B):
    return [B*(8/120), B*(6/120), B]


"""
Function which returns parameters for pv_schedule for given subject or patient 
based on patient index and artificially generated lambda, if applicable

Inputs: 
index:          Patient index as string 'F..' or 'P..' as given above
lambda_version: 'F' subjects are generated with 5 different pv_lambdas.
                 This is number 0 to 4 to specify one of those lambdas
                 
Outputs:
gamma:          Patient parameter gamma
beta:           Patient parameter beta
B:              Steady state value B of x3
bvol:           Estimated total blood volume of patient
x0:             Steady state value of x
pv_lambda:      Patient parameter pv_lambda
"""
# =============================================================================
def return_parameters(index, lambda_version=None):

    lambda_dict = {'F01': [0.512578730891727, 0.40498120621217, 0.418333992907552, 0.521260759828438, 0.512377581556603],
    'F02': [0.385240151041006, 0.55848996431784, 0.706279751615796, 0.708757007422721, 0.497865051237947],
    'F03': [0.549030304326684, 0.412796531167805, 0.479513379466922, 0.393353003493433, 0.326410638146832],
    'F04': [0.877678463480118, 0.603778991724502, 0.522472234037192, 0.888048716840321, 0.544069470266772],
    'F05': [0.321030679108264, 0.334118176201712, 0.191811852796399, 0.36691717633693, 0.418743123463179],
    'F06': [0.290027700482221, 0.343217974988489, 0.414918246241808, 0.354305674583423, 0.330829164959082],
    'F07': [0.423571527754871, 0.342717356797054, 0.479772413084511, 0.370252222974169, 0.349496893333417],
    'F08': [0.370570797064878, 0.339705322872325, 0.432992666555345, 0.215528648117589, 0.243275263179487],
    'F09': [0.604957230986552, 0.601683855811865, 0.555161220440021, 0.558678810821486, 0.366282719729957],
    'F10': [0.199171591549389, 0.396397654599178, 0.391467264060761, 0.207055695418981, 0.298446163554796],
    'F11': [0.385133730979017, 0.383873233825338, 0.532897656163912, 0.615139582523103, 0.651585525659608],
    'F12': [0.581086555604184, 0.528473402749661, 0.540769551234917, 0.631293612463677, 0.566534741402981],
    'F13': [0.583411004244528, 0.427785566320457, 0.50567483769101, 0.197864903321417, 0.368726673404813],
    'F14': [0.767107067454641, 0.639191512186174, 0.845035242247026, 0.787492363802794, 0.74286982574285],
    'F15': [0.333668696101623, 0.387229958295863, 0.407832227597826, 0.288989142652047, 0.396683879409587],
    'F16': [0.339124515667045, 0.379286476373994, 0.226292279605389, 0.350138465828293, 0.348863131548444],
    'F17': [0.541173939123054, 0.690811600374238, 0.705248966034992, 0.757539680322406, 0.513658288347558],
    'F18': [0.661177818695674, 0.689214635813903, 0.695490671649802, 0.529427828044276, 0.847362875119706],
    'F19': [0.393191552667039, 0.539861381797223, 0.400533253386759, 0.450952525912197, 0.511916912869378],
    'F20': [0.334419971548197, 0.327756849035277, 0.341559395898565, 0.341115792565011, 0.359533383673401],
    'F21': [0.389189764747967, 0.463349563574065, 0.404331648221089, 0.408079122468663, 0.344943797935412],
    'F23': [0.336531172668368, 0.302567508871834, 0.304972455755748, 0.277803767930832, 0.341864686366354],
    'F24': [0.317578712817824, 0.502451211579867, 0.517682024980399, 0.563407276272935, 0.532926268243279],
    'F25': [0.434938721639482, 0.453896502273388, 0.478822883900443, 0.475571052313771, 0.469234893520139],
    'F26': [0.624878960721628, 0.575124838050053, 0.599810917230367, 0.559344088111987, 0.407811132788059],
    'F27': [0.680557244471105, 0.745904282006242, 0.681821903608136, 0.659197796918977, 0.708113081323333],
    'F28': [0.440131308664239, 0.466034475240221, 0.443867051476406, 0.48971991233144, 0.55512123229294],
    'F29': [0.471857908301574, 0.496790602957417, 0.506635602870428, 0.52671552906319, 0.551452253402328]}

    x0 = None
    if lambda_version != None:
        pv_lambda = lambda_dict[index][lambda_version]
    else:
        # print('No lambda version was given, using default of pv_lambda = 0.5')
        pv_lambda = 0.5
    if index == 'F01':
        gamma = 0.769
        beta = 1.65
        B = 865.4478782713
        bvol = calc_blood_volume(B, 787.1981165782)
    elif index == 'F02':
        gamma = 0.388
        beta  = 0.867
        B = (887.9670710908 + 882.8658195301) / 2
        bvol = calc_blood_volume(B, 790.5385640564)
    elif index == 'F03':
        gamma = 0.51
        beta  = 1.617
        B = 863.9664324528
        bvol = calc_blood_volume(B, 781.9327847027)
    elif index == 'F04':
        gamma = 0.323
        beta  = 0.424
        B = 854.1470644573
        bvol = calc_blood_volume(B, 782.7861500105)
    elif index == 'F05':
        gamma = 0.061
        beta  = 1.381
        B = 971.6744459046
        bvol = calc_blood_volume(B, 908.8573962829)
    elif index == 'F06':
        gamma = 0.59
        beta  = 2.615
        B = 1001.4244181644
        bvol = calc_blood_volume(B, 903.1810500387)
    elif index == 'F07':
        gamma = 0.262
        beta  = 1.518
        B = 964.5942160412
        bvol = calc_blood_volume(B, 898.254997679)
    elif index == 'F08':
        gamma = 0.324
        beta  = 2.676
        B = 704.4165719154
        bvol = calc_blood_volume(B, 618.3272187083)
    elif index == 'F09':
        gamma = 0.356
        beta  = 0.891
        B = 958.5515600646
        bvol = calc_blood_volume(B, 906.9209507009)
    elif index == 'F10':
        gamma = 0.089
        beta  = 2.557
        B = 851.702312115
        bvol = calc_blood_volume(B, 758.8964352941)    
    elif index == 'F11':
        gamma = 0.243
        beta  = 0.925
        B = 1006.4477644294
        bvol = calc_blood_volume(B, 897.2949939085)        
    elif index == 'F12':
        gamma = 1.003
        beta  = 1.409
        B = 932.5110925716
        bvol = calc_blood_volume(B, 856.4186547737)     
    elif index == 'F13':
        gamma = 0.057
        beta  = 0.879
        B = 647.976181699
        bvol = calc_blood_volume(B, 567.3358695386)  
    elif index == 'F14':
        gamma = 0.762
        beta  = 0.46
        B = (1073.0035700144 + 1089.6781694937) / 2
        bvol = calc_blood_volume(B, 1015.8921337506)
    elif index == 'F15':
        gamma = 0.344
        beta  = 2.132
        B = 939.6110449136
        bvol = calc_blood_volume(B, 870.3018643068) 
    elif index == 'F16':
        gamma = 0.141
        beta  = 1.661
        B = (750.4670973864 + 756.0191012968) / 2
        bvol = calc_blood_volume(B, 700.2177604608)        
    elif index == 'F17':
        gamma = 0.47
        beta  = 0.544
        B = 900.5306496112
        bvol = calc_blood_volume(B, 823.3312912543)         
    elif index == 'F18':
        gamma = 0.525
        beta  = 0.631
        B = (850.5661000629 + 832.6554649928) / 2
        bvol = calc_blood_volume(B, 755.2416887191)         
    elif index == 'F19':
        gamma = 0.423
        beta  = 1.525
        B = (778.6938305164 + 794.2542354547) / 2
        bvol = calc_blood_volume(B, 709.5150208333)           
    elif index == 'F20':
        gamma = 0.661
        beta  = 2.798
        B = 765.994349054
        bvol = calc_blood_volume(B, 720.8625844018)           
    elif index == 'F21':
        gamma = 0.686
        beta  = 1.943
        B = (895.623697123 + 921.5732097495) / 2
        bvol = calc_blood_volume(B, 829.2582912455)         
    elif index == 'F23':
        gamma = 0.613
        beta  = 3.142
        B = 893.0599199985
        bvol = calc_blood_volume(B, 810.9540053454)          
    elif index == 'F24':
        gamma = 0.421
        beta  = 1.528
        B = 695.0543535551
        bvol = calc_blood_volume(B, 625.3963209704)            
    elif index == 'F25':
        gamma = 0.863
        beta  = 2.078
        B = (757.8750037449 + 779.78560612) / 2
        bvol = calc_blood_volume(B, 706.6521049009)         
    elif index == 'F26':
        gamma = 0.414
        beta  = 1.172
        B = 687.849173192
        bvol = calc_blood_volume(B, 627.8652949085)          
    elif index == 'F27':
        gamma = 0.635
        beta  = 0.836
        B = 925.6179043476
        bvol = calc_blood_volume(B, 850.5903319175)         
    elif index == 'F28':
        gamma = 0.952
        beta  = 1.596
        B = (867.2672679748 + 870.723751137) / 2
        bvol = calc_blood_volume(B, 800.5847782191)         
    elif index == 'F29':
        gamma = 0.805
        beta  = 1.486
        B = (796.6758676854 + 821.6792069444) / 2
        bvol = calc_blood_volume(B, 741.6011011007)         
    elif index == 'P01':
        gamma = 0.294
        beta  = 0.588
        B = 767
        bvol = 5081
        x0   = [54.8403, 72.5159, 881.86]
        pv_lambda = 0.905
    elif index == 'P02':
        gamma = 0.1
        beta  = 0.2
        B = 501.04
        bvol = 4277
        x0   = [37.5485, 40.9281, 570.566]
        pv_lambda = 0.62
    elif index == 'P03':
        gamma = 0.1
        beta  = 0.2
        B = 540.17
        bvol = 4532
        x0   = [80.68, 74.64, 515.432]
        pv_lambda = 0.61
    else:
        print('Index', index, 'was not found! Error!')
    if x0 is None:
        x0 = calc_x0(B)
    return gamma, beta, B, bvol, x0, pv_lambda



