# -*- coding: utf-8 -*-
"""
Created on Sun Feb 9 14:57:09 2017

@author: jnzzp5
"""

"""
from v1.x, it will handle with multiple layer deposition.
run the code in the IPython Console. run "%matplotlib qt" at first, then run "%run C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer\main_program.py"
main_program_v9_vect_stopgrowth.py is to employed a decentered polygon algorithm considering cell growth stop.(unnecessary to iterate every cell.) it is for fine(5um) cell size.
Integral_coefficient is needed.
finer temperature file.
Growth conditions： 8 neighbors are sold, the center cell stops growing.
including epitaxial growth in Nucleation() function.
including remelt part in main program.
02272017 revise
03022017 revise (add remelt and air part.)
03172017 revise
03222017revise for pickle
04262017revise for deleting "with open('critical_length_index', 'w') as handle1:" and "with open('nucleation_index', 'w') as handle1:" because of error "IOError: [Errno 13] Permission denied: 'critical_length_index'"
05192017 revised
"""
"""
copied from v1.1 on 04282017
This file is built because there are some bad points in the temperature field. So these bad temperature points will be eliminated this time.
In the Nucleation function, some conditions for nucleation are removed.
05012017 revised. Line565 - Line572
Integral_coefficient_s = 100.0, Integral_coefficient_v = 250.0
"""
"""
06132017 created. copied from v1.2 because this file will input Yanlei's temperature file.
09212017 revised from v1.3. All the grey comments are deleted compared to v1.3.
09262017: the problem may exist in for j, four_idx in enumerate(neighvertices_global_idx) part. It feels like that some grains donot stop and always keep enclosing more cells.
10012017 created from v1.4
10012017 created from v1.5 to redesign remelt.
10172017 created from v1.6 to redesign remelt and use git to control. Only one file is needed for different versions. hh
"""

import numpy as np
import scipy.integrate as integrate
from io import StringIO
import os
from Parameters import *
import random
import math
import timeit
from datetime import datetime
from sympy import Point, Polygon
import pandas as pd
import DataProcess as dp
import DisplayData as dd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import path
import Graingrowth as gg
import cPickle as pickle

DTNUCL = np.zeros((RowNum,ColumnNum),dtype=np.int)
SIDX = np.zeros((RowNum,ColumnNum),dtype=np.int)
NUCIDX = np.zeros((RowNum,ColumnNum),dtype=np.int)
LIDXN1 = np.zeros((RowNum,ColumnNum),dtype=np.int) # only represents the initial(first temperature file) temperature
DEACTIVATION = np.zeros((RowNum,ColumnNum),dtype=np.int)
TempInterpolateT = np.zeros((RowNum,ColumnNum),dtype='float16')
STOP = np.zeros((RowNum,ColumnNum),dtype=np.int)
TIPLEN = np.zeros((RowNum,ColumnNum),dtype='float16')
CLASS2 = np.zeros((RowNum,ColumnNum), dtype='int32') # 10102017 added
"""make it invalid to reduce the computational cost! This is temporary."""
TIPLEN1 = np.zeros((RowNum,ColumnNum),dtype='float16')
TIPLEN2 = np.zeros((RowNum,ColumnNum),dtype='float16')
TIPLEN3 = np.zeros((RowNum,ColumnNum),dtype='float16')
TIPLEN4 = np.zeros((RowNum,ColumnNum),dtype='float16')
EnclosePoint = np.full((RowNum,ColumnNum), False, dtype=bool)
LenCritical = np.zeros((RowNum,ColumnNum), dtype='int32')
re_center = np.zeros((RowNum,ColumnNum), dtype='int32')
REMELT = np.zeros((RowNum,ColumnNum), dtype='int32')
AIR = np.zeros((RowNum,ColumnNum), dtype='int32') # 1 represents it is air.
"""********************************"""
TempInterpolate_Old = np.full((RowNum,ColumnNum), TempLiquid, dtype='float16')
NUCINC, NUC, growth_count, nucleat_count = 0, 0, 0, 0
NUCFLG = 0
TempADJ = TempLiquid - TempAdjustor # adjusted TempLiquid
x_grid, y_grid = 10, 10 # control the size of substrate grain in the base.

Cell = dp.Substrate(RowNum, ColumnNum, x_grid, y_grid) # Cell is a class object.
CLASS = Cell.substrate(.2, .25)
ANGLE = (CLASS*90/48.0 - 45) * PI / 180.0
Destination_path = r"C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer\\" # make sure not odd number blackslash end. so \\" is applied, not \"
nuc_dict = dict()
crit_len_idx = dict()

def get_temp_yanlei(_Row_Num, _Column_Num, _TempFileNumber): # revised on 07172017
    temp1 = np.genfromtxt(r"C:/Users/jnzzp5/OneDrive/study/research/CA05122016/CA_2D_multiple_layer/LeiYanTemperature/15layers/Interpolate2d_"+str(_TempFileNumber), delimiter=",", skip_header = 1)
    print "temp1 shape is {0}".format(np.shape(temp1))
    return temp1

def get_temp(_Row_Num, _Column_Num, _TempFileNumber):
    Temp=np.zeros((_Row_Num,_Column_Num),dtype='float64') # this is created to store temperature data
    try:
        c = 4
        if c == 0: # old temperature 2014
            f=file("C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate_finer"+str(_TempFileNumber),'r')
        elif c == 1: # artificial corrected temperature
            f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_correct2014error_02142017\Nodeset-SET-FACE2_1layerTempTi64_Interpolate_finer"+str(_TempFileNumber),'r') 
        elif c == 2:
            temp = np.genfromtxt(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\Nodeset-SET-CS_CLAD9_25layerTempTi64_Interpolate_finer"+str(_TempFileNumber), delimiter=",")
        elif c == 3: # 04282017 added. extract part temperature from original data. 
            temp1 = np.genfromtxt(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\Nodeset-SET-CS_CLAD9_25layerTempTi64_Interpolate_finer"+str(_TempFileNumber), delimiter=",")
            temp = temp1[0:200, 100:500]
        elif c == 4:
            temp = get_temp_yanlei(_Row_Num, _Column_Num,  _TempFileNumber)
    except IOError:
        print "This temperature file cannot be found"

    return temp
    
def update_index(_Row_Num,_Column_Num,_TempADJ):
    Temp = get_temp(_Row_Num,_Column_Num,TempFileNumber_start) # input 400th temperature as the starting temperature file for solidification. update_index only update "SIDX", "NUCIDX", "LIDXN1" of the 1st temperature file.
    for i in xrange(1,_Row_Num-1):
        for j in xrange(1,_Column_Num-1):
            if Temp[i,j] > _TempADJ:
                SIDX[i,j] = 0 #liquid from the beginning
                NUCIDX[i,j] = 0 #liquid at the beginning
            elif Temp[i,j] <= _TempADJ:
                SIDX[i,j] = -1 #solid from beginning
                NUCIDX[i,j] = -1 #solid fron beginning
            if Temp[i,j] > _TempADJ and Temp[i + 1,j] > _TempADJ and Temp[i - 1,j] > _TempADJ and Temp[i,j + 1] > _TempADJ and Temp[i,j - 1] > _TempADJ:
                LIDXN1[i,j] = 1 #LIDXN1 is location index, 1 represents in the volume liquid
            if Temp[i,j] > _TempADJ and (Temp[i + 1,j] < _TempADJ or Temp[i - 1,j] < _TempADJ or Temp[i,j + 1] < _TempADJ or Temp[i,j - 1] < _TempADJ):
                LIDXN1[i,j] = 0 # 0 represents on the surface at the beginning
            if Temp[i,j] <= _TempADJ:
                LIDXN1[i,j] = -1 # represents in the solid at the beginning
                """
                four boundary conditions
                """

    for k in xrange(_Column_Num):
        if Temp[0,k] <= _TempADJ:        
            LIDXN1[0,k] = -1
            NUCIDX[0,k] = -1
            SIDX[0,k] = -1
        if Temp[0,k] > _TempADJ:
            LIDXN1[0,k] = 0
            NUCIDX[0,k] = 0
            SIDX[0,k] = 0
        if Temp[_Row_Num-1,k] <= _TempADJ:
            LIDXN1[_Row_Num-1,k] = -1
            NUCIDX[_Row_Num-1,k] = -1
            SIDX[_Row_Num-1,k] = -1
        if Temp[_Row_Num-1,k] > _TempADJ:
            LIDXN1[_Row_Num-1,k] = 0
            NUCIDX[_Row_Num-1,k] = 0
            SIDX[_Row_Num-1,k] = 0
    for l in xrange(_Row_Num):
        if Temp[l,0] <= _TempADJ:
            LIDXN1[l,0] = -1
            NUCIDX[l,0] = -1
            SIDX[l,0] = -1
        if Temp[l,0] > _TempADJ:
            LIDXN1[l,0] = 0
            NUCIDX[l,0] = 0
            SIDX[l,0] = 0
        if Temp[l,_Column_Num-1] <= _TempADJ:
            LIDXN1[l,_Column_Num-1] = -1
            NUCIDX[l,_Column_Num-1] = -1
            SIDX[l,_Column_Num-1] = -1
        if Temp[l,_Column_Num-1] > _TempADJ:
            LIDXN1[l,_Column_Num-1] = 0    
            NUCIDX[l,_Column_Num-1] = 0
            SIDX[l,_Column_Num-1] = 0
    return Temp
    
"""
initialization should not consider nucleation and growth
"""        


def initialization2(_Row_Num,_Column_Num,_Temp_Liquid,_NUCINC):
    Temp = update_index(RowNum,ColumnNum,TempADJ) # Temp, not TempInterpolateT, stores temperature. Temp is a local variable. 
    solid_idx = np.argwhere(LIDXN1 == -1)
#    no_solid_idx = np.argwhere(LIDXN1 != -1)
    DTNUCL[solid_idx[:, 0], solid_idx[:, 1]] = 0
    bulk_idx = np.argwhere((Temp < _Temp_Liquid) & (LIDXN1 == 1))
    r_arr = np.sin(np.random.rand(_Row_Num,_Column_Num))
    r_arr_abs = np.absolute(r_arr)
    CLASS[bulk_idx[:, 0], bulk_idx[:, 1]] = -1
    PV_arr = Fv_gaussian_ini * gaussian_int(0, abs(_Temp_Liquid - Temp[i,j]), 8.0, .1, 1)[0] / _Row_Num / _Column_Num
    

def initialization(_Row_Num, _Column_Num, _Temp_Liquid, _NUCINC):
    Temp = update_index(RowNum,ColumnNum,TempADJ) # Temp, not TempInterpolateT, stores temperature. Temp is a local variable. 
    for i in xrange(_Row_Num):
        for j in xrange(_Column_Num):
            if LIDXN1[i,j] == -1: #-1 represents in the solid for the first temperature file
                DTNUCL[i,j] = 0 # it is solid originally.
                continue
            if _Temp_Liquid - Temp[i,j] > 0.0:# Temp is not too higher than liquidus temperature
                r = abs(math.sin(random.random()))
                if LIDXN1[i,j] == 1: #1 represents four neighbours and itself are all above TempADJ ORIGINALLY!
                    CLASS[i, j] = -1 # determine the meltpool CLASS
                    PV = Fv_gaussian_ini * gaussian_int(0, abs(_Temp_Liquid - Temp[i,j]), 8.0, .1, 1)[0] / _Row_Num / _Column_Num # nucleation is selected
                    CLASS_block_ini = CLASS[i - 3:i + 3, j - 3:j + 3]
                    if r <= PV and np.all(CLASS_block_ini == -1):
                        DTNUCL[i,j] = 10 # it nucleats in the bulk liquid
                        _NUCINC += 1
                        NUCIDX[i,j] = _NUCINC
                        SIDX[i,j] = 1
                        CLASS[i,j] = random.randint(0,49)
                        ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0
                        print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the initialization bulk."
                        DTSQ = min((_Temp_Liquid - Temp[i,j])**2, DTSQLIM) # TempInterpolateT should not be called here since TempInterpolateT has bot been offered values before Nucleation() function. It is totally different from Nucleation() and Growth(). So Temp[i,j] is referenced here.
                        Vtip = .729*DTSQ + .0103 * DTSQ**2
                        TIPLEN[i, j] = -1*Vtip * DTIME / math.sqrt(2.0)
                    else:
                        DTNUCL[i,j] = 0 #  it is liquid but it doesn't nucleat.
                if LIDXN1[i,j] == 0: # 0 represents at the interface
                    CLASS[i, j] = -1 # determine the meltpool CLASS
                    PS = Fs_gaussian_ini * gaussian_int(0, abs(_Temp_Liquid - Temp[i,j]),.5,.1,1)[0] /_Row_Num / _Column_Num # nucleation is selected
                    if r <= PS:
                        DTNUCL[i,j] = 5 # it nucleats at the S/L surface
                        _NUCINC += 1
                        NUCIDX[i,j] = _NUCINC
                        SIDX[i,j] = 1
                        CLASS[i,j] = random.randint(0,49)
                        ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0                        
                        print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the initialization S/L surface."
                        DTSQ = min((_Temp_Liquid - Temp[i,j])**2, DTSQLIM) #TempInterpolateT should not be called here since TempInterpolateT has bot been offered values before Nucleation() function. It is totally different from Nucleation() and Growth(). So Temp[i,j] is referenced here.
                        Vtip = .729*DTSQ + .0103 * DTSQ**2
                        TIPLEN[i,j] = -1*Vtip * DTIME / math.sqrt(2.0)
                    else:
                        DTNUCL[i,j] = 0 # #  it is liquid but it doesn't nucleat.
            else:
                CLASS[i, j] = -1 # Temperature is higher than  _Temp_Liquid + 10
    return _NUCINC, Temp
    
def Nucleation(_Row_Num,_Column_Num,_Temp_Liquid, _Temp_Solid, _UnderCoolMean, _UnderCoolDelta, _NUCINC, _TempADJ, File_seq):
    start_time = timeit.default_timer()
    global nucleat_count, NUCFLG
    nucleat_count += 1

    """
    determine the nucleation region
    """
    list_m = list()
    list_n = list()
    if not np.any(TempInterpolateT > _TempADJ): #All temperature are below solidus temp.
        return 
    for m in xrange(Invalid_rows, _Row_Num): # replace (_Row_Num) with (Invalid_rows, _Row_Num) because some rows far away from the melt pool won't melt and solidify. this can improve efficiency.
        for n in xrange(_Column_Num):
            if TempInterpolateT[m, n] > _TempADJ: # temperature is higher than solidus temperature.
                list_m.append(m)
                list_n.append(n)
    m_min, m_max, n_min, n_max  = min(list_m), max(list_m), min(list_n), max(list_n) # determine the possible region for nucleation.
    
    liquid_nuc_idx = np.argwhere(TempInterpolateT >= TempADJ)
    CLASS[liquid_nuc_idx[:, 0], liquid_nuc_idx[:, 1]] = -1
    SIDX[liquid_nuc_idx[:, 0], liquid_nuc_idx[:, 1]] = 0 # 05022017 added
    
    """Remelt grain nucleation. 10182017 added"""
#==============================================================================
#     if File_seq > first_layer:
#         haz_idx = np.argwhere(REMELT == 2)
#         for ii in haz_idx:
#             if ii[1] >= _Column_Num - 1 or ii[1] <= 1 or ii[0] >= _Row_Num - 1 or ii[0] <= 1: # eliminate the boundary
#                 continue
#             r_remelt = abs(math.sin(random.random()))
#             if CLASS[ii[0], ii[1] + 1] == -1 and CLASS[ii[0], ii[1] - 1] == -1 and CLASS[ii[0] + 1, ii[1]] == -1 and CLASS[ii[0] - 1, ii[1]] == -1 and CLASS[ii[0] + 1, ii[1] + 1] == -1 and CLASS[ii[0] - 1, ii[1] + 1] == -1 and CLASS[ii[0] + 1, ii[1] - 1] == -1 and CLASS[ii[0] - 1, ii[1] - 1] == -1:
#                 PV_remelt = Fv_gaussian_nuc * gaussian_int((_Temp_Liquid - TempInterpolate_Old[ii[0], ii[1]])/Integral_coefficient_v, (_Temp_Liquid - TempInterpolateT[ii[0], ii[1]])/Integral_coefficient_v, 1.32, .1, 1)[0] / _Row_Num / _Column_Num # it nucleats
#                 CLASS_block = CLASS[ii[0] - 3:ii[0] + 3, ii[1] - 3:ii[1] + 3] # two nucleation center cannot be too close. The distance is at least 6 cell.
#                 if r_remelt <= PV_remelt and np.all(CLASS_block == -1):
#                     DTNUCL[ii[0], ii[1]] = 10 # it nucleats in the bulk liquid.8
#                     NUCFLG = 1 # 1 represents there is at least one nucleation for this time step.
#                     _NUCINC = int(0 if _NUCINC is None else _NUCINC)
#                     _NUCINC += 1 # Record how many 
#                     CLASS[ii[0], ii[1]] = random.randint(0,49)
#     #                CLASS[i, j] = CLASS2[i, j] if  CLASS2[i, j] > 0 else CLASS[i,j] # 10102017 added
#                     SIDX[ii[0], ii[1]] = 1 # become solid because of nucleation
#                     ANGLE[ii[0], ii[1]] = (CLASS[ii[0], ii[1]]*90/48.0 - 45) * PI / 180.0
#                     NUCIDX[ii[0], ii[1]] = _NUCINC # the NUCINCth grain
#                     DTSQ = min((_Temp_Liquid - TempInterpolateT[ii[0], ii[1]])**2, DTSQLIM)
#                     Vtip = .729*DTSQ + .0103 * DTSQ**2
#                     TIPLEN[ii[0], ii[1]] = -1*Vtip * DTIME / math.sqrt(2.0)
#             else:
#                 PS_remelt = Fs_gaussian_nuc * gaussian_int((_Temp_Liquid - TempInterpolate_Old[ii[0], ii[1]])/Integral_coefficient_s, (_Temp_Liquid - TempInterpolateT[ii[0],ii[1]])/Integral_coefficient_s, .5, .1, 1)[0] / _Row_Num / _Column_Num # it nucleats
#                 if r_remelt <= PS_remelt:
#                     DTNUCL[ii[0], ii[1]] = 5 # it nucleats at the S/L surface
#                     NUCFLG = 1 # 1 represents there is at least one nucleation for this time step.
#                     _NUCINC = int(0 if _NUCINC is None else _NUCINC)
#                     _NUCINC += 1 # Record how many nucleations
#                     CLASS[ii[0], ii[1]] = random.randint(0,49)
#     #                    CLASS[i, j] = CLASS2[i, j] if  CLASS2[i, j] > 0 else CLASS[i,j] # 10102017 added
#                     a = [CLASS[ii[0] + 1, ii[1]], CLASS[ii[0] - 1, ii[1]], CLASS[ii[0], ii[1] + 1], CLASS[ii[0], ii[1] - 1]]               
#                     CLASS_interface = a[min(range(len(a)), key=lambda i: abs(a[i] - Class_mean))] # prescribe the CLASS of interface cell when nucleation. This is to look for which neighbouring cell is solid. find i where abs(a[i] - Class_mean) is minimum. min() returns a number within range(len(a)). simulate epitaxial growth.
#                     if CLASS_interface > 0:
#                         CLASS[ii[0], ii[1]] =  CLASS_interface
#     #                    print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the Nucleation surface."
#                     ANGLE[ii[0], ii[1]] = (CLASS[ii[0], ii[1]]*90/48.0 - 45) * PI / 180.0
#                     SIDX[ii[0], ii[1]] = 1 # become solid because of nucleation
#     #                    print "ANGLE = ", ANGLE[i,j]
#                     NUCIDX[ii[0], ii[1]] = _NUCINC # the NUCINCth grain
#                     DTSQ = min((_Temp_Liquid - TempInterpolateT[ii[0], ii[1]])**2, DTSQLIM)
#                     Vtip = .729*DTSQ + .0103 * DTSQ**2
#                     TIPLEN[ii[0], ii[1]] = -1*Vtip * DTIME / math.sqrt(2.0)
#==============================================================================
                    
        
    
    for i in xrange(m_min + 1, m_max - 1):
        for j in xrange(n_min + 1, n_max - 1):
            if AIR[i, j] == 1:
                continue

            """if SIDX[i, j] == 1 or SIDX[i, j] == 2: # nucleat and grow before"""
            #if SIDX[i, j] != 0: # not liquid. added on 05022017
            if SIDX[i, j] == 1 or SIDX[i, j] == 2 or SIDX[i, j] == 3: # not liquid. added on 10062017
                continue # if it has nucleated or grown, skip this cell.
            if _Temp_Liquid - TempInterpolateT[i,j] < 0.0 or _Temp_Liquid - TempInterpolate_Old[i,j] < 0.0: # cell temperature is higher than adjusted liquidus temperature.
                CLASS[i, j] = -1
                SIDX[i, j] = 0 # 05022017 added
                continue
            """in order to increase nucleation. 04282017"""
            if TempInterpolateT[i,j] < TempADJ:
                continue
            if TempInterpolate_Old[i,j] < TempInterpolateT[i,j]: # temperature is rising
                continue
            # Temperature should be declining
            """
            Gaussian Distribution
            scipy.stats.norm.pdf(, _UnderCoolMean, _UnderCoolDelta)
            for i in xrange (a.size):
            b.append(gaussian_int(d[i],.5, .1, 1))
            """
            r = abs(math.sin(random.random()))
            if TempInterpolateT[i,j] > _TempADJ and TempInterpolateT[i + 1,j] > _TempADJ and TempInterpolateT[i - 1,j] > _TempADJ and TempInterpolateT[i,j + 1] > _TempADJ and TempInterpolateT[i,j - 1] > _TempADJ: # represents four neighbours and itself are all above TempLQD
                PV = Fv_gaussian_nuc * gaussian_int((_Temp_Liquid - TempInterpolate_Old[i,j])/Integral_coefficient_v, (_Temp_Liquid - TempInterpolateT[i,j])/Integral_coefficient_v, 1.32, .1, 1)[0] / _Row_Num / _Column_Num # it nucleats
                CLASS_block = CLASS[i - 3:i + 3, j - 3:j + 3] # two nucleation center cannot be too close. The distance is at least 6 cell. 
                if r <= PV and np.all(CLASS_block == -1):
                    DTNUCL[i,j] = 10 # it nucleats in the bulk liquid.
                    print "**********************************************"
                    print "nucleat in the bulk"
                    print "**********************************************"
                    NUCFLG = 1 # 1 represents there is at least one nucleation for this time step.
                    print "_NUCINC = ", _NUCINC
                    _NUCINC = int(0 if _NUCINC is None else _NUCINC)
                    _NUCINC += 1 # Record how many 
                    CLASS[i,j] = random.randint(0,49)
                    CLASS[i, j] = CLASS2[i, j] if  CLASS2[i, j] > 0 else CLASS[i,j] # 10102017 added
                    SIDX[i,j] = 1 # become solid because of nucleation
                    print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the Nucleation bulk."
                    ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0
                    NUCIDX[i,j] = _NUCINC # the NUCINCth grain
                    DTSQ = min((_Temp_Liquid - TempInterpolateT[i,j])**2, DTSQLIM)
                    Vtip = .729*DTSQ + .0103 * DTSQ**2
                    TIPLEN[i,j] = -1*Vtip * DTIME / math.sqrt(2.0)
                else:
                    DTNUCL[i,j] = 0
            if TempInterpolateT[i,j] > _TempADJ and (TempInterpolateT[i + 1,j] < _TempADJ or TempInterpolateT[i - 1,j] < _TempADJ or TempInterpolateT[i,j + 1] < _TempADJ or TempInterpolateT[i,j - 1] < _TempADJ): # represents at least one neighbour is lower than TempLQD
                PS = Fs_gaussian_nuc * gaussian_int((_Temp_Liquid - TempInterpolate_Old[i,j])/Integral_coefficient_s, (_Temp_Liquid - TempInterpolateT[i,j])/Integral_coefficient_s, .5, .1, 1)[0] / _Row_Num / _Column_Num # it nucleats
                if r <= PS:
                    DTNUCL[i,j] = 5 # it nucleats at the S/L surface
                    print "**********************************************"
                    print "nucleat at S/L surface"
                    print "**********************************************"
                    NUCFLG = 1 # 1 represents there is at least one nucleation for this time step.
                    print "_NUCINC = ", _NUCINC
                    _NUCINC = int(0 if _NUCINC is None else _NUCINC)
                    _NUCINC += 1 # Record how many nucleations
                    CLASS[i,j] = random.randint(0,49)
                    CLASS[i, j] = CLASS2[i, j] if  CLASS2[i, j] > 0 else CLASS[i,j] # 10102017 added
                    a = [CLASS[i + 1, j], CLASS[i - 1, j], CLASS[i, j + 1], CLASS[i, j - 1]]               
                    CLASS_interface = a[min(range(len(a)), key=lambda i: abs(a[i] - Class_mean))] # prescribe the CLASS of interface cell when nucleation. This is to look for which neighbouring cell is solid. find i where abs(a[i] - Class_mean) is minimum. min() returns a number within range(len(a)). simulate epitaxial growth.
                    if CLASS_interface > 0:
                        CLASS[i,j] =  CLASS_interface
                    print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the Nucleation surface."
                    ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0
                    SIDX[i, j] = 1 # become solid because of nucleation
#                    print "ANGLE = ", ANGLE[i,j]
                    NUCIDX[i,j] = _NUCINC # the NUCINCth grain
                    DTSQ = min((_Temp_Liquid - TempInterpolateT[i,j])**2, DTSQLIM)
                    Vtip = .729*DTSQ + .0103 * DTSQ**2
                    TIPLEN[i,j] = -1*Vtip * DTIME / math.sqrt(2.0)
                else:
                    DTNUCL[i,j] = 0 #The cell won't nucleat

            """created 04282017 to solve no nucleation problem."""         
    # record the time for each nucleation
    for j in xrange(_Column_Num): # determine the CLASS of the toppest cell.
        if _Temp_Liquid - TempInterpolateT[_Row_Num-1,j] < 0.0:
            CLASS[_Row_Num-1, j] = -1
            SIDX[_Row_Num-1, j] = 0 # 05022017 added
        if _Temp_Liquid - TempInterpolateT[0,j] < 0.0:
            CLASS[0, j] = -1
            SIDX[0, j] = 0 # 05022017 added

    elapsed_time = timeit.default_timer() - start_time
    f = open(os.path.join(Destination_path, 'MonitorNucleation.txt'), 'a')
    time = str(datetime.now()) + '   ' + '?th run nucleation: ' + str(nucleat_count) + '   ' + 'time to run this nucleation: ' + str(elapsed_time) + '\n'
    f.writelines(time)
    f.close()
    return _NUCINC

def Growth(_Row_Num,_Column_Num, _Temp_Liquid,_Temp_Solid,temp_seq):
    start_time = timeit.default_timer()
    global growth_count, NUCFLG
    growth_count += 1
    growth_times = 0
    if NUCFLG == 0: #there is no nucleation in the history and at this step.
        return
    STOP.fill(0) #All cells can grow at the new time step.
    count, FLG = 0, 0 # the number of liquid CA cell in the computational domain. FLG == 0 represents that there is no liquid CA cell in the domain. FLG = 1 represents that there is still liquid CA cells in the domain. (FLG ONLY appears within Growth function)
    list_i = list()
    list_j = list()
    """
    determine the nucleation region
    """
    for i in xrange(Invalid_rows, _Row_Num): # replace _Row_Num with Invalid_rows because some rows far away from the melt pool won't melt and solidify. 
        for j in xrange(_Column_Num):
            if SIDX[i,j] == 1 or SIDX[i,j] == 2: # become solid because of nucleation or growth
                list_i.append(i)
                list_j.append(j)
    imin = min(list_i)
    imax = max(list_i)
    jmin = min(list_j)
    jmax = max(list_j)
    print "imin: ", imin
    print "imax: ", imax
    print "jmin: ", jmin
    print "jmax: ", jmax

    
    
    """
    detenmine which grains are growing(capturing the neighbor cells), and which are not. LenCritical array represents if the grain length is larger than a critical value Len0. 
    """
    nucleation_idx1 = np.argwhere(SIDX == 1) # index of all nucleated grain
    print "nucleation_idx1 = ", nucleation_idx1
    print "CLASS = ", CLASS[nucleation_idx1[:, 0], nucleation_idx1[:, 1]]
    nucleation_idx2 = np.argwhere(re_center == 1) # index of polygon center because of growing capture, not original nucleation.
    print "CLASS2 = ", CLASS[nucleation_idx2[:, 0], nucleation_idx2[:, 1]]
    nucleation_idx3 = np.concatenate((nucleation_idx1, nucleation_idx2), axis=0) # there may be some repeated indices.
    # Perform lex sort and get sorted data
    sorted_idx = np.lexsort(nucleation_idx3.T)
    print "sorted_idx = ", sorted_idx
    sorted_data =  nucleation_idx3[sorted_idx,:]    
    # Get unique row mask
    row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))
    print row_mask
    # Get unique rows
    if len(sorted_data) == 0:
        row_mask = []
    nucleation_idx = sorted_data[row_mask] # obtain the unionset of nucleation_idx1 and nucleation_idx2. there is no repeated indices.
    print "nucleation_idx is {0}, shape is {1}, dtype is {2}".format(nucleation_idx, nucleation_idx.shape, nucleation_idx.dtype)
    print "CLASS[nucleation_idx[:, 0], nucleation_idx[:, 1]] = ", CLASS[nucleation_idx[:, 0], nucleation_idx[:, 1]]
    """store nucleation_idx into nucleation_index through pickle. 04272017"""
    with open('nucleation_index', 'w') as handle1:
        nuc_dict[temp_seq] = nucleation_idx
        pickle.dump(nuc_dict, handle1)
        
    """figure out which cells are capturing, which stop capturing."""    
    for i in nucleation_idx: # iterate every nucleated(including original nucleated and capture nucleated) cell.
        if STOP[i[0], i[1]] == 1:
            continue
        else:
#        if NUCIDX[i[0], i[1]] > 0 and STOP[i[0], i[1]] == 0: deleted on 05232017
            Len0 = SDX*(abs(math.cos(ANGLE[i[0], i[1]]))+abs(math.sin(ANGLE[i[0], i[1]])))/math.sqrt(2.0) # critical length
            if TIPLEN[i[0], i[1]] >= Len0:
                LenCritical[i[0], i[1]] = 1 # length index array. 1 represents TIPLEN is beyond the critical length and neighboring cells are captured.
                #SIDX[i[0], i[1]] = 1 # 05262017 added. ???
                
                if  i[0] < _Row_Num - 1 and i[0] > 0 and i[1] > 0 and i[1] < _Column_Num - 1: # within the computation domain. avoid overflow the boundary.
                    neighbor_idx =np.array([[i[0]-1, i[1]], [i[0]+1, i[1]], [i[0], i[1]-1], [i[0], i[1]+1], [i[0]-1, i[1]-1], [i[0]-1, i[1]+1], [i[0]+1, i[1]+1], [i[0]+1, i[1]-1]]) # all eight neighboring cells index
                    if np.all(CLASS[neighbor_idx[:,0],neighbor_idx[:,1]] >= 0): # All eight neighbors have been captured. The center cell began to stop growing.
                        SIDX[i[0], i[1]] = 3 # All eight neighbors have been captured. The center cell state is 3.
                        LenCritical[i[0], i[1]] = 2 # 2 represents all 8 neighboring cells have been captured.
                        re_center[i[0], i[1]] = 2 # 2 represents that this cell cannot be the center of growing grain because 8 neighbors have been solid state.
                """boundary condition."""
                if i[1] == 0 and i[0] < _Row_Num - 1:
                    neighbor_idx_l = np.array([[i[0] - 1, i[1]], [i[0] - 1, i[1] + 1], [i[0], i[1] + 1], [i[0] + 1, i[1] + 1], [i[0] + 1, i[1]]])
                    if np.all(CLASS[neighbor_idx_l[:,0],neighbor_idx_l[:,1]] >= 0):
                        SIDX[i[0], i[1]] = 3 
                        LenCritical[i[0], i[1]] = 2 # 2 represents all 8 neighboring cells have been captured.
                        re_center[i[0], i[1]] = 2
                if i[1] == _Column_Num - 1 and i[0] < _Row_Num - 1:
                    neighbor_idx_r = np.array([[i[0] - 1, i[1]], [i[0] - 1, i[1] - 1], [i[0], i[1] - 1], [i[0] + 1, i[1] - 1], [i[0] + 1, i[1]]])
                    if np.all(CLASS[neighbor_idx_r[:,0],neighbor_idx_r[:,1]] >= 0):
                        SIDX[i[0], i[1]] = 3
                        LenCritical[i[0], i[1]] = 2
                        re_center[i[0], i[1]] = 2
            else: # TIPLEN is shorter than critical length, so it continues to grow. it cannot capture neighbor at this step.
                DTSQ = min((_Temp_Liquid - TempInterpolateT[i[0], i[1]])**2, DTSQLIM)
                Vtip = .729*DTSQ + .0103 * DTSQ**2
                TIPLEN[i[0], i[1]] = TIPLEN[i[0], i[1]] + Vtip * DTIME / math.sqrt(2.0)
                """in order to reduce computational cost, make it invalid temporarily."""
                TIPLEN1[i[0], i[1]] = TIPLEN2[i[0], i[1]] = TIPLEN3[i[0], i[1]] = TIPLEN4[i[0], i[1]] = TIPLEN[i[0], i[1]]


 
    """
    detenmine the surrounding polygons. Cover the cells within the fixed polygon.
    """ 
    
    """obtain all polygon vertices index whose centers have stopped growing because of capturing 8 neighbors. neighvertices_global_idx is what we need."""
    if np.count_nonzero(LenCritical == 2) != 0:# there is grains who have captured 8 neighbors. Only assigned the value of 2 at above sentence in Growth function. i.e. np.all(CLASS[neighbor_idx[:,0],neighbor_idx[:,1]] >= 0) is True. 05192917
        neighbor_solid_idx = np.argwhere(LenCritical == 2) # global index of center cell whose eight neighbors are all captured and solid. 
        GrowthClass_neigh = gg.Growth(ANGLE[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]], neighbor_solid_idx[:,0], neighbor_solid_idx[:,1], len(neighbor_solid_idx[:,0]))
        # calculate local coordinates(X1_neigh, Y1_neigh) and global index(SDX1_neigh, SDY1_neigh) of 4 polygon vertices  
        X1_neigh, Y1_neigh, SDX1_neigh, SDY1_neigh = GrowthClass_neigh.growth1(TIPLEN1[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
        X2_neigh, Y2_neigh, SDX2_neigh, SDY2_neigh = GrowthClass_neigh.growth2(TIPLEN2[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
        X3_neigh, Y3_neigh, SDX3_neigh, SDY3_neigh = GrowthClass_neigh.growth3(TIPLEN3[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
        X4_neigh, Y4_neigh, SDX4_neigh, SDY4_neigh = GrowthClass_neigh.growth4(TIPLEN4[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
        # neighvertices_global_idx are vertices global index.
        neighvertices_global_idx1 = np.dstack([SDX1_neigh,SDY1_neigh])[0] # first vertex coordinates of all selected grains
        neighvertices_global_idx2 = np.dstack([SDX2_neigh,SDY2_neigh])[0] # second vertex coordinates of all selected grains
        neighvertices_global_idx3 = np.dstack([SDX3_neigh,SDY3_neigh])[0] # third vertex coordinates of all selected grains
        neighvertices_global_idx4 = np.dstack([SDX4_neigh,SDY4_neigh])[0] # fourth vertex coordinates of all selected grains
        neighvertices_global_idxv = np.vstack([neighvertices_global_idx1, neighvertices_global_idx2, neighvertices_global_idx3, neighvertices_global_idx4]) # [[idx1], [idx2], [idx3], [idx4]], every idx1 contains the first vertex index of all grains.
        re_center[neighvertices_global_idxv[:,0], neighvertices_global_idxv[:,1]] = 1 # mark these cells will grow because of being captured.
        re_center[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] = 2 # 09252017 added.
        LenCritical[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] = 0 # 09252017 added. It is to prevent polygon growing if this polygon have covered 8 neighboring cells.
#        LenCritical[neighvertices_global_idxv[:,0], neighvertices_global_idxv[:,1]] = 2 # 10032017 added. 10092017 commented. 
#==============================================================================
#         STOP[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] = 1 # 09252017 revised. It is to prevent polygon growing if this polygon have covered 8 neighboring cells.
#==============================================================================
#==============================================================================
# ### 10092017 added after it is compared with v9 version.
#         for i in neighvertices_global_idxv: # there are other grains nearby. 
#             if np.unique(CLASS[i[0] - 2:i[0] + 2, i[1] - 2:i[1] + 2]).size > 2:
#                 re_center[i[0], i[1]] = 0
#         """10092017 commented because this part may not be useful."""
#==============================================================================

        """there are other grains nearby or all neighbors are solid."""
        for i in neighvertices_global_idxv: # For neighvertices_global_idxv array, the first vertices of each grain are together. 
            # 05192017 revised
            ###########################################
            if  i[0] < _Row_Num - 1 and i[0] > 0 and i[1] > 0 and i[1] < _Column_Num - 1:
                neighvertices_global_idxv_nei = np.array([[i[0]-1, i[1]], [i[0]+1, i[1]], [i[0], i[1]-1], [i[0], i[1]+1], [i[0]-1, i[1]-1], [i[0]-1, i[1]+1], [i[0]+1, i[1]+1], [i[0]+1, i[1]-1]])
                if np.all(CLASS[neighvertices_global_idxv_nei[:, 0], neighvertices_global_idxv_nei[:, 1]] >= 0):
                    re_center[i[0], i[1]] = 2 # eight neighbors are solid.
#==============================================================================
#                     LenCritical[i[0], i[1]] = 0 # 09252017 revised. it is to prevent polygon growing if this polygon have covered 8 neighboring cells.
#==============================================================================
                    
            #########################################
            if np.unique(CLASS[i[0] - nucleat_spacing:i[0] + nucleat_spacing, i[1] - nucleat_spacing:i[1] + nucleat_spacing]).size > 2: # 
                re_center[i[0], i[1]] = 0
                print "growing cell centers are too close"
            if np.unique(CLASS[i[0] - 1:i[0] + 1, i[1] - 1:i[1] + 1]).size > 1 and np.all(CLASS[i[0] - 1:i[0] + 1, i[1] - 1:i[1] + 1] != -1): # CLASS == -1 represents it is liquid. forbid that one grain cell nucleat at another grain domain. 05032017 added.
                re_center[i[0], i[1]] = 0
        
        neighvertices_global_idx = np.hstack([neighvertices_global_idx1, neighvertices_global_idx2, neighvertices_global_idx3, neighvertices_global_idx4]).reshape(-1, 4, 2) # 每一个grain的四个顶点的index在一起. neighvertices_global_idx stores all polygon vertices index whose centers have stopped growing because of capturing 8 neighbors.
        print neighvertices_global_idx
        
        """only select some polygon to cover enclosed cells randomly. 10102017added."""
#==============================================================================
#         if len(neighvertices_global_idx) >= 4: # 10102017 added.
#             neighvertices_global_idx = extract_array(neighvertices_global_idx, 0.95) # 10102017 added because only 60% polygon are in effect.
#==============================================================================

        """Each edge of polygon is offseted outward offset unit of SDX in order to enclose more cells. calculate local coordinates and global index of enlarged offset*SDX polygon vertices. This offset is to obtain larger polygon envelope."""
        X1_neigh_offset, Y1_neigh_offset, SDX1_neigh_offset, SDY1_neigh_offset = GrowthClass_neigh.growth1(TIPLEN1[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + offset*SDX)
        X2_neigh_offset, Y2_neigh_offset, SDX2_neigh_offset, SDY2_neigh_offset = GrowthClass_neigh.growth2(TIPLEN2[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + offset*SDX)
        X3_neigh_offset, Y3_neigh_offset, SDX3_neigh_offset, SDY3_neigh_offset = GrowthClass_neigh.growth3(TIPLEN3[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + offset*SDX)
        X4_neigh_offset, Y4_neigh_offset, SDX4_neigh_offset, SDY4_neigh_offset = GrowthClass_neigh.growth4(TIPLEN4[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + offset*SDX)       
        neighvertices1_global_off = np.dstack([X1_neigh_offset,Y1_neigh_offset])[0]
        neighvertices2_global_off = np.dstack([X2_neigh_offset,Y2_neigh_offset])[0]
        neighvertices3_global_off = np.dstack([X3_neigh_offset,Y3_neigh_offset])[0]
        neighvertices4_global_off = np.dstack([X4_neigh_offset,Y4_neigh_offset])[0]
        neighvertices_global_off = np.hstack([neighvertices1_global_off, neighvertices2_global_off, neighvertices3_global_off, neighvertices4_global_off]).reshape(-1, 4, 2) # stack four vertices index of all enlarged polygon. four vertices of each polygon are stacked together.
        
        """Prevent the very large grains appearing. The large polygon will be forbidden to grow again. 10052017 added."""
#==============================================================================
#         print neighvertices_global_idx[1]
#         print neighvertices_global_idx[1, 0, 1]
#         print neighvertices_global_idx.shape
#==============================================================================
        squaresum = np.sum((np.amax(neighvertices_global_idx, axis = 1) - np.amin(neighvertices_global_idx, axis = 1))**2, axis = 1)
        neighvertices_global_idx_index = np.argwhere(squaresum <= polygon_dia)
        neighvertices_global_idx = np.squeeze(neighvertices_global_idx[neighvertices_global_idx_index], axis = 1)
#==============================================================================
#         print neighvertices_global_idx[1]
#         print neighvertices_global_idx[1, 0, 1]
#         print neighvertices_global_idx.shape
#==============================================================================
#        simcard()

        
#==============================================================================
#         """enclose all cells within the polygon, which has covered eight neighbors. change these cells' CLASS state to a certain integer."""
#         for j, four_idx in enumerate(neighvertices_global_idx): # following sentence "print neighvertices_global_idx".
# #   enclose some cells into the polygon.
#             Quadrilateral_neigh = path.Path([(neighvertices_global_off[j, 0, 0], neighvertices_global_off[j, 0, 1]), (neighvertices_global_off[j, 1, 0], neighvertices_global_off[j, 1, 1]), (neighvertices_global_off[j, 2, 0], neighvertices_global_off[j, 2, 1]), (neighvertices_global_off[j, 3, 0], neighvertices_global_off[j, 3, 1])]) # a path determined by 4 vertices of the enlarged 2SDX polygon.
#             neigh_cells = CLASS[np.amin(four_idx, axis = 0)[0]:np.amax(four_idx, axis = 0)[0], np.amin(four_idx, axis = 0)[1]:np.amax(four_idx, axis = 0)[1]] # neigh_cells is part of CLASS, which determine the rectangular region covering 4 vertices location.
#             neigh_cells_global = np.argwhere(neigh_cells == -1) # all liquid cell local index with the above rectangular region.
#             neigh_cells_global[:, 0] =  neigh_cells_global[:, 0] + np.amin(four_idx, axis = 0)[0] # all liquid cell global index of first column with the above rectangular region.
#             neigh_cells_global[:, 1] =  neigh_cells_global[:, 1] + np.amin(four_idx, axis = 0)[1] # all liquid cell global index of second column with the above rectangular region.
#             neigh_cells_local_coord = SDX * (neigh_cells_global - np.tile(neighbor_solid_idx[j], (len(neigh_cells_global), 1))) # all liquid cell local coordinate  with the above rectangular region.
#             liquid_neigh_coord_tuple = totuple(neigh_cells_local_coord) # type conversion.
#             if liquid_neigh_coord_tuple != (): # there is liquid cell within the above rectangular region.
#                 EnclosePoint_neigh = Quadrilateral_neigh.contains_points(liquid_neigh_coord_tuple) # check which liquid cells are enclosed within the enlarged polygon. return boolean values. EnclosePoint_neigh is a boolean array.
#                 enclose_neigh_idx = np.argwhere(EnclosePoint_neigh == True) # local (NOT global) index of enclosed liquid cells.
#                 print enclose_neigh_idx
#                 CLASS[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = CLASS[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]] # assign same CLASS value to capturing neighbor cells
#                 NUCIDX[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = NUCIDX[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]] # assign same NUCIDX value to capturing neighbor cells
#                 ANGLE[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = ANGLE[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]] # assign same ANGLE value to capturing neighbor cells
# 
#             """give the initial state of vertices cell(will grow in future)."""
#             idx_valid = np.argwhere(CLASS[four_idx[:, 0], four_idx[:, 1]] == -1) # idx_valid represents all vertices index whose state is liquid. the purpose is to avoid nucleate on the substrate or solidified grain.
#             CLASS[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = CLASS[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
#             NUCIDX[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = NUCIDX[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
#             ANGLE[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = ANGLE[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
#==============================================================================

        """enclose all cells within the polygon, which has covered eight neighbors. change these cells' CLASS state to a certain integer. Revised on 10092017"""
        for j, four_idx in enumerate(neighvertices_global_idx): # following sentence "print neighvertices_global_idx".
#   enclose some cells into the polygon.
            Quadrilateral_neigh = path.Path([(neighvertices_global_off[j, 0, 0], neighvertices_global_off[j, 0, 1]), (neighvertices_global_off[j, 1, 0], neighvertices_global_off[j, 1, 1]), (neighvertices_global_off[j, 2, 0], neighvertices_global_off[j, 2, 1]), (neighvertices_global_off[j, 3, 0], neighvertices_global_off[j, 3, 1])]) # a path determined by 4 vertices of the enlarged 2SDX polygon.
            neigh_cells = CLASS[np.amin(four_idx, axis = 0)[0]:np.amax(four_idx, axis = 0)[0], np.amin(four_idx, axis = 0)[1]:np.amax(four_idx, axis = 0)[1]] # neigh_cells is part of CLASS, which determine the rectangular region covering 4 vertices location.
            neigh_cells_global = np.argwhere(neigh_cells == -1) # all liquid cell local index with the above rectangular region.
            neigh_cells_global[:, 0] =  neigh_cells_global[:, 0] + np.amin(four_idx, axis = 0)[0] # all liquid cell global index of first column with the above rectangular region.
            neigh_cells_global[:, 1] =  neigh_cells_global[:, 1] + np.amin(four_idx, axis = 0)[1] # all liquid cell global index of second column with the above rectangular region.
            neigh_cells_local_coord = SDX * (neigh_cells_global - np.tile(neighbor_solid_idx[j], (len(neigh_cells_global), 1))) # all liquid cell local coordinate  with the above rectangular region.
            liquid_neigh_coord_tuple = totuple(neigh_cells_local_coord) # type conversion.
            if liquid_neigh_coord_tuple != (): # there is liquid cell within the above rectangular region.
                EnclosePoint_neigh = Quadrilateral_neigh.contains_points(liquid_neigh_coord_tuple) # check which liquid cells are enclosed within the enlarged polygon. return boolean values. EnclosePoint_neigh is a boolean array.
                enclose_neigh_idx = np.argwhere(EnclosePoint_neigh == True) # local (NOT global) index of enclosed liquid cells.
                print enclose_neigh_idx
                CLASS[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = CLASS[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]] # assign same CLASS value to capturing neighbor cells
                REMELT_index = np.argwhere(CLASS2 > 0) # 10102017 added. make sure remelted grain keeps its original orientation
                CLASS[REMELT_index[:, 0], REMELT_index[:, 1]] = CLASS2[REMELT_index[:, 0], REMELT_index[:, 1]] # 10102017 added. make sure remelted grain keeps its original orientation
                NUCIDX[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = NUCIDX[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]] # assign same NUCIDX value to capturing neighbor cells
                ANGLE[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = ANGLE[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]] # assign same ANGLE value to capturing neighbor cells

            """give the initial state of vertices cell(will grow in future)."""
            idx_valid = np.argwhere(CLASS[four_idx[:, 0], four_idx[:, 1]] == -1) # idx_valid represents all vertices index whose state is liquid. the purpose is to avoid nucleate on the substrate or solidified grain.
            CLASS[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = CLASS[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
            REMELT_index = np.argwhere(CLASS2 > 0) # 10102017 added. make sure remelted grain keeps its original orientation
            CLASS[REMELT_index[:, 0], REMELT_index[:, 1]] = CLASS2[REMELT_index[:, 0], REMELT_index[:, 1]] # 10102017 added. make sure remelted grain keeps its original orientation
            NUCIDX[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = NUCIDX[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
            ANGLE[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = ANGLE[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]



    """
    grain growth (larger than Len0 and don't enclose eight neighbors completely)
    """
    
    """four vertices global index"""
    critical_len_idx = np.argwhere(LenCritical == 1) # global index of nucleated grains(a grain may include multiple growing centers) larger than critical length Len0.
    GrowthClass = gg.Growth(ANGLE[critical_len_idx[:,0], critical_len_idx[:,1]], critical_len_idx[:,0], critical_len_idx[:,1], len(critical_len_idx[:,0]))
    # SDX1, SDY1 are both global index. X1, Y1 are both local coordinates. 
    X1, Y1, SDX1, SDY1 = GrowthClass.growth1(TIPLEN1[critical_len_idx[:,0], critical_len_idx[:,1]])
    X2, Y2, SDX2, SDY2 = GrowthClass.growth2(TIPLEN2[critical_len_idx[:,0], critical_len_idx[:,1]])
    X3, Y3, SDX3, SDY3 = GrowthClass.growth3(TIPLEN3[critical_len_idx[:,0], critical_len_idx[:,1]])
    X4, Y4, SDX4, SDY4 = GrowthClass.growth4(TIPLEN4[critical_len_idx[:,0], critical_len_idx[:,1]])
    
    # Grain length tip growth.
    TIPLEN1[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN1[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX1, SDY1])
    TIPLEN2[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN2[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX2, SDY2])
    TIPLEN3[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN3[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX3, SDY3])
    TIPLEN4[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN4[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX4, SDY4])
    """comment on 10112017"""
#==============================================================================
#     # in order to determine the rectangular region covering the growing polygon.
#     X_MIN = np.minimum.reduce([SDX1, SDX2, SDX3, SDX4])
#     Y_MIN = np.minimum.reduce([SDY1, SDY2, SDY3, SDY4])
#     X_MAX = np.maximum.reduce([SDX1, SDX2, SDX3, SDX4])
#     Y_MAX = np.maximum.reduce([SDY1, SDY2, SDY3, SDY4])
#==============================================================================
    
    """Prevent the very large grains appearing. The large polygon will be forbidden to grow again. 10112017 added."""
    neighvertices_global_growing_idx1 = np.dstack([SDX1,SDY1])[0]
    neighvertices_global_growing_idx2 = np.dstack([SDX2,SDY2])[0]
    neighvertices_global_growing_idx3 = np.dstack([SDX3,SDY3])[0]
    neighvertices_global_growing_idx4 = np.dstack([SDX4,SDY4])[0]
    neighvertices_global_growing_idxv = np.vstack([neighvertices_global_growing_idx1, neighvertices_global_growing_idx2, neighvertices_global_growing_idx3, neighvertices_global_growing_idx4])

    neighvertices_global_growing_idx = np.hstack([neighvertices_global_growing_idx1, neighvertices_global_growing_idx2, neighvertices_global_growing_idx3, neighvertices_global_growing_idx4]).reshape(-1, 4, 2) # 每一个grain的四个顶点的index在一起. neighvertices_global_idx stores all polygon vertices index whose centers have stopped growing because of capturing 8 neighbors.
    squaresum = np.sum((np.amax(neighvertices_global_growing_idx, axis = 1) - np.amin(neighvertices_global_growing_idx, axis = 1))**2, axis = 1)
    neighvertices_global_growing_idxidx = np.argwhere(squaresum <= polygon_dia)
    neighvertices_global_growing_idx = np.squeeze(neighvertices_global_growing_idx[neighvertices_global_growing_idxidx], axis = 1)
    X_MIN = np.minimum.reduce([neighvertices_global_growing_idx[:, 0, 0], neighvertices_global_growing_idx[:, 1, 0], neighvertices_global_growing_idx[:, 2, 0], neighvertices_global_growing_idx[:, 3, 0]])
    Y_MIN = np.minimum.reduce([neighvertices_global_growing_idx[:, 0, 1], neighvertices_global_growing_idx[:, 1, 1], neighvertices_global_growing_idx[:, 2, 1], neighvertices_global_growing_idx[:, 3, 1]])
    X_MAX = np.maximum.reduce([neighvertices_global_growing_idx[:, 0, 0], neighvertices_global_growing_idx[:, 1, 0], neighvertices_global_growing_idx[:, 2, 0], neighvertices_global_growing_idx[:, 3, 0]])
    Y_MAX = np.maximum.reduce([neighvertices_global_growing_idx[:, 0, 1], neighvertices_global_growing_idx[:, 1, 1], neighvertices_global_growing_idx[:, 2, 1], neighvertices_global_growing_idx[:, 3, 1]])
    
    """there are other grains nearby or all neighbors are solid. 10122017 added."""
    for i in neighvertices_global_growing_idxv: # For neighvertices_global_growing_idxv array, the first vertices of each grain are together. 
#==============================================================================
#         if  i[0] < _Row_Num - 1 and i[0] > 0 and i[1] > 0 and i[1] < _Column_Num - 1:
#             neighvertices_global_growing_idxv_nei = np.array([[i[0]-1, i[1]], [i[0]+1, i[1]], [i[0], i[1]-1], [i[0], i[1]+1], [i[0]-1, i[1]-1], [i[0]-1, i[1]+1], [i[0]+1, i[1]+1], [i[0]+1, i[1]-1]])
#             if np.all(CLASS[neighvertices_global_growing_idxv_nei[:, 0], neighvertices_global_growing_idxv_nei[:, 1]] >= 0):
#                 re_center[i[0], i[1]] = 2 # eight neighbors are solid.
#==============================================================================
        if np.unique(CLASS[i[0] - nucleat_spacing:i[0] + nucleat_spacing, i[1] - nucleat_spacing:i[1] + nucleat_spacing]).size > 2: # 
            re_center[i[0], i[1]] = 0
            print "growing cell centers are too close"
        if np.unique(CLASS[i[0] - 1:i[0] + 1, i[1] - 1:i[1] + 1]).size > 1 and np.all(CLASS[i[0] - 1:i[0] + 1, i[1] - 1:i[1] + 1] != -1): # CLASS == -1 represents it is liquid. forbid that one grain cell nucleat at another grain domain. 05032017 added.
            re_center[i[0], i[1]] = 0

    """enclose the neighbor cells and change their states."""
    for i in range(len(X_MIN)): # ith nucleated grain.
        part_cells_class = CLASS[X_MIN[i]:X_MAX[i] + 1, Y_MIN[i]:Y_MAX[i] + 1] # a rectangle determined by polygon vertices, x_min, x_max, y_min, y_max.
        liquid_idx = np.argwhere(part_cells_class == -1) # index of all liquid cells within part_cells_class.
        Quadrilateral = path.Path([(X1[i], Y1[i]), (X2[i], Y2[i]), (X3[i], Y3[i]), (X4[i], Y4[i])])
        
        # Quadrilateral = Polygon(Point(X1[i], Y1[i]), Point(X2[i], Y2[i]), Point(X3[i], Y3[i]), Point(X4[i], Y4[i]))
        liquid_idx_global = liquid_idx + np.tile(np.dstack([X_MIN[i],Y_MIN[i]]), (len(liquid_idx), 1)) # liquid_idx_global is all global index of liquid cells within part_cells_class domain.
        liquid_coord_local = SDX * (liquid_idx_global - np.tile(critical_len_idx[i], (len(liquid_idx), 1))) # liquid_coord_local is local x and y coordinates of all liquid cells within part_cells_class domain.
        liquid_coord_tuple = totuple(liquid_coord_local)
        if liquid_coord_tuple[0] == ():
            continue
        EnclosePoint = Quadrilateral.contains_points(liquid_coord_tuple[0]) # boolean numpy . when it is true, it means it is enclosed by this grain polygon.
        enclose_idx = np.argwhere(EnclosePoint == True) # index of "True" value in  EnclosePoint. However, keep in mind that this index is not global index!!!
        print 'CLASS[{0}, {1}] is {2}'.format(critical_len_idx[i, 0], critical_len_idx[i, 1], CLASS[critical_len_idx[i, 0], critical_len_idx[i, 1]])
        CLASS[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = CLASS[critical_len_idx[i, 0], critical_len_idx[i, 1]] # this index is pretty important. liquid_idx_global is (1L, len(critical_len_idx), 2L) size. enclose_idx[:, 0] is used to select part of len(critical_len_idx).
        REMELT_index = np.argwhere(CLASS2 > 0) # 10102017 added. make sure remelted grain keeps its original orientation
        CLASS[REMELT_index[:, 0], REMELT_index[:, 1]] = CLASS2[REMELT_index[:, 0], REMELT_index[:, 1]] # 10102017 added. make sure remelted grain keeps its original orientation
        NUCIDX[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = NUCIDX[critical_len_idx[i, 0], critical_len_idx[i, 1]]
        ANGLE[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = ANGLE[critical_len_idx[i, 0], critical_len_idx[i, 1]]
        SIDX[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = 2 # become solid because of growth
        STOP[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = 1 # it won't grow again within this time step. It may be useless under current code.
        EnclosePoint.fill(False) # avoid all the grains have the same CLASS value.
   
     

    """
    Postprocess and save some results
    """
    
    """store critical_len_idx into critical_length_index file through pickle. 04272017"""
    with open('critical_length_index', 'w') as handle1:
        crit_len_idx[temp_seq] = critical_len_idx
        pickle.dump(crit_len_idx, handle1)
    # LenCritical.fill(0) #??? 05262017 delete
    growth_times += 1
    print '{0}th growth'.format(growth_times)
    
    """Make sure all the cells greater than liquidus temperature are -1 for CLASS variable. 09152017 revised"""
    AllLiquidIdx = np.argwhere(TempInterpolateT > _Temp_Liquid)
    CLASS[AllLiquidIdx[:, 0], AllLiquidIdx[:, 1]] = -1

    elapsed_time = timeit.default_timer() - start_time
    f = open(os.path.join(Destination_path, "MonitorGrowth.txt"),'a')
#    f = open(Destination_path + "MonitorGrowth.txt",'a')
    time = str(datetime.now()) + '   ' + '?th run nucleation: ' + str(growth_count) + '   ' + 'time to run this nucleation: ' + str(elapsed_time) + '   growth times: ' + str(growth_times) + '\n'
    f.writelines(time)
    f.close()
    return
    
def Output(_Row_Num,_Column_Num,i):
    os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer')
    FileName = "Class"
    f = open(FileName + str(i) + '.txt','w')
    f.writelines("Microstructure\n")
    f.writelines("Variables= 'I' 'J' 'CLASS'\n")
    f.writelines("ZONE I=681 J=401 F=POINT\n")
    for j in xrange(_Row_Num):
        for k in xrange(_Column_Num):
            seq = [str(j),' ',str(k),' ',str(CLASS[j,k]),'\n']
            f.writelines(seq)
    f.close()

    
def integrand(Delta_T, Delta_TN, Delta, N_max):
    return N_max/(np.sqrt(2*np.pi)*Delta)*np.exp(-1*np.power(Delta_TN - Delta_T, 2)/(2*Delta**2)) #from this paper "Three-Dimensional Probabilistic Simulation of Solidification Grain Structures : Application to Superalloy Precision Castings, 1993"

def gaussian_int(Delta_T_Old, Delta_T, Delta_TN, Delta, N_max):
    return integrate.quad(integrand, Delta_T_Old, Delta_T, args = (Delta_TN, Delta, N_max))

def CenterNuc(_NUCINC, _Temp_Liquid, *args): # invalid, 04272017
    global nucleat_count, NUCFLG
    nucleat_count += 1
#    r = abs(math.sin(random.random()))
    for i in xrange(args[2], args[3]):
        for j in xrange(args[0], args[1]):
            r = random.random()
            if r < Pc_nuc and CLASS[i, j] == -1:
                CLASS[i, j] = random.randint(0,48)
                DTNUCL[i,j] = 10 # it nucleats in the bulk liquid
                print "**********************************************"
                print "nucleat in the bulk"
                print "**********************************************"
                NUCFLG = 1 # 1 represents there is at least one nucleation for this time step.
                _NUCINC += 1 # Record how many.
                print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the CenterNuc."
                ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0
#                    print "ANGLE = ", ANGLE[i,j]
                SIDX[i,j] = 1 # become solid because of nucleation
                NUCIDX[i,j] = _NUCINC # the NUCINCth grain
                DTSQ = min((_Temp_Liquid - TempInterpolateT[i,j])**2, DTSQLIM)
                Vtip = .729*DTSQ + .0103 * DTSQ**2
                TIPLEN[i,j] = -1*Vtip * DTIME / math.sqrt(2.0)
            else:
                continue
    return _NUCINC

def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a
        
def air(): # called two times, 04272017
    air_idx1 = np.argwhere(TempInterpolateT < room_temp)
    air_idx = air_idx1[(air_idx1[:, 0] > air_sub_BDY)]
    CLASS[air_idx[:, 0], air_idx[:, 1]] = -10 # -10 represents it is air.
    AIR[air_idx[:, 0], air_idx[:, 1]] = 1 # 1 represents it is air.
    SIDX[air_idx[:, 0], air_idx[:, 1]] = -2
    re_center[air_idx[:, 0], air_idx[:, 1]] = -2 # ???

"""air cells are determined by FEA material addition time step. 09182017"""
def air2(col_count, CurrentTempFile, LayerEleThick):
    air_sub = air_sub_BDY + LayerEleThick + LayerEleThick * (CurrentTempFile // col_count)
    air_idx1 = np.argwhere(TempInterpolateT)
    air_liquid_idx = air_idx1[(air_idx1[:, 0] > air_sub - LayerEleThick)]
    CLASS[air_liquid_idx[:, 0], air_liquid_idx[:, 1]], AIR[air_liquid_idx[:, 0], air_liquid_idx[:, 1]], SIDX[air_liquid_idx[:, 0], air_liquid_idx[:, 1]], re_center[air_liquid_idx[:, 0], air_liquid_idx[:, 1]] = -1, 0, -1, 0
    air_idx2 = air_idx1[(air_idx1[:, 0] > air_sub)]
    
    #laser scan along forward and backward direction.
    if (CurrentTempFile // col_count) % 2 == 0:
        air_idx3 = air_idx1[((air_idx1[:, 0] > air_sub - LayerEleThick) & (air_idx1[:, 0] <= air_sub) & (air_idx1[:, 1] > InterpolationTimes * (CurrentTempFile % col_count) + LaserEleNum))]
    else:
        air_idx3 = air_idx1[((air_idx1[:, 0] > air_sub - LayerEleThick) & (air_idx1[:, 0] <= air_sub) & (air_idx1[:, 1] < ColumnNum - InterpolationTimes * (CurrentTempFile % col_count) - LaserEleNum))]
    air_idx = np.concatenate([air_idx2, air_idx3])
    CLASS[air_idx[:, 0], air_idx[:, 1]] = -10 # -10 represents it is air.
    AIR[air_idx[:, 0], air_idx[:, 1]] = 1 # 1 represents it is air.
    SIDX[air_idx[:, 0], air_idx[:, 1]] = -2
    re_center[air_idx[:, 0], air_idx[:, 1]] = -2 # ???
### 10102017 added ###
def extract_array(arr1, percentage): # extract an array according to random ordering.
    arr_conv = np.random.permutation(arr1)
    arr2 = arr_conv[0: int(percentage * len(arr_conv))]
    return arr2

### 10102017 added ###      
def extract_array2(arr_temp, percentage):# extract an array according to minimum/maximum ordering.
    pass
#==============================================================================
#     index_ordering = np.dstack(np.unravel_index(np.argsort(arr_temp.ravel()), arr_temp.shape))
#     index_ordering = np.squeeze(index_ordering, axis = 0)
#==============================================================================

    
"""
Main program
"""
NUCINC, TempInterpolate_Old = initialization(RowNum,ColumnNum,TempLiquid,0) #NUCINC value after initialization, DTNUCL, SIDX,NUCIDX will be changed in initialization function.
print "The number of nucleated grains in the initialization is ", NUCINC
for i in xrange(TempFileNumber_start, TempFileNumber): # provide temperature file range.
    if i == 0: # The first temperature file
        continue

    if i > 0:
        print "**********************************************"
        print i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i,i
        print "**********************************************"
        TempNew = get_temp(RowNum,ColumnNum,i)
        TempOld = get_temp(RowNum,ColumnNum,i-1)
        for x in xrange(MicroTime):
            print "MicroTime = ", x
            TempInterpolateT = TempOld + (x + 1) * 1.0 / MicroTime * (TempNew - TempOld)
            if __name__ == '__main__': # BUG test
                if i == remelt_flag: # begin remelting. 
                    liquid_remelt_idx = np.argwhere(TempInterpolateT > TempADJ) # find all indices whose values are greater than TempADJ
                    CLASS[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = -1 # liquid again because of remelting
                    SIDX[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = 0
                    DTNUCL[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = 0
                    NUCIDX[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = 0
                    air() # assign all air cells CLASS -10
                    NUCINC = Nucleation(RowNum, ColumnNum, TempLiquid, TempSolid, UnderCoolMean, UnderCoolDelta, NUCINC, TempADJ, i)
                    Growth(RowNum,ColumnNum,TempLiquid,TempSolid)
                    Output(RowNum,ColumnNum,(i-1)*MicroTime+x)
                    ClassExtract = dp.DataProcess('Extracting the data ' + str((i-1)*MicroTime+x), 'Class', (i-1)*MicroTime+x)
                    ClassExtract.ExtractDataMulti()
                    TempInterpolate_Old = TempInterpolateT
                    continue
            """remelt trial 2. 04272017 revise"""
#==============================================================================
#             liquid_remelt_idx = np.argwhere(TempInterpolateT > TempLiquid) # find all indices whose values are greater than TempLiquid
#             CLASS[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = -1 # liquid again because of remelting
#             SIDX[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = 0
#             DTNUCL[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = 0
#             NUCIDX[liquid_remelt_idx[:, 0], liquid_remelt_idx[:, 1]] = 0
#==============================================================================
            """remelt trial 1. doesn't work well. March 2017 revise."""
            if i > first_layer: # 10102017 added.
                remelt_grain_idx1 = np.argwhere(TempInterpolateT > TempLiquid) # obtain remelt grain index.
                remelt_grain_idx2 = np.argwhere((TempInterpolateT <= TempLiquid) & (CLASS == -1)) # obtain un-solidified cell index.
    #            remelt_grain_idx1 = np.argwhere(((SIDX == 1) | (SIDX == 2) | (SIDX == 3)) & (TempInterpolateT > TempADJ))
    #            remelt_grain_idx1 = np.argwhere(((SIDX == 1) | (SIDX == 2) | (SIDX == 3)) & (TempInterpolateT > TempADJ) & (TempInterpolateT > TempInterpolate_Old)) # obtain remelt grain index.
                REMELT[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]], REMELT[remelt_grain_idx2[:, 0], remelt_grain_idx2[:, 1]] = 1, 2 # to solve remelting cells donot solidify again. 10182017 changed.
                CLASS2[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] = CLASS[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] # 10102017 added
                CLASS[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] = -1
                SIDX[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] = 0
                SIDX[remelt_grain_idx2[:, 0], remelt_grain_idx2[:, 1]] = 0  # to solve remelting cells donot solidify again. 10182017 changed.
                NUCIDX[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] = 0
                re_center[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] = -1 #??? 04282017 revise
                LenCritical[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] = -1 #??? 04292017 revise. if it is liquid, it is -1.\
                del remelt_grain_idx1, remelt_grain_idx2 # 10182017 changed.
#==============================================================================
#             STOP[remelt_grain_idx1[:, 0], remelt_grain_idx1[:, 1]] = 0 # 09252017 revised. 0 indicates that this cell can grow because of no overgrowth.
#==============================================================================
            air2(ColCount, i, LayerEleThick) # assign all air cells CLASS -10. Revised @09182017
            NUCINC = Nucleation(RowNum, ColumnNum, TempLiquid, TempSolid, UnderCoolMean, UnderCoolDelta, NUCINC, TempADJ, i)
            print "The number of nucleated grains after nucleation is ", NUCINC
            Growth(RowNum,ColumnNum,TempLiquid,TempSolid, i)
            Output(RowNum,ColumnNum,(i-1)*MicroTime+x)
            ClassExtract = dp.DataProcess('Extracting the data ' + str((i-1)*MicroTime+x), 'Class', (i-1)*MicroTime+x)
            ClassExtract.ExtractDataMulti()
            # vectorized code. 09292017 added.
            TempInterpolate_Old = TempInterpolateT
#==============================================================================
                # non-vectorized code
#             for y in xrange(RowNum):
#                 for z in xrange(ColumnNum):
#                     TempInterpolate_Old[y, z] = TempInterpolateT[y, z]
#==============================================================================
                

for k in xrange((TempFileNumber_start - 1) * MicroTime, (TempFileNumber - 1) * MicroTime):
    dd.ImageToPNGMulti('Class', k, ColumnNum) # then www.gifmaker.me can convert pngs to gif animation.
"""show the animation, but it doesn't work using Python Console, works using IPython Console sometimes."""
            
print "LIDXN1", LIDXN1
print  "NUCIDX", NUCIDX
print "DTNUCL", DTNUCL
print "SIDX", SIDX
print "LenCritical", LenCritical
os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer')

np.savetxt("DTNUCL.txt",DTNUCL,fmt='%i',delimiter=',')
np.savetxt("SIDX.txt",SIDX,fmt='%i',delimiter=',')
np.savetxt("NUCIDX.txt",NUCIDX,fmt='%i',delimiter=',')
np.savetxt("LIDXN1.txt",LIDXN1,fmt='%i',delimiter=',')
np.savetxt("LenCritical.txt",LenCritical,fmt='%i',delimiter=',')

"""it is useful when not for debugging."""
with open('Var_pickle', 'w') as f1:
    pickle.dump([DTNUCL, SIDX, NUCIDX, LIDXN1, LenCritical, re_center, CLASS, REMELT, AIR], f1)
with open( 'Var_pickle2', 'w') as f2:
    pickle.dump([TIPLEN, TIPLEN1, TIPLEN2, TIPLEN3, TIPLEN4], f2)

