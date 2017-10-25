# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:05:02 2017

@author: jnzzp5
"""

"""
run the code in the IPython Console. run "%matplotlib qt" at first, then run "%run C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer\main_program.py"
main_program_v5 is to employ a decentered polygon algorithm.(unnecessary to iterate every cell.) it is for fine(5um) cell size.
Integral_coefficient is needed.
finer temperature file.
Growth conditions： 8 neighbors are solid, it continue growing.
try vectorizing code.
02212017 modified. copy from main_program_v9_vect_invalidcolumn.py
considering stoping grain growth
02272017 revised.
10062017 created. just delete all the comment compared to v9_vect_stopgrowth 
"""
import numpy as np
#import scipy.stats
import scipy.integrate as integrate
from io import StringIO
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

DTNUCL = np.zeros((RowNum,ColumnNum),dtype=np.int)
SIDX = np.zeros((RowNum,ColumnNum),dtype=np.int)
NUCIDX = np.zeros((RowNum,ColumnNum),dtype=np.int)
LIDXN1 = np.zeros((RowNum,ColumnNum),dtype=np.int) # only represents the initial(first temperature file) temperature
DEACTIVATION = np.zeros((RowNum,ColumnNum),dtype=np.int)
TempInterpolateT = np.zeros((RowNum,ColumnNum),dtype='float16')
STOP = np.zeros((RowNum,ColumnNum),dtype=np.int)
TIPLEN = np.zeros((RowNum,ColumnNum),dtype='float16')
"""make it invalid to reduce the computational cost! This is temporary."""
TIPLEN1 = np.zeros((RowNum,ColumnNum),dtype='float16')
TIPLEN2 = np.zeros((RowNum,ColumnNum),dtype='float16')
TIPLEN3 = np.zeros((RowNum,ColumnNum),dtype='float16')
TIPLEN4 = np.zeros((RowNum,ColumnNum),dtype='float16')
EnclosePoint = np.full((RowNum,ColumnNum), False, dtype=bool)
LenCritical = np.zeros((RowNum,ColumnNum), dtype='int32')
re_center = np.zeros((RowNum,ColumnNum), dtype='int32')
"""********************************"""
TempInterpolate_Old = np.full((RowNum,ColumnNum), TempLiquid, dtype='float16')
#TIPLEN_pd = pd.DataFrame({'TIPLEN0': np.zeros((RowNum,ColumnNum),dtype='float16'), 'TIPLEN1': np.zeros((RowNum,ColumnNum),dtype='float16'), 'TIPLEN2': np.zeros((RowNum,ColumnNum),dtype='float16'), 'TIPLEN3': np.zeros((RowNum,ColumnNum),dtype='float16'), 'TIPLEN4': np.zeros((RowNum,ColumnNum),dtype='float16')})
NUCINC, NUC, growth_count, nucleat_count = 0, 0, 0, 0
NUCFLG = 0
TempADJ = TempLiquid - TempAdjustor # adjusted TempLiquid
x_grid, y_grid = 3, 3 # control the size of substrate grain in the base.

Cell = dp.Substrate(RowNum, ColumnNum, x_grid, y_grid) # Cell is a class object.
CLASS = Cell.substrate(.2, .25)
ANGLE = (CLASS*90/48.0 - 45) * PI / 180.0 
def get_temp(_Row_Num, _Column_Num, _TempFileNumber):
    Temp=np.zeros((_Row_Num,_Column_Num),dtype='float64') # this is created to store temperature data
    try:
        c = 2
        if c == 0: # old temperature 2014
            f=file("C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate_finer"+str(_TempFileNumber),'r')
        elif c == 1: # artificial corrected temperature
            f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_correct2014error_02142017\Nodeset-SET-FACE2_1layerTempTi64_Interpolate_finer"+str(_TempFileNumber),'r') 
        elif c == 2:
            temp = np.genfromtxt(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_correct2014error_02142017\Nodeset-SET-FACE2_1layerTempTi64_Interpolate_finer"+str(_TempFileNumber), delimiter=",")

    except IOError:
        print "This temperature file cannot be found"
    else:
        temp_flip = np.fliplr(temp) # flip left over right(between columns).
    return temp_flip
    
    
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

def initialization(_Row_Num,_Column_Num,_Temp_Liquid,_NUCINC):
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
    
    for i in xrange(m_min + 1, m_max - 1):
        for j in xrange(n_min + 1, n_max - 1):
            if DTNUCL[i,j] != 0: #have nucleated at the interface or in the bulk liquid at previous time steps
                NUCFLG = 1 # nucleation at this cell has occurred in the history
                continue
            if SIDX[i, j] == 1 or SIDX[i, j] == 2: # nucleat and grow before
                continue # if it has nucleated or grown, skip this cell.
            if _Temp_Liquid - TempInterpolateT[i,j] < 0.0 or _Temp_Liquid - TempInterpolate_Old[i,j] < 0.0: # cell temperature is higher than adjusted liquidus temperature.
                CLASS[i, j] = -1
                continue
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
                PV = 10 * Fv_gaussian_nuc * gaussian_int((_Temp_Liquid - TempInterpolate_Old[i,j])/Integral_coefficient, (_Temp_Liquid - TempInterpolateT[i,j])/Integral_coefficient, 1.32, .1, 1)[0] / _Row_Num / _Column_Num # it nucleats
                CLASS_block = CLASS[i - 3:i + 3, j - 3:j + 3] # two nucleation center cannot be too close. The distance is at least 6 cell. 
                if r <= PV and np.all(CLASS_block == -1):
                    DTNUCL[i,j] = 10 # it nucleats in the bulk liquid.
                    print "**********************************************"
                    print "nucleat in the bulk"
                    print "**********************************************"
                    NUCFLG = 1 # 1 represents there is at least one nucleation for this time step.
                    _NUCINC += 1 # Record how many 
                    CLASS[i,j] = random.randint(0,49)
                    print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the Nucleation."
                    ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0
                    SIDX[i,j] = 1 # become solid because of nucleation
                    NUCIDX[i,j] = _NUCINC # the NUCINCth grain
                    DTSQ = min((_Temp_Liquid - TempInterpolateT[i,j])**2, DTSQLIM)
                    Vtip = .729*DTSQ + .0103 * DTSQ**2
                    TIPLEN[i,j] = -1*Vtip * DTIME / math.sqrt(2.0)
                else:
                    DTNUCL[i,j] = 0
            if TempInterpolateT[i,j] > _TempADJ and (TempInterpolateT[i + 1,j] < _TempADJ or TempInterpolateT[i - 1,j] < _TempADJ or TempInterpolateT[i,j + 1] < _TempADJ or TempInterpolateT[i,j - 1] < _TempADJ): # represents at least one neighbour is lower than TempLQD
                PS = 20 * Fs_gaussian_nuc * gaussian_int((_Temp_Liquid - TempInterpolate_Old[i,j])/Integral_coefficient, (_Temp_Liquid - TempInterpolateT[i,j])/Integral_coefficient, .5, .1, 1)[0] / _Row_Num / _Column_Num # it nucleats
                if r <= PS:
                    DTNUCL[i,j] = 5 # it nucleats at the S/L surface
                    print "**********************************************"
                    print "nucleat at S/L surface"
                    print "**********************************************"
                    NUCFLG = 1 # 1 represents there is at least one nucleation for this time step.
                    _NUCINC += 1 # Record how many nucleations
                    CLASS[i,j] = random.randint(0,49)
                    a = [CLASS[i + 1, j], CLASS[i - 1, j], CLASS[i, j + 1], CLASS[i, j - 1]]               
                    CLASS[i, j] = a[min(range(len(a)), key=lambda i: abs(a[i] - Class_mean))] # prescribe the CLASS of interface cell when nucleation. This is to look for which neighbouring cell is solid. find i where abs(a[i] - Class_mean) is minimum. min() returns a number within range(len(a)).
                    print "CLASS[{0},{1}] = ".format(i, j), CLASS[i, j], ", which is formed in the Nucleation."
                    ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0
                    SIDX[i,j] = 1 # become solid because of nucleation
                    NUCIDX[i,j] = _NUCINC # the NUCINCth grain
                    DTSQ = min((_Temp_Liquid - TempInterpolateT[i,j])**2, DTSQLIM)
                    Vtip = .729*DTSQ + .0103 * DTSQ**2
                    TIPLEN[i,j] = -1*Vtip * DTIME / math.sqrt(2.0)
                else:
                    DTNUCL[i,j] = 0 #The cell won't nucleat

    # record the time for each nucleation
    for j in xrange(_Column_Num): # determine the CLASS of the toppest cell.
        if _Temp_Liquid - TempInterpolateT[_Row_Num-1,j] < 0.0:
            CLASS[_Row_Num-1, j] = -1
        if _Temp_Liquid - TempInterpolateT[0,j] < 0.0:
            CLASS[0, j] = -1

    elapsed_time = timeit.default_timer() - start_time
    f = open('MonitorNucleation.txt','a')
    time = str(datetime.now()) + '   ' + '?th run nucleation: ' + str(nucleat_count) + '   ' + 'time to run this nucleation: ' + str(elapsed_time) + '\n'
    f.writelines(time)
    f.close()
    return _NUCINC

def Growth(_Row_Num,_Column_Num, _Temp_Liquid,_Temp_Solid):
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
    detenmine which grains are growing(capturing the neighbor cells) 
    """
    nucleation_idx1 = np.argwhere(SIDX == 1) # index of all nucleated grain
    nucleation_idx2 = np.argwhere(re_center == 1) # index of polygon center because of growing capture, not original nucleation.
    nucleation_idx3 = np.concatenate((nucleation_idx1, nucleation_idx2), axis=0) # there may be some repeated indices.
    print "nucleation_idx1 is {0}, shape is {1}, dtype is {2}".format(nucleation_idx1, nucleation_idx1.shape, nucleation_idx1.dtype)
    print "nucleation_idx2 is {0}, shape is {1}, dtype is {2}".format(nucleation_idx2, nucleation_idx2.shape, nucleation_idx2.dtype)
    print "nucleation_idx3 is {0}, shape is {1}".format(nucleation_idx3, nucleation_idx3.shape)
        # Perform lex sort and get sorted data
    sorted_idx = np.lexsort(nucleation_idx3.T)
    sorted_data =  nucleation_idx3[sorted_idx,:]    
    # Get unique row mask
    row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))
    # Get unique rows
    nucleation_idx = sorted_data[row_mask] # obtain the unionset of nucleation_idx1 and nucleation_idx2
    print "nucleation_idx is {0}, shape is {1}, dtype is {2}".format(nucleation_idx, nucleation_idx.shape, nucleation_idx.dtype)
    
    """give the initial state of each cell"""
    for i in nucleation_idx:
        if NUCIDX[i[0], i[1]] > 0 and STOP[i[0], i[1]] == 0:
            Len0 = SDX*(abs(math.cos(ANGLE[i[0], i[1]]))+abs(math.sin(ANGLE[i[0], i[1]])))/math.sqrt(2.0)
            if TIPLEN[i[0], i[1]] >= Len0:
                LenCritical[i[0], i[1]] = 1 # length index array. 1 represents neighboring cells are captured.
                if  i[0] < _Row_Num - 1 and i[0] > 0 and i[1] > 0 and i[1] < _Column_Num - 1: # within the computation domain. avoid overflow the boundary.
                    neighbor_idx =np.array([[i[0]-1, i[1]], [i[0]+1, i[1]], [i[0], i[1]-1], [i[0], i[1]+1], [i[0]-1, i[1]-1], [i[0]-1, i[1]+1], [i[0]+1, i[1]+1], [i[0]+1, i[1]-1]]) # all eight neighboring cells index
                    if np.all(CLASS[neighbor_idx[:,0],neighbor_idx[:,1]] >= 0): # All eight neighbors have been captured.
                        SIDX[i[0], i[1]] = 3 #???
                        LenCritical[i[0], i[1]] = 2 # 2 represents all 8 neighboring cells have been captured.
                        re_center[i[0], i[1]] = 0 # 0 represents that this cell cannot be the center of growing grain.
                # boundary condition.
                if i[1] == 0 and i[0] < _Row_Num - 1:
                    neighbor_idx_l = np.array([[i[0] - 1, i[1]], [i[0] - 1, i[1] + 1], [i[0], i[1] + 1], [i[0] + 1, i[1] + 1], [i[0] + 1, i[1]]])
                    if np.all(CLASS[neighbor_idx_l[:,0],neighbor_idx_l[:,1]] >= 0):
                        SIDX[i[0], i[1]] = 3 #???
                        LenCritical[i[0], i[1]] = 2 # 2 represents all 8 neighboring cells have been captured.
                        re_center[i[0], i[1]] = 0
                if i[1] == _Column_Num - 1 and i[0] < _Row_Num - 1:
                    neighbor_idx_r = np.array([[i[0] - 1, i[1]], [i[0] - 1, i[1] - 1], [i[0], i[1] - 1], [i[0] + 1, i[1] - 1], [i[0] + 1, i[1]]])
                    if np.all(CLASS[neighbor_idx_r[:,0],neighbor_idx_r[:,1]] >= 0):
                        SIDX[i[0], i[1]] = 3
                        LenCritical[i[0], i[1]] = 2
                        re_center[i[0], i[1]] = 0
            else: # control the growth of nucleated cells whose length is shorter than critical length Len0. doesn't control the growth of nucleated cells whose length is longer than critical length Len0. 05192017
                DTSQ = min((_Temp_Liquid - TempInterpolateT[i[0], i[1]])**2, DTSQLIM)
                Vtip = .729*DTSQ + .0103 * DTSQ**2
                TIPLEN[i[0], i[1]] = TIPLEN[i[0], i[1]] + .45 * Vtip * DTIME / math.sqrt(2.0)
                """in order to reduce computational cost, make it invalid temporarily."""
                TIPLEN1[i[0], i[1]] = TIPLEN2[i[0], i[1]] = TIPLEN3[i[0], i[1]] = TIPLEN4[i[0], i[1]] = TIPLEN[i[0], i[1]]
        else:
            continue

    """
    detenmine the surrounding polygons 
    """        
    if np.count_nonzero(LenCritical == 2) != 0:# there is grain capturing 8 neighbors.
        neighbor_solid_idx = np.argwhere(LenCritical == 2) # global index of whose eight neighbors are all captured and solid.    
        GrowthClass_neigh = gg.Growth(ANGLE[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]], neighbor_solid_idx[:,0], neighbor_solid_idx[:,1], len(neighbor_solid_idx[:,0]))
        X1_neigh, Y1_neigh, SDX1_neigh, SDY1_neigh = GrowthClass_neigh.growth1(TIPLEN1[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
        X2_neigh, Y2_neigh, SDX2_neigh, SDY2_neigh = GrowthClass_neigh.growth2(TIPLEN2[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
        X3_neigh, Y3_neigh, SDX3_neigh, SDY3_neigh = GrowthClass_neigh.growth3(TIPLEN3[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
        X4_neigh, Y4_neigh, SDX4_neigh, SDY4_neigh = GrowthClass_neigh.growth4(TIPLEN4[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]])
    # neighvertices_global_idx are vertices global index.
        neighvertices_global_idx1 = np.dstack([SDX1_neigh,SDY1_neigh])[0] # first vertex coordinates of all selected grains
        neighvertices_global_idx2 = np.dstack([SDX2_neigh,SDY2_neigh])[0] # second vertex coordinates of all selected grains
        neighvertices_global_idx3 = np.dstack([SDX3_neigh,SDY3_neigh])[0] # third vertex coordinates of all selected grains
        neighvertices_global_idx4 = np.dstack([SDX4_neigh,SDY4_neigh])[0] # fourth vertex coordinates of all selected grains
        neighvertices_global_idxv = np.vstack([neighvertices_global_idx1, neighvertices_global_idx2, neighvertices_global_idx3, neighvertices_global_idx4]) # [[idx1], [idx2], [idx3], [idx4]], every idx1 contains all grains.
        re_center[neighvertices_global_idxv[:,0], neighvertices_global_idxv[:,1]] = 1 # mark these cells will grow because of being captured.
        for i in neighvertices_global_idxv: # there are other grains nearby. 
            if np.unique(CLASS[i[0] - 2:i[0] + 2, i[1] - 2:i[1] + 2]).size > 2:
                re_center[i[0], i[1]] = 0

        neighvertices_global_idx = np.hstack([neighvertices_global_idx1, neighvertices_global_idx2, neighvertices_global_idx3, neighvertices_global_idx4]).reshape(-1, 4, 2) # 每一个grain的四个顶点的index在一起
        print neighvertices_global_idx1, neighvertices_global_idx2, neighvertices_global_idx3, neighvertices_global_idx4
        print neighvertices_global_idx

        # in order to enclose more cells.
        X1_neigh_offset, Y1_neigh_offset, SDX1_neigh_offset, SDY1_neigh_offset = GrowthClass_neigh.growth1(TIPLEN1[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + 2*SDX)
        X2_neigh_offset, Y2_neigh_offset, SDX2_neigh_offset, SDY2_neigh_offset = GrowthClass_neigh.growth2(TIPLEN2[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + 2*SDX)
        X3_neigh_offset, Y3_neigh_offset, SDX3_neigh_offset, SDY3_neigh_offset = GrowthClass_neigh.growth3(TIPLEN3[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + 2*SDX)
        X4_neigh_offset, Y4_neigh_offset, SDX4_neigh_offset, SDY4_neigh_offset = GrowthClass_neigh.growth4(TIPLEN4[neighbor_solid_idx[:,0], neighbor_solid_idx[:,1]] + 2*SDX)       
        neighvertices1_global_off = np.dstack([X1_neigh_offset,Y1_neigh_offset])[0]
        neighvertices2_global_off = np.dstack([X2_neigh_offset,Y2_neigh_offset])[0]
        neighvertices3_global_off = np.dstack([X3_neigh_offset,Y3_neigh_offset])[0]
        neighvertices4_global_off = np.dstack([X4_neigh_offset,Y4_neigh_offset])[0]
        neighvertices_global_off = np.hstack([neighvertices1_global_off, neighvertices2_global_off, neighvertices3_global_off, neighvertices4_global_off]).reshape(-1, 4, 2)
        """enclose all cells within the polygon, which has covered eight neighbors"""
        for j, four_idx in enumerate(neighvertices_global_idx):
#   enclose some cells into the polygon.
            Quadrilateral_neigh = path.Path([(neighvertices_global_off[j, 0, 0], neighvertices_global_off[j, 0, 1]), (neighvertices_global_off[j, 1, 0], neighvertices_global_off[j, 1, 1]), (neighvertices_global_off[j, 2, 0], neighvertices_global_off[j, 2, 1]), (neighvertices_global_off[j, 3, 0], neighvertices_global_off[j, 3, 1])])
            neigh_cells = CLASS[np.amin(four_idx, axis = 0)[0]:np.amax(four_idx, axis = 0)[0], np.amin(four_idx, axis = 0)[1]:np.amax(four_idx, axis = 0)[1]]
            neigh_cells_global = np.argwhere(neigh_cells == -1)
            neigh_cells_global[:, 0] =  neigh_cells_global[:, 0] + np.amin(four_idx, axis = 0)[0]
            neigh_cells_global[:, 1] =  neigh_cells_global[:, 1] + np.amin(four_idx, axis = 0)[1]
            neigh_cells_local_coord = SDX * (neigh_cells_global - np.tile(neighbor_solid_idx[j], (len(neigh_cells_global), 1)))
            liquid_neigh_coord_tuple = totuple(neigh_cells_local_coord)
            if liquid_neigh_coord_tuple != ():
                EnclosePoint_neigh = Quadrilateral_neigh.contains_points(liquid_neigh_coord_tuple)
                enclose_neigh_idx = np.argwhere(EnclosePoint_neigh == True)
                print enclose_neigh_idx
                CLASS[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = CLASS[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
                NUCIDX[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = NUCIDX[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
                ANGLE[neigh_cells_global[enclose_neigh_idx[:, 0],0], neigh_cells_global[enclose_neigh_idx[:, 0],1]] = ANGLE[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]

# give the initial state of vertices cell(will grow in foture).
            idx_valid = np.argwhere(CLASS[four_idx[:, 0], four_idx[:, 1]] == -1)
            CLASS[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = CLASS[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
            NUCIDX[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = NUCIDX[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]
            ANGLE[four_idx[idx_valid, 0], four_idx[idx_valid, 1]] = ANGLE[neighbor_solid_idx[j, 0], neighbor_solid_idx[j, 1]]

    """grain growth (larger than Len0 and don't enclose eight neighbors completely)"""
    critical_len_idx = np.argwhere(LenCritical == 1) # global index of nucleated grains larger then critical length Len0.
    GrowthClass = gg.Growth(ANGLE[critical_len_idx[:,0], critical_len_idx[:,1]], critical_len_idx[:,0], critical_len_idx[:,1], len(critical_len_idx[:,0]))
# SDX1, SDY1 are both global index. X1, Y1 are both local coordinates. 
    X1, Y1, SDX1, SDY1 = GrowthClass.growth1(TIPLEN1[critical_len_idx[:,0], critical_len_idx[:,1]])
    X2, Y2, SDX2, SDY2 = GrowthClass.growth2(TIPLEN2[critical_len_idx[:,0], critical_len_idx[:,1]])
    X3, Y3, SDX3, SDY3 = GrowthClass.growth3(TIPLEN3[critical_len_idx[:,0], critical_len_idx[:,1]])
    X4, Y4, SDX4, SDY4 = GrowthClass.growth4(TIPLEN4[critical_len_idx[:,0], critical_len_idx[:,1]])
    
    TIPLEN1[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN1[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX1, SDY1])
    TIPLEN2[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN2[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX2, SDY2])
    TIPLEN3[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN3[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX3, SDY3])
    TIPLEN4[critical_len_idx[:,0], critical_len_idx[:,1]] = GrowthClass.growth_polygon(TIPLEN4[critical_len_idx[:,0], critical_len_idx[:,1]], TempInterpolateT[SDX4, SDY4])
    X_MIN = np.minimum.reduce([SDX1, SDX2, SDX3, SDX4])
    Y_MIN = np.minimum.reduce([SDY1, SDY2, SDY3, SDY4])
    X_MAX = np.maximum.reduce([SDX1, SDX2, SDX3, SDX4])
    Y_MAX = np.maximum.reduce([SDY1, SDY2, SDY3, SDY4])
    for i in range(len(X_MIN)): # ith nucleated grain (including capturing).
        part_cells_class = CLASS[X_MIN[i]:X_MAX[i] + 1, Y_MIN[i]:Y_MAX[i] + 1] # a rectangle determined by polygon vertices, x_min, x_max, y_min, y_max.
        liquid_idx = np.argwhere(part_cells_class == -1) # index of all liquid cells within part_cells_class.
        Quadrilateral = path.Path([(X1[i], Y1[i]), (X2[i], Y2[i]), (X3[i], Y3[i]), (X4[i], Y4[i])])
        liquid_idx_global = liquid_idx + np.tile(np.dstack([X_MIN[i],Y_MIN[i]]), (len(liquid_idx), 1)) # liquid_idx_global is all global index of liquid cells within part_cells_class domain.
        liquid_coord_local = SDX * (liquid_idx_global - np.tile(critical_len_idx[i], (len(liquid_idx), 1))) # liquid_coord_local is local x and y coordinates of all liquid cells within part_cells_class domain.
        liquid_coord_tuple = totuple(liquid_coord_local)
        if liquid_coord_tuple[0] == ():
            continue
        EnclosePoint = Quadrilateral.contains_points(liquid_coord_tuple[0]) # boolean numpy . when it is true, it means it is enclosed by this grain polygon.
        enclose_idx = np.argwhere(EnclosePoint == True) # index of "True" value in  EnclosePoint. However, keep in mind that this index is not global index!!!
        print 'CLASS[{0}, {1}] is {2}'.format(critical_len_idx[i, 0], critical_len_idx[i, 1], CLASS[critical_len_idx[i, 0], critical_len_idx[i, 1]])
        CLASS[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = CLASS[critical_len_idx[i, 0], critical_len_idx[i, 1]] # this index is pretty important. liquid_idx_global is (1L, len(critical_len_idx), 2L) size. enclose_idx[:, 0] is used to select part of len(critical_len_idx). 
        NUCIDX[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = NUCIDX[critical_len_idx[i, 0], critical_len_idx[i, 1]]
        ANGLE[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = ANGLE[critical_len_idx[i, 0], critical_len_idx[i, 1]]
        SIDX[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = 2
        STOP[liquid_idx_global[:,enclose_idx[:, 0],0], liquid_idx_global[:,enclose_idx[:, 0],1]] = 1
        EnclosePoint.fill(False)

    LenCritical.fill(0)
    growth_times += 1
    print '{0}th growth'.format(growth_times)
    elapsed_time = timeit.default_timer() - start_time
    f = open('MonitorGrowth.txt','a')
    time = str(datetime.now()) + '   ' + '?th run nucleation: ' + str(growth_count) + '   ' + 'time to run this nucleation: ' + str(elapsed_time) + '   growth times: ' + str(growth_times) + '\n'
    f.writelines(time)
    f.close()
    return
    
def Output(_Row_Num,_Column_Num,i):
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
"""
This is test for temperature files. NOT necessary!!!
"""
#for i in xrange(170, TempFileNumber):
#    new = get_temp(RowNum, ColumnNum, i)[166, 50]
#    old = get_temp(RowNum, ColumnNum, i - 1)[166, 50]
#    if new >= old:
#        print "{0} is not valid.".format(i)
#        continue
#    print "{0} file is the starting file".format(i)
#print "***********pause!*****************"
def CenterNuc(_NUCINC, _Temp_Liquid, *args):
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
                _NUCINC += 1 # Record how many 
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
            for y in xrange(RowNum):
                for z in xrange(ColumnNum):
                    TempInterpolateT[y,z] = TempOld[y,z] + (x + 1)*1.0/MicroTime*(TempNew[y,z] - TempOld[y,z]) # interpolate temperature file along time dimension.
            if i - TempFileNumber_start > critical_CET: # control CET grain size
                Fv_gaussian_ini = 6e02 #assumption
                Fs_gaussian_ini = 2.4e04 #assumption.  PS is higher than PV.
                Fv_gaussian_nuc = 3e07 #assumption
                Fs_gaussian_nuc = 1e09 #assumption. PS is higher than PV.
                DTIME = 4e-10
                NUCINC = Nucleation(RowNum, ColumnNum, TempLiquid, TempSolid, UnderCoolMean, UnderCoolDelta, NUCINC, TempADJ, i)
                Growth(RowNum,ColumnNum,TempLiquid,TempSolid)
                Output(RowNum,ColumnNum,(i-1)*MicroTime+x)
                ClassExtract = dp.DataProcess('Extracting the data ' + str((i-1)*MicroTime+x), 'Class', (i-1)*MicroTime+x)
                ClassExtract.ExtractData()
                for y in xrange(RowNum):
                    for z in xrange(ColumnNum):
                        TempInterpolate_Old[y, z] = TempInterpolateT[y, z]
                continue
            
            NUCINC = Nucleation(RowNum, ColumnNum, TempLiquid, TempSolid, UnderCoolMean, UnderCoolDelta, NUCINC, TempADJ, i)
            print "The number of nucleated grains after nucleation is ", NUCINC
            Growth(RowNum,ColumnNum,TempLiquid,TempSolid)
            Output(RowNum,ColumnNum,(i-1)*MicroTime+x)
#==============================================================================
#             np.savetxt("TIPLEN_" + str(i) + ".txt", TIPLEN, fmt='%16.8e', delimiter=',')
#             np.savetxt("TIPLEN1_" + str(i) + ".txt", TIPLEN1, fmt = '%16.8e', delimiter=',')
#             np.savetxt("TIPLEN2_" + str(i) + ".txt", TIPLEN2, fmt = '%16.8e', delimiter=',')
#             np.savetxt("TIPLEN3_" + str(i) + ".txt", TIPLEN3, fmt = '%16.8e', delimiter=',')
#             np.savetxt("TIPLEN4_" + str(i) + ".txt", TIPLEN4, fmt = '%16.8e', delimiter=',')
#==============================================================================
            
            ClassExtract = dp.DataProcess('Extracting the data ' + str((i-1)*MicroTime+x), 'Class', (i-1)*MicroTime+x)
            ClassExtract.ExtractData()
            for y in xrange(RowNum):
                for z in xrange(ColumnNum):
                    TempInterpolate_Old[y, z] = TempInterpolateT[y, z]

#try:
#    import IPython
#    shell = IPython.get_ipython()
#    shell.enable_matplotlib(gui='qt')
#except:
#    pass
for k in xrange((TempFileNumber_start - 1) * MicroTime, (TempFileNumber - 1) * MicroTime):
    dd.ImageToPNG('Class', k, ColumnNum) # then www.gifmaker.me can convert pngs to gif animation.

"""show the animation, but it doesn't work using Python Console, works using IPython Console sometimes."""
#aa = dd.Ani((TempFileNumber_start - 1) * MicroTime, (TempFileNumber - 1) * MicroTime, ColumnNum)
#    """show the animation"""
#aa.ani() # run the code in the IPython Console. run "%matplotlib qt" at first, then run "%run C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer\main_program.py"
#    """show the images"""
##aa.ani_static()
            
print "LIDXN1", LIDXN1
print  "NUCIDX", NUCIDX
print "DTNUCL", DTNUCL
print "SIDX", SIDX
print "LenCritical", LenCritical
np.savetxt("DTNUCL.txt",DTNUCL,fmt='%i',delimiter=',')
np.savetxt("SIDX.txt",SIDX,fmt='%i',delimiter=',')
np.savetxt("NUCIDX.txt",NUCIDX,fmt='%i',delimiter=',')
np.savetxt("LIDXN1.txt",LIDXN1,fmt='%i',delimiter=',')
np.savetxt("LenCritical.txt",LenCritical,fmt='%i',delimiter=',')


    
#    TIMEFLG = 0
#    for j in xrange(RowNum):
#        for k in xrange(ColumnNum):
#            if DTNUCL[j,k] == 0: # the cell is solid originally, or it is liquid but not selected to nucleat.
#                continue
#            if TempNew[j,k] < TempLiquid + 70.0: #70.0 is questionable, don't allow temperature too high? 
#                TIMEFLG = 1
#    if TIMEFLG == 0: # at this time step, no nucleation occurred
#        continue
#    if i == 0: # The first temperature file
#        TempInterpolateT = TempNew
#    if i > 0:    
#        for x in xrange(MicroTime):
#            for y in xrange(RowNum):
#                for z in xrange(ColumnNum):
#                    TempInterpolateT[y,z] = TempOld[y,z] + x*1.0/MicroTime*(TempNew[y,z] - TempOld[y,z])
#            NUCINC = Nucleation(RowNum, ColumnNum, TempLiquid, TempSolid, UnderCoolMean, UnderCoolDelta, NUCINC)[0]     
##    print "NUCIDX", NUCIDX
#            Growth(RowNum,ColumnNum,TempLiquid,TempSolid)
#            print TIPLEN
#            Output(RowNum,ColumnNum,i*MicroTime+x)

#        print Nucleation(RowNum,ColumnNum,TempLiquid, TempSolid, UnderCoolMean, UnderCoolDelta, NUCINC)
#            Growth()
#        Output()
#            print ".............................."                    
#            print TempInterpolateT
#            print TempInterpolateT.shape
    
    
"""
Comment:
Gaussian distribution is not contained.
cell capture rules are not contained.
growth velocity is not based on physics, instead randomly given
nucleation probability is not based on physics, instead randomly given

the following cides are meaningless: 
#    TIMEFLG = 0
#    for j in xrange(RowNum):
#        for k in xrange(ColumnNum):
#            if DTNUCL[j,k] == 0: # the cell is solid originally, or it is liquid but not selected to nucleat.
#                continue
#            if TempNew[j,k] < TempLiquid + 70.0: #70.0 is questionable, don't allow temperature too high? 
#                TIMEFLG = 1
#    if TIMEFLG == 0: # at this time step, no nucleation occurred
#        continue

Vtip equation is questionable
Delta_TN, Delta, N_max are just uncertified assumption
DTIME = 4e-9 is questionable, it can be an adjustor
Fs_gaussian, Fv_gaussian can be adjusted. 
consider using pandas dataframe structure to replace ndarray in TIPLEN
instead of defining global variable, pack all the global variables in a seperate file.

"""  

""" Important Comment"""
#the locations that refer to starting temperature file:
#    def update_index(_Row_Num,_Column_Num,_TempADJ):, Temp = get_temp(_Row_Num,_Column_Num,TempFileNumber_start)
#    for i in xrange(400, TempFileNumber): # main program
#    In "DisplayData.py" file, class Ani(): in def showclass(self, k): method, for i in xrange(self._ClassFileNumber):
#   In "DisplayData.py" file, class Ani(): in def showclass(self, k): method, classimagelist = [t.GetClass(30, self._Column_Num, i) for i in xrange(self.start, self._ClassFileNumber)] 
#   In "DisplayData.py" file, class Ani(): in def ani() method, ani = animation.FuncAnimation(self.fig, self.showclass, frames=xrange(self.start, self._ClassFileNumber), interval=50, blit=True)

# 01/20/2017 changes:
    # In Initialization(), line 173, 
    # if _Temp_Liquid - Temp[i,j] > -10.0:# Temp is not too higher than liquidus temperature, -10.0 is unnecessary. replace -10.0 with 0.0
    
    # In Nucleation(), line 234, if _Temp_Liquid - TempInterpolateT[i,j] < 0.0 or _Temp_Liquid - TempInterpolate_Old[i,j] < 0.0:
    # delete _Temp_Liquid - TempInterpolate_Old[i,j] < 0.0!!!

    # Nucleation(), line230                
    # if SIDX[i, j] != 1 and SIDX[i, j] != 2: # not nucleat and grow before
        # NUCIDX[i, j] = 0
    # this is not correct. NUCIDX[i, j] can also be -1(solid beginning)