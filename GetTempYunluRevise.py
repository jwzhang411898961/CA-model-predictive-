# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 10:03:41 2017

@author: jnzzp5
"""

"""
This is for Yunlu's interpolation on linux system.
"""

import numpy as np
from Parameters import *
from io import StringIO
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import pandas as pad
from scipy import interpolate
from scipy.interpolate import interpn
import matplotlib.pyplot as plt
#from mayavi.mlab import *
import timeit
import os

"""....................create reading, writing and interpolating temperature class........................"""
class GetTemp:
    def __init__(self, _Row_Num, _Column_Num, _TempFileNumber):
        self.RowNum = _Row_Num
        self.ColumnNum = _Column_Num
        self.TempFileNumber = _TempFileNumber
        
#    def NodeTemp(self, TempFile_start): # One node temperature; pandas.series is called.
#        TempList = list()
#        IndexList = list()
#        for i in range(TempFile_start, self.TempFileNumber):
#            TempND = np.loadtxt(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_correct2014error_02142017\Nodeset-SET-FACE2_1layerTempTi64_Interpolate"+str(i),delimiter = ',')[155, 50]
#            TempList.append(TempND)
#            IndexList.append(str(i))
#        TempSer = pad.Series(TempList, index = IndexList)
#        TempSer.plot(kind = 'line', title = '[{0},{1}] node temperature history'.format(61, 160))
#        print("TempSer is:\n", TempSer)
#        print("TempSer shape is ", TempSer.shape)
#    
#    def NodesTemp(self): #Several nodes temprature. pandas.dataframe is called.
#        TempList = list()
#        IndexList = list()
#        TempDict = dict()
#        for i in range(self.TempFileNumber):
#            TempND = list(np.loadtxt(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate" + str(i),delimiter = ',')[141:171, 50]) #???????is list necessary?
#            TempList.append(TempND)
#            IndexList.append(str(i))
#            TempArray = np.array(TempList) #?????????
##        print TempArray.shape
##        print type(TempList)
##        names = ['160','161','162','163','164']
#        names = range(141, 171) # The columns
#        for j, name in enumerate(names):
#            TempDict[name] =  TempArray[:, j]
#        Tempdf = pad.DataFrame(TempDict, columns=names)
#        Tempdf.to_csv(path_or_buf = r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer\NodesTempHistory')
#        writer = pad.ExcelWriter(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer\NodesTempHistory.xlsx', engine='xlsxwriter') # Create a Pandas Excel writer using XlsxWriter as the engine.
#        Tempdf.to_excel(writer, sheet_name='Sheet1') # Convert the dataframe to an XlsxWriter Excel object.
#        writer.save() # Close the Pandas Excel writer and output the Excel file.
#        print(Tempdf)
#
##        TempSer = pad.Series(TempList, index = IndexList)
##        TempSer.plot(kind = 'line', title = '[{0},{1}] node temperature history'.format(166, 50))
##        print "TempSer is:\n", TempSer
##        print "TempSer shape is ", TempSer.shape
#
#    def NodesTempMultiLayer(self): #Several nodes temprature. pandas.dataframe is called. 02092017. 03012017 revise.
#        TempList = list()
#        IndexList = list()
#        TempDict = dict()
#        for i in range(self.TempFileNumber):
#            TempND = list(np.loadtxt(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\Nodeset-SET-CS_CLAD9_25layerTempTi64_Interpolate_finer" + str(i),delimiter = ',')[0:233:8, 300]) #???????is list necessary? In raw data, there are 233 rows and 601 columns. Cell size is 10um.
#            TempList.append(TempND)
#            IndexList.append(str(i))
#            TempArray = np.array(TempList) #????????? stores temperature.
#            print i
##        print TempArray.shape
##        s = raw_input('-->')
##        print type(TempList)
##        names = ['160','161','162','163','164']
#        names = range(0, 233, 8) # The columns
#        for j, name in enumerate(names):
#            TempDict[name] =  TempArray[:, j]
#        Tempdf = pad.DataFrame(TempDict, columns=names)
#        Tempdf.to_csv(path_or_buf = r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer\25LayerNodesTempHistory')
#        writer = pad.ExcelWriter(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer\25LayerNodesTempHistory.xlsx', engine='xlsxwriter') # Create a Pandas Excel writer using XlsxWriter as the engine.
#        Tempdf.to_excel(writer, sheet_name='Sheet1') # Convert the dataframe to an XlsxWriter Excel object.
#        
#
#
#        Tempdf.plot(kind = 'line')
#        writer.save() # Close the Pandas Excel writer and output the Excel file.
#        print Tempdf
        
    def get_temp_yanlei(self): # 06132017 created. get subset node temperature from yanlei's temperature files.
        selector = 2 # determine which part temperature will be extracted and sorted.
        datalist = list()
        ExtractTemp = np.array(datalist) # create an empty array.
        
        """select the temperature domain based on the selector value."""
        if selector == 0: # only deposit part
            for i in range(self.TempFileNumber):
                temp_all = np.genfromtxt(r"/home/yunlu/YunluInterpolation/" + str(i + 1) + ".txt")
                idx = np.argwhere(temp_all[:, 2] >= 0.006) # all row number of deposit location.
                temp_deposit = temp_all[idx[:, 0], :] # extract deposit part node temperature
                sort_deposit_idx = np.lexsort((temp_deposit[:, 0], temp_deposit[:, 1], temp_deposit[:, 2]))
                temp_middle_sort = temp_deposit[sort_deposit_idx] # temperature after sorting
                np.savetxt(r"/home/yunlu/YunluInterpolation/deposit" + str(i + 1) + ".txt", temp_middle_sort, fmt = '%1.6f', delimiter = ',')
                
        if selector == 1: # deposit and underneath substrate part for three layers
            for i in range(self.TempFileNumber):
                temp_all = np.genfromtxt(r"/home/yunlu/YunluInterpolation/" + str(i + 1) + ".txt")
                middle_idx = np.argwhere((temp_all[:, 0] >= 0.005) & (temp_all[:, 0] <= 0.02) & (temp_all[:, 1] >= 0.0055) & (temp_all[:, 1] <= 0.0075))
                temp_middle = temp_all[middle_idx[:, 0], :] # extract middle part node temperature
                sort_middle_idx = np.lexsort((temp_middle[:, 0], temp_middle[:, 1], temp_middle[:, 2]))
                temp_middle_sort = temp_middle[sort_middle_idx]
                np.savetxt(r"/home/yunlu/YunluInterpolation/middle" + str(i + 1) + ".txt", temp_middle_sort, fmt = '%1.6f', delimiter = ',')
                ExtractTemp = temp_middle_sort[:, 5]
                print(ExtractTemp.shape)
                print(temp_middle_sort.shape)
                
        if selector == 2: # deposit and underneath substrate part for 15 layer result.
            for i in range(self.TempFileNumber):
                temp_all = np.genfromtxt(r"/home/yunlu/YunluInterpolation/" + str(i + 1) + ".txt")
                middle_idx = np.argwhere((temp_all[:, 0] >= 0.005) & (temp_all[:, 0] <= 0.02) & (temp_all[:, 1] >= 0.0055) & (temp_all[:, 1] <= 0.0075))
                temp_middle = temp_all[middle_idx[:, 0], :] # extract middle part node temperature
                sort_middle_idx = np.lexsort((temp_middle[:, 0], temp_middle[:, 1], temp_middle[:, 2]))
                temp_middle_sort = temp_middle[sort_middle_idx]
                file_name = r"/home/yunlu/YunluInterpolation/middle" + str(i + 1) + ".txt"
                print(i)
                if not os.path.exists(file_name):
                    np.savetxt(file_name, temp_middle_sort, fmt = '%1.6f', delimiter = ',')
                ExtractTemp = np.append(ExtractTemp, temp_middle_sort[:, 5])
            ExtractTemp = ExtractTemp.reshape(self.TempFileNumber, -1)
            print("ExtractTemp shape is {0}".format(np.shape(ExtractTemp)))
                
        """interpolate 3D temperature"""
        z, y, x = np.arange(1, 43), np.arange(1, 11), np.arange(1, 70) #original index
        znew, ynew, xnew = np.arange(1, 42.2, 0.2), np.arange(1, 10.2, 0.2), np.arange(1, 69.2, 0.2) #interpolated index
        points = (z, y, x) # three 1D grid along 3 direction.       
        inter_mesh1 = np.asarray(np.meshgrid(znew, ynew, xnew)) # generate 3D grid. keep consistent with points sequence.
        interp_points1 = tuple(map(tuple, inter_mesh1))
        print("interp_points1 shape is {0}".format(np.shape(interp_points1)))
        points2 = (znew, ynew, xnew)
        znew2, ynew2, xnew2 = np.arange(1.05, 42, 0.05), np.arange(1.05, 10, 0.05), np.arange(1.05, 69, 0.05)
        inter_mesh2 = np.asarray(np.meshgrid(znew2, ynew2, xnew2))
        interp_points2 = tuple(map(tuple, inter_mesh2))
        print("interp_points2 shape is {0}".format(np.shape(interp_points2)))
        for i in range(self.TempFileNumber):# i represents the time step
            ExtractTempArr=ExtractTemp[i, :].reshape(len(z),len(y),len(x)) # shape is (42, 10, 69)!!! This ordering is very important!
            print("ExtractTempArr shape is {0}".format(np.shape(ExtractTempArr)))
            
            # use mayavi to show 3D temperature.
#            self.DisplayTemp3D(ExtractTempArr)
            
            plt.imshow(ExtractTempArr[:, 5, :], cmap='hot', interpolation='nearest') # show 2D temperature before interpolation. "5" is chosen randomly. ("5" is middle plane.)
            plt.savefig(r'/home/yunlu/YunluInterpolation/Temp_15layers_NoInter_2D' + str(i) + '.png')
            start_time1 = timeit.default_timer()
            InterpTemp1 = interpn(points, ExtractTempArr, interp_points2, method = "nearest").reshape(len(znew2), len(ynew2), len(xnew2)) # core method to return interpolated points. points are sampled(original) grids. ExtractTempArr are the values of sampled grid. xi is new interpolated grid. returns are the values of interpolated grids.
            elapsed1 = timeit.default_timer() - start_time1
            print("InterpTemp1 execution time is {0}".format(elapsed1))
#==============================================================================
#             start_time2 = timeit.default_timer()
#             InterpTemp2 = interpn(points2, InterpTemp1, interp_points2, method = "nearest")
#             elapsed2 = timeit.default_timer() - start_time2
#             print "InterpTemp1 execution time is {0}, \nInterpTemp2 execution time is {1}".format(elapsed1, elapsed2)
#==============================================================================
            
            """save 3D and 2D interpolated results to files"""
            with open(r"/home/yunlu/YunluInterpolation/Interpolate"+str(i), 'w') as outfile, open(r"/home/yunlu/YunluInterpolation/Interpolate2D_XZPlane"+str(i), 'w') as outfile2d: # open multiple files using "with...as..." structure
                outfile.write('# Array shape: {0}\n'.format(InterpTemp1.shape))
                for index, data_slice in enumerate(InterpTemp1):
                    
                    """3D interpolated storage. in order to save computation time, it is commented here."""
#==============================================================================
#                     outfile.write('#{0} New slice\n'.format(index))
#                     np.savetxt(outfile, data_slice, fmt='%1.3e', delimiter=',')
#==============================================================================

                    if index == int(round(len(ynew2) / 2)):
                        print("index is {0}".format(index))
                        start_time3 = timeit.default_timer()
                        np.savetxt(outfile2d, data_slice, fmt = '%1.3e', delimiter = ',')
                        elapsed3 = timeit.default_timer() - start_time3
                        print("np.savetxt execution time is {0}".format(elapsed3))
                        """display the interpolated temperature"""
                        plt.imshow(data_slice, cmap='hot', interpolation='nearest')
                        plt.savefig(r'/home/yunlu/YunluInterpolation/Temp_15layers_Interpolate2D_XZPlane' + str(i) + '.png')
        print("done!")
        
    def DisplayTemp3D(self, temp): # show 3d result.
        obj = contour3d(temp, contours=4, transparent=True)
        return obj

#    def get_temp(self):
#        Temp=np.zeros((self.RowNum,self.ColumnNum),dtype='float16')# this is created to store temperature data
#        try:
#            f=file("C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate"+str(self.TempFileNumber),'r')
#        except IOError:
#            print "This temperature file cannot be found"
#        else:
#            j=0
#            for line in f.readlines():
#                line = line.replace("\n","") #delete '\n' in the line string
#                line = line.split(',') #split string data with comma
#                for i in xrange(self.ColumnNum):
#                    Temp[j,i] = float(line[i]) # convert string data to float data, delimiter is comma; store every data to ndarray
#                j+=1
#            f.close()
#        return Temp

"""....................operation........................"""
c = 3 # select operation
if c == 0: # 25 multiple layers
    TempClass = GetTemp(233, 601, 200) # RowNum is 233, ColumnNum is 601, Time step number is 827.
    TempClass.NodesTempMultiLayer()
if c == 1: # Onelayer one node temperature history
    TempClass = GetTemp(171, 101, 519) # RowNum is 233, ColumnNum is 601, Time step number is 827.
    TempClass.NodeTemp(400)
if c == 2:
    TempClass = GetTemp(171, 101, 2) # 180 temperature files for 3 layers 06152017 revised. 
    TempClass.get_temp_yanlei()
if c == 3: # leiyan 15layers result
    TempClass = GetTemp(206, 341, 959) # 960 temperature files for 15 layers 07172017 revised. 
    TempClass.get_temp_yanlei()
