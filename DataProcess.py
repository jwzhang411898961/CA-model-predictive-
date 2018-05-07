# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 16:13:07 2016

@author: jnzzp5
"""

import sys
import numpy as np
import random
import math
import os
import matplotlib.pyplot as plt
from Parameters import *

class DataProcess:# Extract certain lines CLASS data.
    FileCount = 0
    def __init__(self, Description, Filename, Sequence):
        self.Description = Description
        self.Filename = Filename
        self.Sequence = Sequence
        print self.Description
        DataProcess.FileCount += 1
    
    def ExtractData(self):
        with open(self.Filename + str(self.Sequence) + '.txt') as handle1: # only read, not write.
            with open('PostProcess' + self.Filename + str(self.Sequence) + '.txt','w') as handle2: # extract the rows that close to melt pool.
                lines = handle1.readlines()
                handle2.writelines("Microstructure\n")
                handle2.writelines("Variables= 'I' 'J' 'CLASS'\n")
                handle2.writelines("ZONE I=171 J=101 F=POINT\n")
                for line in lines:
                    if line.startswith(("M", "V", "Z")):
                        continue
                    words = line.strip().split()
                    if int(words[0]) > Invalid_rows: # greater than Invalid_rows row. close to meltpool.
                        handle2.write(line)
                        
    def ExtractDataMulti(self): # multiple layer. 02102017
        os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer')
        with open(self.Filename + str(self.Sequence) + '.txt') as handle1: # only read, not write.
            with open('PostProcess' + self.Filename + str(self.Sequence) + '.txt','w') as handle2: # extract the rows that close to melt pool.
                lines = handle1.readlines()
                handle2.writelines("Microstructure\n")
                handle2.writelines("Variables= 'I' 'J' 'CLASS'\n")
                handle2.writelines("ZONE I=171 J=101 F=POINT\n")
                for line in lines:
                    if line.startswith(("M", "V", "Z")):
                        continue
                    words = line.strip().split()
                    if int(words[0]) > Invalid_rows: # greater than Invalid_rows row. close to meltpool.
                        handle2.write(line)
                        
class Remelt:
    def __init__(self, CLASS):
        self.CLASS = CLASS
    
    def remelt(self, _TempADJ):
        for i in xrange(_Row_Num):
            for j in xrange(_Column_Num):
                if Temp[i, j] > _TempADJ:
                    self.CLASS[i, j] = -1

    def epitaxial(self, sidx):
        indices = np.transpose(np.nonzero(self.CLASS == -1)) # find the indices of all CLASS cell with -1.
        for i, j in indices:
            if self.CLASS[i + 1, j] > 0 or self.CLASS[i - 1, j] > 0 or self.CLASS[i, j + 1] > 0 or self.CLASS[i, j - 1] > 0 or self.CLASS[i + 1, j + 1] > 0 or self.CLASS[i + 1, j - 1] > 0 or self.CLASS[i - 1, j - 1] > 0 or self.CLASS[i - 1, j + 1] > 0:
                sidx[i, j] == 3 # 3 represents this cell is at the S/L interface.
                a = [self.CLASS[i + 1, j], self.CLASS[i - 1, j], self.CLASS[i, j + 1], self.CLASS[i, j - 1], self.CLASS[i + 1, j + 1], self.CLASS[i + 1, j - 1], self.CLASS[i - 1, j - 1], self.CLASS[i - 1, j + 1]]
                self.CLASS[i, j] = a[min(range(len(a)), key=lambda i: abs(a[i] - Class_mean))] # prescribe the CLASS of interface cell when nucleation. This is to look for which neighbouring cell is solid.
                ANGLE[i,j] = (CLASS[i,j]*90/48.0 - 45) * PI / 180.0
        return sidx, CLASS

class Substrate:
    count, num = 0, 0 # count is to record how many cells have been defined. num is the number of iteration.
    def __init__(self, _Row_Num, _Column_Num, grid_x, grid_y):
        self._Row_Num =  _Row_Num
        self._Column_Num = _Column_Num
        self.grid_x = grid_x # grid_x defines the seed density along one direction
        self.grid_y = grid_y # grid_y defines the seed density along the other direction
    
    def air(self):
        pass
        
    def substrate(self, r1, r2): # r1 and r2 represent the probability of substrate grain growth along two directions.
        substrate_row = int(self._Row_Num / substrate_ratio) # how many substrate rows   
        air_row = self._Row_Num - substrate_row # how many air rows
        base = np.zeros((self._Row_Num,self._Column_Num), dtype = 'int32') # base[] represents the CLASS(grain orientation) for each cell.
        indicator = np.zeros((self._Row_Num,self._Column_Num), dtype = 'int32') # indicator[] represents whether the substrate cell is defined. when it is 1, the cell is defined. when it is 0, the cell has not been defined.
        stop = np.zeros((self._Row_Num,self._Column_Num), dtype=np.int) 
        for i in xrange(substrate_row, self._Row_Num):
            for j in xrange(0, self._Column_Num):
                base[i, j] = -1
        for i in xrange(0, substrate_row, self.grid_x):
            for j in xrange(0, self._Column_Num, self.grid_y):
                base[i, j] = random.randint(1,49)
                self.count += 1
                indicator[i, j] = 1
        while True:
            self.num += 1
            stop.fill(0)
            for i in xrange(1, substrate_row - 1):
                for j in xrange(1, self._Column_Num - 1):
                    if random.random() < r1 and indicator[i, j] == 1 and stop[i, j] == 0:
                        if indicator[i + 1, j] == 0:
                            base[i + 1, j] = base[i, j]
                            indicator[i + 1, j] = 1
                            self.count += 1
                            stop[i + 1, j] = -1
                        if indicator[i + 1, j + 1] == 0:
                            base[i + 1, j + 1] = base[i, j]
                            indicator[i + 1, j + 1] = 1
                            self.count += 1
                            stop[i + 1, j + 1] = -1
                        if indicator[i, j + 1] == 0:
                            base[i, j + 1] = base[i, j]
                            indicator[i, j + 1] = 1
                            self.count += 1
                            stop[i, j + 1] = -1
                        if indicator[i - 1, j + 1] == 0:
                            base[i - 1, j + 1] = base[i, j]
                            indicator[i - 1, j + 1] = 1
                            self.count += 1
                            stop[i - 1, j + 1] = -1
                    if random.random() < r2 and indicator[i, j] == 1 and stop[i, j] == 0:
                        if indicator[i - 1, j] == 0:
                            base[i - 1, j] = base[i, j]
                            indicator[i - 1, j] = 1
                            self.count += 1
                            stop[i - 1, j] = -1
                        if indicator[i - 1, j - 1] == 0:
                            base[i - 1, j - 1] = base[i, j]
                            indicator[i - 1, j - 1] = 1
                            self.count += 1
                            stop[i - 1, j - 1] = -1
                        if indicator[i, j - 1] == 0:
                            base[i, j - 1] = base[i, j]
                            indicator[i, j - 1] = 1
                            self.count += 1
                            stop[i, j - 1] = -1
                        if indicator[i + 1, j - 1] == 0:
                            base[i + 1, j - 1] = base[i, j]
                            indicator[i + 1, j - 1] = 1
                            self.count += 1
                            stop[i + 1, j - 1] = -1
#                    if random.random() < r1 and indicator[i + 1, j] == 0 and indicator[i, j] == 1:
#                        base[i + 1, j] = base[i, j]
#                        indicator[i + 1, j] = 1
#                        self.count += 1
#                    if random.random() < r3 and indicator[i - 1, j] == 0 and indicator[i, j] == 1:
#                        base[i - 1, j] = base[i, j]
#                        indicator[i - 1, j] = 1
#                        self.count += 1
#                    if random.random() < r2 and indicator[i, j] == 1 and indicator[i, j + 1] == 0:
#                        base[i, j + 1] = base[i, j]
#                        indicator[i, j + 1] = 1
#                        self.count += 1
#                    if random.random() < r4 and indicator[i, j] == 1 and indicator[i, j - 1] == 0:
#                        base[i, j - 1] = base[i, j]
#                        indicator[i, j - 1] = 1
#                        self.count += 1
#                    if random.random() < r1 and indicator[i, j] == 1 and indicator[i + 1, j + 1] == 0:
#                        base[i + 1, j + 1] = base[i, j]
#                        indicator[i + 1, j + 1] = 1
#                        self.count += 1
#                    if random.random() < r2 and indicator[i, j] == 1 and indicator[i - 1, j + 1] == 0:
#                        base[i - 1, j + 1] = base[i, j]
#                        indicator[i - 1, j + 1] = 1
#                        self.count += 1
#                    if random.random() < r3 and indicator[i, j] == 1 and indicator[i - 1, j - 1] == 0:
#                        base[i - 1, j - 1] = base[i, j]
#                        indicator[i - 1, j - 1] = 1
#                        self.count += 1
#                    if random.random() < r4 and indicator[i, j] == 1 and indicator[i + 1, j - 1] == 0:
#                        base[i + 1, j - 1] = base[i, j]
#                        indicator[i + 1, j - 1] = 1
#                        self.count += 1
            print "the number of iteration is {0}".format(self.num)
            print "the count is {0}".format(self.count)
            if self.count == substrate_row * self._Column_Num:
                break
        #print base, indicator
        return base
#    def air(self, Air_Row, Air_Column):
#        air = np.zeros((Air_Row,Air_Column), dtype = 'int32')
#        air.fill(-1)
        
    def showimage(self, r1, r2):
        base = self.substrate(r1, r2)
        plt.imshow(base, cmap='hot', interpolation='nearest')
        plt.colorbar()
        plt.show()
        return base
"""operation"""
if __name__ == '__main__':
    test = Substrate(171, 101, 6, 6)
    a = test.showimage(0.2, 0.25)
    print a
                    
#==============================================================================
# ClassFileNumber = 3
# for i in xrange(ClassFileNumber):
#     t = DataProcess('Postprocessing the data ' + str(i), 'Class', i)
#     t.method()
#==============================================================================
