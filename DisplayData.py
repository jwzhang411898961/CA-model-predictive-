# -*- coding: utf-8 -*-
"""
Created on Wed Nov 02 10:09:18 2016

@author: jnzzp5
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import os
from Parameters import *
from mayavi.mlab import *

"""Show temparature files not consecutively. 11022016"""
class ShowTemp:
    def __init__(self, _Row_Num, _Column_Num):
        self._Row_Num = _Row_Num
        self._Column_Num = _Column_Num
        self.Temp = np.zeros((self._Row_Num, self._Column_Num), dtype='float16') # this is created to store temperature data
        #RowNum,ColumnNum = 171, 101
        #Temp=np.zeros((RowNum, ColumnNum), dtype='float16')# this is created to store temperature data
        
    def get_temp(self, _TempFileNumber):
        try:
            f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate" + str(_TempFileNumber),'r')
            #f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Sort_extract"+str(_TempFileNumber),'r')
        except IOError:
            print "This temperature file cannot be found"
        else:
            j=0
            for line in f.readlines():
                line = line.replace("\n","") #delete '\n' in the line string
                line = line.split(',') #split string data with comma
                for i in xrange(self._Column_Num):
                    self.Temp[j,i]=float(line[i]) # convert string data to float data, delimiter is comma; store every data to ndarray
                j+=1
            f.close()
        return
        
    def show_temp(self): # show temperature result not in the animation.
        for i in xrange(10):        
            self.get_temp(i)
            plt.imshow(self.Temp, cmap='hot', interpolation='nearest')
            #plt.colorbar()
            plt.show()
"""Operation"""
#==============================================================================
# if __name__ == '__main__':
#     p = ShowTemp(171, 101)
#     p.show_temp()
# 
#==============================================================================
"""***********************************************************************"""


"""Show the temperature files in animation consecutively. 11022016/ 02132017"""
class ShowTempAnimate:
    def __init__(self, _Row_Num, _Column_Num, _TempFileNumber, _StartNum):
        self._StartNum = _StartNum
        self._Row_Num = _Row_Num
        self._Column_Num = _Column_Num
        self.fig = plt.figure()
        self._TempFileNumber = _TempFileNumber
        #from Parameters import *
        #RowNum,ColumnNum = 171, 101
    def get_temp(self, FileNumber):
        Temp=np.zeros((self._Row_Num, self._Column_Num), dtype='float16')# this is created to store temperature data
        try:
            i = 7 # control which temperature file to display.
            if i == 0:
                f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate"+str(self._StartNum + FileNumber),'r') # 400 may change according to starting temperature file.
            if i == 1: # SET-FACE2
                f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_02132017\Nodeset-SET-FACE2_1layerTempTi64_gau_Interpolate"+str(self._StartNum + FileNumber),'r') # FACE2! 50 may change according to starting temperature file.
            if i == 2: # SET-FACE3
                f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_02132017\Nodeset-SET-FACE3_1layerTempTi64_gau_Interpolate"+str(self._StartNum + FileNumber),'r') # FACE3! 50 may change according to starting temperature file.
            if i == 3: # SET-FACE2. 2014 original temperature field odb.
                f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_correct2014error_02142017\Nodeset-SET-FACE2_1layerTempTi64_Interpolate"+str(self._StartNum + FileNumber),'r')
            if i == 4: # SET-FACE2. 400W result for Final Defense.
                f = file(r"C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\Nodeset-SET-FACE2_1layerTempTi64_gau_400W_Interpolate_finer" + str(self._StartNum + FileNumber),'r')
            if i == 5: # SET-FACE3. 400W result for Final Defense.
                f = file(r"C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\Nodeset-SET-FACE3_1layerTempTi64_gau_400W_Interpolate_finer" + str(self._StartNum + FileNumber),'r')
            if i == 6: # SET-FACE3. 1000W result for Final Defense.
                f = file(r"C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\Nodeset-SET-FACE3_1layerTempTi64_gau_1000W_Interpolate_finer" + str(self._StartNum + FileNumber),'r')
            if i == 7: # SET-FACE2. 1000W result for Final Defense.
                f = file(r"C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\Nodeset-SET-FACE2_1layerTempTi64_gau_1000W_Interpolate_finer" + str(self._StartNum + FileNumber),'r')
            #f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Sort_extract"+str(_TempFileNumber),'r')
        except IOError:
            print "This temperature file cannot be found"
        #f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate"+str(400 + FileNumber),'r')
        else:
            j=0
            for line in f.readlines():
                line = line.replace("\n","") #delete '\n' in the line string
                line = line.split(',') #split string data with comma
                for i in xrange(self._Column_Num):
                    Temp[j,i]=float(line[i]) # convert string data to float data, delimiter is comma; store every data to ndarray
                j+=1
            f.close()
        return Temp
        
    def show_temp_animate(self, k):
        imagelist = [self.get_temp(i) for i in xrange(self._TempFileNumber)]
        #print imagelist
        im = plt.imshow(imagelist[0], cmap='hot', interpolation='nearest')
        im.set_array(imagelist[k])
        return [im]

    def ImageToPNG(self):
        for i in xrange(self._TempFileNumber):
            image_temp = self.get_temp(i)
            plt.imshow(image_temp, interpolation='bilinear')
            #plt.legend()
            pos = 5 # select position
            if pos == 1:
                plt.savefig(r'C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\ti64_1layer_temp_correct2014error_02142017\OneLayerTemp_SET-FACE2_' + str(i + self._StartNum) + '.png')
            if pos == 2: # SET-FACE 2, 400W
                plt.savefig(r'C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\OneLayerTemp_SET-FACE2_gau_400W' + str(i + self._StartNum) + '.png')
            if pos == 3: # SET-FACE 3, 400W
                plt.savefig(r'C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\OneLayerTemp_SET-FACE3_gau_400W' + str(i + self._StartNum) + '.png')
            if pos == 4: # SET-FACE 3, 1000W
                plt.savefig(r'C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\OneLayerTemp_SET-FACE3_gau_1000W' + str(i + self._StartNum) + '.png')
            if pos == 5: # SET-FACE 2, 1000W
                plt.savefig(r'C:\Program Files (x86)\SIMULIA\Abaqus614\Temp\FinalDefense\OneLayerTemp_SET-FACE2_gau_1000W' + str(i + self._StartNum) + '.png')
#                raw_input("_____________")
                
    def ani(self):
        ani = animation.FuncAnimation(self.fig, self.show_temp_animate, frames=xrange(self._TempFileNumber), interval=50, blit=True)
        plt.show()
        
"""Operation"""
#==============================================================================
# if __name__ == '__main__':
#     #q = ShowTempAnimate(171, 101, 50) # 161 is for SET-FACE3, 171 is for SET-FACE2.
#     #q.ani()
#     j = 5
#     if j == 1: # 750W temperature result
#         q = ShowTempAnimate(171, 101, 118, 400) # 161 is for SET-FACE3, 171 is for SET-FACE2.
#         q.ImageToPNG()
#     if j == 2: # 400W SET-FACE2 temperature result. 11132017
#         q = ShowTempAnimate(681, 401, 107, 0) # for finer temperature file.
#         q.ImageToPNG()
#     if j == 3: # 400W SET-FACE3 temperature result. 11132017
#         q = ShowTempAnimate(641, 401, 100, 0) # for finer temperature file.
#         q.ImageToPNG()
#     if j == 4: # 1000W SET-FACE3 temperature result. 11132017
#         q = ShowTempAnimate(641, 401, 384, 0) # for finer temperature file.
#         q.ImageToPNG()
#     if j == 5: # 1000W SET-FACE2 temperature result. 11142017
#         q = ShowTempAnimate(681, 401, 384, 0) # for finer temperature file.
#         q.ImageToPNG()
#==============================================================================
"""***********************************************************************"""



"""Show the non-interpolated temperature files in animation consecutively. 02132017"""
class ShowTempAnimate2:
    def __init__(self, _Row_Num, _Column_Num, _TempFileNumber):
        self._Row_Num = _Row_Num
        self._Column_Num = _Column_Num
        self.fig = plt.figure()
        self._TempFileNumber = _TempFileNumber
        #from Parameters import *
        #RowNum,ColumnNum = 171, 101
    def get_temp(self, FileNumber):
        Temp=np.zeros((self._Row_Num, self._Column_Num), dtype='float16')# this is created to store temperature data
        try:
            f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Sort_extract"+str(400 + FileNumber),'r') # 400 may change according to starting temperature file.
            #f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Sort_extract"+str(_TempFileNumber),'r')
        except IOError:
            print "This temperature file cannot be found"
        #f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate"+str(400 + FileNumber),'r')
        else:
            j=0
            for line in f.readlines():
                line = line.replace("\n","") #delete '\n' in the line string
                line = line.split(',') #split string data with comma
                for i in xrange(self._Column_Num):
                    Temp[j,i]=float(line[i]) # convert string data to float data, delimiter is comma; store every data to ndarray
                j+=1
            f.close()
        return Temp
    def show_temp_animate(self, k):
        imagelist = [self.get_temp(i) for i in xrange(self._TempFileNumber)]
        #print imagelist
        im = plt.imshow(imagelist[0], cmap='hot', interpolation='nearest')
        im.set_array(imagelist[k])
        return [im]

    def ani(self):
        ani = animation.FuncAnimation(self.fig, self.show_temp_animate, frames=xrange(self._TempFileNumber), interval=50, blit=True)
        plt.show()
"""Operation"""
#==============================================================================
# if __name__ == '__main__':
#     q = ShowTempAnimate2(18, 11, 50)
#     q.ani()
#==============================================================================
    

"""***********************************************************************"""



"""Show the multiple-layer temperature files in animation consecutively. 02092017"""
class ShowTempAnimateMulti:
    def __init__(self, _Row_Num, _Column_Num, _TempFileNumber):
        self._Row_Num = _Row_Num
        self._Column_Num = _Column_Num
        self.fig = plt.figure()
        self._TempFileNumber = _TempFileNumber
        #from Parameters import *
        #RowNum,ColumnNum = 171, 101
        
    def get_temp(self, FileNumber):
        Temp=np.zeros((self._Row_Num, self._Column_Num), dtype='float16')# this is created to store temperature data
        start_file = 0
        try:
            f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\Nodeset-SET-CS_CLAD9_25layerTempTi64_Interpolate_finer"+str(start_file + FileNumber),'r') # 400 may change according to starting temperature file.
            #f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Sort_extract"+str(_TempFileNumber),'r')
        except IOError:
            print "This temperature file cannot be found"
        #f=file(r"C:\Program Files (x86)\SIMULIA\Abaqus\Temp\JingweiZhang\PostprocessData\NodesetSET-FACE2OnelayerTempTi6Al4VEff04_Interpolate"+str(400 + FileNumber),'r')
        else:
            j=0
            for line in f.readlines():
                line = line.replace("\n","") #delete '\n' in the line string
                line = line.split(',') #split string data with comma
                for i in xrange(self._Column_Num):
                    Temp[j,i]=float(line[i]) # convert string data to float data, delimiter is comma; store every data to ndarray
                j+=1
            f.close()
        return Temp
        
    def ImageToPNGMulti(self): # show temperature PNGs for multiple layers. 02092017
        for i in xrange(self._TempFileNumber):
            image_temp = self.get_temp(i)
            plt.imshow(image_temp, interpolation='nearest')
            plt.savefig(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer\MultiTemp' + str(i) + '.png')
            
    def show_temp_animate(self, k):
        imagelist = [self.get_temp(i) for i in xrange(self._TempFileNumber)]
        #print imagelist
        im = plt.imshow(imagelist[0], cmap='hot', interpolation='nearest')
        im.set_array(imagelist[k])
        return [im]

    def ani(self):
        ani = animation.FuncAnimation(self.fig, self.show_temp_animate, frames=xrange(self._TempFileNumber), interval=50, blit=True)
        plt.show()
"""Operation"""
#==============================================================================
# if __name__ == '__main__':
#     q = ShowTempAnimateMulti(233, 601, 300)
#     #q.ani()
#     q.ImageToPNGMulti()
#==============================================================================
    


"""***********************************************************************"""

"""Show the class files in animation consecutively. 11032016"""
#x, y = np.linspace(141, 171, 31), np.linspace(1, 101, 101)
#X, Y = np.meshgrid(x, y)
class ShowClass:
    FileCount = 0
    def __init__(self, Description, Filename, Sequence):
        self.Description = Description
        self.Filename = Filename
        self.Sequence = Sequence
        print self.Description
        ShowClass.FileCount += 1
    
    def SaveClass(self, SecondColumn, TitleLines = -1): #save CLASS according to temperature files format. 
        newline = ''
        self.SecondColumn = SecondColumn
        os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer')
        with open(r'PostProcess' + self.Filename + str(self.Sequence) + '.txt') as handle1: # 'PostProcess' + self.Filename + str(self.Sequence) is from "def ExtractData(self):" method in DataProcess.py file.
            with open('WrapPostProcess' + self.Filename + str(self.Sequence) + '.txt','w') as handle2: # just change the pattern of CLASS data. Let CLASS of each row display in the same row. So there will be 30*101.
                lines = handle1.readlines()
                for j, line in enumerate(lines):
                    if line.startswith(("M", "V", "Z")):
                        TitleLines += 1 # record how many title lines in the file.
                        continue
                    words = line.strip().split(' ')
                    #newline.append(words[2])
                    newline = newline + words[2] + ','
                    if (j - TitleLines) != 0 and (j - TitleLines) % SecondColumn == 0:
                        handle2.writelines(newline.strip(","))
                        handle2.write('\n')
                        newline = ''
                        
    def SaveClassMulti(self, SecondColumn, TitleLines = -1): #save CLASS according to temperature files format. 
        newline = ''
        self.SecondColumn = SecondColumn
        os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer')
        with open(r'PostProcess' + self.Filename + str(self.Sequence) + '.txt') as handle1: # 'PostProcess' + self.Filename + str(self.Sequence) is from "def ExtractData(self):" method in DataProcess.py file.
            with open('WrapPostProcess' + self.Filename + str(self.Sequence) + '.txt','w') as handle2: # just change the pattern of CLASS data. Let CLASS of each row display in the same row. So there will be 30*101.
                lines = handle1.readlines()
                for j, line in enumerate(lines):
                    if line.startswith(("M", "V", "Z")):
                        TitleLines += 1 # record how many title lines in the file.
                        continue
                    words = line.strip().split(' ')
                    #newline.append(words[2])
                    newline = newline + words[2] + ','
                    if (j - TitleLines) != 0 and (j - TitleLines) % SecondColumn == 0:
                        handle2.writelines(newline.strip(","))
                        handle2.write('\n')
                        newline = ''
    
    def GetClass(self, _Row, _Column, Seq):# read class into memory
        self._Row = _Row
        self._Column = _Column
        self.Seq = Seq
        Class = np.zeros((_Row, _Column), dtype='int16')# this is created to store temperature data
        try:
            f=file(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer\WrapPostProcess'+self.Filename + str(self.Seq) + '.txt', 'r')
        except IOError:
            print "This class file cannot be found"
        else:
            j=0
            for line in f.readlines():
                line = line.replace("\n","").split(',') #delete '\n' in the line string
                #line = line.split(',') #split string data with comma
                for i in xrange(_Column):
                    Class[j,i]=int(line[i]) # convert string data to float data, delimiter is comma; store every data to ndarray
                j+=1
            f.close()
        return Class

class Ani:
    def __init__(self, start, _ClassFileNumber, _Column_Num):  
        self._Column_Num = _Column_Num
        self.start = start
        self._ClassFileNumber = _ClassFileNumber
        self.fig = plt.figure()
        
    def showclass(self, k):
        for i in xrange(self.start, self._ClassFileNumber): # because the starting temperature file is 400th file, 
            t = ShowClass('Postprocessing the data ' + str(i), 'Class', i) # ShowClass is a class, t is an object of this class.
            t.SaveClass(self._Column_Num)
        classimagelist = [t.GetClass(Valid_rows, self._Column_Num, i) for i in xrange(self.start, self._ClassFileNumber)] # return a list, not a generator
        #fig = plt.figure()
        im = plt.imshow(classimagelist[0], cmap='hot', interpolation='nearest')
        im.set_array(classimagelist[k])
        return [im]

    def ani_static(self): # Show CLASS files not consecutively.
        for i in xrange(self.start, self._ClassFileNumber): # because the starting temperature file is 400th file, 
            t = ShowClass('Postprocessing the data ' + str(i), 'Class', i) # ShowClass is a class, t is an object of this class.
            t.SaveClass(self._Column_Num)
            plt.imshow(t.GetClass(Valid_rows, self._Column_Num, i), cmap='hot', interpolation='nearest')
            plt.show()
#        classimagelist = [t.GetClass(30, self._Column_Num, i) for i in xrange(self.start, self._ClassFileNumber)] # 30 represents there are 30 rows.
#        for i in xrange((self._ClassFileNumber - self.start) * MicroTime):        
##            self.get_temp(i)
#            plt.imshow(classimagelist[i], cmap='hot', interpolation='nearest')
#            plt.show()

    #function to update figure    
    #def updatefig(k):
    #    # set the data in the axesimage object
    #    im.set_array(classimagelist[k])
    #    # return the artists set
    #    return [im]
    #kick off the animation
    
    def ani(self):
#        inlinemode = False
#        if inlinemode:   
#            print "plot will be inline..."
#        else:
#            print "plot will be qt..."
#            plt.switch_backend()
        ani = animation.FuncAnimation(self.fig, self.showclass, frames=xrange((TempFileNumber - TempFileNumber_start) * MicroTime), interval=50, blit = True) # keep in mind that frames should range from 0 to a certain number. 
        plt.show()
"""Operation"""
#==============================================================================
# aa = Ani(10, 101)
# aa.ani()
#==============================================================================

"""***************************************************************"""
class ShowTemp3D:
    def __init__(self):
        pass
    
    def temp3d(self, temp):
        
        obj = contour3d(temp, contours=4, transparent=True)




"""***************************************************************"""

"""Save the class images in to png files. 01202017"""
def ImageToPNG(Filename, Sequence, _Column_Num): 
    t = ShowClass('Postprocessing the data ' + str(Sequence), 'Class', Sequence) # ShowClass is a class, t is an object of this class.
    t.SaveClass(_Column_Num) # 'WrapPostProcessClass' files can be obtained.
    os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer')
    class_wrap = np.genfromtxt(r'WrapPostProcess' + Filename + str(Sequence) + '.txt', delimiter = ',')
    class_wrap_flip = np.flipud(class_wrap)
    plt.imshow(class_wrap_flip, interpolation='nearest')
    plt.savefig(r'WrapPostProcess' + Filename + str(Sequence) + '.png')
    
"""Save the class images in to png files. 11202017. laser path direction for final defense"""
def ImageToPNG_laserpath(Filename, Sequence, _Column_Num): 
    t = ShowClass('Postprocessing the data ' + str(Sequence), 'Class', Sequence) # ShowClass is a class, t is an object of this class.
    t.SaveClass(_Column_Num) # 'WrapPostProcessClass' files can be obtained.
    os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer')
    class_wrap = np.genfromtxt(r'WrapPostProcess' + Filename + str(Sequence) + '.txt', delimiter = ',')
    class_wrap_flip = np.flipud(class_wrap)
    plt.imshow(class_wrap_flip, interpolation='nearest')
    plt.savefig(r'WrapPostProcess' + Filename + str(Sequence) + '.png')
    
"""Save the class images in to png files. 12072017. laser path direction for final defense"""
def ImageToPNG_laserpath_1207(Filename, Sequence, _Column_Num): 
    t = ShowClass('Postprocessing the data ' + str(Sequence), 'Class', Sequence) # ShowClass is a class, t is an object of this class.
    t.SaveClassMulti(_Column_Num) # 'WrapPostProcessClass' files can be obtained.
    os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer')
    class_wrap = np.genfromtxt(r'WrapPostProcess' + Filename + str(Sequence) + '.txt', delimiter = ',')
    class_wrap_flip = np.flipud(class_wrap)
    plt.imshow(class_wrap_flip, interpolation='nearest')
    plt.savefig(r'WrapPostProcess' + Filename + str(Sequence) + '.png')
"""Operation"""
#==============================================================================
# for k in xrange((TempFileNumber_start - 1) * MicroTime, (TempFileNumber - 1) * MicroTime):
#     ImageToPNG('Class', k, ColumnNum)
#==============================================================================

#==============================================================================
# if __name__ == '__main__': # plot grain growth
#     for k in xrange(4004 * MicroTime, 4994 * MicroTime, 5):
#         ImageToPNG_laserpath('Class', k, Valid_columns)
#==============================================================================
        
if __name__ == '__main__': # plot grain growth
    for k in xrange(0* MicroTime, 199 * MicroTime, 1):
        ImageToPNG_laserpath_1207('Class', k, ColumnNum)

#==============================================================================
#     for k in xrange(TempFileNumber_start * MicroTime, TempFileNumber * MicroTime, 5):
# #        ImageToPNG_laserpath('Class', k, ColumnNum)
#         # select some columns to show in the final figures.
#         ImageToPNG_laserpath('Class', k, Valid_columns)
#==============================================================================

def ImageToPNGMulti(Filename, Sequence, _Column_Num): # 02102017. multiple layers.
    t = ShowClass('Postprocessing the data ' + str(Sequence), 'Class', Sequence) # ShowClass is a class, t is an object of this class.
    t.SaveClassMulti(_Column_Num) # 'WrapPostProcessClass' files can be obtained.
    os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_multiple_layer')
    class_wrap = np.genfromtxt(r'WrapPostProcess' + Filename + str(Sequence) + '.txt', delimiter = ',')
    class_wrap_flip = np.flipud(class_wrap)
#    plt.imshow(class_wrap_flip, interpolation='nearest')
    plt.imshow(class_wrap_flip)
    plt.savefig(r'WrapPostProcess' + Filename + str(Sequence) + '.png')
"""Operation"""
#==============================================================================
# if __name__ == '__main__': # plot grain growth
#     TempFileNumber_start = 100
#     TempFileNumber = 250
#     for k in xrange((TempFileNumber_start - 1) * MicroTime, (TempFileNumber - 1) * MicroTime):
#         ImageToPNGMulti('Class', k, ColumnNum)
#==============================================================================
        
"""This is for comprehensive exam presentation. The purpose is to display grain growth within a smaller region. 04212017"""
def ImagetoPNGComp():
    t = ShowClass('Postprocessing the data ' + str(Sequence), 'Class', Sequence) # ShowClass is a class, t is an object of this class.
    t.SaveClassMulti(_Column_Num) # 'WrapPostProcessClass' files can be obtained.
    os.chdir(r'C:\Users\jnzzp5\OneDrive\study\research\CA05122016\CA_2D_onelayer')