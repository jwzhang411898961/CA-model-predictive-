# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:37:06 2016

@author: jnzzp5
"""

"""coarse cell"""
#TempFileNumber = 410 #the number of interpolated temperature files: 519
#TempFileNumber_start = 400 # determined by where the temperature will drop.
#ColumnNum = 101 #there are 101 columns in the interpolated temperature file
#RowNum = 171 #there are 171 rows in the interpolated temperature file
#TempLiquid = 1923.0
#TempSolid = 1877.0
#TempMelt = 1900.0
##TempLiquid = 1323.0 # it should be 1423.0
##TempSolid = 1277.0 # it should be 1377.0
##TempMelt = 1300.0 # it should be 1400
#DTSQLIM = (TempLiquid - TempSolid)**2 
#MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
##PV = .1
##PS = .05
##NUCFLG = 0
#UnderCoolMean = 19.5 # assumption
#UnderCoolDelta = 5 # assumption
#PI = 3.1415926
#SDX = 2e-5
#SDY = 2e-5
#TempAdjustor = 100 #assumption
#Fv_gaussian_ini = 5e02 #assumption
#Fs_gaussian_ini = 7e02 #assumption
#Fv_gaussian_nuc = 1e04 #assumption
#Fs_gaussian_nuc = 2e03 #assumption
#DTIME = 9e-9 # determine the growth velocity directly. DTIME higher, grain velocity higher. 
#Invalid_rows = 110
#Valid_rows = RowNum - Invalid_rows - 1
#Pc_nuc = 4e-3
#critical_value1 = 7
#critical_value2 = 14
#    """2D range"""
#row_min1 = 150
#row_max1 = 160
#row_max2 = 171
#column_min1 = 10
#column_max1 = 90



"""fine cell, 01-23-2017 for main_program_v3 and main_program_v4 file"""
#==============================================================================
# TempFileNumber = 403 #the number of interpolated temperature files: 519
# TempFileNumber_start = 400 # determined by where the temperature will drop.
# ColumnNum = 401 #there are 401 columns in the interpolated temperature file
# RowNum = 681 #there are 681 rows in the interpolated temperature file
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
# #TempLiquid = 1323.0 # it should be 1423.0
# #TempSolid = 1277.0 # it should be 1377.0
# #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2 
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
# #PV = .1
# #PS = .05
# #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# DTIME = 9e-9 # determine the growth velocity directly. DTIME higher, grain velocity higher. 
# Invalid_rows = 440 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# critical_value1 = 10
# critical_value2 = 20
# critical_value3 = 30
# critical_value4 = 40
# """2D range"""
# row_min1 = 550
# row_max1 = 583
# row_max2 = 616
# row_max3 = 650
# row_max4 = 681
# column_min1 = 40
# column_max1 = 360
# Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# substrate_ratio = 1.133 # referenced in DataProcess.py file, class Substrate, def substrate().
# Class_mean = 24.0
#==============================================================================


"""fine cell, 01-25-2017 for main_program_v5 file,main_program_v6 file. also works for v7 and v8. It works for v9_vect"""
#==============================================================================
# TempFileNumber = 430 #the number of interpolated temperature files: 519
# TempFileNumber_start = 400 # determined by where the temperature will drop.
# #ColumnNum = 401 #there are 401 columns in the interpolated temperature file
# ColumnNum = 321 # in order to keep symmetry, 80 culumns are deleted.
# RowNum = 681 #there are 681 rows in the interpolated temperature file
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
# #TempLiquid = 1323.0 # it should be 1423.0
# #TempSolid = 1277.0 # it should be 1377.0
# #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2 
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
# #PV = .1
# #PS = .05
# #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 5e02 #assumption
# Fs_gaussian_ini = 6e03 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 1e08 #assumption
# Fs_gaussian_nuc = 2e08 #assumption. PS is higher than PV.
# #==============================================================================
# # # it works for small columnar grains
# # Fv_gaussian_ini = 1e03 #assumption
# # Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# # Fv_gaussian_nuc = 2e08 #assumption
# # Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# #==============================================================================
# DTIME = 2e-10 # determine the growth velocity directly. DTIME higher, grain velocity higher. for instance, DTIME = 9e-9, Vtip * DTIME / math.sqrt(2.0) is approximately 3e-04.
# Invalid_rows = 440 # it may change for different case
# Invalid_columns = 80 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# critical_value1 = 60
# critical_value2 = 80
# critical_value3 = 100
# critical_value4 = 110
# """2D range"""
# row_min1 = 550
# row_max1 = 583
# row_max2 = 616
# row_max3 = 650
# row_max4 = 681
# column_min1 = 40
# column_max1 = 360
# Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# substrate_ratio = 1.133 # referenced in DataProcess.py file, class Substrate, def substrate().
# Class_mean = 24.0 # referenced in def epitaxial() in DataProcess.py file.
# remelt_flag = 800 # determine when the solidified grain remelt.
# file_low = 25
# file_up = 60
#==============================================================================




"""fine cell, 02-22-2017 for main_program_v9_vect_stopgrowth. revised on 02272017. works well for one layer!!!"""
TempFileNumber = 519 #the number of interpolated temperature files: 519
TempFileNumber_start = 400 # determined by where the temperature will drop.
#ColumnNum = 401 #there are 401 columns in the interpolated temperature file
ColumnNum = 321 # in order to keep symmetry, 80 culumns are deleted.
RowNum = 681 #there are 681 rows in the interpolated temperature file
TempLiquid = 1923.0
TempSolid = 1877.0
TempMelt = 1900.0
#TempLiquid = 1323.0 # it should be 1423.0
#TempSolid = 1277.0 # it should be 1377.0
#TempMelt = 1300.0 # it should be 1400
DTSQLIM = (TempLiquid - TempSolid)**2 
MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
#PV = .1
#PS = .05
#NUCFLG = 0
UnderCoolMean = 19.5 # assumption
UnderCoolDelta = 5 # assumption
PI = 3.1415926
SDX = 5e-6
SDY = 5e-6
TempAdjustor = 100 #assumption
Fv_gaussian_ini = 6e02 #assumption
Fs_gaussian_ini = 2.4e04 #assumption.  PS is higher than PV.
Fv_gaussian_nuc = 9e07 #assumption
Fs_gaussian_nuc = 2e09 #assumption. PS is higher than PV.
#==============================================================================
# # it works for small columnar grains
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
#==============================================================================
DTIME = 2.7e-10 # determine the growth velocity directly. DTIME higher, grain velocity higher. for instance, DTIME = 9e-9, Vtip * DTIME / math.sqrt(2.0) is approximately 3e-04.
Invalid_rows = 440 # it may change for different case
Invalid_columns = 80 # it may change for different case
Valid_rows = RowNum - Invalid_rows - 1
Pc_nuc = 1e-3
critical_value1 = 60
critical_value2 = 80
critical_value3 = 100
critical_value4 = 110
critical_CET = 100 # Control CET size.
"""2D range"""
row_min1 = 550
row_max1 = 583
row_max2 = 616
row_max3 = 650
row_max4 = 681
column_min1 = 40
column_max1 = 360
Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
substrate_ratio = 1.133 # referenced in DataProcess.py file, class Substrate, def substrate().
Class_mean = 24.0 # referenced in def epitaxial() in DataProcess.py file.
remelt_flag = 800 # determine when the solidified grain remelt.
file_low = 25
file_up = 60



"""fine cell, 02-09-2017 for main_program_v1.0 file for multiple layer deposition. revised 02282017"""
#==============================================================================
# TempFileNumber =170 #the number of interpolated temperature files: 828
# TempFileNumber_start = 100 # determined by where the temperature will drop.
# ColumnNum = 601 #there are 601 columns in the interpolated temperature file
# RowNum = 233 #there are 233 rows in the interpolated temperature file
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
# #TempLiquid = 1323.0 # it should be 1423.0
# #TempSolid = 1277.0 # it should be 1377.0
# #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2 
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
# #PV = .1
# #PS = .05
# #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# DTIME = 1e-10 # determine the growth velocity directly. DTIME higher, grain velocity higher. for instance, DTIME = 9e-9, Vtip * DTIME / math.sqrt(2.0) is approximately 3e-04.
# Invalid_rows = 0 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# critical_value1 = 10 #?
# critical_value2 = 20
# critical_value3 = 30
# critical_value4 = 40
# """2D range"""
# row_min1 = 1
# row_max1 = 59
# row_max2 = 117
# row_max3 = 175
# row_max4 = 233
# column_min1 = 0
# column_max1 = 600
# Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# substrate_ratio = 7.25 # referenced in DataProcess.py file, class Substrate, def substrate(). 232/32 == 7.25. only 32 rows are substrate.
# Class_mean = 24.0 # referenced in def epitaxial() in DataProcess.py file.
# remelt_flag = 800 # determine when the solidified grain remelt.
# room_temp = 300
# air_sub_BDY = 30 # 30 is specifically for this case. It represents 0-30 rows are substrate. rows higher than 30 are air originally.
#==============================================================================


"""fine cell, 02-09-2017 for main_program_v1.1 file for multiple layer deposition. revised 03172017"""

#==============================================================================
# TempFileNumber =450 #the number of interpolated temperature files: 828
# TempFileNumber_start = 100 # determined by where the temperature will drop.
# ColumnNum = 601 #there are 601 columns in the interpolated temperature file
# RowNum = 233 #there are 233 rows in the interpolated temperature file
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
#  #TempLiquid = 1323.0 # it should be 1423.0
#  #TempSolid = 1277.0 # it should be 1377.0
#  #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2 
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
# #PV = .1
# #PS = .05
# #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# DTIME = 1e-10 # determine the growth velocity directly. DTIME higher, grain velocity higher. for instance, DTIME = 9e-9, Vtip * DTIME / math.sqrt(2.0) is approximately 3e-04.
# Invalid_rows = 0 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# critical_value1 = 10 #?
# critical_value2 = 20
# critical_value3 = 30
# critical_value4 = 40
# """2D range"""
# row_min1 = 1
# row_max1 = 59
# row_max2 = 117
# row_max3 = 175
# row_max4 = 233
# column_min1 = 0
# column_max1 = 600
# Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# substrate_ratio = 7.25 # referenced in DataProcess.py file, class Substrate, def substrate(). 232/32 == 7.25. only 32 rows are substrate.
# Class_mean = 24.0 # referenced in def epitaxial() in DataProcess.py file.
# remelt_flag = 800 # determine when the solidified grain remelt.
# room_temp = 300
# air_sub_BDY = 30 # 30 is specifically for this case. It represents 0-30 rows are substrate. rows higher than 30 are air originally.
#==============================================================================


"""fine cell, 04-28-2017 for main_program_v1.2 file for multiple layer deposition. revised 04-28-2017"""
#==============================================================================
# 
# TempFileNumber =207 #the number of interpolated temperature files: 828
# TempFileNumber_start = 100 # determined by where the temperature will drop.
# ColumnNum = 400 #there are 400 columns in the extracted temperature file
# RowNum = 200 #there are 200 rows in the extracted temperature file
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
#  #TempLiquid = 1323.0 # it should be 1423.0
#  #TempSolid = 1277.0 # it should be 1377.0
#  #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2 
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
# #PV = .1
# #PS = .05
# #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# Fm_gaussian_nuc = 5e01 #multiple layer re-nucleation. 04282017
# DTIME = 2e-10 # determine the growth velocity directly. DTIME higher, grain velocity higher. for instance, DTIME = 9e-9, Vtip * DTIME / math.sqrt(2.0) is approximately 3e-04.
# Invalid_rows = 0 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# #==============================================================================
# # critical_value1 = 10 #?
# # critical_value2 = 20
# # critical_value3 = 30
# # critical_value4 = 40
# #==============================================================================
# """2D range"""
# row_min1 = 1
# row_max1 = 59
# row_max2 = 117
# row_max3 = 175
# row_max4 = 233
# column_min1 = 0
# column_max1 = 600
# #Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# Integral_coefficient2 = 200.0 # in order to prevent PM from zero.
# Integral_coefficient_s = 100.0
# Integral_coefficient_v = 200.0
# substrate_ratio = 7.25 # referenced in DataProcess.py file, class Substrate, def substrate(). 232/32 == 7.25. only 32 rows are substrate.
# Class_mean = 24.0 # referenced in def epitaxial() in DataProcess.py file.
# remelt_flag = 1000 # determine when the solidified grain remelt.
# room_temp = 300
# air_sub_BDY = 30 # 30 is specifically for this case. It represents 0-30 rows are substrate. rows higher than 30 are air originally.
# nucleat_spacing = 3 # 04292017 create. nucleation cell minimum spacing. prevent nucleation cells are too close.
# offset = 2
#==============================================================================


"""fine cell, 04-21-2017 for main_program_v3_comprehensiveexam and main_program_v3_comprehensiveexam files."""
#==============================================================================
# TempFileNumber = 518 #the number of interpolated temperature files: 519
# TempFileNumber_start = 400 # determined by where the temperature will drop.
# #ColumnNum = 401 #there are 401 columns in the interpolated temperature file
# ColumnNum = 321 # in order to keep symmetry, 80 culumns are deleted.
# RowNum = 681 #there are 681 rows in the interpolated temperature file
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
#  #TempLiquid = 1323.0 # it should be 1423.0
#  #TempSolid = 1277.0 # it should be 1377.0
#  #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2 
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
#  #PV = .1
#  #PS = .05
#  #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# DTIME = 9e-9 # determine the growth velocity directly. DTIME higher, grain velocity higher. 
# Invalid_rows = 440 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# critical_value1 = 10
# critical_value2 = 20
# critical_value3 = 30
# critical_value4 = 40
# """2D range"""
# row_min1 = 550
# row_max1 = 583
# row_max2 = 616
# row_max3 = 650
# row_max4 = 681
# #column_min1 = 40
# #column_max1 = 360
# column_min1 = 20
# column_max1 = 300 # because the maximum ColumnNum is 321.
# Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# substrate_ratio = 1.133 # referenced in DataProcess.py file, class Substrate, def substrate().
# Class_mean = 24.0
#==============================================================================

"""06-13-2017 for main_program_v1.3 file for multiple layer deposition. revised 06-13-2017. this is for yanlei's ansys temperature three layer result."""
#==============================================================================
# TempFileNumber =207 #the number of interpolated temperature files: 828
# TempFileNumber_start = 100 # determined by where the temperature will drop.
# ColumnNum = 400 #there are 400 columns in the extracted temperature file in yanlei's temperature result
# RowNum = 200 #there are 200 rows in the extracted temperature file in yanlei's temperature result
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
#  #TempLiquid = 1323.0 # it should be 1423.0
#  #TempSolid = 1277.0 # it should be 1377.0
#  #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2 
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
# #PV = .1
# #PS = .05
# #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# Fm_gaussian_nuc = 5e01 #multiple layer re-nucleation. 04282017
# DTIME = 2e-10 # determine the growth velocity directly. DTIME higher, grain velocity higher. for instance, DTIME = 9e-9, Vtip * DTIME / math.sqrt(2.0) is approximately 3e-04.
# Invalid_rows = 0 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# #==============================================================================
# # critical_value1 = 10 #?
# # critical_value2 = 20
# # critical_value3 = 30
# # critical_value4 = 40
# #==============================================================================
# """2D range"""
# row_min1 = 1
# row_max1 = 59
# row_max2 = 117
# row_max3 = 175
# row_max4 = 233
# column_min1 = 0
# column_max1 = 600
# #Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# Integral_coefficient2 = 200.0 # in order to prevent PM from zero.
# Integral_coefficient_s = 100.0
# Integral_coefficient_v = 200.0
# substrate_ratio = 7.25 # referenced in DataProcess.py file, class Substrate, def substrate(). 232/32 == 7.25. only 32 rows are substrate.
# Class_mean = 24.0 # referenced in def epitaxial() in DataProcess.py file.
# remelt_flag = 1000 # determine when the solidified grain remelt.
# room_temp = 300
# air_sub_BDY = 30 # 30 is specifically for this case. It represents 0-30 rows are substrate. rows higher than 30 are air originally.
# nucleat_spacing = 3 # 04292017 create. nucleation cell minimum spacing. prevent nucleation cells are too close.
# offset = 2
#==============================================================================



"""07-18-2017 for main_program_v1.3 file for multiple layer deposition. revised 07-18-2017. this is for yanlei's ansys temperature fifteen layer result. 09112017 revision again. It also works for v1.4 version"""
#==============================================================================
# TempFileTotal = 960
# TempFileNumber = 130 #the number of interpolated temperature files: 
# TempFileNumber_start = 1 # determined by where the temperature will drop.
# ColumnNum = 1359 #there are 681 columns in the extracted temperature file in yanlei's temperature result
# RowNum = 819 #there are 411 rows in the extracted temperature file in yanlei's temperature result
# TempLiquid = 1923.0
# TempSolid = 1877.0
# TempMelt = 1900.0
#  #TempLiquid = 1323.0 # it should be 1423.0
#  #TempSolid = 1277.0 # it should be 1377.0
#  #TempMelt = 1300.0 # it should be 1400
# DTSQLIM = (TempLiquid - TempSolid)**2
# MicroTime = 1 # one temperature file at temporal dimension in FEA corresponds to MicroTime temperature files in CA 
# #PV = .1
# #PS = .05
# #NUCFLG = 0
# UnderCoolMean = 19.5 # assumption
# UnderCoolDelta = 5 # assumption
# PI = 3.1415926
# SDX = 5e-6
# SDY = 5e-6
# TempAdjustor = 100 #assumption
# Fv_gaussian_ini = 1e03 #assumption
# Fs_gaussian_ini = 1.2e04 #assumption.  PS is higher than PV.
# Fv_gaussian_nuc = 2e08 #assumption
# Fs_gaussian_nuc = 4e08 #assumption. PS is higher than PV.
# Fm_gaussian_nuc = 5e01 #multiple layer re-nucleation. 04282017
# DTIME = 2e-10 # determine the growth velocity directly. DTIME higher, grain velocity higher. for instance, DTIME = 9e-9, Vtip * DTIME / math.sqrt(2.0) is approximately 3e-04.
# Invalid_rows = 0 # it may change for different case
# Valid_rows = RowNum - Invalid_rows - 1
# Pc_nuc = 1e-3
# #==============================================================================
# # critical_value1 = 10 #?
# # critical_value2 = 20
# # critical_value3 = 30
# # critical_value4 = 40
# #==============================================================================
# """2D range"""
# row_min1 = 1
# row_max1 = 59
# row_max2 = 117
# row_max3 = 175
# row_max4 = 233
# column_min1 = 0
# column_max1 = 600
# #Integral_coefficient = 100.0 # in order to prevent PS and PV from zero.
# Integral_coefficient2 = 200.0 # in order to prevent PM from zero.
# Integral_coefficient_s = 100.0
# Integral_coefficient_v = 200.0
# substrate_ratio = 3.7 # referenced in DataProcess.py file, class Substrate, def substrate(). 232/32 == 7.25. only 32 rows are substrate.
# Class_mean = 24.0 # referenced in def epitaxial() in DataProcess.py file.
# remelt_flag = TempFileNumber + 1 # determine when the solidified grain remelt. make sure the code part under remelt_flag won't be execuated.
# room_temp = 300 # not 300 because some starting point is a little bit lower than 300, e.g. 298
# air_sub_BDY = int(RowNum / substrate_ratio) # RowNum / substrate_ratio is specifically for this case. It represents 0-air_sub_BDY rows are substrate. rows higher than air_sub_BDY are air originally.
# nucleat_spacing = 4 # 04292017 create. nucleation cell minimum spacing. prevent nucleation cells are too close.
# offset = 2
# 
# #09182017 add
# InterpolationTimes = 20
# ColCount = 60
# LayerEleThick = InterpolationTimes * 2 # 20 denotes 20 times interpolation, and 2 is original element thickness.
# LaserEleNum = 200
# Multi_path1 = r"C:/Users/jnzzp5/OneDrive/study/research/CA05122016/CA_2D_multiple_layer/"
# polygon_dia = 2*50000 # the largest polygon nominal area that are allowed
# first_layer = 60 # the number of time step of the first layer.
#==============================================================================
