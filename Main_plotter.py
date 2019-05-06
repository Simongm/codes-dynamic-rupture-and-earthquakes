# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 14:46:33 2017

@author: Snhh26
"""
import os,sys
import matplotlib.pyplot as plt
import numpy as np


## ADD THE PATH TO FIND AND IMPORT THE READ FUNCTION ####
sys.path.append("C:\\Users\\Snhh26\\Desktop\\DU965_9June2017\\functions")
from read_function import read_function

##SET THE WORKING DIRECTORY: ######################################
os.chdir('C:\\Users\\Snhh26\\Desktop\\DU965_9June2017')
print(os.getcwd()) # to check working directory where the files are

################################parameters to change#################################
#richter='richter-m' #name of richter 
first=5 #number of 1st richter file
last=5 #number of last richter-file
### if we want to display  only a few datasets:
dataset1=1 # 1st dataset
dataset2='all' #last dataset, write dataset2='all' if you want to display everything

## if you want to take only few points in the dataset, otherwise leave default values for 1MHz
firstindex=0
lastindex=4980736
cuto=10000000 # low pass cutoff frequency

ratio=1 ## decimating ration, chose a multiple of 2 : [1,2,4,7,8,14,16,28,32,56,64,112,128,224,256,448,512,896,1024,1792,2048,3584,4096,7168,8192,14336,16384,28672,32768,57344,65536,114688,229376,458752]
filtercondition=0 # 0= no filter; 1=lowpass at cuto freq
indochannels=[0] #indices of channels to display, =[0,1,2,3] if 4 channels



[data_m,time_m]=read_function('richter-m',first,last,dataset1,dataset2,firstindex,lastindex,cuto,ratio,filtercondition,indochannels)
[data_s1,time_s1]=read_function('richter-s1',first,last,dataset1,dataset2,firstindex,lastindex,cuto,ratio,filtercondition,indochannels)
#
#example of plot
plt.plot(time_m,data_m)
plt.xlabel('time in s')
plt.ylabel('voltage in V')

## example of event counter on 1 channel..........................
threshold=0.5
Window=0.01

cond=data_m>0.5 ## find all the points with amplitude > 0.5V
indall=(np.where(cond))[0]

event=data_m[indall]
WindowL=200 # window length to adjust, not to pick the same event

## loop to remove points in a same event
t0=0
for i in range (0,len(indall)):
    t1=time_m[indall[i]]
    if t1-t0<0.005: 
        indall[i]=0
    t0=t1
    
plt.plot(time_m[indall],data_m[indall],'r.')
