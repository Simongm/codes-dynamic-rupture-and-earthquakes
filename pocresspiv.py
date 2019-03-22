# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:28:07 2019

@author: Snhh26
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
import scipy
plt.close('all')
cmtopx=33;
pxtom=1/cmtopx*0.01;

mat=np.loadtxt('pivImJ.txt')
mat=mat*pxtom

x=mat[:,0]
y=mat[:,1]
ux=mat[:,2]
#ux = np.expand_dims(ux, axis=1)
uy=mat[:,3]
#uy = np.expand_dims(uy, axis=1)
#X, Y = np.meshgrid(x,y)
#grid_ux = griddata(X,Y, (x,y), method='nearest')


def fit_params(x, y, a, b, c):
    return a*x+b*y+c

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


############_select points for regression_#######################
r=0.04

exx=np.zeros(x.shape)
exy=np.zeros(x.shape)
eyy=np.zeros(x.shape)

for i, (xp, yp) in enumerate(zip(x,y)):
    ind_reg=np.where(((x-xp)**2+(y-yp)**2)**0.5 < r)[0]
    UX = np.c_[x[ind_reg],y[ind_reg],ux[ind_reg]]
    UY = np.c_[x[ind_reg],y[ind_reg],uy[ind_reg]]
    # best-fit linear plane (1st-order)
    AX = np.c_[UX[:,0], UX[:,1], np.ones(UX.shape[0])]
    CX,_,_,_ = scipy.linalg.lstsq(AX, UX[:,2])    # coefficients
    AY = np.c_[UY[:,0], UY[:,1], np.ones(UY.shape[0])]
    CY,_,_,_ = scipy.linalg.lstsq(AY, UY[:,2])    # coefficients
    exx[i]=CX[0]
    exy[i]=(CX[1]+CY[0])/2
    eyy[i]=CY[1]


fig = plt.figure(1)

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')
ax.scatter(x,y,exy,c='y')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('eyy')


## regular grid covering the domain of the data
f2 = plt.figure(2)

ax2 = f2.add_subplot(111)
X=np.unique(x)
Y=np.unique(y)
XX,YY=np.meshgrid(X,Y)
exy.reshape(XX.shape)
ax2.tricontourf(x,y,exy)



