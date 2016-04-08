# -*- coding: utf-8 -*-
"""
Created on Tue Aug 04 11:42:52 2015

@author: BRENNICH
"""

import numpy as np
import os.path as op
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import os
import h5py
#import autorg
#import rti
from scipy.optimize import curve_fit


if os.name =='nt':
    syspath = op.normpath("Y:/inhouse/opd29/20150417/ION_X")

bufferpath = op.normpath(op.join(syspath,'bufferAwo_step'))
bufferpath1D = op.normpath(op.join(bufferpath,'1d'))


framenr = 1
filename=op.normpath(op.join(bufferpath1D,"buffer_007_%(frame_number)05i.dat" % {"frame_number":framenr}))
data_in_file = np.genfromtxt(filename)
buffers = data_in_file[:,1]
q = data_in_file[:,0]
buffstd = data_in_file[:,2]
framenr = framenr+1
filename=op.normpath(op.join(bufferpath1D,"buffer_007_%(frame_number)05i.dat" % {"frame_number":framenr}))
while op.exists(filename):
    print(filename)
    data_in_file = np.genfromtxt(filename)
    buffers = np.vstack((buffers,data_in_file[:,1]))
    print(framenr, buffstd.shape)
    buffstd = np.vstack((buffstd,data_in_file[:,2]))
    framenr = framenr+1
    filename=op.normpath(op.join(bufferpath1D,"buffer_007_%(frame_number)05i.dat" % {"frame_number":framenr}))
    



buffers_normed = buffers/buffers[0,:] 
    
for i in range(0,1000,100):
    plt.plot(buffers_normed.transpose()[i,0:3000]+i/1000)

plt.show()

datapath = op.normpath(op.join(syspath,'bsaA_Wo_Fewsteps'))
datapath1D = op.normpath(op.join(datapath,'1d'))

framenr = 1
filename=op.normpath(op.join(datapath1D,"BSAfewsteps_012_%(frame_number)05i.dat" % {"frame_number":framenr}))
data_in_file = np.genfromtxt(filename)
data = data_in_file[:,1]
datastd = data_in_file[:,2]
framenr = framenr+1
filename=op.normpath(op.join(datapath1D,"BSAfewsteps_012_%(frame_number)05i.dat" % {"frame_number":framenr}))
while op.exists(filename):
    print(filename)
    data_in_file = np.genfromtxt(filename)
    data = np.vstack((data,data_in_file[:,1]))
    datastd = np.vstack((datastd,data_in_file[:,2]))
    framenr = framenr+1
    filename=op.normpath(op.join(datapath1D,"BSAfewsteps_012_%(frame_number)05i.dat" % {"frame_number":framenr}))
    
    
data_normed = data/data[0,:]    

for i in range(0,1100,100):
    plt.plot(buffers_normed.transpose()[i,0:3000]+i/1000)
    plt.plot(data_normed.transpose()[i,0:3000]+i/1000)

plt.show()



plt.plot(buffers[:,898:1004].sum(axis=1))
plt.plot(data[:,898:1004].sum(axis=1))

plt.show()




plt.plot(buffers[:,1004:1034].sum(axis=1))
plt.plot(data[:,1004:1034].sum(axis=1))


plt.show()



plt.plot(data[:,:].sum(axis=1))


plt.show()


plt.plot(data[:,:].sum(axis=1)/data[:,:].sum())
plt.plot(data[:,950:].sum(axis=1)/data[:,950:].sum())
plt.show()

bufferA = data[1100:1150,:].mean(axis=0)
bufferB = data[1520:1570,:].mean(axis=0)

stdA = datastd[1100:1150,:].mean(axis=0)
stdB = datastd[1520:1570,:].mean(axis=0)


meanA = bufferA[950:].mean()
meanB = bufferB[950:].mean()

def func(N, N0, rate):
    return meanB - (meanB-meanA)*(np.exp(-((N-N0)/rate)))
        
    
xdata = np.arange(1290,1360)
ydata = data[1290:1360,950:].mean(axis=1)
    
popt, pcov = curve_fit(func, xdata, ydata, p0 = (1200,20))    

plt.plot(xdata,func(xdata,popt[0],popt[1]))
plt.plot(xdata,ydata)

xcontrol = np.arange(1290,1500)
ycontrol = data[1290:1500,950:].mean(axis=1)

plt.plot(xcontrol,func(xcontrol,popt[0],popt[1]))
plt.plot(xcontrol,ycontrol)

def bufferFit(N):
    N0 = popt[0]
    rate = popt[1]
    return bufferB - (bufferB-bufferA)*(np.exp(-((N-N0)/rate)))    
    
def bufferStd(N):    
    N0 = popt[0]
    rate = popt[1]
    return stdB - (stdB-stdA)*(np.exp(-((N-N0)/rate)))    

BC = data[1289,:] - bufferFit(1289)
BCstd = np.sqrt((datastd[1289,:]*datastd[1289,:] + bufferStd(1289)**2))/2
for N in xcontrol:
    BC = np.vstack((BC,data[N,:] - bufferFit(N)))
    BCstd = np.vstack((BCstd,np.sqrt((datastd[N,:]*datastd[N,:] + bufferStd(N)**2))/2))

rg0 = []
RTI0 = []
    
for i in xcontrol:
    filename = "bsa_012_%(frame_number)05i_exp.dat" % {"frame_number":i}
    scat = np.array([q,BC[i-1289,:],BCstd[i-1289,:]])
   # np.savetxt(filename,scat.transpose(),fmt='%.4e')
    #rgt = autorg.runautorg(filename,'temprg.dat')
    #if rgt is None:
    #    rgt = tuple([0] * 6)
   #     rtt = (0)
    #else: 
   #     rtt= (rti.RamboTainerInvariant(scat.transpose(), rgt[0], rgt[1], rgt[2], rgt[3], rgt[4], qmax=2)['mass'])
    #RTI0.append(rtt)   
    #rg0.append(rgt)
   # print(rg0[i-1290])
    
def calculate_ratio(BC):
    ratio = BC[:,16:101].mean(axis=1)/BC[:,310:530].mean(axis=1)
    return ratio   
    
plt.plot(calculate_ratio(BC)/calculate_ratio(BC[50:52,:]).mean(), c = 'k')
plt.plot(calculate_ratio(data[1290:1500,:] - bufferA)/calculate_ratio(data[1340:1342,:] - bufferA).mean(), c = 'r')
plt.plot(calculate_ratio(data[1290:1500,:] - bufferB)/calculate_ratio(data[1340:1342,:] - bufferB).mean(), c = 'b')