# -*- coding: utf-8 -*-
"""
Created on Sat Jun 13 10:15:56 2015

@author: BRENNICH
"""

import numpy as np
import os.path as op
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import os

if os.name =='nt':
    syspath = op.normpath("Y:/inhouse/opd29/20150611")

bufferpath = op.normpath(op.join(syspath,'IECBuffer'))
bufferpath1D = op.normpath(op.join(bufferpath,'1d'))


framenr = 1
filename=op.normpath(op.join(bufferpath1D,"buffer_007_%(frame_number)05i.dat" % {"frame_number":framenr}))
data_in_file = np.genfromtxt(filename)
buffers = data_in_file[:,1]
framenr = framenr+1
filename=op.normpath(op.join(bufferpath1D,"buffer_007_%(frame_number)05i.dat" % {"frame_number":framenr}))
while op.exists(filename):
    print(filename)
    data_in_file = np.genfromtxt(filename)
    buffers = np.vstack((buffers,data_in_file[:,1]))
    framenr = framenr+1
    filename=op.normpath(op.join(bufferpath1D,"buffer_007_%(frame_number)05i.dat" % {"frame_number":framenr}))
#    
#fig = pl.figure("buffers")
#ax = fig.add_subplot(1,1,1)  
#ax.plot(np.arange(0,10,0.1) ) 
#ax.set_xlim((0.04,1.7))
##ax.set_yscale('log')
#fig.show()
    
    
buffers_normed = buffers/buffers[0,:]    

plt.plot(buffers[:,898:1004].sum(axis=1))

plt.show()

datapath = op.normpath(op.join(syspath,'D5GNK'))
datapath1D = op.normpath(op.join(datapath,'1d'))

framenr = 1
filename=op.normpath(op.join(datapath1D,"D5GNK_008_%(frame_number)05i.dat" % {"frame_number":framenr}))
data_in_file = np.genfromtxt(filename)
data = data_in_file[:,1]
framenr = framenr+1
filename=op.normpath(op.join(datapath1D,"D5GNK_008_%(frame_number)05i.dat" % {"frame_number":framenr}))
while op.exists(filename):
    print(filename)
    data_in_file = np.genfromtxt(filename)
    data = np.vstack((data,data_in_file[:,1]))
    framenr = framenr+1
    filename=op.normpath(op.join(datapath1D,"D5GNK_008_%(frame_number)05i.dat" % {"frame_number":framenr}))

data_normed = data/data[0,:]  

buffers_normed = buffers/buffers[0,:]     
    
plt.plot(buffers[:,898:1004].sum(axis=1))
plt.plot(data[:,898:1004].sum(axis=1))

plt.show()

plt.xlim(750,1500)
plt.ylim(3050,3170)


plt.plot(buffers[:,1004:1034].sum(axis=1))
plt.plot(data[:,1004:1034].sum(axis=1))

plt.show()

plt.xlim(750,1500)

plt.plot(data[:,950:1034].sum(axis=1))
plt.plot(buffers[:,950:1034].sum(axis=1))
plt.xlim(750,2000)

def calculate_ratio(data,buff1,buff2,perc):
    buff = perc/100.0 * buff1 + (100.0-perc)/100.0 * buff2
    BC = data-buff
    ratio = BC[:,16:101].mean(axis=1)/BC[:,310:530].mean(axis=1)
    return ratio/ratio[832:841].mean()
    
    
buffer2100_2200 = buffers[2100:2200,:].mean(axis=0)

buffer1300_1400 = buffers[1300:1400,:].mean(axis=0)

buffer1850_1950 = buffers[1850:1950,:].mean(axis=0)

pts = np.arange(0,105,5)
ratios = [calculate_ratio(data,buffer2100_2200,buffer1850_1950,i)[832:841].var(axis=0) for i in pts]

plt.plot(calculate_ratio(data,buffer1300_1400,buffer1850_1950,90))

plt.xlim(800, 860)
plt.ylim(-2,2)
plt.plot(calculate_ratio(data,buffer1300_1400,buffer1850_1950,10))
