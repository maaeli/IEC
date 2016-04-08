# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 17:01:26 2015

@author: BRENNICH
"""


import numpy as np
import os.path as op
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import os
import h5py
import autorg
import rti


if os.name =='nt':
    syspath = op.normpath("Y:/inhouse/opd29/20150417/ION_X")

bufferpath = op.normpath(op.join(syspath,'bufferAwo_linera_try2'))
bufferpath1D = op.normpath(op.join(bufferpath,'1d'))


framenr = 1
filename=op.normpath(op.join(bufferpath1D,"buffer_006_%(frame_number)05i.dat" % {"frame_number":framenr}))
data_in_file = np.genfromtxt(filename)
buffers = data_in_file[:,1]
q = data_in_file[:,0]
buffstd = data_in_file[:,2]
framenr = framenr+1
filename=op.normpath(op.join(bufferpath1D,"buffer_006_%(frame_number)05i.dat" % {"frame_number":framenr}))
while op.exists(filename):
    print(filename)
    data_in_file = np.genfromtxt(filename)
    buffers = np.vstack((buffers,data_in_file[:,1]))
    print(framenr, buffstd.shape)
    buffstd = np.vstack((buffstd,data_in_file[:,2]))
    framenr = framenr+1
    filename=op.normpath(op.join(bufferpath1D,"buffer_006_%(frame_number)05i.dat" % {"frame_number":framenr}))
#    

#plt.plot(buffers[:,898:1004].sum(axis=1))



buffers_normed = buffers/buffers[0,:] 

for i in range(0,1000,100):
    plt.plot(buffers_normed.transpose()[i,0:3000]+i/1000)

plt.show()



datapath = op.normpath(op.join(syspath,'bsaAwo_linear2'))
datapath1D = op.normpath(op.join(datapath,'1d'))

framenr = 1
filename=op.normpath(op.join(datapath1D,"bsa_010_%(frame_number)05i.dat" % {"frame_number":framenr}))
data_in_file = np.genfromtxt(filename)
data = data_in_file[:,1]
datastd = data_in_file[:,2]
framenr = framenr+1
filename=op.normpath(op.join(datapath1D,"bsa_010_%(frame_number)05i.dat" % {"frame_number":framenr}))
while op.exists(filename):
    print(filename)
    data_in_file = np.genfromtxt(filename)
    data = np.vstack((data,data_in_file[:,1]))
    datastd = np.vstack((datastd,data_in_file[:,2]))
    framenr = framenr+1
    filename=op.normpath(op.join(datapath1D,"bsa_010_%(frame_number)05i.dat" % {"frame_number":framenr}))
    
    
data_normed = data/data[0,:]    

for i in range(0,1100,100):
    plt.plot(buffers_normed.transpose()[i,0:3000]+i/1000)
    plt.plot(data_normed.transpose()[i,0:3000]+i/1000)

plt.show()

def calculate_ratio(data,buff,shift):
    BC = data-np.roll(buff,shift,axis=0)
    ratio = BC[:,16:101].mean(axis=1)/BC[:,310:530].mean(axis=1)
    return ratio


ratioshift = calculate_ratio(data,buffers,80)
for i in range(90,185,10):
    plt.plot(calculate_ratio(data,buffers,i),label=str(i))
    ratioshift = np.vstack((ratioshift,calculate_ratio(data,buffers,i)))

plt.ylim(60,200)
plt.xlim(1000,1400)
plt.legend()

for i in range(0,1100,100):
    plt.plot(np.roll(buffers_normed,130,axis=0).transpose()[i,0:3000]+i/1000)
    plt.plot(data_normed.transpose()[i,0:3000]+i/1000)


#BCstd130 = np.sqrt(datastd*datastd  + np.roll(buffstd*buffstd,130,axis=0))
#BC130 = data-np.roll(buffers,130,axis=0)
#
#rg130 = []
#RTI130 = []
#
#for i in range(2900):
#    filename = "bsa_010_%(frame_number)05i_shift130.dat" % {"frame_number":i}
#    scat = np.array([q,BC130[i,:],BCstd130[i,:]])
#    np.savetxt(filename,scat.transpose(),fmt='%.4e')
#    rgt = autorg.runautorg(filename,'temprg.dat')
#    if rgt is None:
#        rgt = tuple([0] * 6)
#        rtt = (0)
#    else: 
#        rtt= (rti.RamboTainerInvariant(scat.transpose(), rgt[0], rgt[1], rgt[2], rgt[3], rgt[4], qmax=2)['mass'])
#    RTI130.append(rtt)   
#    rg130.append(rgt)
#    print(rg130[i])
#    
#    
#for i in range(140,185,10):
#    plt.plot(calculate_ratio(data,buffers,i),label=str(i))
#
#plt.ylim(50,150)
#plt.xlim(1350,1650)
#plt.legend()    
#
#for i in range(0,1100,100):
#    plt.plot(np.roll(buffers_normed,170,axis=0).transpose()[i,0:3000]+i/1000)
#    plt.plot(data_normed.transpose()[i,0:3000]+i/1000)
#
#BCstd170 = np.sqrt(datastd*datastd  + np.roll(buffstd*buffstd,170,axis=0))
#BC170 = data-np.roll(buffers,170,axis=0)
#
#rg170 = []
#RTI170 = []
#
#for i in range(2900):
#    filename = "bsa_010_%(frame_number)05i_shift170.dat" % {"frame_number":i}
#    scat = np.array([q,BC170[i,:],BCstd170[i,:]])
#    np.savetxt(filename,scat.transpose(),fmt='%.4e')
#    rgt = autorg.runautorg(filename,'temprg.dat')
#    if rgt is None:
#        rgt = tuple([0] * 6)
#        rtt = (0)
#    else: 
#        rtt= (rti.RamboTainerInvariant(scat.transpose(), rgt[0], rgt[1], rgt[2], rgt[3], rgt[4], qmax=2)['mass'])
#    RTI170.append(rtt)   
#    rg170.append(rgt)
#    print(rg170[i])
#
#BCstd0 = np.sqrt(datastd*datastd  + np.roll(buffstd*buffstd,0,axis=0))
#BC0 = data-np.roll(buffers,0,axis=0)
#
#rg0 = []
#RTI0 = []
#
#for i in range(2900):
#    filename = "bsa_010_%(frame_number)05i_shift000.dat" % {"frame_number":i}
#    scat = np.array([q,BC0[i,:],BCstd0[i,:]])
#    np.savetxt(filename,scat.transpose(),fmt='%.4e')
#    rgt = autorg.runautorg(filename,'temprg.dat')
#    if rgt is None:
#        rgt = tuple([0] * 6)
#        rtt = (0)
#    else: 
#        rtt= (rti.RamboTainerInvariant(scat.transpose(), rgt[0], rgt[1], rgt[2], rgt[3], rgt[4], qmax=2)['mass'])
#    RTI0.append(rtt)   
#    rg0.append(rgt)
#    print(rg0[i])
#        
BCstd120 = np.sqrt(datastd*datastd  + np.roll(buffstd*buffstd,120,axis=0))
BC120 = data-np.roll(buffers,120,axis=0)

rg120 = []
RTI120 = []

for i in range(2900):
    filename = "bsa_010_%(frame_number)05i_shift120.dat" % {"frame_number":i}
    scat = np.array([q,BC120[i,:],BCstd120[i,:]])
    np.savetxt(filename,scat.transpose(),fmt='%.4e')
    rgt = autorg.runautorg(filename,'temprg.dat')
    if rgt is None:
        rgt = tuple([0] * 6)
        rtt = (0)
    else: 
        rtt= (rti.RamboTainerInvariant(scat.transpose(), rgt[0], rgt[1], rgt[2], rgt[3], rgt[4], qmax=2)['mass'])
    RTI120.append(rtt)   
    rg120.append(rgt)
    print(rg120[i])
    
BCstd140 = np.sqrt(datastd*datastd  + np.roll(buffstd*buffstd,140,axis=0))
BC140 = data-np.roll(buffers,140,axis=0)

rg140 = []
RTI140 = []

for i in range(2900):
    filename = "bsa_010_%(frame_number)05i_shift140.dat" % {"frame_number":i}
    scat = np.array([q,BC140[i,:],BCstd140[i,:]])
    np.savetxt(filename,scat.transpose(),fmt='%.4e')
    rgt = autorg.runautorg(filename,'temprg.dat')
    if rgt is None:
        rgt = tuple([0] * 6)
        rtt = (0)
    else: 
        rtt= (rti.RamboTainerInvariant(scat.transpose(), rgt[0], rgt[1], rgt[2], rgt[3], rgt[4], qmax=2)['mass'])
    RTI140.append(rtt)   
    rg140.append(rgt)
    print(rg140[i])