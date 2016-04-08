# -*- coding: utf-8 -*-
#"""
#Created on Thu Aug 01 15:17:31 2013

#@author: BRENNICH
#"""

import numpy as np
import os
import os.path as op
import glob
import subprocess as sp
import matplotlib.pylab as pl
from math import floor,sqrt
import scipy.optimize as scop


class FitError(Exception):
    def __init__(self):
        self.value = "No good fit found"
    def __str__(self):
        return repr(self.value)
        
        
 
def chisq(x,data):
    """Note dataformat: x**2, y, 1/sigma**2"""
    return (data[:,2]*(data[:,1]-x[0] + x[1]*data[:,0])**2).sum()

def blocklin(data):
    sigscale = data[:,2].mean()
    data1 = np.array(data)
    data1[:,0] = data[:,0]**2
    data1[:,2] = (sigscale/data[:,2])**2
    chimax = 10000
    chimin = chimax
    for index1 in range(1,490):
        sigma = (data1[index1:index1+9,2]).sum()
        sigmay = (data1[index1:index1+9,1]*data1[index1:index1+9,2]).sum()
        sigmax = (data1[index1:index1+9,0]*data1[index1:index1+9,2]).sum()
        sigmaxy = (data1[index1:index1+9,0]*data1[index1:index1+9,1]*data1[index1:index1+9,2]).sum() 
        sigmaxx = (data1[index1:index1+9,0]**2*data1[index1:index1+9,2]).sum()
        for index2 in range(index1+10,500):
           dat = data1[index1:index2,:]
           sigma += data1[index2,2]
           sigmay += (data1[index2,1]*data1[index2,2])
           sigmax += (data1[index2,0]*data1[index2,2])
           sigmaxy += (data1[index2,0]*data1[index2,1]*data1[index2,2])
           sigmaxx += (data1[index2,0]**2*data1[index2,2])
           detA = sigmax**2-sigma*sigmaxx 
           if not(detA==0):
               x = [-sigmaxx*sigmay + sigmax*sigmaxy,-sigmax*sigmay+sigma*sigmaxy]/detA
               chicurr = chisq(x,dat)**1.1/((index2-index1)**2*dat[:,1].sum())
               if (chicurr < chimin) and (chicurr > 0) and x[0] > dat[:,1].mean() and x[1]/x[0] > 1.5e-3:
                   xret = x
                   chimin = chicurr
                   imin = index1
                   imax = index2
    if chimin < chimax:
        n = (imax-imin)
        xmean = data1[imin:imax,0].mean() 
        ymean = data1[imin:imax,1].mean() 
        ssxx = (data1[imin:imax,0]**2).sum() -n*xmean**2
        ssyy = (data1[imin:imax,1]**2).sum() -n*ymean**2
        ssxy = (data1[imin:imax,1]*data1[imin:imax,0]).sum() -n*ymean*xmean
        s = ((ssyy + xret[0]*ssxy)/(n-2))**0.5
        da = s*(1/n + (xmean**2)/ssxx)**0.5
        db = s/ssxx**0.5
        a = xret[0]
        b = xret[1]
        xopt = (a,(3*b/a)**0.5)
        dx = (da,0.5*(3*(b/a**(3)*da**2 + db**2/(a*b)))**0.5)
        return xopt, dx, chimin/sigscale**2, imin, imax
    else: 
        raise FitError()
 
 
def runautorg(filename,outputname):
    if os.name =='nt':
        commandname = [op.normpath("C:/atsas/bin/autorg")]
    else:
        commandname = ["autorg", "--PXSOFT_VERSION", "v2.5.1" ]  
    try:   
        sp.call(commandname + ["-o",outputname, "-f","ssv",filename])
    except:
        print("Error running autorg")
    try:    
        f = open(outputname,'r')
        index = range(6)
        for line in f:
            words = line.split(None, 8)
            if len(words) < 8:
                break
            Rg = float(words[index[0]])
            dRg = float(words[index[1]])*Rg
            I0 = float(words[index[2]])
            dI0 = float(words[index[3]])*I0
            imin = int(words[index[4]])
            imax = int(words[index[5]])
            return (Rg,dRg,I0,dI0,imin,imax)
    except:
        print("Error reading autorg results")
        
def funcdevpres2(x,data):
    a = x[0]
    b = x[1]
    dat2 = data[:,0]
    ep = np.exp(-b*dat2)
    deva = -2*(data[:,2])*(data[:,1]-a*ep)*ep
    return np.array([deva.sum(),-a*np.einsum('i,i',deva,dat2)])

def funcpres2(x,data):
    a = x[0]
    b = x[1]
    dat2 = data[:,0]
    ep = np.exp(-b*dat2)
    arr2 = data[:,1]-a*ep
    return np.einsum('i,i,i',data[:,2],arr2,arr2)
    
        
def blockmin2(data):
    datloc = data[data[:,1]>0.1*data[:,1].max()]
    indmax = datloc.shape[0]
    minint = 10
    amax = 10*data[:,1].max()
    dq = 0.005
    dind = max(1,int(floor(dq/((data[1,0])-(data[0,0])))))
    chimin = 10000
    chimax = chimin
    datpresloc = np.array([data[0:indmax,0]**2,data[0:indmax,1],1/data[0:indmax,2]**2]).transpose()
    for index1 in range(1,indmax-minint,dind):
        for index2 in range(index1+9,indmax,dind):
           datpres = datpresloc[index1:index2,:]
           x0 = (datpres[:,1].max(),(datpres[0,1]-datpres[-1,1])/(datpres[0,0]-datpres[-1,0])) 
           amin = datpres[:,1].min()
           bmin = 0
           bmax = 10/(datpres[-1,0]*amin)
           bmax = None
           x = scop.minimize(funcpres2,x0,args = [datpres], method = 'L-BFGS-B',  jac = funcdevpres2, tol = 0.001, bounds = ((amin,amax),(bmin,bmax)))
           chicurr = x.fun
           if (chicurr < chimin) and (chicurr > 0):
               xret = [x.x[0],sqrt(3*x.x[1])]
               chimin = chicurr
               imin = index1
               imax = index2
    if chimin < chimax:
        return xret, chimin, imin, imax
    else: 
        raise FitError()        

if __name__ == "__main__":
        
    if os.name =='nt':
        syspath = op.normpath("N:/")
    else:
        syspath = op.normpath('/users/brennich/')
    datapath = op.normpath(op.join(syspath,'noisydata'))
    
    os.chdir(datapath)
    
    filenames = glob.iglob("*.dat")
    
    
    for filename in filenames:
        atsasfit =  runautorg(filename,'temp')
        print(atsasfit)
        try:
            data = np.genfromtxt(filename)
            #parabolafit = blocklin(data)
            #print parabolafit
            fullfit = blockmin2(data)
            print( fullfit)
            fig1 = pl.figure(filename)
            ax1 = fig1.add_subplot(1,1,1)
            pl1 = ax1.plot(data[0:4*fullfit[3],0],data[0:4*fullfit[3],1],'ob')
            #pl2 = ax1.plot(data[0:4*parabolafit[4],0],parabolafit[0][0]-parabolafit[0][1]**2*parabolafit[0][0]/3*data[0:4*parabolafit[4],0]**2,'-g')
            pl3 = ax1.plot(data[0:4*fullfit[3],0],atsasfit[2]*np.exp(-((atsasfit[0]*data[0:4*fullfit[3],0])**2)/3),'-r')
            pl4 = ax1.plot(data[0:4*fullfit[3],0],fullfit[0][0]*np.exp(-((fullfit[0][1]*data[0:4*fullfit[3],0])**2)/3),'-g')
            #ax1.set_yscale('log')
            #fig1.show()
        except:
            pass