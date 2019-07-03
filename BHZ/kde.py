#!/usr/bin/env python
from numpy import *
from scipy.integrate import simps
import numpy as np
import sys
if len(sys.argv) != 2:  # the program name and the two arguments
    # stop the program and print an error message
    sys.exit("Range argument is not defined")
Range=sys.argv[1]
a,b,N=Range.split(":")
a=float(a)
b=float(b)
N=int(N)
x=linspace(a,b,N)
DataFile=sys.stdin
data=abs(np.asarray(loadtxt(DataFile,unpack=True,comments='#',dtype=np.float32)))
Ldata=len(data)
hdata=(4.0/3.0/float(Ldata))**(1.0/5.0)*np.std(data,dtype=np.float32)
def gaussian_kernel(x,mean,sigma):
    return np.exp(-0.5*((x-mean)/sigma)**2)/np.sqrt(2.0*np.pi*sigma**2)
pdf = np.zeros(shape=(N))
for i in range(Ldata):
    pdf+=gaussian_kernel(x,data[i],hdata)/float(Ldata)
pdf=pdf/simps(pdf,x)
savetxt(sys.stdout,np.transpose([x,pdf]),fmt='%2.6f')
