
# coding: utf-8

# In[164]:

#!/usr/bin/env python
# Plots data in binary dumps.
# Usage: python plot.py nl-shwave.0024.bin

from numpy import *
import sys
from pylab import figure, axes, pie, title, show
from matplotlib import use
use('Agg')
from matplotlib import pyplot

#get_ipython().magic('matplotlib inline')


# 

# In[165]:

#
# Read binary file
#
istep=234
sstep=str(istep)
imdir='../../bin/'
imname='rt'
inname=imdir+imname+sstep+'.bin'

if istep>9:
	if istep>99:
		if istep>999:
			ssstep=str(istep)
		else:
			ssstep='0'+str(istep)
	else:
		ssstep='00'+str(istep)
else:
	ssstep='000'+str(istep)
sinname=imdir+imname+ssstep+'.bin'
print(sinname)
try:
  file = open('../../bin/rt.0234.bin','rb')
#  file = open(sinname,'rb')
except:
  print('Usage: ./read.py <binary_dump>')
  sys.exit()

file.seek(0,2)
eof = file.tell()
file.seek(0,0)

coordsys = fromfile(file,dtype=int32,count=1)[0]

#nx,ny,nz = fromfile(file,dtype=int32,count=7)[:3]

ndata = fromfile(file,dtype=int32,count=7)[:7]
nx= ndata[0]
ny= ndata[1]
nz= ndata[2]


gamma1,cs = fromfile(file,dtype=float,count=2)

t,dt = fromfile(file,dtype=float,count=2)

x = fromfile(file,dtype=float,count=nx)
y = fromfile(file,dtype=float,count=ny)
z = fromfile(file,dtype=float,count=nz)

shape = (nz,ny,nx)
count = prod(shape)

rho = fromfile(file,dtype=float,count=count).reshape(shape)
m1 = fromfile(file,dtype=float,count=count).reshape(shape)
m2 = fromfile(file,dtype=float,count=count).reshape(shape)
m3 = fromfile(file,dtype=float,count=count).reshape(shape)
e = fromfile(file,dtype=float,count=count).reshape(shape)
b1 = fromfile(file,dtype=float,count=count).reshape(shape)
b2 = fromfile(file,dtype=float,count=count).reshape(shape)
b3 = fromfile(file,dtype=float,count=count).reshape(shape)

if file.tell() != eof: print('Error: Too few bytes read.')

file.close()


# In[166]:

print(shape)
print(nx,ny,nz)
print(ndata)
print(t,dt)
print(gamma1,cs)
print(coordsys)
print(count)


# In[167]:
lay=59
print(z[0],z[1],z[63],z[119])
print(y[0],y[1],y[63],y[119])
print(x[0],x[1],x[63],x[119])
rhosec=rho[:,:,63]
m3sec=m3[:,:,63]
print(rhosec.shape)
print(rhosec[lay][0],rhosec[lay][30],rhosec[lay][63],rhosec[lay][118],rhosec[lay][119])
m7=m3sec/rhosec

v2sec=m3sec/rhosec
print(v2sec[lay][0],v2sec[lay][30],v2sec[lay][63],v2sec[lay][118],v2sec[lay][119])
#print(rho[5][0],rho[5][30],rho[5][63],rho[5][119])
#print(rho[lay][0],rho[lay][1],rho[lay][30],rho[lay][63],rho[lay][117],rho[lay][118],rho[lay][119])
#print(m7[lay][0],m7[lay][1],m7[lay][30],m7[lay][63],m7[lay][117],m7[lay][118],m7[lay][119])
#print(m7[:,59])
#print(rho[:,59])
#print(v2sec[:,59])


# In[168]:

#
# Plot density field
#
fig1,ax1 = pyplot.subplots(1,2,num=1)

ax1[0].imshow(rhosec)

ax1[1].imshow(v2sec)


#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg

#imgplot = plt.imshow(rho)
#imgplot.set_cmap('spectral')

#plt.plot(rho[96,:])
#print rho[96][:]

ffig=imdir+imname+sstep+'.png'
fig1.savefig(ffig)

#fig1.savefig('fig1.png',dpi=fig1.get_dpi())

