
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
#matplotlib.use()
#get_ipython().magic('matplotlib inline')


# 

# In[165]:

#
# Read binary file
#
istep=654
sstep=str(istep)
imdir='../../bin/'
imname='rt.'
inname=imdir+imname+sstep+'.bin'
for istep in range(654,678):
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
	sinname=str(imdir+imname+ssstep+'.bin')
	print(sinname)
	try:
	  file = open(sinname,'rb')
	except:
	  print('Usage: ./read.py <binary_dump>')
	  sys.exit()
	file.seek(0,2)
	eof = file.tell()
	file.seek(0,0)
        ffig=imdir+imname+ssstep+'.png'
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
	rhosec=rho[:,:,63]
	m3sec=m3[:,:,63]
	v2sec=m3sec/rhosec
	#
	# Plot density field
	#
	fig1,ax1 = pyplot.subplots()
	ax1.imshow(v2sec)
	print(ffig)
	fig1.savefig(ffig)

