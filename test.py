# !/usr/bin/python

import numpy
import matplotlib.pylab as plt
import pylab
from pylab import *

n=input('enter number of simulation')
f=input('enter frequency of simulation')
foldername=("Sim%d_%dhz" % (n,f))

print foldername

outfolderpath=('exp%dhz/Sim%d' %(f,n))
print outfolderpath

#Projlist = list()
#for projLoop in range(7):
#    print projLoop
#    
#    Projlist.append(projLoop)
#    
#for i in Projlist:
#	print i+1
#	print '\n'
    
#    array = ("Proj" + str(projLoop))
#    print array
