import os
import numpy as np
import pylab as pl

steps = 500000
fname = 'data/tmp2.bin'
if(1):
    os.system('rm '+fname)
    os.system('./Karbowski_model '+str(steps)+' '+fname)

k = 14 # number of arrays stored
x = np.fromfile(fname)
x = np.reshape(x,[k,len(x)/k]).T
A = x[:,:5]
C = x[:,5:10]
P = x[:,10:13]
H = x[:,13]

L = np.tile(np.linspace(0,40,(40+0.25)/0.25),[5,1]).T

F = pl.figure(figsize=(10,8))
f = F.add_subplot(211)
f.plot(L,A)
f.set_ylabel('A(x)')
f.set_xlabel('x')
f.set_xlim(0,40)
f = F.add_subplot(212)
f.plot(L,C)
f.set_ylabel('C(x)')
f.set_xlabel('x')
f.set_xlim(0,40)
f.set_ylim(0,1)

F.savefig('data/figure2.pdf')

os.system('evince data/figure2.pdf')
#pl.show()
