import os
import numpy as np
import pylab as pl

steps=100000
k = 14 # number of arrays stored
fname = 'data/tmp2.bin'
for ii in range(1000,steps,1000):
    x = np.fromfile(fname+'.'+str(ii).zfill(6))
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
    F.savefig(fname+'.'+str(ii).zfill(6)+'.jpg')
    pl.close(F)

x = np.fromfile(fname+'.G')
x = np.reshape(x,[5,161]).T
g0 = x[:,0]
g1 = x[:,1]

L = np.tile(np.linspace(0,40,(40+0.25)/0.25),[5,1]).T

F = pl.figure(figsize=(10,8))
f = F.add_subplot(211)
f.plot(L,g0)
f.set_ylabel('g[0](x)')
f.set_xlabel('x')
f.set_xlim(0,40)
f = F.add_subplot(212)
f.plot(L,g1)
f.set_ylabel('g[1](x)')
f.set_xlabel('x')
f.set_xlim(0,40)
f.set_ylim(0,1)
F.savefig(fname+'.G.jpg')
pl.close(F)
