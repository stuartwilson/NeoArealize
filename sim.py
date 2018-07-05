import sys
import time
import processes as P
from math import *
import numpy as np

basePort = 8000
worlds = [P.pTemp('build/sim','world00','logs/log00',1,basePort+0),
          P.pTemp('build/sim','world01','logs/log01',2,basePort+1)]

# These are the diffusion constants
#Dn = np.array([100.,5.])
#Dc = np.array([100.*0.3,5.*0.3])
Dn = np.array([200.,10.])
Dc = Dn*0.3

# Set simulation params
nField = 5
for i,w in enumerate(worlds):
    for j in range(nField):
        w.stream('4,0,'+str(j)+','+str((1.*(j+1.)/nField)*Dn[i]))
        w.stream('4,1,'+str(j)+','+str((1.*(j+1.)/nField)*Dc[i]))

''' # Uncomment to load
for i,w in enumerate(worlds):
    w.stream('6,'+'logs/data'+str(i)+'.bin')
'''

# MAIN SIMULATION LOOP
for t in range(100000):

    for i,w in enumerate(worlds):
        w.stream('1,')

    for i,w in enumerate(worlds):
        if(t%1==0):
            w.stream('8,') # display image
            # w.stream('3,') #Save image to file

''' # Uncomment to save
for i,w in enumerate(worlds):
    w.stream('5,'+'logs/data'+str(i)+'.bin')
'''

w.quit()
