#
# Testing process2.cpp
#

import sys
import time
import processes as P
from math import *
import numpy as np

basePort = 8000
#worlds = [P.pTemp('build/sim/process2','world00','logs/log00',1,basePort+0),
#          P.pTemp('build/sim/process2','world01','logs/log01',2,basePort+1)]
worlds = [P.pTemp('build/sim/process2','world00','logs/log00',1,basePort+0)]

# These are the diffusion constants
Dn = np.array([200.,10.])
Dc = Dn*0.3

# Set simulation params
for i,w in enumerate(worlds):

    # Send 2 commands to each world:
    strn = '4,0,'+str(0)+','+str((1.*(0+1.))*Dn[i])
    print (strn)
    w.stream(strn)
    strn = '4,1,'+str(0)+','+str((1.*(0+1.))*Dc[i])
    w.stream(strn)

''' # Uncomment to load
for i,w in enumerate(worlds):
    w.stream('6,'+'logs/data'+str(i)+'.bin')
'''

# MAIN SIMULATION LOOP
for t in range(1000):

    for i,w in enumerate(worlds):
        w.stream('1,')

    #for i,w in enumerate(worlds):
    #    if(t%1==0):
    #        w.stream('2,') # display image
    #        # w.stream('3,') #Save image to file

''' # Uncomment to save
for i,w in enumerate(worlds):
    w.stream('5,'+'logs/data'+str(i)+'.bin')
'''

w.quit()
