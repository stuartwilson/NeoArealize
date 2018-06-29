import sys
import time
import processes as P
from math import *
import numpy as np

basePort = 8000
worlds = [P.pTemp('processes/sim/build','world00','logs/log00',1,basePort+0),
          P.pTemp('processes/sim/build','world01','logs/log01',2,basePort+1)]

Dn = np.array([100.,5.])
Dc = np.array([100.*0.3,5.*0.3])

# Set simulation params
for i,w in enumerate(worlds):
    w.stream('4,0,'+str(Dn[i]))
    w.stream('4,1,'+str(Dc[i]))

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
            w.stream('2,') # display image
            # w.stream('3,') #Save image to file

''' # Uncomment to save
for i,w in enumerate(worlds):
    w.stream('5,'+'logs/data'+str(i)+'.bin')
'''

w.quit()
