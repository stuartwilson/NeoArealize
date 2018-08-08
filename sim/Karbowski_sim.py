import os

steps = 100000 # Seems enough to reach steady state
fname = 'data/tmp2.bin'
os.system ('./Karbowski_model ' + str(steps) + ' ' + fname)
