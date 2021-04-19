#!/usr/bin/env python

#
# tuplcomplex.py: More complex tUPL framework.
#
#     The more complex tUPL framework is suitable for modeling whilelem
#     loops with multiple conditions and serial codes. In this case, when
#     multiple conditions are True, one of the serial codes is to be
#     chosen arbitrarily. The loop terminates when no tuple exists for
#     which one of the conditions is true.
#
# Copyright (C) Harry Wijshoff, Kristian Rietveld
#

# For our Python3 friends.
from __future__ import print_function, division, absolute_import, unicode_literals

import collections

import os
import random


# We define our Tuple type here with fields "u" and "v".
Tuple2 = collections.namedtuple('Tuple2', 'u v')
Tuple3 = collections.namedtuple('Tuple3', 'm x y')

# gravitational constant
G = 6.67408/100000000000 # 6.67408e-11
# timestep size (arbitrary for now)
dt = 0.01
sim_years = 1
n_timesteps = (60*60*24*365.25*sim_years)/dt
cur_timestep = 1

# calculate force j exerts on i for each dimension
def F_gravitational(i, j):
	global G
	r = i.x-j.x
	d = sqrt(x.x**2 + x.y**2 + x.z**2) #TODO this line should not be compatible yet
	u = r/d
	f = (-1*G*i.m*j.m)/(d**2)
	return u*f





def init():
    '''Initialize some tuple reservoirs (e.g. T) and necessary shared
    spaces (e.g. S). Shared spaces are usually modeled with a dictionary
    that is indexed with tuples.'''
    global T, S, N, Kx, Kv, K_mask, T_mask

    #TODO delete this, from example
    T = []
    T.append(Tuple2(1, 1))
    S = {}

    #TODO load data from some file into tuples
    N = []
    N.append(Tuple3(1,(10,10,10),(20,1,1)) #N.append(Tuple3(m,X,V))

    # K_mask is a mask for n,k,t combinations, counting the number of modifications comparing to len(N).
    Kx={}
    Kv={}
    K_mask={}

    #TODO is this a valid shared space?
    T_mask=0


###
# Serial code 1
##

def cond1(t):
    '''Implement the condition for serial code 1, return a boolean'''
    if t.i != t.j and t.kn == 1:
        return True:
    else:
        return False

def SC1body(t):
    '''The actual serial code to execute if cond1 is True'''
    a = F_gravitational(t.i,t.j)/t.i.m 	# vector3 a contains x, y and z components
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = t.i.v * dt
    K_mask[t.i,t.kn,cur_timestep] = K_mask[t.i,t.kn,cur_timestep] + 1


###
# Serial code 2
#

def cond2(t):
    '''Implement the condition for serial code 2, return a boolean'''
    if t.i != t.j and t.kn == 2 and K_mask[t.i,t.kn-1] == len(N) and K_mask[t.j,t.kn-1] == len(N):
        return True
    else:
        return False

def SC2body(t):
    '''The actual serial code to execute if cond2 is True'''
    temp_i = t.i
    temp_j = t.j
    temp_i.x = temp_i.x + Kx[t.i,t.kn-1]/2
    temp_j.x = temp_j.x + Kx[t.j,t.kn-1]/2
    a = F_gravitational(temp_i,temp_j)/t.i.m
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = (t.i.v + Kv[t.i,t.kn-1]/2) * dt
    K_mask[t.i,t.kn,cur_timestep] = K_mask[t.i,t.kn,cur_timestep] + 1


###
# Serial code 3
#

def cond3(t):
    '''Implement the condition for serial code 3, return a boolean'''
    if t.i != t.j and t.kn == 3 and K_mask[t.i,t.kn-1] == len(N) and K_mask[t.j,t.kn-1] == len(N):
        return True
    else:
        return False

def SC3body(t):
    '''The actual serial code to execute if cond3 is True'''
    temp_i = t.i
    temp_j = t.j
    temp_i.x = temp_i.x + Kx[t.i,t.kn-1]/2
    temp_j.x = temp_j.x + Kx[t.j,t.kn-1]/2
    a = F_gravitational(temp_i,temp_j)/t.i.m
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = (t.i.v + Kv[t.i,t.kn-1]/2) * dt
    K_mask[t.i,t.kn,cur_timestep] = K_mask[t.i,t.kn,cur_timestep] + 1


###
# Serial code 4
#

def cond4(t):
    '''Implement the condition for serial code 4, return a boolean'''
    if t.i != t.j and t.kn == 4 and K_mask[t.i,t.kn-1] == len(N) and K_mask[t.j,t.kn-1] == len(N):
        return True
    else:
        return False

def SC4body(t):
    '''The actual serial code to execute if cond4 is True'''
    temp_i = t.i
    temp_j = t.j
    temp_i.x = temp_i.x + Kx[t.i,t.kn-1]
    temp_j.x = temp_j.x + Kx[t.j,t.kn-1]
    a = F_gravitational(temp_i,temp_j)/t.i.m
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = (t.i.v + Kv[t.i,t.kn-1]) * dt
    K_mask[t.i,t.kn,cur_timestep] = K_mask[t.i,t.kn,cur_timestep] + 1


###
# Serial code 5
#

def cond5(t):
    '''Implement the condition for serial code 5, return a boolean'''
    if t.kn == 5 and K_mask[t.i,t.kn,cur_timestep] == 0 and K_mask[t.i,t.kn-1] == len(N) and K_mask[t.j,t.kn-1] == len(N) and T_mask < len(N)
        return True
    else:
        return False

def SC5body(t):
    '''The actual serial code to execute if cond5 is True'''
    t.i.x = t.i.x + (Kx[t.i,1] + 2*Kx[t.i,2] + 2*Kx[t.i,3] + Kx[t.i,4])/6
    t.i.v = t.i.v + (Kv[t.i,1] + 2*Kv[t.i,2] + 2*Kv[t.i,3] + Kv[t.i,4])/6 
    Kx[t.i,:] = 0
    Kv[t.i,:] = 0
    K_mask[t.i,t.kn,cur_timestep] = 1
    T_mask = T_mask + 1 #TODO become global?


###
# Serial code 6
#

def cond6(t):
    '''Implement the condition for serial code 6, return a boolean'''
    if T_mask == len(N):
        return True
    else
        return False

def SC6body(t):
    '''The actual serial code to execute if cond6 is True'''
    cur_timestep = cur_timestep + 1 #TODO become global?
    T_mask = 0 #TODO become global?



###
# Loop mechanics
##

def body(t):
    ''' Return True when modification has been made '''
    #assert t in T 
    #TODO above does not work, but list comprehension could be done in init, such that this will work

    # FIXME: of course this can be made more generic in case you add more
    # than 2 serial codes.

    c1 = cond1(t)
    c2 = cond2(t)
    c3 = cond3(t)
    c4 = cond4(t)
    c5 = cond5(t)
    c6 = cond6(t)

    if c1 and c2:
        # Randomly select a serial code to run.
        choice = random.randint(1,2)
        if choice == 1:
            SC1body(t)
            return
        if choice == 2:
            SC2body(t)
            return
    elif c1:
        SC1body(t)
        return
    elif c2:
        SC2body(t)
        return
    assert "Shouldn't be reached"

def checkConditions():
    '''Check if there is any tuple left for which at least one condition
    is True.'''
    tmp = False
    for t in T:
        tmp |= cond1(t)
        tmp |= cond2(t)

    return tmp

def checkConditionsForelem(t):
    '''Check if there is any tuple left for which at least one condition
    is True.'''
    tmp = False
    tmp |= cond1(t)
    tmp |= cond2(t)
    tmp |= cond3(t)
    tmp |= cond4(t)
    tmp |= cond5(t)
    tmp |= cond6(t)
    return tmp


def loop():
    '''Implementation of the actual forelem loop.'''
    #TODO at this moment, in at least the python implementation, possibly also tUPL, randomly assigning a timestep t and then checking it against cur_timestep does not 
    # make any sense. As such, timestep t has been deleted from the tuplereservoir


    #TODO This is not officially a forelem, as no guarantees are given on execution at least once
    tmp = list(N)
    tmp2 = List(N)
    random.shuffle(tmp)
    random.shuffle(tmp2)
    tuples = [(i,j,kn) for i in tmp for j in tmp2 for kn in [1,5]] #List comprehension
    for t in tuples:
        if checkConditionsForelem(t):
            body(t)

def initRandom():
    seed = 0
    for i, c in enumerate(bytearray(os.urandom(4))):
        seed += c << i * 8

    random.seed(seed)

    print("Initialized random number generator with seed: 0x{:x}\n".format(seed))


if __name__ == '__main__':
    initRandom()

    init()
    loop()
