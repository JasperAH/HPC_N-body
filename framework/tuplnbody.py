#!/usr/bin/env python

#
# tuplnbody.py: N-body implementation based on tuplcomplex.py
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

import sys
import numpy as np
import math


# gravitational constant
G = 6.67408/100000000000 # 6.67408e-11
# timestep size dt (arbitrary for now)
dt = 0.01
sim_years = 1
n_timesteps = 1 #(60*60*24*365.25*sim_years)/dt
cur_timestep = 1

# We define our Tuple type here with fields m:'mass' x:'x' 'y' 'z' v:'vx' 'vy' 'vz'
Tuple = collections.namedtuple('Tuple', 'i j kn t')


# calculate force j exerts on i for each dimension
def f_gravitational(i, j): # bodies contain m:[0] x:[1] v:[2]
    global G
    r = np.subtract(i[1],j[1])
    d = math.sqrt(np.dot(r,r))
    u = np.true_divide(r,d)
    f = (-1*G*i[0]*j[0])/(d**2)
    return u*f



def init():
    '''Initialize some tuple reservoirs (e.g. T) and necessary shared
    spaces (e.g. S). Shared spaces are usually modeled with a dictionary
    that is indexed with tuples.'''
    global T, Kx, Kv, K_mask, T_mask
    global N, n_bodies
    global file 
    global M, X, V
    n_timesteps = sim_years #should be more like: (60*60*24*365.25*sim_years)/dt

    inputdata = np.loadtxt(file) #  n rows, and columns that are 'mass' 'x' 'y' 'z' 'vx' 'vy' 'vz'
    N = [[row[0],np.array([row[1],row[2],row[3]]),np.array([row[4],row[5],row[6]])] for row in inputdata]
    n_bodies = len(N)
    T = []

    Kx = {}
    Kv = {}
    K_mask = {}
    T_mask = {}

    M = {}
    X = {}
    V = {}

    T_mask = [0] * (n_timesteps+1)
    T_mask[0] = len(N)

    for i in range(len(N)):
        M[i] = N[i][0]
        X[i,cur_timestep] = N[i][1]
        V[i,cur_timestep] = N[i][2]
        for j in range(len(N)):
            if i == j:
                continue
            for kn in range(1,6): # [1,5]
                for t in range(1,n_timesteps+1): # [1,n_timesteps]
                    T.append(Tuple(i,j,kn,t))
                    K_mask[i,kn,t] = 0
                Kv[i,kn] = 0
                Kx[i,kn] = 0

######
#   SERIAL CODE CONDITIONS AND BODIES
###
sc_cond = []
sc_body = []

###
# Serial code 1
##

def cond1(t):
    '''Implement the condition for serial code 1, return a boolean'''
    return (t.i is not t.j and t.kn == 1 and t.t == cur_timestep)
            # and K_mask[t.i,t.kn,t.t] < n_bodies-1 
            # and K_mask[t.j,t.kn,t.t] < n_bodies-1)

def SC1body(t):
    '''The actual serial code to execute if cond1 is True'''
    temp_i = [M[t.i],X[t.i,t.t],V[t.i,t.t]]
    temp_j = [M[t.j],X[t.j,t.t],V[t.j,t.t]]
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = V[t.i,t.t] * dt
    K_mask[t.i,t.kn,t.t] = K_mask[t.i,t.kn,t.t] + 1

    print(cur_timestep, t.i, 1)

sc_cond.append(cond1)
sc_body.append(SC1body)

###
# Serial code 2
#

def cond2(t):
    '''Implement the condition for serial code 2, return a boolean'''
    return (t.i is not t.j and t.kn == 2 and t.t == cur_timestep 
            and K_mask[t.i,(t.kn-1),t.t] == n_bodies-1 
            and K_mask[t.j,(t.kn-1),t.t] == n_bodies-1)

def SC2body(t):
    '''The actual serial code to execute if cond2 is True'''
    temp_i = [M[t.i],X[t.i,t.t],V[t.i,t.t]]
    temp_j = [M[t.j],X[t.j,t.t],V[t.j,t.t]]
    temp_i[1] = temp_i[1] + Kx[t.i,(t.kn-1)]/2
    temp_j[1] = temp_j[1] + Kx[t.j,(t.kn-1)]/2
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = (V[t.i,t.t] + Kv[t.i,(t.kn-1)]/2) * dt
    K_mask[t.i,t.kn,t.t] = K_mask[t.i,t.kn,t.t] + 1

    print(cur_timestep, t.i, 2)

sc_cond.append(cond2)
sc_body.append(SC2body)

###
# Serial code 3
#

def cond3(t):
    return (t.i is not t.j and t.kn == 3 and t.t == cur_timestep 
            and K_mask[t.i,(t.kn-1),t.t] == n_bodies-1 
            and K_mask[t.j,(t.kn-1),t.t] == n_bodies-1)

def SC3body(t):
    temp_i = [M[t.i],X[t.i,t.t],V[t.i,t.t]]
    temp_j = [M[t.j],X[t.j,t.t],V[t.j,t.t]]
    temp_i[1] = temp_i[1] + Kx[t.i,(t.kn-1)]/2
    temp_j[1] = temp_j[1] + Kx[t.j,(t.kn-1)]/2
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = (V[t.i,t.t] + Kv[t.i,(t.kn-1)]/2) * dt
    K_mask[t.i,t.kn,t.t] = K_mask[t.i,t.kn,t.t] + 1

    print(cur_timestep, t.i, 3)

sc_cond.append(cond3)
sc_body.append(SC3body)

###
# Serial code 4
#

def cond4(t):
    return (t.i is not t.j and t.kn == 4 and t.t == cur_timestep 
            and K_mask[t.i,(t.kn-1),t.t] == n_bodies-1 
            and K_mask[t.j,(t.kn-1),t.t] == n_bodies-1)

def SC4body(t):
    temp_i = [M[t.i],X[t.i,t.t],V[t.i,t.t]]
    temp_j = [M[t.j],X[t.j,t.t],V[t.j,t.t]]
    temp_i[1] = temp_i[1] + Kx[t.i,(t.kn-1)]
    temp_j[1] = temp_j[1] + Kx[t.j,(t.kn-1)]
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    Kv[t.i,t.kn] = Kv[t.i,t.kn] + a*dt
    Kx[t.i,t.kn] = (V[t.i,t.t] + Kv[t.i,(t.kn-1)]) * dt
    K_mask[t.i,t.kn,t.t] = K_mask[t.i,t.kn,t.t] + 1

    print(cur_timestep, t.i, 4)

sc_cond.append(cond4)
sc_body.append(SC4body)

###
# Serial code 5
#

def cond5(t):
    return (t.kn == 5 and t.t == cur_timestep and K_mask[t.i,t.kn,t.t] == 0
            and K_mask[t.i,(t.kn-1),t.t] == n_bodies-1 
            and K_mask[t.j,(t.kn-1),t.t] == n_bodies-1)

def SC5body(t):
    global T_mask, cur_timestep
    # FIXME: x and v before = need to be from next timestep somehow
    # maybe fixable using another shared space to store values somehow?
    X[t.i,t.t+1] = X[t.i,t.t] + (Kx[t.i,1] + 2*Kx[t.i,2] + 2*Kx[t.i,3] + Kx[t.i,4])/6
    V[t.i,t.t+1] = V[t.i,t.t] + (Kv[t.i,1] + 2*Kv[t.i,2] + 2*Kv[t.i,3] + Kv[t.i,4])/6
    for index in range(1,5):
        Kx[t.i,index] = 0
        Kv[t.i,index] = 0
    K_mask[t.i,t.kn,t.t] = 1
    T_mask[cur_timestep] = T_mask[cur_timestep] + 1
    if T_mask[cur_timestep] == len(N):
        cur_timestep = cur_timestep + 1
    print(cur_timestep,t.i,X[t.i,t.t+1],V[t.i,t.t+1])
    
sc_cond.append(cond5)
sc_body.append(SC5body)

###
# Serial code 6
#

#def cond6(t):
#    return (T_mask == n_bodies)

#def SC6body(t):
#    global T_mask, cur_timestep
#    cur_timestep = cur_timestep + 1
#    T_mask = 0
#    print("SC6body")
#    print(t)

#sc_cond.append(cond6)
#sc_body.append(SC6body)

###
# Loop mechanics
##

def body(t):
    global T
    ''' Return True when modification has been made '''
    #assert t in T

    # FIXME: of course this can be made more generic in case you add more
    # than 2 serial codes.

    # c1 = cond1(t)
    # c2 = cond2(t)
    bodies_with_true_cond = []
    for i in range(len(sc_cond)):
        if sc_cond[i](t):
            bodies_with_true_cond.append(sc_body[i])
    
    if(len(bodies_with_true_cond) > 0): # len is always 1 or 0 I believe
        choice = random.randint(0,len(bodies_with_true_cond)-1)
        bodies_with_true_cond[choice](t)
        print("removing")
        T.remove(t) # attempt to convert whilelem to forelem
    return

    # if c1 and c2:
    #     # Randomly select a serial code to run.
    #     choice = random.randint(1,2)
    #     if choice == 1:
    #         SC1body(t)
    #         return
    #     if choice == 2:
    #         SC2body(t)
    #         return
    # elif c1:
    #     SC1body(t)
    #     return
    # elif c2:
    #     SC2body(t)
    #     return
    assert "Shouldn't be reached"

def checkConditions():
    '''Check if there is any tuple left for which at least one condition
    is True.'''
    tmp = False
    for t in T:
        for c in sc_cond:
            tmp |= c(t)
        # tmp |= cond1(t)
        # tmp |= cond2(t)

    return tmp

def loop():
    '''attempt at forelem loop. Keeps executing a shrinking T in 
    random order until size doesn't change anymore'''
    done = False
    while not done:
        tmp = list(T)
        random.shuffle(tmp)
        size = len(T)
        print("pre ",size)
        for t in tmp:
            body(t)
        print("post",len(T))
        if len(T) == size:
            print("Stuck, bad exit")
            print(T)
            return
        done = len(T) == 0 #len(T) == size

def loop_whilelem():
    '''Implementation of the actual whilelem loop.'''

    # The whilelem loop is run until there is no tuple left for which at
    # least one of the conditions is True.
    changed = True
    while changed:
        changed = False

        # Execute loop body for each tuple in random order.
        # For testing it is sometimes useful to temporarily disable the
        # call to "shuffle".
        tmp = list(T)
        random.shuffle(tmp)
        for t in tmp:
            body(t)

        # Check if there's a tuple left for which at least one of the
        # conditions is True.
        changed |= checkConditions()

def initRandom():
    seed = 0
    for i, c in enumerate(bytearray(os.urandom(4))):
        seed += c << i * 8

    seed = 0x3f93ba96 # TODO temp
    random.seed(seed)

    print("Initialized random number generator with seed: 0x{:x}\n".format(seed))


if __name__ == '__main__':
    global file
    if(len(sys.argv) < 2):
        print("Usage:",sys.argv[0],"<input.dat> <duration (years)>")
        sys.exit()
    
    file = sys.argv[1] # n rows, and columns that are 'mass' 'x' 'y' 'z' 'vx' 'vy' 'vz'
    sim_years = int(sys.argv[2])
    #output = sys.argv[3]

    initRandom()

    init()
    print(T)
    loop()
    print(T)
