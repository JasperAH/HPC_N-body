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
n_outputs = 1000
sim_years = 1
n_timesteps = 1 #(60*60*24*365.25*sim_years)/dt
cur_timestep = 1

k_means_interval = 1
k = 3
k_means_Tuple = collections.namedtuple('kTuple', 'i k')

# We define our Tuple type here with fields m:'mass' x:'x' 'y' 'z' v:'vx' 'vy' 'vz'
Tuple = collections.namedtuple('Tuple', 'i j kn')


# START functions used for error checking / verification
def verify(isInit = False):
    global system_am,system_e,Time
    global cur_timestep
    am = angular_momentum(Time[0,1])
    am_mag = np.sqrt(am.dot(am))
    e = total_energy(Time[0,1])
    if(isInit):
        system_am = am_mag
        system_e = e
    else:
        am_err = abs(-1+am_mag/system_am)
        e_err = abs(-1+e/system_e)
        print("initial AM:",system_am,"initial E:",system_e)
        print("current AM:",am_mag,"current E:",e)
        print("Angular Momentum error:",am_err,"Energy error:",e_err)

def total_energy(t):
    global n_bodies, G
    kin = 0
    pot = 0 

    for i in range(n_bodies): 
        kin = kin + M[i]*((V[i,t][0]**2)+(V[i,t][1]**2)+(V[i,t][2]**2))/2

    for i in range(n_bodies):
        for j in range(i):
            r = np.subtract(X[i,t],X[j,t])
            d = math.sqrt(np.dot(r,r))
            pot = pot + -1*G*M[i]*M[j]*d

    return kin+pot

def angular_momentum(t):
    global M,X,V
    global n_bodies
    am = np.array([0,0,0])

    for i in range(n_bodies):
        cross = np.cross(X[i,t],V[i,t])
        am = am + M[i]*cross
    
    return am

# END error checking


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
    global T, Kx, Kv, Done, Time, Printed
    global N, n_bodies
    global file 
    global M, X, V
    global K, k
    global n_timesteps, dt, n_outputs
    dt = sim_years/(n_outputs) # with 10 years and 1000 outputs this remains 0.01
    print("dt",dt)
    n_timesteps = int(2*math.pi*sim_years/dt) #should be more like: (60*60*24*365.25*sim_years)/dt?
    print("n_timesteps",n_timesteps)

    inputdata = np.loadtxt(file) #  n rows, and columns that are 'mass' 'x' 'y' 'z' 'vx' 'vy' 'vz'
    N = [[row[0],np.array([row[1],row[2],row[3]]),np.array([row[4],row[5],row[6]])] for row in inputdata]
    n_bodies = len(N)
    T = []

    Kx = {}
    Kv = {}

    M = {}
    X = {}
    V = {}

    Done = {}
    Time = {}
    Printed = {}

    for i in range(len(N)):
        for j in range(len(N)):
            if i is not j:
                for kn in range(1,6): # [1,5]
                    for t in range(1,n_timesteps+1): # [1,n_timesteps]
                        Kv[i,j,kn,t] = 0
                        Kx[i,kn,t] = 0
                        Done[i,j,kn,t] = False
                        Printed[i,t] = False
                        X[i,t] = 0
                        V[i,t] = 0
                    T.append(Tuple(i,j,kn))
                Time[i,j] = 1
        filedrop(0, i, N[i][1], N[i][2])
        Printed[i,0] = True 
        M[i] = N[i][0]*(15000000000) 
        X[i,cur_timestep] = N[i][1]
        V[i,cur_timestep] = N[i][2]
    print("bodies initialized")
    for i in range(-k, 0):
        V[i,cur_timestep] = np.zeros(3)
        for j in range(len(N)):
            if i is not j:
                for kn in range(1,6): # [1,5]
                    for t in range(1,n_timesteps+1): # [1,n_timesteps]
                        Kv[i,j,kn,t] = 0
                        Kx[i,kn,t] = 0
                        Done[i,j,kn,t] = False
                        Printed[i,t] = False
                    T.append(Tuple(i,j,kn))
                Time[i,j] = 1
        Printed[i,0] = True 

#dataplotting function 
def filedrop(t, i, x, v): #cur_timestep,t.i,X[t.i,t.t+1],V[t.i,t.t+1]
    file_name = "particle_" + str(i+1) + ".dat"
    location = "output_data/"
    white = "        "
    data = "{:10.10f}".format(t) + white + "{:10.10f}".format(x[0]) + white + "{:10.10f}".format(x[1]) + white + "{:10.10f}".format(x[2]) + white + "{:10.10f}".format(v[0]) + white + "{:10.10f}".format(v[1]) + white + "{:10.10f}".format(v[2])
    if t >= 1:
        outF = open(location+file_name, "a")
        # write line to output file
        outF.write(data)
        outF.write("\n")
        outF.close()
    else:
        outF = open(location+file_name, "w")
        # write line to output file
        outF.write(data)
        outF.write("\n")
        outF.close()

def sumJs(i,kn,t):
    global n_bodies, k, K
    global Kx, Kv
    return_value = 0

    if i < 0:
        for j in range(n_bodies):
            if (i is not j) and K[i, Time[i,j]] != K[j, Time[i,j]]:
                return_value = return_value + Kv[i,j,kn,t]
    else:
        for j in range(-k, 0):
            if (i is not j) and K[i, Time[i,j]] != K[j, Time[i,j]]:
                return_value = return_value + Kv[i,j,kn,t]
        for j in range(0, n_bodies):
            if (i is not j) and K[i, Time[i,j]] == K[j, Time[i,j]]:
                return_value = return_value + Kv[i,j,kn,t]

    return return_value

def doneJs(i,kn,t):
    global Done, n_bodies, k, K
    if i<0:
        for j in range(n_bodies):
            if (i is not j) and K[i, Time[i,j]] != K[j, Time[i,j]]:
                if not Done[i,j,kn,t]:
                    return False
    else:
        for j in range(-k, 0):
            if (i is not j) and K[i, Time[i,j]] != K[j, Time[i,j]]:
                if not Done[i,j,kn,t]:
                    return False
        for j in range(0, n_bodies):
            if (i is not j) and K[i, Time[i,j]] == K[j, Time[i,j]]:
                if not Done[i,j,kn,t]:
                    return False
    return True

def printedJs(t):
    global Printed, n_bodies
    for j in range(0,n_bodies): 
        if not Printed[j,t]:
            return False
    return True

######
#   SERIAL CODE CONDITIONS AND BODIES
###
sc_cond = []
sc_body = []

###
# Serial code 1
##

def rcond1(t):
    global Done, Time, K
    '''Implement the condition for serial code 1, return a boolean
        Only tuples that are not yet done and were printed in the
        previous timestep can be executed as a baseline. Next to
        this, this enforces body-body interactions within a cluster
        and body-cluster interactions when a body is in a different
        cluster.'''
    res = ((not Done[t.i,t.j,t.kn,Time[t.i,t.j]]) and (
            ((t.j < 0 or t.i < 0) and not (t.j < 0 and t.i < 0) and not (K[t.i,Time[t.i,t.j]] == K[t.j,Time[t.i,t.j]])) 
            or ((t.j >= 0 and t.i >= 0) and K[t.i,Time[t.i,t.j]] == K[t.j,Time[t.i,t.j]])
        ) and Printed[t.j,Time[t.i,t.j]-1])
    return res

def rSC1body(t):
    global T, Time, Done, Printed
    global cur_timestep, n_timesteps, N
    global kmeans_to_print
    global K_size, K_m, K_x
    global Kx, Kv
    global doneCounter
    '''The actual serial code to execute if cond1 is True
        Executes RK4 steps, aggregates result (as kn=5) for 
        next timestep. Determines when tuple is actually done
        and when results are final such that they can be
        output to a file. Initiates K-means step at interval.'''
    temp_i = [M[t.i],X[t.i,Time[t.i,t.j]],V[t.i,Time[t.i,t.j]]]
    temp_j = [M[t.j],X[t.j,Time[t.i,t.j]],0]
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    
    Kv[t.i,t.j,1,Time[t.i,t.j]] = a*dt
    Kx[t.i,1,Time[t.i,t.j]] = V[t.i,Time[t.i,t.j]] * dt
    

    temp_i = [M[t.i],X[t.i,Time[t.i,t.j]],V[t.i,Time[t.i,t.j]]]
    temp_j = [M[t.j],X[t.j,Time[t.i,t.j]],0]
    temp_i[1] = temp_i[1] + Kx[t.i,1,Time[t.i,t.j]]/2
    temp_j[1] = temp_j[1] + Kx[t.j,1,Time[t.i,t.j]]/2
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    
    Kv[t.i,t.j,2,Time[t.i,t.j]] = a*dt
    kv1 = sumJs(t.i,1,Time[t.i,t.j])
    Kx[t.i,2,Time[t.i,t.j]] = (V[t.i,Time[t.i,t.j]] + kv1/2) * dt


    temp_i = [M[t.i],X[t.i,Time[t.i,t.j]],V[t.i,Time[t.i,t.j]]]
    temp_j = [M[t.j],X[t.j,Time[t.i,t.j]],0]
    temp_i[1] = temp_i[1] + Kx[t.i,2,Time[t.i,t.j]]/2
    temp_j[1] = temp_j[1] + Kx[t.j,2,Time[t.i,t.j]]/2
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    
    Kv[t.i,t.j,3,Time[t.i,t.j]] = a*dt
    kv2 = sumJs(t.i, 2, Time[t.i,t.j])
    Kx[t.i,3,Time[t.i,t.j]] = (V[t.i,Time[t.i,t.j]] + kv2/2) * dt


    temp_i = [M[t.i],X[t.i,Time[t.i,t.j]],V[t.i,Time[t.i,t.j]]]
    temp_j = [M[t.j],X[t.j,Time[t.i,t.j]],0]
    temp_i[1] = temp_i[1] + Kx[t.i,3,Time[t.i,t.j]]
    temp_j[1] = temp_j[1] + Kx[t.j,3,Time[t.i,t.j]]
    a = f_gravitational(temp_i,temp_j)/M[t.i]
    
    Kv[t.i,t.j,4,Time[t.i,t.j]] = a*dt
    kv3 = sumJs(t.i, 3, Time[t.i,t.j])
    Kx[t.i,4,Time[t.i,t.j]] = (V[t.i,Time[t.i,t.j]] + kv3) * dt
    
    
    X[t.i,Time[t.i,t.j]+1] = X[t.i,Time[t.i,t.j]] + (Kx[t.i,1,Time[t.i,t.j]] + 2*Kx[t.i,2,Time[t.i,t.j]] + 2*Kx[t.i,3,Time[t.i,t.j]] + Kx[t.i,4,Time[t.i,t.j]])/6
    kv4 = sumJs(t.i,4,Time[t.i,t.j])
    V[t.i,Time[t.i,t.j]+1] = V[t.i,Time[t.i,t.j]] + (kv1 + 2*kv2 + 2*kv3 + kv4)/6

    if((t.kn == 1 and (Time[t.i,t.j] == 1 or printedJs(Time[t.i,t.j]-1))) 
        or (t.kn > 1 and doneJs(t.i,t.kn-1,Time[t.i,t.j]))):

        Done[t.i,t.j,t.kn,Time[t.i,t.j]] = True
        doneCounter += 1 # count successfully executed tuples for experiments
        
        if(t.kn == 5 and not Printed[t.i,Time[t.i,t.j]]):
            Printed[t.i,Time[t.i,t.j]] = True
            if(t.i == 0 and Time[t.i,t.j] % (n_timesteps/n_outputs) < 1.0): print("{:0.1f}%".format((Time[t.i,t.j]/n_timesteps)*100))
            if(t.i == 0 and Time[t.i,t.j] % k_means_interval == 0): kmeans_doTimestep(Time[t.i,t.j]) # do kmeans timestep on results of previous timestep because this definitely has all bodies done
            if Time[t.i,t.j] % (n_timesteps/n_outputs) < 1.0: #division should produce float, no cast required (from __future__)
                if t.i >= 0:
                    filedrop(Time[t.i,t.j],t.i,X[t.i,Time[t.i,t.j]+1],V[t.i,Time[t.i,t.j]+1])
                if(t.i == 0): kmeans_to_print += 1
            if(Time[t.i,t.j] < n_timesteps):
                if t.i < 0:
                    for j in range(n_bodies):
                        Time[t.i,j] = Time[t.i,j] + 1
                else:
                    for j in range(-k,n_bodies):
                        if t.i is not j:
                            Time[t.i,j] = Time[t.i,j] + 1

def kmeans_doTimestep(time):
    global K, K_size, k, k_means_interval, k_means_Tuple, T_K, K_x, K_m
    global N, n_bodies
    global cur_timestep
    # To update clusters, new body values need to be used
    # to make sure cluster ID stays consistent, 
    # use existing body/cluster membership to initialize before executing 
    for _k in range(1,k+1):
        K_size[_k,time] = 0
        K_m[_k,time] = 0
        K_x[_k,time] = np.array([0,0,0])
    for i in range(n_bodies):
        cluster = K[i,cur_timestep]
        K[i,time] = cluster
        K_x[cluster,time] = (K_x[cluster,time]*K_m[cluster,time] + X[i,time]*M[i])/(K_m[cluster,time]+M[i])
        K_m[cluster,time] = K_m[cluster,time] + M[i]
        K_size[cluster,time] = K_size[cluster,time] + 1

    cur_timestep = time
    kmeans_loop()


    for i in range(n_bodies):
        for _interval in range(k_means_interval+2): #some headroom for async execution
            cluster = K[i,cur_timestep]
            K[i,cur_timestep + _interval] = cluster
            K_m[cluster,cur_timestep + _interval] = K_m[cluster,cur_timestep]
            K_x[cluster,cur_timestep + _interval] = K_x[cluster,cur_timestep]
            K_size[cluster,cur_timestep + _interval] = K_size[cluster,cur_timestep]



def kmeans_dist(a, b):
    r = np.subtract(a,b)
    d = math.sqrt(np.dot(r,r))
    return d

def kmeans_init():
    global K, K_size, k, k_means_interval, k_means_Tuple, T_K, K_x, K_m, kmeans_to_print
    global N, n_bodies
    global M, X, V, T
    global Kx, Kv, cur_timestep
    # k means:T
    T_K = []
    K = {} # address using K[i,Time[i,j]] or something like that
    K_size = {} # address using K_size[k,K[i,Time[i,j]]] or K_size[k,0] where number is time
    K_x = {} #idem
    K_m = {} #idem

    kmeans_to_print = 0

    for _k in range(1,k+1):
        K_size[_k,cur_timestep] = 0
        K_m[_k,cur_timestep] = 0
        K_x[_k,cur_timestep] = np.array([0,0,0])

    cluster = 1
    for i in range(n_bodies):
        for t in range(1,n_timesteps+1):
            K[i,t] = 0
        for _k in range(1,k+1):
            T_K.append(k_means_Tuple(i,_k))

        if cluster > k: # initialize remaining bodies to closest cluster
            clus = random.randint(1,k)
            K[i,cur_timestep] = clus
            K_x[clus,cur_timestep] = (K_x[clus,cur_timestep]*K_m[clus,cur_timestep] + X[i,cur_timestep]*M[i])/(K_m[clus,cur_timestep]+M[i])
            K_m[clus,cur_timestep] = K_m[clus,cur_timestep] + M[i]
            K_size[clus,cur_timestep] = K_size[clus,cur_timestep] + 1  
        else: # initialize first k bodies to the k clusters such that clusters are not empty
            K[i,cur_timestep] = cluster 
            K_x[cluster,cur_timestep] = (K_x[cluster,cur_timestep]*K_m[cluster,cur_timestep] + X[i,cur_timestep]*M[i])/(K_m[cluster,cur_timestep]+M[i])
            K_m[cluster,cur_timestep] = K_m[cluster,cur_timestep] + M[i]
            K_size[cluster,cur_timestep] = K_size[cluster,cur_timestep] + 1
            cluster = cluster + 1


    for i in range(n_bodies):
        for _k in range(1,k+1):
            for kn in range(1, 6):
                Time[i,-_k] = 1
                for t in range(1,n_timesteps+1): 
                    Kv[i,-_k,kn,t] = 0
                    Kx[-_k,kn,t] = 0
                    Done[i,-_k,kn,t] = False
                    K[-_k, t] = _k
                T.append(Tuple(i,-_k,kn))
    



def kmeans_cond(t):
    global K, K_size, k, k_means_interval, k_means_Tuple, T_K, K_x, K_m
    global M, X, V
    # return true if:
    # t.k is not the current cluster of i
    # and
    # the distance to t.k from i is less than the distance between i and its current cluster
    return (K[t.i,cur_timestep] != t.k and 
        (kmeans_dist(X[t.i,cur_timestep],K_x[t.k,cur_timestep]) < kmeans_dist(X[t.i,cur_timestep],K_x[K[t.i,cur_timestep],cur_timestep])) and
        K_size[K[t.i,cur_timestep],cur_timestep] > 1) #divide by zero protection

def kmeans_body(t):
    global K, K_size, k, k_means_interval, k_means_Tuple, T_K, K_x, K_m
    global M, X, V
    K_x[K[t.i,cur_timestep],cur_timestep] = (K_x[K[t.i,cur_timestep],cur_timestep]*K_m[K[t.i,cur_timestep],cur_timestep] - X[t.i,cur_timestep]*M[t.i])/(K_m[K[t.i,cur_timestep],cur_timestep] - M[t.i])
    K_m[K[t.i,cur_timestep],cur_timestep] = K_m[K[t.i,cur_timestep],cur_timestep] - M[t.i]
    K_size[K[t.i,cur_timestep],cur_timestep] = K_size[K[t.i,cur_timestep],cur_timestep]-1

    K_x[t.k,cur_timestep] = (K_x[t.k,cur_timestep]*K_m[t.k,cur_timestep] + X[t.i,cur_timestep]*M[t.i])/(K_m[t.k,cur_timestep] + M[t.i])
    K_m[t.k,cur_timestep] = K_m[t.k,cur_timestep] + M[t.i]
    K_size[t.k,cur_timestep] = K_size[t.k,cur_timestep]+1

    K[t.i,cur_timestep] = t.k


def kmeans_loop():
    global K, K_size, k, k_means_interval, k_means_Tuple, T_K, K_x, K_m, kmeans_to_print
    global N, n_bodies
    global M, X, V
    done = False
    while not done:
        tmp = list(T_K)
        random.shuffle(tmp)
        size = len(T_K)
        for t in tmp:
            if kmeans_cond(t):
                kmeans_body(t)

        valid_left = False
        for t in tmp:
            valid_left |= kmeans_cond(t)
        done = not valid_left
    if(cur_timestep == 1): kmeans_to_print = 1
    for _k in range(1,k+1):
        for _ in range(kmeans_to_print):
            filedrop(cur_timestep-1,n_bodies+_k-1,K_x[_k,cur_timestep],[0, 0, 0])
    kmeans_to_print = 0

    for _k in range(1,k+1):
        M[-_k] = K_m[_k,cur_timestep]
        X[-_k,cur_timestep] = K_x[_k, cur_timestep]


    if(cur_timestep == 1): 
        for _interval in range(k_means_interval+2): #some headroom for async execution
            for i in range(n_bodies):
                cluster = K[i,cur_timestep]
                K[i,cur_timestep + _interval] = cluster
                K_m[cluster,cur_timestep + _interval] = K_m[cluster,cur_timestep]
                K_x[cluster,cur_timestep + _interval] = K_x[cluster,cur_timestep]
                K_size[cluster,cur_timestep + _interval] = K_size[cluster,cur_timestep]     


###
# Loop mechanics
##

def body(t):
    global T
    ''' Return True when modification has been made '''
    if rcond1(t) == True:
        rSC1body(t)
    return


def checkConditions():
    '''Check if there is any tuple left for which at least one condition
    is True.'''
    tmp = False
    for t in T:
        tmp |= rcond1(t)

    return tmp

def loop():
    global T
    '''Attempt at whilelem loop.'''
    done = True
    while done:
        tmp = list(T)
        random.shuffle(tmp)
        size = len(T)
        for t in tmp:
            body(t)
        done = checkConditions()


def initRandom():
    seed = 0
    for i, c in enumerate(bytearray(os.urandom(4))):
        seed += c << i * 8
    random.seed(seed)

    print("Initialized random number generator with seed: 0x{:x}\n".format(seed))


if __name__ == '__main__':
    global file, doneCounter
    if(len(sys.argv) < 2):
        print("Usage:",sys.argv[0],"<input.dat> <duration (years)>")
        sys.exit()
    
    file = sys.argv[1] # n rows, and columns that are 'mass' 'x' 'y' 'z' 'vx' 'vy' 'vz'
    sim_years = int(sys.argv[2]) 

    doneCounter = 0

    initRandom()

    init()

    kmeans_init()
    kmeans_loop()

    verify(isInit=True)

    loop()

    kmeans_loop()

    verify()

    print("number of successfully executed tuples",doneCounter)

