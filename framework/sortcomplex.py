#!/usr/bin/env python

#
# sortcomplex.py: Complex tUPL sort example.
#
# Copyright (C) Harry Wijshoff, Kristian Rietveld
#

# For our Python3 friends.
from __future__ import print_function, division, absolute_import, unicode_literals

import collections

import os
import random


N = 128
numbers = [44, 56, 31, 20, 23, 64, 1, 98, 54, 99, 50, 22, 47, 52, 46, 80, 72, 8, 100, 49, 28, 77, 88, 34, 25, 88, 37, 55, 83, 64, 49, 44, 55, 49, 79, 41, 59, 11, 90, 37, 51, 22, 4, 26, 2, 38, 61, 88, 93, 46, 43, 53, 41, 61, 73, 52, 64, 50, 99, 53, 80, 26, 80, 8, 77, 71, 91, 72, 30, 37, 67, 96, 22, 97, 92, 71, 16, 18, 1, 78, 25, 85, 14, 80, 76, 62, 6, 90, 19, 30, 55, 18, 36, 85, 57, 83, 47, 70, 48, 17, 15, 64, 62, 15, 75, 15, 32, 63, 78, 49, 88, 51, 69, 76, 19, 44, 1, 20, 34, 19, 13, 85, 20, 24, 46, 42, 53, 16 ]


Tuple = collections.namedtuple('Tuple', 'i j')


def init():
    '''Initialize some tuple reservoirs (e.g. T) and necessary shared
    spaces (e.g. S). Shared spaces are usually modeled with a dictionary
    that is indexed with tuples.'''
    global T, A

    T = []
    for i in range(N-1):
        T.append(Tuple(i, i+1))

    A = list(numbers)


###
# Serial code 1
##

def cond1(t):
    '''Implement the condition for serial code 1, return a boolean'''
    return A[t.i] > A[t.j]

def SC1body(t):
    '''The actual serial code to execute if cond1 is True'''
    A[t.i], A[t.j] = A[t.j], A[t.i]


###
# Loop mechanics
##

def body(t):
    ''' Return True when modification has been made '''
    assert t in T

    c1 = cond1(t)

    if c1:
        SC1body(t)
        return
    assert "Shouldn't be reached"

def checkConditions():
    '''Check if there is any tuple left for which at least one condition
    is True.'''
    tmp = False
    for t in T:
        tmp |= cond1(t)

    return tmp


def loop():
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

    random.seed(seed)

    print("Initialized random number generator with seed: 0x{:x}\n".format(seed))


if __name__ == '__main__':
    initRandom()

    init()
    loop()

    # Output resulting state
    print("Result:")
    for i in range(N):
        print(A[i], end=" ")
    print()
