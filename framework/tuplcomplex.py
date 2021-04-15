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
Tuple = collections.namedtuple('Tuple', 'u v')


def init():
    '''Initialize some tuple reservoirs (e.g. T) and necessary shared
    spaces (e.g. S). Shared spaces are usually modeled with a dictionary
    that is indexed with tuples.'''
    global T, S

    T = []
    T.append(Tuple(1, 1))

    S = {}


###
# Serial code 1
##

def cond1(t):
    '''Implement the condition for serial code 1, return a boolean'''
    return False

def SC1body(t):
    '''The actual serial code to execute if cond1 is True'''
    pass


###
# Serial code 2
#

def cond2(t):
    '''Implement the condition for serial code 2, return a boolean'''
    return False

def SC2body(t):
    '''The actual serial code to execute if cond2 is True'''
    pass


###
# Loop mechanics
##

def body(t):
    ''' Return True when modification has been made '''
    assert t in T

    # FIXME: of course this can be made more generic in case you add more
    # than 2 serial codes.

    c1 = cond1(t)
    c2 = cond2(t)

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
