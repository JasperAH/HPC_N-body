#!/usr/bin/env python

#
# tuplsimple.py: Simple tUPL framework.
#
#     The simple tUPL framework is suitable for modeling simple whilelem
#     loops with no, or a single, condition.
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


def algo():
    '''Here the algorithm is implemented by means of a simple, single
    whilelem loop.'''
    global T, S

    changed = True

    # We make a copy to preserve the order of the original tuple reservoir.
    tmp = list(T)

    # A whilelem loop is run until the state of the program no longer changes.
    while changed:
        changed = False

        # Execute loop body for each tuple in random order.
        # For testing it is sometimes useful to temporarily disable the
        # call to "shuffle".
        random.shuffle(tmp)
        for t in tmp:
            # TODO: write loop body that operates on "t"
            # If the state of the shared spaces changes: set changed to
            # True.
            print(t)

            pass


def initRandom():
    seed = 0
    for i, c in enumerate(bytearray(os.urandom(4))):
        seed += c << i * 8

    random.seed(seed)

    print("Initialized random number generator with seed: 0x{:x}\n".format(seed))


if __name__ == '__main__':
    initRandom()

    init()
    algo()
