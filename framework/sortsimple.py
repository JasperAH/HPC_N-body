#!/usr/bin/env python

#
# sortsimple.py: Simple tUPL sort example
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
            if A[t.i] > A[t.j]:
                A[t.i], A[t.j] = A[t.j], A[t.i]
                changed = True


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

    # Output resulting state
    print("Result:")
    for i in range(N):
        print(A[i], end=" ")
    print()
