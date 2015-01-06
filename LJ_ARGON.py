__author__ = 'mhumbert'
# This file models 864 argon molecules initially in an FCC configuration with distributed velocities

import array
import random
import math
import time

class LJ_ARGON:
    # constants used in calculation
    N = 864 # number of atoms
    sigma = 3.4e-10 # Lennard Jones parameter, m
    kb = 1.38e-23 # Boltzmann constant, J/K
    M = 39.95*1.6747e-24 # mass per atom (Argon)
    epsilon = 120*kb # depth of potential well, J
    L = 10.229*sigma # length of box
    R = 2.25*sigma # maximum radius of interactions
    R2 = R**2 # maximum radius of interaction squared
    n = 10 # number of atoms per side of box in initial conditions
    nstep = 500 # number of time steps
    temp = 90 # initial temperature

    atoms = array.array('f')


    def __init__(self):
        pass
