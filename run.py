# -*- coding: utf-8 -*-

from LJ_ARGON import LJ_ARGON
from fcc_positions import FCC_Positions

simulation = LJ_ARGON()
#fcc = FCC_Positions()

nsteps = 100

simulation.timestep()
#fcc.initialpositions()