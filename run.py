# -*- coding: utf-8 -*-

import time

from LJ_ARGON import LJ_ARGON
from fcc_positions import FCC_Positions

tic = time.time()

simulation = LJ_ARGON()


#simulation.timestep()
#simulation.pairdistrobutionfunction()

toc = time.time()

#print(str(toc-tic))