# -*- coding: utf-8 -*-

import time

from LJ_ARGON import LJ_ARGON
from fcc_positions import FCC_Positions

tic = time.time()

simulation = LJ_ARGON()


simulation.timestep()
simulation.pairdistrobutionfunction()
simulation.writetemp()
simulation.writevacf()
simulation.writetime()
simulation.writeradius()
simulation.writepairdistrobution()
toc = time.time()

print("length of run: " + str((toc-tic)/60) + " minutes")