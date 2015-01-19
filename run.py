# -*- coding: utf-8 -*-

import time
import math

from LJ_ARGON import LJ_ARGON

tic = time.time()

simulation = LJ_ARGON()


simulation.timestep()
simulation.pairdistributionfunction()
simulation.writetemp()
simulation.writevacf()
simulation.writetime()
simulation.writeradius()
simulation.writepairdistribution()
toc = time.time()

hours = math.floor((toc - tic)/3600)
minutes = math.floor((toc-tic-3600*hours)/60)
seconds = (toc-tic-3600*hours - 60*minutes)

print("length of run: " + str(hours) + " hours " + str(minutes) + " minutes " + str(seconds) + " seconds")