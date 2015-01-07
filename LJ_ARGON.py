__author__ = 'mhumbert'
# This file models 864 argon molecules initially in an FCC configuration with distributed velocities

import array
import random
import math
import time


class LJ_ARGON:
    # constants used in calculation
    N = 2 # number of atoms
    sigma = 3.4e-10 # Lennard Jones parameter, m
    kb = 1.38e-23 # Boltzmann constant, J/K
    M = 39.95*1.6747e-27 # mass per atom (Argon)
    epsilon = 120*kb # depth of potential well, J
    L = 10.229*sigma # length of box
    R = 2.25*sigma # maximum radius of interactions
    R2 = R**2 # maximum radius of interaction squared
    a = L/6 # length of unit cell for fcc
    a2 = a/2 # half of the unit cell length
    nstep = 500000 # number of time steps
    temp = 90 # initial temperature
    dt = 1e-14 # time step, seconds
    count = 1
    simtemp = 0 # simulation temperature
    vacf = 0 #velocity autocorrelation function
    vacf1 = 0 #velocity autocorrelation function of first step to normalize
    
    xpositions = array.array('f')
    ypositions = array.array('f')
    zpositions = array.array('f')

    xvelocities = array.array('f')
    yvelocities = array.array('f')
    zvelocities = array.array('f')
    
    xforces = array.array('f')
    yforces = array.array('f')
    zforces = array.array('f')
    

    def __init__(self):
        for i in range(0, self.N):
            self.xpositions.append(0)
            self.ypositions.append(0)
            self.zpositions.append(0)
            self.xvelocities.append(0)
            self.yvelocities.append(0)
            self.zvelocities.append(0)
            self.xforces.append(0)
            self.yforces.append(0)
            self.zforces.append(0)
        #self.initialposition()
        self.initialvelocities()
        self.xpositions[1]= self.sigma
        
            
    def initialposition(self):
        #assigns initial postions of atoms
        particle = 0
        for x in range(0,6):
            for y in range(0,6):
                for z in range(0,6):
                    self.xpositions[particle] = x*self.a
                    self.ypositions[particle] = y*self.a
                    self.zpositions[particle] = z*self.a
                    particle += 1
            
            for y in range(0,6):
                for z in range(0,6):
                    self.xpositions[particle] = x*self.a
                    self.ypositions[particle] = y*self.a + self.a2
                    self.zpositions[particle] = z*self.a + self.a2
                    particle += 1
                    
        for x in range(0,6):
            for y in range(0,6):
                for z in range(0,6):
                    self.xpositions[particle] = x*self.a + self.a2
                    self.ypositions[particle] = y*self.a
                    self.zpositions[particle] = z*self.a + self.a2
                    particle += 1
                    
            for y in range(0,6):
                for z in range(0,6):
                    self.xpositions[particle] = x*self.a + self.a2
                    self.ypositions[particle] = y*self.a + self.a2
                    self.zpositions[particle] = z*self.a
                    particle += 1
        
                        
    def initialvelocities(self):
        # assigns initial velocities to atoms according to a Boltzmann distrobution
        normdist = array.array('f')
        mean_velocity = math.sqrt(self.kb*self.temp/self.M)
        
        for i in range(0,3*self.N):
            normdist.append(random.gauss(0,1))
            normdist[i] *= mean_velocity
            
        for atom in range(0, self.N):
            self.xvelocities[atom] = normdist[3*atom]
            self.yvelocities[atom] = normdist[3*atom+1]
            self.zvelocities[atom] = normdist[3*atom+2]
            
        self.initialxvelocities = self.xvelocities
        self.initialyvelocities = self.yvelocities
        self.initialzvelocities = self.zvelocities
            
    def timestep(self):
        for step in range(0, self.nstep):
            self.updateforces()
            self.updatevelocities()
            self.updatepositions()
            print("Step Number " + str(self.count))            
            #print('position ' + str(self.xpositions[238]))
            #print('velocity ' +str(self.xvelocities[238]))
            #print('force ' +str(self.xforces[238]))
            self.temperature()
            self.temprecalibration()
            self.velocityautocorrelation()
            self.count += 1

    def updateforces(self):
        for atom in range(0,self.N):
            self.xforces[atom] = 0
            self.yforces[atom] = 0
            self.zforces[atom] = 0
        
        for atom1 in range(0,self.N-1):
            for atom2 in range(atom1+1,self.N):
                dx = self.xpositions[atom1]-self.xpositions[atom2]
                dy = self.ypositions[atom1]-self.ypositions[atom2]
                dz = self.zpositions[atom1]-self.zpositions[atom2]
                
                dx -= self.L*round(dx/self.L)
                dy -= self.L*round(dy/self.L)
                dz -= self.L*round(dz/self.L)
                
                r2 = dx**2 + dy**2 + dz**2
                print("radius: " +str(math.sqrt(r2)))
                
                if r2 < self.R2:
                    sr2 = (self.sigma**2)/r2
                    sr6 = sr2**3
                    force = 48*self.epsilon*sr6*(sr6-0.5)/r2
                    
                    self.xforces[atom1] += force*dx
                    self.xforces[atom2] -= force*dx
                    self.yforces[atom1] += force*dy
                    self.yforces[atom2] -= force*dy
                    self.zforces[atom1] += force*dz
                    self.zforces[atom2] -= force*dz
    def updatevelocities(self):
        for atom in range(0,self.N):
            self.xvelocities[atom] += self.xforces[atom]/self.M*self.dt
            self.yvelocities[atom] += self.yforces[atom]/self.M*self.dt
            self.zvelocities[atom] += self.zforces[atom]/self.M*self.dt
            
    def updatepositions(self):
        for atom in range(0,self.N):
            self.xpositions[atom] += self.xvelocities[atom]*self.dt
            self.ypositions[atom] += self.yvelocities[atom]*self.dt
            self.zpositions[atom] += self.zvelocities[atom]*self.dt
            
            if self.xpositions[atom] < 0:
                self.xpositions[atom] += self.L
            elif self.xpositions[atom] > self.L:
                self.xpositions[atom] -= self.L
                
            if self.ypositions[atom] < 0:
                self.ypositions[atom] += self.L
            elif self.ypositions[atom] > self.L:
                self.ypositions[atom] -= self.L
                
            if self.zpositions[atom] < 0:
                self.zpositions[atom] += self.L
            elif self.zpositions[atom] > self.L:
                self.zpositions[atom] -= self.L
            
    def temperature(self):
       sumv2 = 0
       for atom in range(0,self.N):
           sumv2 += self.xvelocities[atom]**2 + self.yvelocities[atom]**2 + self.zvelocities[atom]**2
       
       self.simtemp = self.M/3/self.N/self.kb*sumv2
       print("TEMP: " + str(self.simtemp))
       
    def temprecalibration(self):
       if self.count > 125:
           if self.simtemp > 100.0 or self.simtemp<80:
               print("temperature recalibration")
               for atom in range(0,self.N):
                   self.xvelocities[atom] *=math.sqrt(self.temp/self.simtemp)
                   self.yvelocities[atom] *=math.sqrt(self.temp/self.simtemp)
                   self.zvelocities[atom] *=math.sqrt(self.temp/self.simtemp)
                   
    def velocityautocorrelation(self):
        if self.count == 1:
            sumvdot = 0
            for atom in range(0,self.N):
                sumvdot += self.xvelocities[atom]*self.initialxvelocities[atom]
                sumvdot += self.yvelocities[atom]*self.initialyvelocities[atom]
                sumvdot += self.zvelocities[atom]*self.initialzvelocities[atom]
            self.vacf1 = sumvdot/self.N
            
            
        sumvdot = 0
        for atom in range(0,self.N):
            sumvdot += self.xvelocities[atom]*self.initialxvelocities[atom]
            sumvdot += self.yvelocities[atom]*self.initialyvelocities[atom]
            sumvdot += self.zvelocities[atom]*self.initialzvelocities[atom]
        self.vacf = (sumvdot/self.N)/self.vacf1
        print("Velocity Autocorrelation Function: " + str(self.vacf))
        cvacf = self.vacf *math.sqrt(self.temp/self.simtemp)
        print("Corrected Autocorrelation Function: " + str(cvacf))