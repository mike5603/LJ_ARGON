__author__ = 'mhumbert'
# This file models 864 argon molecules initially in an FCC configuration with distributed velocities

import array
import random
import math
import time
import os

class LJ_ARGON:
    # constants used in calculation
    N = 864 # number of atoms
    sigma = 3.4e-10 # Lennard Jones parameter, m
    kb = 1.38e-23 # Boltzmann constant, J/K
    M = 39.95*1.6747e-27 # mass per atom (Argon)
    epsilon = 120*kb # depth of potential well, J
    L = 10.229*sigma # length of box
    R = 2.25*sigma # maximum radius of interactions
    R2 = R**2 # maximum radius of interaction squared
    ncell = int(math.ceil((N/4)**(1.0/3.0))) # number of fc unit cells in box
    a = L/ncell # length of unit cell for fcc
    a2 = a/2 # half of the unit cell length
    nstep = 2000 # number of time steps
    temp = 90 # initial temperature
    dt = 1e-14 # time step, seconds
    count = 1 # count of timesteps
    simtemp = 0 # simulation temperature
    vacf = 0 # velocity autocorrelation function
    vacf1 = 0 # velocity autocorrelation function of first step to normalize
    sumvx = 0 # placeholder for later calculation
    sumvy = 0 # placeholder for later calculation
    sumvz = 0 # placeholder for later calculation
    dr = sigma/10 # thickness of shell in pair distrobution function
    maxr = 5*sigma # the maximum radius for pair distrobution function
    npair = int(math.ceil(maxr/dr)) # number of shells in pair distrobution function
    
    # creates position arrays
    xpositions = array.array('f')
    ypositions = array.array('f')
    zpositions = array.array('f')

    # creates velocity arrays    
    xvelocities = array.array('f')
    yvelocities = array.array('f')
    zvelocities = array.array('f')
    
    # creates force arrays
    xforces = array.array('f')
    yforces = array.array('f')
    zforces = array.array('f')
    
    # creates velocity arrays
    initialxvelocities = array.array('f')
    initialyvelocities = array.array('f')
    initialzvelocities = array.array('f')
    
    # creates arrays for pair-distrobution function
    n = array.array('f')
    g = array.array('f')
    
    # creates arrays for file writing
    temperatures = array.array('f')
    velacf = array.array('f')

    def __init__(self): # initialize method including initial positions and velocities
        
        try:
            os.remove("argon.xyz")
        except OSError:
            pass        
        
        try:
            os.remove("temp.csv")
        except OSError:
            pass
        
        try:
            os.remove("time.csv")
        except OSError:
            pass
        
        try:
            os.remove("vacf.csv")
        except OSError:
            pass
        
        try:
            os.remove("radius.csv")
        except OSError:
            pass
        
        try:
            os.remove("pairdistrobution.csv")
        except OSError:
            pass
        
        
        
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
            self.initialxvelocities.append(0)
            self.initialyvelocities.append(0)
            self.initialzvelocities.append(0)
        self.initialposition()
        self.initialvelocities()
        
        for i in range(0,self.npair):
            self.n.append(0)
            self.g.append(0)
        
    def initialposition(self):
        #assigns initial postions of atoms to FCC structure
        particle = 0
        
        # assigns simple cubic positions         
        for x in range(0,self.ncell):
            for y in range(0,self.ncell):
                for z in range(0,self.ncell):
                    self.xpositions[particle] = x*self.a
                    self.ypositions[particle] = y*self.a
                    self.zpositions[particle] = z*self.a
                    particle += 1
            
            # assigns atoms on x face
            for y in range(0,self.ncell):
                for z in range(0,self.ncell):
                    self.xpositions[particle] = x*self.a
                    self.ypositions[particle] = y*self.a + self.a2
                    self.zpositions[particle] = z*self.a + self.a2
                    particle += 1
                    
        # assigns atoms on y face
        for x in range(0,self.ncell):
            for y in range(0,self.ncell):
                for z in range(0,self.ncell):
                    self.xpositions[particle] = x*self.a + self.a2
                    self.ypositions[particle] = y*self.a
                    self.zpositions[particle] = z*self.a + self.a2
                    particle += 1
                    
            # assigns atoms on z face
            for y in range(0,self.ncell):
                for z in range(0,self.ncell):
                    self.xpositions[particle] = x*self.a + self.a2
                    self.ypositions[particle] = y*self.a + self.a2
                    self.zpositions[particle] = z*self.a
                    particle += 1
                    
        self.writetoxyz()
        
    def initialvelocities(self):
        # assigns initial velocities to atoms according to a Boltzmann distrobution
        normdist = array.array('f')
        mean_velocity = math.sqrt(self.kb*self.temp/self.M)
        
        # creating a gaussian distrobution around 0 with std. dev. of mean velocity   
        for i in range(0,3*self.N):
            normdist.append(random.gauss(0,1))
            normdist[i] *= mean_velocity
        
        # assigning atoms velocities from the gaussian distrobution
        for atom in range(0, self.N):
            self.xvelocities[atom] = normdist[3*atom]
            self.yvelocities[atom] = normdist[3*atom+1]
            self.zvelocities[atom] = normdist[3*atom+2]
            self.sumvx += self.xvelocities[atom]
            self.sumvy += self.yvelocities[atom]
            self.sumvz += self.zvelocities[atom]
        
        # correcting overall momentum (should be 0)
        for atom in range(0,self.N):
            self.xvelocities[atom] -= self.sumvx/self.N
            self.yvelocities[atom] -= self.sumvy/self.N
            self.zvelocities[atom] -= self.sumvz/self.N
        
        #setting initial velocity to be same as atom velocity
        for atom in range(0, self.N):
            self.initialxvelocities[atom] = normdist[3*atom]
            self.initialyvelocities[atom] = normdist[3*atom+1]
            self.initialzvelocities[atom] = normdist[3*atom+2]
            
        for atom in range(0,self.N):
            self.initialxvelocities[atom] -= self.sumvx/self.N
            self.initialyvelocities[atom] -= self.sumvy/self.N
            self.initialzvelocities[atom] -= self.sumvz/self.N
            
    def timestep(self):
        # main time step that calls functions to perform posistion adjustments
        # for each time step
        for step in range(0, self.nstep):
            self.updateforces()
            self.updatevelocities()
            self.updatepositions()
            self.temperature()
            self.velocityautocorrelation()
            self.temprecalibration()
            self.writetoxyz()
            print("--------------------Completed Step Number " + str(self.count) + "--------------------")       
            self.count += 1

    def updateforces(self):
        # calculates the forces acting on each atom over one timestep
        for atom in range(0,self.N):
            self.xforces[atom] = 0
            self.yforces[atom] = 0
            self.zforces[atom] = 0
        
        # calculating the distance between each pair of atoms        
        for atom1 in range(0,self.N-1):
            for atom2 in range(atom1+1,self.N):
                dx = self.xpositions[atom1]-self.xpositions[atom2]
                dy = self.ypositions[atom1]-self.ypositions[atom2]
                dz = self.zpositions[atom1]-self.zpositions[atom2] 
                
                # making sure we use the closest image                
                dx -= self.L*round(dx/self.L)
                dy -= self.L*round(dy/self.L)
                dz -= self.L*round(dz/self.L)
                
                r2 = dx**2 + dy**2 + dz**2
                
                # if atom is within range, calculating Lennard-Jones force                
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
        # updates the current velocity based on the previous velocity and the 
        # forces acting on the atoms
        for atom in range(0,self.N):
            self.xvelocities[atom] += self.xforces[atom]/self.M*self.dt
            self.yvelocities[atom] += self.yforces[atom]/self.M*self.dt
            self.zvelocities[atom] += self.zforces[atom]/self.M*self.dt
            
            
    def updatepositions(self):
        # updates the position of each atom based on the previous position
        # and the current velocity of the atom
        for atom in range(0,self.N):
            self.xpositions[atom] += self.xvelocities[atom]*self.dt
            self.ypositions[atom] += self.yvelocities[atom]*self.dt
            self.zpositions[atom] += self.zvelocities[atom]*self.dt
            
            # implementing periodic boundary conditions            
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
        # finds the current temperature of the system based on velocities
       sumv2 = 0
       for atom in range(0,self.N):
           sumv2 += self.xvelocities[atom]**2 + self.yvelocities[atom]**2 + self.zvelocities[atom]**2
       
       self.simtemp = self.M/3/self.N/self.kb*sumv2
       print("TEMP: " + str(self.simtemp))
       self.temperatures.append(self.simtemp)
       
    def temprecalibration(self):
        # adjusts the velocities of the atoms to adjust the temperature to 
        # match the wanted temperature
       if self.count > 15:
           if self.simtemp > self.temp + 10 or self.simtemp < self.temp-10:
               print("temperature recalibration")
               for atom in range(0,self.N):
                   self.xvelocities[atom] *=math.sqrt(self.temp/self.simtemp)
                   self.yvelocities[atom] *=math.sqrt(self.temp/self.simtemp)
                   self.zvelocities[atom] *=math.sqrt(self.temp/self.simtemp)
                   
    def velocityautocorrelation(self):
        # computes the average dot product of velocity, should tend toward 0
        # normalized by the correlation for the first time step
        if self.count == 1: # calculating for first step to normalize
            sumvdot = 0
            for atom in range(0,self.N):
                sumvdot += self.xvelocities[atom]*self.initialxvelocities[atom]
                sumvdot += self.yvelocities[atom]*self.initialyvelocities[atom]
                sumvdot += self.zvelocities[atom]*self.initialzvelocities[atom]
            self.vacf1 = sumvdot/self.N
            
        sumvdot = 0 # reset each timestep
        for atom in range(0,self.N): # calculate function for timestep
            sumvdot += self.xvelocities[atom]*self.initialxvelocities[atom]
            sumvdot += self.yvelocities[atom]*self.initialyvelocities[atom]
            sumvdot += self.zvelocities[atom]*self.initialzvelocities[atom]
        self.vacf = (sumvdot/self.N)/self.vacf1
        print("Velocity Autocorrelation Function: " + str(self.vacf))
        self.velacf.append(self.vacf)
    
    def pairdistrobutionfunction(self):
        # calculates the number of atoms in a shell around the central atom
        # averaged over all atoms
              
        for atom1 in range(0,self.N-1):
            for atom2 in range(atom1,self.N):
                dx = self.xpositions[atom1]-self.xpositions[atom2]
                dy = self.ypositions[atom1]-self.ypositions[atom2]
                dz = self.zpositions[atom1]-self.zpositions[atom2]
                
                dx -= self.L*round(dx/self.L)
                dy -= self.L*round(dy/self.L)
                dz -= self.L*round(dz/self.L)
                
                r2 = dx**2 + dy**2 + dz**2
            
                for radius in range(0,self.npair):
                    if r2 > (radius*self.dr)**2 and r2 < ((radius+1)*(self.dr))**2:
                        self.n[radius] += 1
                    
        for radius in range(1,self.npair):
            self.g[radius] = 2*self.L**3/self.N**2*self.n[radius]/4/math.pi/(radius*self.dr)**2/self.dr
            
        print("r                        g(r)")
        
        for radius in range(0,self.npair):
            print(str(radius*self.dr) + "                " + str(self.g[radius]))
            
    def writetoxyz(self):
        xyz = open("argon.xyz", "a")
        xyz.write(str(self.N) + "\n")
        xyz.write("positions of argon atom for timestep " + str(self.count) + "\n")
        for atom in range(0,self.N):
            xyz.write("Ar " + str(self.xpositions[atom]) + " " + str(self.ypositions[atom]) + " " + str(self.zpositions[atom]) + "\n")
        xyz.close()
        
    def writetemp(self):
        tempfile = open("temp.csv", "a")
        for entry in range(0,self.nstep):
            tempfile.write(str(self.temperatures[entry]) + "\n")
        tempfile.close()
        
    def writevacf(self):
        vacffile = open("vacf.csv", "a")
        for entry in range(0,self.nstep):
            vacffile.write(str(self.velacf[entry]) + "\n")
        vacffile.close()
            
    def writetime(self):
        timefile = open("time.csv", "a")
        for entry in range(0,self.nstep):
            timefile.write(str(entry*self.dt)+ "\n")
        timefile.close()
            
    def writeradius(self):
        radiusfile = open("radius.csv", "a")
        for radius in range(0,self.npair):
            radiusfile.write(str(radius*self.dr) + "\n")
        radiusfile.close()
            
    def writepairdistrobution(self):
        gfile = open("pairdistrobution.csv", "a")
        for radius in range(0,self.npair):
            gfile.write(str(self.g[radius]) + "\n")
        gfile.close()