# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 12:30:32 2015

@author: mhumbert
"""
import array

class FCC_Positions:
    
    N = 864
    sigma = 3.4e-10
    L = 10.229*sigma
    a = L/6
    a2 = a/2
    
    xpositions = array.array('f')
    ypositions = array.array('f')
    zpositions = array.array('f')
    
    def __init__(self):
        for i in range(0, self.N):
            self.xpositions.append(0)
            self.ypositions.append(0)
            self.zpositions.append(0)    
    
    def initialpositions(self):
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
        
        print(str(particle))