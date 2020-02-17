# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:51:28 2020

@author: thijs
Centroid calculations

input list with [y,z] locations for stringers
and of course values for other dimensions
"""
import math
import numpy as np

def zcentroid(self):
    #be careful with units. did not look at this yet
    r = self.h / 2 
    
    #now calculate centroids of separate parts
    c_halfcircle = 4/3*r*math.pi 
    c_skin = 0.5*(self.c_a - r)
    c_spar = r
    
    #the length of one of the diagonal skins
    l_skin = math.sqrt((self.c_a-r)^2 + r^2)
    #and the areas
    A_halfcircle = math.pi *r* self.t_sk
    A_skin = 2 * l_skin * self.t_sk   #for two skins
    A_spar = self.t_sp * self.h
    A_stif = self.w_st * t_st + self.h_st *t_st
    
    Q_stiff = 0 #same as z~ * A
    for i in range (self.n_st):
        z_i = stringercoordinates[i,1] #z coordinate of ith stringer
        Q_stiff += z_i *A_stiff    #+= is the same as Q_stiff + .....
        
    
    
    z_centroid = (Q_stiff + c_halfcircle*A_halfcircle + c_skin * A_skin + c_spar * A_spar)\
    /(A_halfcircle + A_skin + A_spar + n_st * A_stif)
    return(z_centroid)