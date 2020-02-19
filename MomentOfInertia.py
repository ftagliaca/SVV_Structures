# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:42:55 2020

@author: Thijs Bolscher and Filippo
Calculation of moments of inertia of aileron
"""

import math
import numpy as np


def momInertia(self, stringercoordinates):
    #be careful with units. did not look at this yet
    r = self.h / 2 
    

    #the length of one of the diagonal skins
    l_skin = math.sqrt((self.c_a-r)^2 + r^2)
    #and the areas
    A_halfcircle = math.pi *r* self.t_sk
    A_skin = 2 * l_skin * self.t_sk          #for both two diagonal skins
    A_spar = self.t_sp * self.h
    A_stif = self.w_st * t_st + self.h_st *t_st #for ONE stiffener
    
    
    #I_yy calculations
    #Important to note that the dz's are chosen in accordance with points we 
    #calculated the individual Moments of Inertia around
    dz_skin = -((self.c_a-r)/2+r)
    dz_spar = -r
    dz_halfcircle = -r  
           
    dz_stif = (Q_stiff)/(self.n_st * A_stif)
    
    I_yy = A_halfcircle * dz_halfcircle**2 + A_skin * dz_skin**2 + A_spar * dz_spar**2
    
    #I_zz calculations
    dy_halfcircle = 0
    dy_spar = 0
    dy_skin = 0.5*r
    
    I_zz = A_halfcircle * dy_halfcircle**2 + A_skin * dy_skin**2 + A_spar * dy_spar**2
    
    #now we add the influence of the stringers on I_yy
    for i in range(self.n_st):
        dz_st = stringercoordinates[i,1]
        dy_st = stringercoordinates[i,0]
        
        I_yy = I_yy + A_st * dz_st**2
        I_zz = I_zz + A_st * dy_st**2        
     
    
    
    #Moments of inertia of separate parts around own centroids.
    Beta = math.atan(r/(self.c_a-r))
    I_zzskin= self.t_sk*(l_skin)**(3)*math.sin(Beta)**(2)/12
    I_yyskin= self.t_sk*(l_skin)**(3)*math.cos(Beta)**(2)/12
    
    I_zzspar = self.t_sp * self.h**3 / 12
    #I_yyspar = self.h * self.t_sp**3 / 12
    #I_yyspar is 0, thinwalled assumption.
    
    I_yyhc = I_zzhc = 0.5*math.pi*self.t_sk*r**3 
    
    I_zztot = I_zz + I_zzskin + I_zzhc + I_zzspar
    I_yytot = I_yy + I_yyskin + I_yyhc
    
    return(I_yytot, I_zztot)
    
    
    
    

    
    
