import math
import numpy as np
from matplotlib import pyplot as plt

class Aileron():
    def __init__(self, C_a, l_a, x_1, x_2, x_3, x_a, h, t_sk, t_sp, t_st, h_st, w_st, n_st, d_1, d_3, theta, P):
        self.C_a = C_a  #m
        self.l_a = l_a #m
        self.x_1 = x_1 #m
        self.x_2 = x_2 #m
        self.x_3 = x_3 #m
        self.x_a = x_a*0.01 #m
        self.h = h*0.01 #m
        self.t_sk = t_sk*0.001 #m
        self.t_sp = t_sp*0.001 #m
        self.t_st = t_st*0.001 #m
        self.h_st = h_st*0.01 #m
        self.w_st = w_st*0.01 #m
        self.n_st = n_st #(number of stringers)
        self.d_1 = d_1*0.01 #m
        self.d_3 = d_3*0.01 #m
        self.theta = theta #deg
        self.P = P*1000 #N

    def stringersPosition(self):

        r = 0.5*self.h
        l_s = math.sqrt((self.C_a-r)**2 + r**2)

        perimeter = math.pi*r + 2*l_s #perimeter of the aileron
        d_st = perimeter/self.n_st #distance between stiffners

        self.st_pos = np.zeros([self.n_st, 2])

        for i in range(self.n_st):

            d = i*d_st
            alpha = d/r
            delta_d = d - math.pi*r/2

            if alpha <= math.pi/2:
                self.st_pos[i,0] = r*math.sin(alpha)
                self.st_pos[i,1] = -r + r*math.cos(alpha)
                print("A")

            elif i <= self.n_st/2:
                self.st_pos[i,0] = (l_s-delta_d)*r/l_s
                self.st_pos[i,1] = -r-delta_d*(self.C_a-r)/l_s
                print("B")

            else:
                self.st_pos[i,0] = - self.st_pos[int(self.n_st) - i,0]
                self.st_pos[i,1] = self.st_pos[int(self.n_st) - i,1]
                print("C")



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
            z_i = self.st_pos[i,1] #z coordinate of ith stringer
            Q_stiff += z_i *A_stiff    #+= is the same as Q_stiff + .....



        z_centroid = (Q_stiff + c_halfcircle*A_halfcircle + c_skin * A_skin + c_spar * A_spar)\
        /(A_halfcircle + A_skin + A_spar + self.n_st * A_stif)
        return(z_centroid)
    def z_i(self, i, N_z = 81):
        '''
        Inputs:
        i = ith row in spanwise direction
        N_z = number of rows (optional)

        Output:
        z_i = z-coordinate of station
        '''

        theta = (i-1)*math.pi/N_z
        theta_1 = i*math.pi/N_z
        z_i = -0.5*(0.5*self.C_a*(1-math.cos(theta))+0.5*self.C_a*(1-math.cos(theta_1)))

        return z_i

    def x_i(self, i, N_x = 41):
        '''
        Inputs:
        i = ith column in chordwise direction
        N_x = number of column (optional)

        Output:
        x_i = x-coordinate of station
        '''

        theta = (i-1)*math.pi/N_x
        theta_1 = i*math.pi/N_x
        x_i = 0.5*(0.5*self.l_a*(1-math.cos(theta))+0.5*self.l_a*(1-math.cos(theta_1)))

        return x_i


    def moi(self):
        pass

A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
A320.stringersPosition()
print(A320.st_pos)
plt.plot(-A320.st_pos[:,1], A320.st_pos[:,0])
plt.show()
