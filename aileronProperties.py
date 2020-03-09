from aero_loads import AerodynamicLoad
from matplotlib import pyplot as plt
import os
import numpy as np
import math

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
        self.theta = math.radians(theta) #rad
        self.P = P*1000 #N

        self.x_I = self.x_2 - 0.5*self.x_a
        self.x_II = self.x_2 + 0.5*self.x_a

        #Material properties obtained from http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T3

        self.E = 73.1e9
        self.G = 28e9

        self.stringersPosition()
        self.zCentroid()
        self.momInertia()
        self.interpolate_grid()

    def stringersPosition(self):
        '''
        Inputs:
        Only the class itself is necessary

        Outputs:
        st_pos = A 2D numpy array containing the locations of
                 the stiffners.
                 1st column is y-coordinate, 2nd column is z-coordinate
        '''
        r = 0.5*self.h #radius of the circular section
        self.r = r
        self.l_s = math.sqrt((self.C_a-r)**2 + r**2) #length of the straight skin section

        perimeter = math.pi*r + 2*self.l_s #perimeter of the aileron
        d_st = perimeter/self.n_st #distance between stiffners

        self.st_pos = np.zeros([self.n_st, 2]) #create list of stiffners locations

        for i in range(self.n_st):

            d = i*d_st
            alpha = d/r
            delta_d = d - math.pi*r/2

            if alpha <= math.pi/2:
                self.st_pos[i,0] = r*math.sin(alpha)
                self.st_pos[i,1] = -r + r*math.cos(alpha)

            elif i <= self.n_st/2:
                self.st_pos[i,0] = (self.l_s-delta_d)*r/self.l_s
                self.st_pos[i,1] = -r-delta_d*(self.C_a-r)/self.l_s

            else:
                self.st_pos[i,0] = - self.st_pos[int(self.n_st) - i,0]
                self.st_pos[i,1] = self.st_pos[int(self.n_st) - i,1]

        return self.st_pos


    def zCentroid(self):
    #be careful with units. did not look at this yet
        r = self.h / 2

        #now calculate centroids of separate parts
        c_halfcircle = - (r - 2*r/math.pi)
        c_skin = - 0.5*(self.C_a - r) - r
        c_spar = - r

        #the length of one of the diagonal skins
        l_skin = math.sqrt((self.C_a-r)**2 + r**2)
        #and the areas
        A_halfcircle = math.pi *r* self.t_sk
        A_skin = 2 * l_skin * self.t_sk   #for two skins
        A_spar = self.t_sp * self.h
        A_stif = (self.w_st + self.h_st)* self.t_st

        self.Q_stiff = 0 #same as z~ * A
        for i in range(self.n_st):
            z_i = self.st_pos[i,1] #z coordinate of ith stringer
            self.Q_stiff += z_i * A_stif    #+= is the same as Q_stiff + .....

        self.z_centroid = (self.Q_stiff + c_halfcircle*A_halfcircle + c_skin * A_skin + c_spar * A_spar)\
        /(A_halfcircle + A_skin + A_spar + self.n_st * A_stif)

        self.A_stif = A_stif
        return self.z_centroid


    def crossArea(self):
        '''
        Inputs:
        Only the class itself is necessary

        Outputs:
        A = the total cross sectional area,
            but the single cross sectional areas are available
        '''

        self.A1 = math.pi*(self.r**2)/2 #Area of semicircle
        self.A2 = 0.5*self.h*(self.C_a - self.r) #Area of triagle
        self.A = self.A1 + self.A2

        return self.A

    def momInertia(self):
        #be careful with units. did not look at this yet
        r = self.h / 2


        #the length of one of the diagonal skins
        l_skin = math.sqrt((self.C_a-r)**2 + r**2)
        #and the areas
        A_halfcircle = math.pi *r* self.t_sk
        A_skin = 2 * l_skin * self.t_sk          #for both two diagonal skins
        A_spar = self.t_sp * self.h
        A_stif = (self.w_st + self.h_st)*self.t_st #for ONE stiffener

        zCentroid = np.abs(self.zCentroid())
        #I_yy calculations
        #Important to note that the dz's are chosen in accordance with points we
        #calculated the individual Moments of Inertia around
        dz_skin = -((self.C_a-r)/2+r)+zCentroid
        dz_spar = -r+zCentroid
        dz_halfcircle = -r+zCentroid

        dz_stif = (self.Q_stiff)/(self.n_st * A_stif)

        I_yy = A_halfcircle * dz_halfcircle**2 + A_skin * dz_skin**2 + A_spar * dz_spar**2

        #I_zz calculations
        dy_halfcircle = 0
        dy_spar = 0
        dy_skin = 0.5*r

        I_zz = A_halfcircle * dy_halfcircle**2 + A_skin * dy_skin**2 + A_spar * dy_spar**2

        #now we add the influence of the stringers on I_yy
        for i in range(self.n_st):
            dz_st = self.st_pos[i,1]+zCentroid
            dy_st = self.st_pos[i,0]

            I_yy += A_stif * dz_st**2
            I_zz += A_stif * dy_st**2



        #Moments of inertia of separate parts around own centroids.
        Beta = math.atan(r/(self.C_a-r))
        I_zzskin= self.t_sk*(l_skin)**3*(math.sin(Beta)**2)/12
        I_yyskin= self.t_sk*(l_skin)**3*(math.cos(Beta)**2)/12

        I_zzspar = self.t_sp * self.h**3 / 12
        #I_yyspar = self.h * self.t_sp**3 / 12
        #I_yyspar is 0, thinwalled assumption.

        I_yyhc = I_zzhc = 0.5*math.pi*self.t_sk*r**3

        I_zztot = I_zz + 2*I_zzskin + I_zzhc + I_zzspar
        I_yytot = I_yy + 2*I_yyskin + I_yyhc

        self.Izz = I_zztot
        self.Iyy = I_yytot

        return(I_yytot, I_zztot)

    def torsionalStiffness(self):
        '''
        Inputs:
        Only the class itself is necessary, requires for crossArea to run first

        Outputs:
        J1, J2 = Torsional constant of the two sections
                 (with 1 being the semicircle and 2 being the triagle)
        '''

        _ = self.crossArea()

        self.J1 = 4*(self.A1**2)/(math.pi*self.h/(2*self.t_sk) + self.h/self.t_sp)
        self.J2 = 4*(self.A2**2)/(2*self.l_s/self.t_sk)

        return self.J1, self.J2


    def z_i(self, i, N_z = 81, i1 = 0):
        '''
        Inputs:
        i = ith row in spanwise direction
        N_z = number of rows (optional)

        Output:
        z_i = z-coordinate of station
        '''

        theta = (i - 1 + i1) * math.pi/N_z
        theta_1 = i * math.pi/N_z
        z_i = -0.5 * (0.5 * self.C_a * (1 - np.cos(theta)) + 0.5 * self.C_a * (1 - np.cos(theta_1)))

        return z_i

    def x_i(self, i, N_x = 41, i1 = 0):
        '''
        Inputs:
        i = ith column in chordwise direction
        N_x = number of column (optional)

        Output:
        x_i = x-coordinate of station
        '''

        theta = (i - 1 + i1) * math.pi/N_x
        theta_1 = i * math.pi/N_x
        x_i = 0.5 * (0.5 * self.l_a * (1 - np.cos(theta)) + 0.5 * self.l_a * (1 - np.cos(theta_1)))

        return x_i

def interpolate_grid(self):
    source_data_file = 'data/aerodynamicloada320.dat'
    interpolated_data_file = 'data/a320_interpolated_loading.npy'

    if os.path.exists(interpolated_data_file):
        self.aero_loads = np.load(interpolated_data_file)
    else:
        n_x = 41 * 10
        n_z = 81 * 10

        x = self.x_i(np.arange(n_x) + 1, N_x=n_x)

        for x_val in [self.x_1, self.x_2, self.x_3, self.x_a, self.x_I, self.x_II]:
            idx = x.searchsorted(x_val)
            x = np.concatenate((x[:idx], [x_val], x[idx:]))
        
        z = self.z_i(np.arange(n_z) + 1, N_x=n_z)

        self.aero_loads = AerodynamicLoad(self, filename=source_data_file).get_values_grid(z, x)
        np.save(interpolated_data_file, self.aero_loads)



A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
A320.stringersPosition()

if __name__ == "__main__":
    print(A320.st_pos)
    print(A320.zCentroid())
    plt.plot(-A320.st_pos[:,1], A320.st_pos[:,0])
    plt.show()
