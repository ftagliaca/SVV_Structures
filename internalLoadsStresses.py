import math

def normalStress(y, z, Aileron, M_z, M_x):
    '''
    Input:

    y = y-coordinate
    z = z-coordinate
    Aileron = aileron class containing geometrical and material properties
    M_z = moment around z-axis
    M_x = moment around x-axis

    Output:

    sigma_z = normal stress along z axis, in Pa
    '''
    sigma_z = (M_z*Aileron.Iyy*y + M_x*Aileron.Izz*z)/(Aileron.Izz*Aileron.Iyy)
    return sigma_z

def vonMises(sigma, tau):
    '''
    Input:

    sigma = list with len = 3 containing [sigma_xx, sigma_yy, sigma_zz]
    tau = list with len = 3 containing [tau_xy, tau_xz, tau_yz]

    Output:

    sigma_vm = Von Mises stress
    '''
    sigma_vm = math.sqrt(((sigma[0]-sigma[1])**2 + (sigma[1]-sigma[2])**2 + (sigma[0]-sigma[2])**2)/2\
    3*(tau[0]**2+tau[1]**2+tau[2]**2))
    return sigma_vm
