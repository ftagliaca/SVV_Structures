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
