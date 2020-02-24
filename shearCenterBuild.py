"""

Author: Sherman Lee, Sadra Mogaddam
It's given that the formula for the shear flow at any given point from a distance is (insert formula)
Sz and Sy are taken to be nonzero so Izz and Iyy are required.
* The geometry is constant throughout
* The reference point can be set to be constant, and then the integration phase can be simpler.
* NOTE: multi-cell analysis is needed for this system.
* Thickness is constant for the skin
* The analysis for the stringers are done in the boom-skin analysis, with the skin analysed separately from the booms
* Filippo says the coordinates of the booms are provided in a 17x2 matrix, with the boom areas to be constant.
* TODO: set up framework to receive geometry values from Aileron class.
* Note that symmetry means that top and bottom parts of the
"""

import math as m
import numpy as np

# import of class in order to use geometrical properties
# Note that this only imports the class, not the geometric values of the aircraft.
from aileronProperties import Aileron


# solve for base shear flow


def get_constants(Sz, Sy, Izz, Iyy):
    '''
    Step for calculating the constants Lambda = - Sz / Iyy and Lambda = - Sy / Izz. Note Lambda is arbitrarily set, not
    equivalent to anything in literature.
    :param Sz: Shear in Z direction
    :param Sy: Shear in Y direction
    :param Izz: MoI around z axis
    :param Iyy: MoI around y axis
    :return: tuple containing both Lambdas (Lambda_z, Lambda_y)
    '''

    Lambda_z = - Sz / Iyy
    Lambda_y = - Sy / Izz

    return np.array([Lambda_z, Lambda_y])


def get_idealised_shear_flow(boom_area_array, x_y_array, num_stiffeners):
    '''
    Gets the sum of base shear flow of the idealised boom sections only. They are referenced in Sherman's derivations
    as Epsilon_z and Epsilon_y.
    :boom_area_array:
    :x_y_array:
    :return:
    '''
    if boom_area_array.size == num_stiffeners * 2:
        # initialise output sum
        Epsilon_z = 0
        Epsilon_y = 0
        for iter in range(0, num_stiffeners):
            Epsilon_z += boom_area_array[iter] * x_y_array[iter][0]
            Epsilon_y += boom_area_array[iter] * x_y_array[iter][1]

        return np.array([Epsilon_z, Epsilon_y])
    else:
        raise Exception("Something's wrong, I can feel it!")


def get_shear_flow_1(Lambda_array, t, h_spar, z_bar, z_spar, Epsilon_array):
    """
    Gets the first shear flow equation for the base shear flows in the semi-circular profile from Sherman's derived
    equation
    :param Lambda_z:
    :param Lambda_y:
    :param t:
    :param h_spar:
    :param z_bar:
    :return:
    """
    radius = h_spar / 2
    q_11 = Lambda_array[0] * (t * (radius * radius + radius * z_bar * m.pi) + Epsilon_array[0]) + Lambda_array[1] * (
            t * radius * radius) + Epsilon_array[1]
    q_12 = Lambda_array[0] * (t * h_spar * h_spar + Epsilon_array[0]) + Lambda_array[1] * (
            t * (z_spar - z_spar) * h_spar + Epsilon_array[1]) + q_11
    q_13
