"""

Author: Sherman Lee
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
#

from aileronProperties import Aileron, A320

analysed_aircraft = A320


# solve for base shear flow

def z_ii(aircraft_class):
    """
    gets z_II, the z length from spar to trailing edge of the aileron (also z-length of section II).
    :param aircraft_class: This uses the given geometric values from the aileron
    :return: length of z_II.
    """
    z_2 = aircraft_class.C_a - aircraft_class.h
    return z_2


def get_constants(Szy_list, aircraft_class):
    '''
    Step for calculating the constants Lambda = - Sz / Iyy and Lambda = - Sy / Izz. Note Lambda is arbitrarily set \
    during the derivation phase, not equivalent to anything in literature.
    :param Szy_list: list of size 2 of the shear forces in z and y direction respectively.
    : param aircraft_class: Using the MoI values calculated from the aileronProperties section.
    :return: tuple containing both Lambdas (Lambda_z, Lambda_y)
    '''
    iyy, izz = aircraft_class.momInertia()
    Lambda_z = - Szy_list[0] / iyy
    Lambda_y = - Szy_list[1] / izz

    return np.array([Lambda_z, Lambda_y])


def get_idealised_shear_flow(aircraft_class):
    '''
    Gets the sum of base shear flow of the IDEALISED BOOM SECTIONS ONLY. They are referenced in Sherman's derivations
    as Epsilon_z and Epsilon_y, but are truncated in the variables used below for convenience.
    :boom_area_array:
    :x_y_array:
    :return:
    '''
    str_pos = aircraft_class.st_pos  # called st_pos in original class file
    # print(str_pos)
    if str_pos.size == aircraft_class.n_st * 2:  # str_pos is a 17 x 2 numpy array
        # initialise output sum
        eps_z = 0
        eps_y = 0
        # a_stif =  aircraft_class.w_st * aircraft_class.t_st + aircraft_class.h_st * aircraft_class.t_st
        a_stiff = aircraft_class.A_stif
        print(a_stiff)
        for iter in range(0, aircraft_class.n_st):
            eps_z += a_stiff * str_pos[iter][0]
            eps_y += a_stiff * str_pos[iter][1]

        return np.array([eps_z, eps_y])
    else:
        raise Exception("Something's wrong, I can feel it!")


def get_geometric_2ndary(aircraft_class):
    '''
    Gets the relevant derived aircraft geometry values for use in subsequent calculations.
    :param aircraft_class:
    :return:
    '''
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    z_tr = z_ii(aircraft_class)  # z length of the trailing edge section in cross-section(section II)
    radius = h_spar / 2  # half the length of the spar
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin

    return


def get_shear_flow_base(lambd_array, aircraft_class, eps_array):
    """
    Gets the first shear flow equation for the base shear flows in the semi-circular profile from Sherman's derived
    equation

    :return:
    """
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    print(z_bar)
    print(eps_array)
    z_tr = z_ii(aircraft_class)  # z length of the trailing edge section in cross-section(section II)
    radius = h_spar / 2  # half the length of the spar
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin
    q_11 = lambd_array[0] * (t_skin * (radius * radius + radius * z_bar * m.pi / 2) + eps_array[0]) + lambd_array[1] * (
            t_skin * radius * radius) + eps_array[1]
    q_12 = lambd_array[0] * (eps_array[0]) + lambd_array[1] * (
            t_spar * (radius - z_bar) * h_spar + eps_array[1]) + q_11
    # q_13 = q_11  # taking the assumption that symmetry = same shear flow
    q_13 = lambd_array[0] * (t_skin * (radius * radius + radius * z_bar * m.pi / 2) + eps_array[0]) + lambd_array[1] *(
        t_skin * radius * radius + eps_array[1]) + q_12

    q_21 = lambd_array[0] * (t_skin * (z_tr * l_sk + (radius - z_bar)) + eps_array[0]) + lambd_array[1] * (
            t_skin * - l_sk / 2 * radius + eps_array[1])
    q_22 = lambd_array[0] * eps_array[0] + lambd_array[1] * (t_spar * (radius + z_bar) * radius + eps_array[1]) + q_21
    q_23 = q_21  # due the symmetry
    q_23 = lambd_array[0] * (t_skin * (l_sk / 2 * z_tr + l_sk * (radius - z_bar)) + eps_array[0]) + lambd_array[1] * (
            t_skin * radius * l_sk / 2) + q_22

    output_array = np.array([[q_11, q_12, q_13], [q_21, q_22, q_23]])

    return output_array


# def get_shear_flow_redundant(shr_base_array, aircraft_class):
#     '''
#
#     :return:
#     '''
#
#     # taking the moment equilibrium with respect to centre of the spar
#     m11 = (shr_base_array[0, 0] + shr_base_array[0, 2]) * radius
def get_shortest_normal(aircraft_class):
    '''
    Solves for the shortest distance from the trailing edge skin to the centre of spar.
    :param aircraft_class:
    :return:
    '''
    z_tr = z_ii(aircraft_class)  # z length of the trailing edge section in cross-section(section II)
    theta = m.atan(aircraft_class.h / 2 / z_tr)
    return z_tr * m.sin(theta)


def get_shear_centre(aircraft_class):
    '''
    Solves for shear centre.
    :param aircraft_class:
    :return:
    '''
    lambdas = get_constants([0, 1], aircraft_class)
    epses = get_idealised_shear_flow(aircraft_class)
    qb_list = get_shear_flow_base(lambdas, aircraft_class, epses)
    print(qb_list)
    # since only the z coordinate of the shear centre is needed,

    # anotherinstance of quality of life variable naming at a cost of absolute dumbness...
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    z_tr = z_ii(aircraft_class)  # z length of the trailing edge section in cross-section(section II)
    radius = h_spar / 2  # half the length of the spar
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin

    def get_SC_twist():
        '''

        :return:
        '''
        # Calculate the terms in the matrix A and initialise matrix A
        a_11 = (2 * m.pi * radius) * (qb_list[0, 0] + qb_list[0, 2]) / t_skin
        a_12 = - h_spar / t_spar
        a_21 = a_12
        a_22 = l_sk / t_skin + h_spar / t_spar
        a_matrix = np.array([[a_11, a_12],
                             [a_21, a_22]])
        # Matrix B is a 2x1 matrix to be solved for
        b_1 = -1 * ((qb_list[0, 0] + qb_list[0, 2]) / t_skin * (2 * m.pi * radius) +
                    qb_list[0, 1] * h_spar / t_spar)
        b_2 = -1 * ((qb_list[1, 0] + qb_list[1, 2]) * l_sk / t_skin +
                    qb_list[1, 1] * h_spar / t_spar)
        b_matrix = np.array([[b_1], [b_2]])

        # Matrix c is a 2x1 matrix containing the base shear flows
        return np.linalg.solve(a_matrix, b_matrix)

    c_matrix = get_SC_twist()
    d = get_shortest_normal(aircraft_class)

    A_i = aircraft_class.A1
    A_ii = aircraft_class.A2
    # solve for shear centre
    eta = (qb_list[0, 0] + qb_list[0, 2]) * (2 * m.pi * radius * radius) + (
            qb_list[1, 0] + qb_list[1, 2]) * l_sk * d + 2 * A_i * c_matrix[0] + A_ii * c_matrix[1] - radius
    print('eta is: ',eta)


A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
# q = aeroLoad.get_value_at
q = A320.crossArea()
_ = A320.stringersPosition()
_ = A320.zCentroid()
_ = A320.momInertia()
get_shear_centre(A320)
