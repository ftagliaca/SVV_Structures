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
    z_2 = aircraft_class.C_a - aircraft_class.h / 2
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
    # print('stringer positions :\n', str_pos)
    if str_pos.size == aircraft_class.n_st * 2:  # str_pos is a 17 x 2 numpy array
        # initialise output sum
        eps_z = 0
        eps_y = 0
        # a_stif =  aircraft_class.w_st * aircraft_class.t_st + aircraft_class.h_st * aircraft_class.t_st
        a_stiff = aircraft_class.A_stif
        # print(a_stiff)
        for iter in range(0, aircraft_class.n_st):
            eps_z += a_stiff * str_pos[iter][0]
            eps_y += a_stiff * str_pos[iter][1]
        return np.array([eps_z, eps_y])
    else:
        raise Exception("Something's wrong, I can feel it!")


def get_idealised_shear_contribution(aircraft_class):
    '''
    Get the summed carried shear flows on the booms per wall. This is used in the subsequent base
    shear flow calculation (get_shear_flow_base).
    :param aircraft_class:
    :return:
    '''
    # expressing the stringer position coordinates relative to centroid z-position
    z_bar = aircraft_class.z_centroid
    str_pos = aircraft_class.st_pos - np.hstack((np.zeros((17, 1)), z_bar * np.ones((17, 1))))
    # print('stringer positions wrt centroid: \n',str_pos)

    # creating a matrix of stringer areas for multiplication. Note that stringer on point 0,0 is halved since it's
    # 'shared'(a cut is made on the y=0 section of that point)
    a_stiff = aircraft_class.A_stif
    a_stiff_mat = np.vstack((1 / 2 * np.ones((1, 2)), np.ones((16, 2)))) * a_stiff

    # multiply the areas matrix with the stringer coordinates matrix before summation. array is A*[ycoord, zcoord]
    ay_az_mat = np.multiply(str_pos, a_stiff_mat)
    # print('stringer areas matrix:\n', ay_az_mat)  # turn on for debug

    # sum up the individual A*y and A*z, and sort them into numpy arrays for individual walls. Note: i = cell 1,
    # ii = cell 2. Remember that end points are n+1 for python indexing.
    eps_iy = np.array([[np.sum(ay_az_mat[0:2, 0]), 0, ay_az_mat[16, 0] + ay_az_mat[0, 0]]])
    eps_iz = np.array([[np.sum(ay_az_mat[0:2, 1]), 0, ay_az_mat[16, 1] + ay_az_mat[0, 1]]])
    eps_izy = np.transpose(np.vstack((eps_iz, eps_iy)))
    eps_iiy = np.array([[np.sum(ay_az_mat[slice(9, 16, 1), 0]), 0, np.sum(ay_az_mat[slice(2, 9, 1), 0])]])
    eps_iiz = np.array([[np.sum(ay_az_mat[slice(9, 16, 1), 1]), 0, np.sum(ay_az_mat[slice(2, 9, 1), 1])]])
    eps_iizy = np.transpose(np.vstack((eps_iiz, eps_iiy)))
    # print('eps_iizy: \n', eps_iiy,eps_iiz,'\n',eps_iizy)
    # print('eps_izy: \n', eps_iy, eps_iz, '\n', eps_izy)

    # note that both the epsilon arrays have the z in first column and the y in 2nd column
    # (inverted from the str_coords)
    return eps_izy, eps_iizy


def get_shear_flow_base(lambd_array, aircraft_class):
    """
    Gets the first shear flow equation for the base shear flows in the semi-circular profile from Sherman's derived
    equation

    :return:
    """
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    # z_bar = -0.21578
    eps_izy, eps_iizy = get_idealised_shear_contribution(aircraft_class)
    # print(z_bar)
    # print(eps_array)
    # print(h_spar)
    z_tr = - z_ii(aircraft_class)  # z length of the trailing edge section in cross-section(section II)
    radius = h_spar / 2  # half the length of the spar
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin
    q_11 = lambd_array[0] * (t_skin * (- radius * radius + radius * z_bar * m.pi / 2) + eps_izy[0, 0]) + lambd_array[
        1] * ((t_skin * radius * radius) + eps_izy[0, 1])
    q_12 = lambd_array[0] * (eps_izy[1, 0]) + lambd_array[1] * (
            t_spar * (- radius - z_bar) * h_spar + eps_izy[1, 1]) + q_11
    # q_13 = q_11  # taking the assumption that symmetry = same shear flow
    q_13 = lambd_array[0] * (t_skin * (- radius * radius + radius * z_bar * m.pi / 2) + eps_izy[2, 0]) + lambd_array[
        1] * (t_skin * - radius * radius + eps_izy[2, 1]) + q_12

    q_21 = lambd_array[0] * (t_skin * (z_tr * l_sk + (- radius - z_bar)) + eps_iizy[0, 0]) + lambd_array[1] * (
            t_skin * - l_sk / 2 * radius + eps_iizy[0, 1])
    q_22 = lambd_array[0] * eps_iizy[1, 0] + lambd_array[1] * (
            t_spar * (radius + z_bar) * radius + eps_iizy[1, 1]) + q_21
    q_23 = q_21  # due the symmetry
    q_23 = lambd_array[0] * (t_skin * (l_sk / 2 * z_tr + l_sk * (- radius - z_bar)) + eps_iizy[2, 0]) + lambd_array[
        1] * (t_skin * radius * l_sk / 2 + eps_iizy[2, 1]) + q_22

    output_array = np.array([[q_11, q_12, q_13], [q_21, q_22, q_23]])
    # print('output array: \n ', output_array)

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


def get_shear_center(aircraft_class):
    '''
    Solves for shear centre, eta. It's measured from the leading edge of the aileron cross-section.
    :param aircraft_class:
    :param szys: list of [Sz, Sy]. Remember to set sc_calc to False if you want flows
    :param sc_calc: Boolean whether to calculate shear center. Default is True.
    :return:
    '''
    # Since aileron is symmetric on the z axis, only a shear of S_y = 1 is applied.

    lambdas = get_constants([0, 1], aircraft_class)
    # epses = get_idealised_shear_flow(aircraft_class)
    qb_list = get_shear_flow_base(lambdas, aircraft_class)
    # print(qb_list)
    # since only the z coordinate of the shear centre is needed,

    # another instance of quality of life variable naming at a cost of absolute dumbness...
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    z_tr = z_ii(aircraft_class)  # z length of the trailing edge section in cross-section(section II)
    radius = h_spar / 2  # half the length of the spar
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin

    def get_sc_twist():
        '''

        :return:
        '''
        # Calculate the terms in the matrix A and initialise matrix A
        a_11 = (1 / 2 * m.pi * radius) / t_skin + h_spar / t_spar
        a_12 = - h_spar / t_spar
        a_21 = a_12
        a_22 = l_sk / t_skin + h_spar / t_spar
        a_matrix = np.array([[a_11, a_12],
                             [a_21, a_22]])
        # Matrix B is a 2x1 matrix to be solved for
        b_1 = -1 * ((qb_list[0, 0] + qb_list[0, 2]) / t_skin * (1 / 2 * m.pi * radius) +
                    qb_list[0, 1] * h_spar / t_spar)
        b_2 = -1 * ((qb_list[1, 0] + qb_list[1, 2]) * l_sk / t_skin +
                    qb_list[1, 1] * h_spar / t_spar)
        b_matrix = np.array([[b_1], [b_2]])

        # Matrix c is a 2x1 matrix containing the base shear flows
        return np.linalg.solve(a_matrix, b_matrix)

    c_matrix = get_sc_twist()
    print('c_matrix: ', c_matrix)
    d = get_shortest_normal(aircraft_class)

    A_i = aircraft_class.A1
    A_ii = aircraft_class.A2
    # solve for shear centre
    eta = - ((qb_list[0, 0] + qb_list[0, 2]) * (1 / 2 * m.pi * radius * radius) + (
            qb_list[1, 0] + qb_list[1, 2]) * l_sk * d + 2 * A_i * c_matrix[0] + A_ii * c_matrix[1])
    print('eta is: ', eta)

    # Output eta as a distance from the leading edge of the aileron
    return eta


def get_shear_distr(aircraft_class, szy_magnit, szy_applied):
    '''
    Gets the shear distribution from a given shear load (and point of application) and spits out the
    redundant shear flow
    :param aircraft_class: type of aircraft.
    :param szy_magnit: magnitude of shear load, [Sz, Sy] format.
    :param szy_applied: line of application of shear load, [z_coord,y_coord] from leading edge.
                        Uses the stated coordinate system. Units in metres.
    :return: redundant shear flows in cells I and II, [q_0,I, q_0,II]
    '''
    lambds = get_constants(szy_magnit, aircraft_class)
    q_base = get_shear_flow_base(lambds, aircraft_class)

    # another instance of quality of life variable naming at a cost of absolute dumbness...
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    z_tr = z_ii(aircraft_class)  # z length of the trailing edge section in cross-section(section II)
    radius = h_spar / 2  # half the length of the spar
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin

    d = get_shortest_normal(aircraft_class)

    # preparing terms from moment equation.
    mom_known = (q_base[0, 0] + q_base[0, 2]) * (1 / 2 * m.pi * radius * radius) + (
            q_base[1, 0] + q_base[1, 2]) * l_sk * d  # contribution of base shears
    mom_unknown_i = 2 * aircraft_class.A1  # unknown redundant shear cell I to be solved
    mom_unknown_ii = 2 * aircraft_class.A2  # unknown redundant shear cell II to be solved
    mom_rhs_y = szy_magnit[1] * (
            szy_applied[0] - z_bar)  # positive Sy and positive z distance from SC gives positive moment
    mom_rhs_z = -szy_magnit[0] * (szy_applied[1])  # negative Sz and positive y distance gives positive moment
    mom_rhs = mom_rhs_y + mom_rhs_z - mom_known

    # preparing terms from twist compatibility equations
    twist_11 = (1 / 2 * m.pi * radius) / t_skin + h_spar / t_spar
    twist_12 = - h_spar / t_spar
    twist_21 = twist_12
    twist_22 = l_sk / t_skin + h_spar / t_spar

    twist_rhs_1 = -1 * ((q_base[0, 0] + q_base[0, 2]) / t_skin * (1 / 2 * m.pi * radius) +
                        q_base[0, 1] * h_spar / t_spar)
    twist_rhs_2 = -1 * ((q_base[1, 0] + q_base[1, 2]) * l_sk / t_skin +
                        q_base[1, 1] * h_spar / t_spar)

    # combine all into large 3x3 matrix and solve using the power of linear algebra
    big_chungus_r1 = np.array([mom_unknown_i, mom_unknown_ii, 0])
    big_chungus_r2 = np.array([twist_11, twist_12, 1])
    big_chungus_r3 = np.array([twist_21, twist_22, 1])
    big_chungus_lhs = np.vstack((big_chungus_r1, big_chungus_r2, big_chungus_r3))
    big_chungus_rhs = np.array([[mom_rhs], [twist_rhs_1], [twist_rhs_2]])
    q0s = np.linalg.solve(big_chungus_lhs, big_chungus_rhs)
    print(q_base)
    print(np.hstack((q0s[0,0] * np.ones((2,1)), q0s[1,0] * np.ones((2,1)))))
    q_tot = q_base + np.vstack((q0s[0,0] * np.ones((1,3)), q0s[1,0] * np.ones((1,3))))
    print(q_tot)
    q0s = np.transpose(q0s)
    print(q0s)
    # output formats:
    # q0s gives a 1x3 matrix, with q0I, q0II and d/dz theta
    # q_base gives a 2x3 matrix with the q_base on each wall, row is cell, column is wall
    # q_tot gives a 2x3 matrix with q_tot = q_base + q0,cell, row is cell, column is wall
    return q0s, q_base, q_tot


A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
# q = aeroLoad.get_value_at
q = A320.crossArea()
_ = A320.stringersPosition()
_ = A320.zCentroid()
_ = A320.momInertia()

get_shear_center(A320)
get_shear_distr(A320, [0,1], [0,0])
# get_idealised_shear_contribution(A320)
