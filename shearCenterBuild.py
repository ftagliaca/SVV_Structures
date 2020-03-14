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
import internalLoadsStresses as ils
from aileronProperties import A320
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from stress_plotting import colorline

# import of class in order to use geometrical properties
#

# A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
# q = aeroLoad.get_value_at
q = A320.crossArea()
_ = A320.stringersPosition()
_ = A320.zCentroid()
_ = A320.momInertia()

get_shear_centre = None
get_torsional_stiffness_ = None

def shear_calc_env(aircraft_class, szy_list, mesh_size=100):
    global get_shear_centre, get_torsional_stiffness_

    # mesh_size = 100
    # another instance of quality of life variable naming at a cost of absolute dumbness...
    print('are you ready for a bad time?')
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    l_tr = aircraft_class.C_a - aircraft_class.h / 2  # z length of the trailing edge section in cross-section(section II) (INIT)
    z_tr = - aircraft_class.C_a  # z coord of trailing edge section in cross section
    radius = h_spar / 2  # half the length of the spar
    z_spar = - radius
    # print('radius',radius)
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin
    n_st = aircraft_class.n_st  # number of stringers
    a_i = aircraft_class.A1
    a_ii = aircraft_class.A2
    # intermediate distance between spar and stringer before spar(cell 1)
    d_ip = (m.pi * radius) / 2 - 2 * aircraft_class.d_st
    d_st = aircraft_class.d_st  # distance between stringers on perimeter
    d_i = d_st - d_ip  # intermediate distance between spar and stringer after spar (in cell 2)
    yz_str = aircraft_class.n_st

    # check for stringers that are fore/aft of spar
    # strs_cell_placement is whether the stringers are in cell 1 (idx 0) or cell 2 (idx 1)

    def get_constants(szys):
        """
        Step for calculating the constants Lambda = - Sz / Iyy and Lambda = - Sy / Izz. Note Lambda is arbitrarily set \
        during the derivation phase, not equivalent to anything in literature.
        :param szys: list of size 2 of the shear forces in z and y direction respectively.
        : param aircraft_class: Using the MoI values calculated from the aileronProperties section.
        :return: tuple containing both Lambdas (Lambda_z, Lambda_y)
        """
        iyy, izz = aircraft_class.momInertia()
        const_z = - szys[0] / iyy
        const_y = - szys[1] / izz

        return np.array([const_z, const_y])

    def get_strs_wall_sort(ayaz_mat):
        # retrieve indices of cells based on fore/aft
        strs_zcoords = aircraft_class.st_pos[:, 1]
        strs_cell_placement = np.where(abs(strs_zcoords) < radius, 0, 1)
        # find number of strs on cell1 and cell2
        strs_c2_indices = np.nonzero(strs_cell_placement)[0]
        strs_c1_indices = np.nonzero(strs_cell_placement - 1)[0]

        # extract the bottom and top stringers
        c1_strs = ayaz_mat[strs_c1_indices]
        c1_strs = np.vstack((c1_strs, c1_strs[0]))
        c1_top = c1_strs[:int(c1_strs.shape[0] / 2), :]
        c1_bot = c1_strs[int(c1_strs.shape[0] / 2):, :]
        c2_strs = ayaz_mat[strs_c2_indices]
        c2_top = c2_strs[:int(c2_strs.shape[0] / 2), :]
        c2_bot = c2_strs[int(c2_strs.shape[0] / 2):, :]
        # all flipped so to be in line with direction of analysis
        c1_top = np.flip(c1_top, 0)
        c1_bot = np.flip(c1_bot, 0)
        c2_top = np.flip(c2_top, 0)
        c2_bot = np.flip(c2_bot, 0)
        #
        # print('c1', c2_top)
        # print(c2_bot)
        # print(strs_c1_indices, strs_c2_indices)
        # NOTE output is in ayaz format, even though it still works for the stringer coords
        return [c1_bot, c1_top, c2_bot, c2_top]

    def get_idealised_shear_contribution(mode=1):
        """
        Get the summed carried shear flows on the booms per wall. This is used in the subsequent base
        shear flow calculation (get_shear_flow_base_v1).
        :return:
        """
        # looking at the stringer positions and the radius of the front part, it means there's 1 more
        # stringer in the front (5 total) and 6 at the trailing edge.
        # expressing the stringer position coordinates relative to centroid z-position
        str_pos_mat = aircraft_class.st_pos - np.hstack((np.zeros((yz_str, 1)), z_bar * np.ones((yz_str, 1))))
        # since leading edge stringer Area is halved, another instance of LE stringer is cloned to the end of array
        str_pos_mat = np.vstack((str_pos_mat, str_pos_mat[0]))
        # print('stringer positions wrt centroid: \n',str_pos_mat)
        # print('converted stringer positions', str_pos_mat)

        # creating a matrix of stringer areas for multiplication. Note that stringer on LE is halved since it's
        # 'shared'(a cut is made on the y=0 section of that point)
        a_stiff = aircraft_class.A_stif
        # print('A_stiffener',a_stiff)
        a_stiff_mat = np.vstack((1 / 2 * np.ones((1, 2)), np.ones((yz_str - 1, 2)), 1 / 2 * np.ones((1, 2)))) * a_stiff
        # multiply the areas matrix with the stringer coordinates matrix before summation. array is A*[ycoord, zcoord]
        ay_az_mat = np.multiply(str_pos_mat, a_stiff_mat)
        # print('stringer areas matrix:\n', ay_az_mat)  # turn on for debug
        if mode != 1:
            ay_az_mat = get_strs_wall_sort(ay_az_mat)

        # note that both the epsilon arrays have the y in first column and the z in 2nd column
        # (inverted from the str_coords)
        return str_pos_mat, ay_az_mat

    def get_shear_flow_base_v1(lambd_array):
        """
        Gets the first shear flow equation for the base shear flows in the semi-circular profile from Sherman's derived
        equation

        :return:
        """

        # eps_izy, eps_yz = get_idealised_shear_contribution(aircraft_class)
        str_pos, ayaz_mat = get_idealised_shear_contribution()
        n_tr = int((n_st - 3) / 2)  # number of stringers on trailing edge, per skin, must be integer
        # print(n_tr)
        print(ayaz_mat)
        eps_yz = np.vstack((np.sum(ayaz_mat[0:2], 0), np.sum(ayaz_mat[2:n_tr + 2], 0),
                            np.sum(ayaz_mat[n_tr + 2: n_st - 1], 0), np.sum(ayaz_mat[n_st - 1:], 0)))
        # print(eps_yz)
        # print(str_pos)
        z_tr_i = - l_tr  # z length of the trailing edge section in cross-section(section II), but instantiated for
        # this function
        q_11 = lambd_array[0] * (t_skin * (- radius * radius + radius * z_bar * m.pi / 2) + eps_yz[0, 1]) + \
               lambd_array[1] * ((t_skin * radius * radius) + eps_yz[0, 0])
        q_12 = lambd_array[1] * (t_spar * (- radius - z_bar) * h_spar) + q_11
        # q_13 = q_11  # taking the assumption that symmetry = same shear flow
        q_13 = lambd_array[0] * (t_skin * (- radius * radius + radius * z_bar * m.pi / 2) + eps_yz[-1, 1]) + \
               lambd_array[1] * (t_skin * - radius * radius + eps_yz[-1, 0]) + q_12

        q_21 = lambd_array[0] * (t_skin * (z_tr_i * l_sk + (- radius - z_bar)) + eps_yz[2, 1]) + lambd_array[1] * (
                t_skin * - l_sk / 2 * radius + eps_yz[2, 0])
        q_22 = lambd_array[1] * (t_spar * (radius + z_bar) * radius) + q_21
        q_23 = q_21  # due the symmetry
        q_23 = lambd_array[0] * (t_skin * (l_sk / 2 * z_tr_i + l_sk * (- radius - z_bar)) + eps_yz[1, 1]) + \
               lambd_array[
                   1] * (t_skin * radius * l_sk / 2 + eps_yz[1, 0]) + q_22

        output_array = [q_11, q_12, q_13, q_21, q_22, q_23]
        # print('output array: \n ', output_array)

        return output_array

    def get_cumulative_strs_influence(lambd_array, ayaz_mat, order=1):

        # lambd array is in zy format, ayaz_mat is in yz format
        # order of -1 means the procedure will be done in reverse (bottom-up), done by flipping the array
        if order == -1:
            ayaz_mat = np.flip(ayaz_mat, 0)

        # print('waddis',ayaz_mat)
        lambd_ayaz = np.multiply(lambd_array, ayaz_mat)
        # print('who ses im gey', lambd_ayaz)
        # usage of epsilon denoting the summation method
        eps_yz = np.empty(ayaz_mat.shape)
        # print(ayaz_mat.shape)
        for item in range(ayaz_mat.shape[0]):
            # sum all items from prior ayaz
            eps_yz[item] = np.sum((lambd_ayaz[:item + 1]), 0)
        # print('eps:', eps_yz)
        eps_tot = np.sum(eps_yz, 1)
        # print('dedly komandos', eps_tot)
        return eps_tot

    def plot_colour_scatter(qb_toplot, qb_y, qb_z):
        fig = plt.figure()
        ax2 = plt.axes()
        p = ax2.scatter(qb_z, qb_y, c=(qb_toplot)*10e3, cmap='jet')  # wing profile
        clb = fig.colorbar(p)
        ax2.set_title('Shear flow distribution')
        ax2.set_ylabel('y [m]')
        ax2.set_xlabel('z [m]')
        ax2.set_xlim(z_tr - 0.2, 0.2)
        ax2.invert_xaxis()
        clb.set_label(' Shear flow [N/m]')
        plt.show()

    def plot_colour_(qb_toplot, qb_y, qb_z):
        fig, ax = plt.subplots()

        lc = colorline(qb_z, qb_y, z=qb_toplot, cmap='jet')
        clb = plt.colorbar(lc, orientation='horizontal')

        x_margin = 0.05
        ax.set_xlim(-x_margin -0.55, x_margin + 0)
        ax.invert_xaxis()
        ax.set_ylim(-0.125, 0.125)
        ax.set_title('Shear flow distribution')
        ax.set_ylabel('y [m]')
        ax.set_xlabel('z [m]')
        clb.set_label(' Shear flow [N/m]')
        plt.show()


    def get_mesh_coords(range11, range12, range21, range22):
        """
        Gets the yz coords of the mesh. Array sizes change with mesh sizes.
        :param range11: array of s radians for circular segments
        :param range12: array of s lengths for spar
        :param range21: array of s lengths for half-spar
        :param range22: array of s lengths for trailing edge skin panels
        :return: list of lists of arrays in yz format
        """
        # gets zy coordinates of the sliced crosssections
        # use anchor points which are already defined
        yz_11 = [z_spar * np.sin(range11), z_spar + radius * np.cos(range11)]
        yz_12 = [z_spar + range12, z_spar * np.ones(range12.size)]
        yz_13 = [radius * np.cos(range11), z_spar + radius * np.sin(range11)]
        yz_21 = [-range21, z_spar * np.ones(range21.size)]
        yz_22 = [z_spar + range22 / l_sk * radius, z_spar - range22 / l_sk * l_tr]
        yz_23 = [range22 / l_sk * radius, z_tr + range22 / l_sk * l_tr]
        yz_24 = [radius - range21, z_spar * np.ones(range21.size)]

        # for flattening everything into plotting
        big_y = np.vstack((yz_11[0], yz_12[0], yz_13[0], yz_21[0], yz_22[0], yz_23[0], yz_24[0])).flatten()
        big_z = np.vstack((yz_11[1], yz_12[1], yz_13[1], yz_21[1], yz_22[1], yz_23[1], yz_24[1])).flatten()
        # print('waddisshit', big_y,big_z)
        # fig, ax = plt.subplots()
        # ax.scatter(big_z, big_y)
        # plt.show()

        # return [yz_11,yz_12,yz_13, yz_21, yz_22, yz_23, yz_24]
        return big_y, big_z

    def get_shear_flow_base_v2(lambd_array):
        # lambd_array = np.reshape(lambd_array, (1,2))
        str_pos, ayaz_mat = get_idealised_shear_contribution(mode=-1)
        # print(ayaz_mat)
        ayaz11_stif = ayaz_mat[0]
        # qb11 -----------
        range_11 = np.linspace(0, np.pi / 2, mesh_size)  # range_11 contains radians
        # print('range_11', range_11)
        qb11_z = lambd_array[0] * (
                (radius * t_skin * (range_11 * (-radius - z_bar) + np.sin(range_11) * radius)) + ayaz11_stif[0, 1])
        qb11_y = lambd_array[1] * (radius * radius * t_skin * (np.cos(range_11) - 1) + ayaz11_stif[0, 0])
        # inputting the stiffener influence (note that there're two stiffeners on the arc, not incl. the LE one)
        theta_stif = aircraft_class.d_st / radius
        # floor divide to get int (number of strs passed)
        strs11_passed = range_11 * radius / d_st
        strs11_passed = strs11_passed.astype(int)
        eps11_tot = get_cumulative_strs_influence(lambd_array, ayaz11_stif[1:])  # not including LE stringer
        qb11_proto = np.add(qb11_z, qb11_y)
        qb11 = np.where(range_11 < theta_stif, qb11_proto, qb11_proto + eps11_tot[strs11_passed - 1])

        # qb12 -----------
        range_12 = np.linspace(0, h_spar, mesh_size)  # range_12 is about distance, not radians
        qb12_z = lambd_array[0] * (t_spar * (-radius - z_bar) * range_12)
        qb12_y = lambd_array[1] * t_spar / 2 * (range_12 * range_12 - h_spar * range_12)
        # print(np.ones(qb11.shape[0]))
        qb12 = qb12_z + qb12_y + qb11[-1] * np.ones(mesh_size)
        # print(qb12)

        # qb13 -----------
        # range_13 is the same as range_11, so range_11 is used here, also in radians
        ayaz13_stif = ayaz_mat[1]
        # print(ayaz13_stif)
        strs13_first_passed = d_ip / radius
        strs13_passed = (range_11 * radius - d_ip) / d_st
        strs13_passed = strs13_passed.astype(int)
        # print(range_11, strs13_first_passed, d_st / radius)
        # print('dadday',strs13_passed)
        eps13_tot = get_cumulative_strs_influence(lambd_array, ayaz13_stif)
        qb13_z = lambd_array[0] * (radius * t_skin * (range_11 * (-radius - z_bar) + radius * (1 - np.cos(range_11))))
        qb13_y = lambd_array[1] * (radius * radius * t_skin * (np.sin(range_11)))
        qb13_proto = np.add(qb13_y, qb13_z)
        # ayaz_stif3 = lambd_array[0] * ayaz_mat[1, 1] + lambd_array[1] * ayaz_mat[1, 0]
        # print('ayaz', ayaz_stif3)
        # print('y r u gay', eps13_tot, qb13_proto)
        qb13_proto = np.where(range_11 < strs13_first_passed, qb13_proto, qb13_proto + eps13_tot[strs13_passed])
        qb13 = qb13_proto + qb12[-1]

        # qb21 -----------
        # range_21 is similar to range_12, the segment is split up to small segments of a defined mesh size.
        range_21 = np.linspace(0, radius, mesh_size)
        qb21_z = lambd_array[0] * (t_spar * (- radius - z_bar) * range_21)
        qb21_y = lambd_array[1] * (t_spar * -range_12 * range_21)
        qb21 = np.add(qb21_y, qb21_z)

        # qb22 ------------
        range_22 = np.linspace(0, l_sk, mesh_size)
        ayaz22_stif = ayaz_mat[2]
        qb22_z = lambd_array[0] * t_skin * ((- radius - z_bar) * range_22 - (range_22 ** 2) / (2 * l_sk) * l_tr)
        qb22_y = lambd_array[1] * t_skin * (-radius * range_22 + (range_22 ** 2) / (2 * l_sk) * radius)
        qb22_proto = np.add(qb22_y, qb22_z)
        strs22_passed = (range_22 - d_i) / d_st
        strs22_passed = strs22_passed.astype(int)
        eps22 = get_cumulative_strs_influence(lambd_array, ayaz22_stif)
        # print('action is coming', range_22, d_i, d_st, strs22_passed)
        # print('ugandan pawer', qb22_proto, qb21[-1])
        qb22 = np.where(range_22 < d_i, qb22_proto, qb22_proto + eps22[strs22_passed]) + qb21[-1]
        # print(qb22[-1])

        # qb23 ------------
        ayaz23_stif = ayaz_mat[3]
        qb23_z = lambd_array[0] * t_skin * ((z_bar - z_tr) * range_22 + (range_22 ** 2) / (2 * l_sk) * l_tr)
        qb23_y = lambd_array[1] * t_skin * (range_22 ** 2) / (2 * l_sk) * radius
        qb23_proto = np.add(qb23_y, qb23_z)
        strs23_passed = (range_22 - d_st / 2) / d_st
        strs23_passed = strs23_passed.astype(int)
        eps23 = get_cumulative_strs_influence(lambd_array, ayaz23_stif)
        # print(range_22, strs23_passed)
        # print(qb23_proto)
        qb23 = np.where(range_22 < d_i / 2, qb23_proto, qb23_proto + eps23[strs23_passed]) + qb22[-1]
        # print(qb23[-1])

        # qb24 -------------
        qb24_z = lambd_array[0] * t_spar * (-radius - z_bar) * range_21
        qb24_y = lambd_array[1] * t_spar * (radius * range_21 - range_21 ** 2 / 2)
        qb24 = np.add(qb24_y, qb24_z) + qb23[-1]

        qb_discrete = [qb11, qb12, qb13, qb21, qb22, qb23, qb24]
        qb_net = [qb11[-1], qb12[-1], qb13[-1], qb21[-1], qb22[-1], qb23[-1], qb24[-1]]
        # qb_yz_coords = get_mesh_coords(range_11, range_12, range_21, range_22)
        qb_y, qb_z = get_mesh_coords(range_11, range_12, range_21, range_22)

        qb_plot = np.empty((0,))
        for item in qb_discrete:
            qb_plot = np.append(qb_plot, item)

        # plot_colour_(qb_plot, qb_y, qb_z)
        # print('ugandans unite', qb_discrete, qb_plot)

        return qb_discrete, qb_net, qb_plot

    def get_shortest_normal():
        '''
        Solves for the shortest distance from the trailing edge skin to the centre of spar.
        :return:
        '''
        theta = m.atan(aircraft_class.h / 2 / l_tr)
        return l_tr * m.sin(theta)

    def get_shear_center():
        """
        Solves for shear centre, eta. It's measured from the leading edge of the aileron cross-section.
        :param aircraft_class:
        :return:
        """
        # Since aileron is symmetric on the z axis, only a downwards shear of S_y = 1 is applied.

        lambdas = get_constants([0, -1])
        # epses = get_idealised_shear_flow(aircraft_class)
        # qb_list = get_shear_flow_base_v1(lambdas)
        _, qb_list, _ = get_shear_flow_base_v2(lambdas)

        # print(qb_list)
        # since only the z coordinate of the shear centre is needed,
        def get_sc_twist():
            '''

            :return:
            '''
            # Calculate the terms in the matrix A and initialise matrix A
            a_11 = (m.pi * radius) / t_skin + h_spar / t_spar
            a_12 = - h_spar / t_spar
            a_21 = a_12
            a_22 = 2 * l_sk / t_skin + h_spar / t_spar
            a_matrix = np.array([[a_11, a_12],
                                 [a_21, a_22]])
            # Matrix B is a 2x1 matrix to be solved for
            b_1 = -1 * ((qb_list[0] + qb_list[2]) / t_skin * (0.5 * m.pi * radius) +
                        qb_list[1] * h_spar / t_spar)
            b_2 = -1 * ((qb_list[4] + qb_list[5]) * l_sk / t_skin +
                        (qb_list[3] + qb_list[6]) * radius / t_spar)
            b_matrix = np.array([[b_1], [b_2]])

            # Matrix c is a 2x1 matrix containing the base shear flows
            return np.linalg.solve(a_matrix, b_matrix)

        c_matrix = get_sc_twist().flatten()
        print('c_matrix: ', c_matrix)
        d = get_shortest_normal()

        # print(a_i, a_ii)
        # solve for shear centre
        eta = - radius - ((qb_list[0] + qb_list[2]) * (1 / 2 * m.pi * radius * radius) + (
                qb_list[4] + qb_list[5]) * l_sk * d + 2 * a_i * c_matrix[0] + 2 * a_ii * c_matrix[1])
        print('eta is: ', eta, radius)

        # Output eta as a distance from the leading edge of the aileron
        return eta

    get_shear_centre = get_shear_center

    def get_shear_distr(szy_magnit, szy_applied):
        '''
        Gets the shear distribution from a given shear load (and point of application) and spits out the
        redundant shear flow
        :param aircraft_class: type of aircraft.
        :param szy_magnit: magnitude of shear load, [Sz, Sy] format.
        :param szy_applied: line of application of shear load, [z_coord,y_coord] from leading edge.
                            Uses the stated coordinate system. Units in metres.
        :return: redundant shear flows in cells I and II, [q_0,I, q_0,II]
        '''
        lambds = get_constants(szy_magnit)
        # q_base = get_shear_flow_base_v1(lambds)
        _, q_base, _ = get_shear_flow_base_v2(lambds)
        # get shortest distance perpendicular to the line
        d = get_shortest_normal()
        eta = get_shear_center()
        # preparing terms from moment equation.
        mom_known = (q_base[0] + q_base[2]) * (1 / 2 * m.pi * radius * radius) + (
                q_base[4] + q_base[5]) * l_sk * d  # contribution of base shears
        mom_unknown_i = 2 * aircraft_class.A1  # unknown redundant shear cell I to be solved
        mom_unknown_ii = 2 * aircraft_class.A2  # unknown redundant shear cell II to be solved
        mom_rhs_y = - szy_magnit[1] * (szy_applied[0] - eta)  # neg Sy and pos z distance from SC gives positive moment
        mom_rhs_z = szy_magnit[0] * szy_applied[1]  # pos Sz and pos y distance gives positive moment
        mom_rhs = mom_rhs_y + mom_rhs_z - mom_known

        # preparing terms from twist compatibility equations
        twist_a1 = 1 / (2 * aircraft_class.A1)
        twist_a2 = 1 / (2 * aircraft_class.A2)
        twist_11 = ((m.pi * radius) / t_skin + h_spar / t_spar) * twist_a1
        twist_12 = - h_spar / t_spar * twist_a1
        twist_21 = - h_spar / t_spar * twist_a2
        twist_22 = (2 * l_sk / t_skin + h_spar / t_spar) * twist_a2

        twist_rhs_1 = -1 * ((q_base[0] + q_base[2]) / t_skin * (1 / 2 * m.pi * radius) +
                            q_base[1] * h_spar / t_spar) * twist_a1
        twist_rhs_2 = -1 * ((q_base[4] + q_base[5]) * l_sk / t_skin +
                            (q_base[3] + q_base[6]) * h_spar / t_spar) * twist_a2

        # combine all into large 3x3 matrix and solve using the power of linear algebra
        big_chungus_r1 = np.array([mom_unknown_i, mom_unknown_ii, 0])
        big_chungus_r2 = np.array([twist_11, twist_12, -1])
        big_chungus_r3 = np.array([twist_21, twist_22, -1])
        big_chungus_lhs = np.vstack((big_chungus_r1, big_chungus_r2, big_chungus_r3))
        big_chungus_rhs = np.array([[mom_rhs], [twist_rhs_1], [twist_rhs_2]])
        q0s = np.linalg.solve(big_chungus_lhs, big_chungus_rhs)
        # print(q_base)
        # print(np.hstack((q0s[0, 0] * np.ones((2, 1)), q0s[1, 0] * np.ones((2, 1)))))
        q_tot = q_base + np.vstack((q0s[0, 0] * np.ones((1, 3)), q0s[1, 0] * np.ones((1, 3))))
        q_tot = q_tot.flatten()
        # print(q_tot)
        q0s = np.transpose(q0s)
        # print(q0s)
        # output formats:
        # q0s gives a 1x3 matrix, with q0I, q0II and d/dz theta
        # q_base gives a 2x3 matrix with the q_base on each wall, row is cell, column is wall
        # q_tot gives a 2x3 matrix with q_tot = q_base + q0,cell, row is cell, column is wall
        return q0s, q_base, q_tot

    def init_get_szy(aircraft_class, x_mesh_size):
        x_max = aircraft_class.l_a  # aileron length
        # creates a list of mesh points in x direction.
        slice_length = x_max / x_mesh_size
        x_mesh_points = [round(slice_length * n, 4) for n in range(0, x_mesh_size + 1, 1)]
        x_mesh_points = np.asarray(x_mesh_points)  # convert to numpy array for speed
        yield x_mesh_points
        for idx in range(x_mesh_points.size):
            x_point = x_mesh_points[idx]
            s_z = float(ils.S_z(x_point))  # / 1e3
            s_y = float(ils.S_y(x_point))  # / 1e3
            print('shear in z and y: ', np.array([s_z, s_y]))
            yield np.array([s_z, s_y])
        # return x_mesh_points

    def get_torsional_stiffness(aileron):
        h_a = aileron.h
        C_a = aileron.C_a
        
        t_spar = aileron.t_sp
        t_skin = aileron.t_sk

        A_cell1 = aileron.A1
        A_cell2 = aileron.A2

        l_side = np.sqrt((h_a / 2) ** 2 + (C_a - h_a/2) ** 2)

        A = np.array([[-A_cell2,   -1/t_spar,        1/t_spar + 2 * l_side / t_skin * 1 / h_a],
                    [0,          2 * A_cell1,                               2 * A_cell2],
                    [-A_cell1,   1/2 * np.pi / t_skin + 1 / t_spar,           -1/t_spar]])

        b = np.array([[0],
                    [1],
                    [0]])

        J = 2 / (np.linalg.solve(A, b)[0] * h_a)
        
        return J

    get_torsional_stiffness_ = get_torsional_stiffness

    # def get_torsional_j(aircraft_class):

    # def gen_szy(x_mesh_points):
    #
    #
    #     # return szy_list

    def get_shear_analysis(aircraft_class, extx, exty, extz):
        z_SC = float(get_shear_center(aircraft_class))
        print('z_SC', z_SC)
        mesh_size = len(extx)
        szy_gen = init_get_szy(A320, mesh_size)
        mesh_points = szy_gen.__next__()
        j = get_torsional_stiffness(aircraft_class)
        print('j', j)
        q_total_big = np.empty((1, 6))
        if len(exty) > 0:
            print("using verification data)")
            mesh_points = mesh_size
            szy_list = [[extz[n], exty[n]] for n in range(0, mesh_size, 1)]
            for x_idx in range(mesh_points):
                # szy_current = [extz[x_idx],exty[x_idx]]
                q0s, qbases, qtotals = get_shear_distr(A320, szy_list[x_idx], [z_SC, 0])
                # print('qtotals: ', extx[x_idx], '\n', qtotals)


        else:
            print("using our data")
            szy_list = np.zeros((mesh_size + 1, 2))
            qtotals = []
            for idx in range(0, mesh_size):
                try:
                    # print('output', szy_gen.__next__())
                    szy_current = szy_gen.__next__()
                    szy_list[idx, 0] = szy_current[0]
                    szy_list[idx, 1] = szy_current[1]
                    # q0s, qbases, qtotals = get_shear_distr(A320, next(szy_gen), [z_SC, 0])
                    q0s, qbases, qtotals = get_shear_distr(A320, szy_current, [z_SC, 0])
                    # print('q0s: \n', q0s)
                    # print('qbases: \n', qbases)
                    print('qtotals: ', extx[idx], '\n', qtotals)
                    q_total_big = np.vstack((q_total_big, qtotals))
                    # idx += 1

                except StopIteration:
                    break

        return mesh_points, szy_list, q_total_big

    # def grapher(x, y, z):
    #     plt.title("Shear in y axis")
    #     plt.xticks(x, rotation=40)
    #     plt.ylabel('Shear force [N]')
    #     plt.xlabel('x-location [m]')
    #     plt.plot(x, y)
    #     plt.show()
    #
    #     plt.title("Shear in z axis")
    #     plt.xticks(x, rotation=40)
    #     plt.ylabel('Shear force [N]')
    #     plt.xlabel('x-location[m]')
    #     plt.plot(x, z)
    #     plt.show()

    # szy_list = list(szy_gen)
    # print('szy_list: \n', szy_list)
    # get_shear_distr(A320, [0, 1], [0, 0])

    # get_idealised_shear_contribution(A320)
    # mesh_points, szy_outputs, qtot = get_shear_analysis(A320)
    # mesh_size = 200
    # x = np.linspace(0, A320.l_a, num=mesh_size)  # Subsequent functions accept numpy-arrays
    # x = [0, 0.5, 0.6, 2.7]
    # y = veri.aileron.Sy(x)
    # z = veri.aileron.Sz(x)
    # y = []
    # z = []
    # z = [-142927.3208521, 50482.49383354, 47337.01848392, 28483.05827346,
    #      81931.44007969, -47482.17520149, -40350.65987382, -50588.59840216,
    #      -52326.56806738, 109083.76089247]
    # y = [-40780.67606372, -24640.70289515, -23940.76485954, -27257.63765062,
    #      21286.85433852, 14219.16544322, 13295.9750672, 8606.43415073,
    #      9175.02022091, 90179.46565151]

    # mesh_points, szy_outputs, qtot = get_shear_analysis(A320, x, y, z)
    # j = get_torsional_stiffness(A320)
    # # print(j)
    # print(qtot)
    # # grapher(mesh_points, szy_outputs[:, 0], szy_outputs[:, 1])
    # # grapher(x, y, z)
    # # convert to shear stress since I don't have time
    # t_sk = A320.t_sk
    # t_sp = A320.t_sp
    # # thicknesses = np.array([t_sk, t_sp, t_sk, t_sk, t_sp, t_sk])
    # # big_boi = np.empty((0,6))
    # # for iter in range(mesh_size):
    # #     big_boi = np.vstack((big_boi,thicknesses))
    #
    # # q_tot = np.multiply(qtot )
    #
    # q_11t = qtot[:, 0] * t_sk
    # q_12t = qtot[:, 1] * t_sp
    # q_13t = qtot[:, 2] * t_sk
    # q_21t = qtot[:, 3] * t_sk
    # q_22t = qtot[:, 4] * t_sp
    # q_23t = qtot[:, 5] * t_sk
    # # print(q_11t.size)
    #
    # fig = plt.figure("Shear stress on cross section")
    # plotted = plt.axes(projection='3d')
    # plotted.set_xlabel("Wall")
    # plotted.set_ylabel("x_position [m]")
    # plotted.set_zlabel("Shear stress [N]")
    # plotted.plot_wireframe(1 * np.ones(q_11t.size),range(q_11t.size),q_11t, color="black")
    # plt.show()

    # print(A320.st_pos.size)
    # fig = plt.figure("Coordinates of stringers on cross-section")

    # z_plot = q_tot.flatten()
    # print(get_idealised_shear_contribution(A320))
    # print("z-coord of centroid:", z_bar)
    # print("Shear center z_coordinate", float(get_shear_center(aircraft_class)))
    get_shear_flow_base_v2([1, 1])
    get_shear_center()


shear_calc_env(A320, [0, 0])
