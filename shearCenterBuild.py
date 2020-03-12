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

# import of class in order to use geometrical properties
#

# A320 = Aileron(0.547, 2.771, 0.153, 1.281, 2.681, 28.0, 22.5, 1.1, 2.9, 1.2, 1.5, 2.0, 17, 1.103, 1.642, 26, 91.7)
# q = aeroLoad.get_value_at
q = A320.crossArea()
_ = A320.stringersPosition()
_ = A320.zCentroid()
_ = A320.momInertia()


def shear_calc_env(aircraft_class, szy_list, mesh_size=10):
    # another instance of quality of life variable naming at a cost of absolute dumbness...
    h_spar = aircraft_class.h  # height of aileron in y direction, also the length of the spar
    l_sk = aircraft_class.l_s  # length of straight skin section
    z_bar = aircraft_class.z_centroid  # z_coord of the centroid
    z_tr = aircraft_class.C_a - aircraft_class.h / 2  # z length of the trailing edge section in cross-section(section II) (INIT)
    radius = h_spar / 2  # half the length of the spar
    t_spar = aircraft_class.t_sp  # thickness of spar
    t_skin = aircraft_class.t_sk  # thickness of skin
    n_st = aircraft_class.n_st  # number of stringers
    d_ip = (m.pi * radius) / 2 - aircraft_class.d_st  # intermediate distance between spar and stringer before spar(cell 1)
    d_st = aircraft_class.d_st  # distance between stringers on perimeter
    d_i = d_st - d_ip   # intermediate distance between spar and stringer after spar (in cell 2)

    # check for stringers that are fore/aft of spar
    # strs_cell_placement is whether the stringers are in cell 1 (idx 0) or cell 2 (idx 1)
    strs_zcoords = aircraft_class.st_pos[:, 1]
    strs_cell_placement = np.where(abs(strs_zcoords) < radius, 0, 1)
    # strs_cell_placement is whether the stringers are above (idx 1) or below (idx -1) z-axis
    strs_topdown_placement = np.where(aircraft_class.st_pos[:, 0] < 0, -1, 1)
    # print(strs_cell_placement)
    # print(strs_zcoords)
    # print(strs_topdown_placement)
    # print(aircraft_class.st_pos)

    # exit()

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

    def get_idealised_shear_contribution():
        """
        Get the summed carried shear flows on the booms per wall. This is used in the subsequent base
        shear flow calculation (get_shear_flow_base_v1).
        :return:
        """
        n_str = aircraft_class.n_st
        # looking at the stringer positions and the radius of the front part, it means there's 1 more
        # stringer in the front (5 total) and 6 at the trailing edge.
        # expressing the stringer position coordinates relative to centroid z-position
        str_pos_mat = aircraft_class.st_pos - np.hstack((np.zeros((n_str, 1)), z_bar * np.ones((n_str, 1))))
        # since leading edge stringer Area is halved, another instance of LE stringer is cloned to the end of array
        str_pos_mat = np.vstack((str_pos_mat, str_pos_mat[0]))
        # print('stringer positions wrt centroid: \n',str_pos_mat)
        print('converted stringer positions', str_pos_mat)

        # creating a matrix of stringer areas for multiplication. Note that stringer on LE is halved since it's
        # 'shared'(a cut is made on the y=0 section of that point)
        a_stiff = aircraft_class.A_stif
        # print('A_stiffener',a_stiff)
        a_stiff_mat = np.vstack((1 / 2 * np.ones((1, 2)), np.ones((n_str - 1, 2)), 1 / 2 * np.ones((1, 2)))) * a_stiff

        # multiply the areas matrix with the stringer coordinates matrix before summation. array is A*[ycoord, zcoord]
        ay_az_mat = np.multiply(str_pos_mat, a_stiff_mat)
        print('stringer areas matrix:\n', ay_az_mat)  # turn on for debug

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
        # print(ayaz_mat)
        eps_yz = np.vstack((np.sum(ayaz_mat[0:2], 0), np.sum(ayaz_mat[2:n_tr + 2], 0),
                            np.sum(ayaz_mat[n_tr + 2: n_st - 1], 0), np.sum(ayaz_mat[n_st - 1:], 0)))
        print(eps_yz)
        # print(str_pos)
        z_tr_i = - z_tr  # z length of the trailing edge section in cross-section(section II), but instantiated for
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

        output_array = np.array([[q_11, q_12, q_13], [q_21, q_22, q_23]])
        # print('output array: \n ', output_array)

        return output_array

    def get_cumulative_strs_influence(lambd_array, ayaz_mat, order=1):

        # lambd array is in zy format, ayaz_mat is in yz format
        # order of -1 means the procedure will be done in reverse (bottom-up), done by flipping the array
        if order == -1:
            ayaz_mat = np.flip(ayaz_mat, 0)

        # print(ayaz_mat)
        lambd_ayaz = np.multiply(lambd_array, ayaz_mat)
        # usage of epsilon denoting the summation method
        eps_yz = np.empty(ayaz_mat.shape)
        # print(ayaz_mat.shape)
        for item in range(ayaz_mat.shape[1]):
            # sum all items from prior ayaz
            eps_yz[item] = np.sum((lambd_ayaz[:item + 1]), 0)
        # print('eps: \n', eps_yz)
        eps_tot = np.sum(eps_yz, 1)
        return eps_tot

    def get_shear_flow_base_v2(lambd_array):

        str_pos, ayaz_mat = get_idealised_shear_contribution()

        # qb11 -----------
        range_11 = np.linspace(0, np.pi / 2, mesh_size)  # range_11 contains radians
        # print('range_11', range_11)
        qb11_z = lambd_array[0] * (
                (radius * t_skin * (range_11 * (radius - z_bar) + np.sin(range_11) * radius)) + ayaz_mat[-1, 1])
        qb11_y = lambd_array[1] * (radius * radius * t_skin * (np.cos(range_11) - 1) + ayaz_mat[-1, 0])
        # inputting the stiffener influence (note that there're two stiffeners on the arc, not incl. the LE one)
        # print(qb11_z)
        theta_stif = aircraft_class.d_st / radius
        # print(theta_stif)
        ayaz_stif11 = ayaz_mat[-3:-1]
        # print(qb11_z)
        # print('ayaz_11', ayaz_stif11)
        lambd_ayaz_stif11 = np.multiply(lambd_array, ayaz_stif11)
        # print(lambd_ayaz_stif11)
        # floor divide to get int (number of strs passed)
        strs_passed = range_11 * radius / aircraft_class.d_st
        strs_passed = strs_passed.astype(int)
        # print('strs_passed', strs_passed)
        eps11_tot = get_cumulative_strs_influence(lambd_array, ayaz_stif11, -1)
        # print('eps_tot', eps11_tot)
        qb11_pre_eps = np.add(qb11_z, qb11_y)
        qb11 = np.where(range_11 < theta_stif, qb11_pre_eps, qb11_pre_eps + eps11_tot[strs_passed - 1])
        # qb11_z = np.where(range_11 < theta_stif, qb11_z, qb11_z + eps11_yz[strs_passed-1, 1])
        # qb11_y = np.where(range_11 < theta_stif, qb11_y, qb11_y + eps11_yz[strs_passed-1, 0])
        # print('qb11_pre', qb11_pre_eps)
        # print('qb11', qb11)
        # exit()

        # qb12 -----------
        range_12 = np.linspace(0, h_spar, mesh_size)  # range_12 is about distances
        qb12_z = lambd_array[0] * (t_spar * (radius - z_bar) * range_12)
        qb12_y = lambd_array[1] * t_spar / 2 * (range_12 * range_12 - h_spar * range_12)
        # print(np.ones(qb11.shape[0]))
        qb12 = qb12_z + qb12_y + qb11[-1] * np.ones(mesh_size)
        # print(qb12)

        # qb13 -----------
        # range_13 is the same as range_11, so range_11 is used here, also in radians
        qb13_z = lambd_array[0] * (radius * t_skin * (range_11 * (radius - z_bar) + radius * (1 - np.cos(range_11))))
        qb13_y = lambd_array[1] * (radius * radius * t_skin * (np.sin(range_11)))
        qb13_proto = np.add(qb13_y, qb13_z)
        ayaz_stif3 = lambd_array[0] * ayaz_mat[1, 1] + lambd_array[1] * ayaz_mat[1, 0]
        # print('ayaz', ayaz_stif3)
        qb13_proto = np.where(range_11 < theta_stif, qb13_proto, qb13_proto + ayaz_stif3)
        qb13 = qb13_proto + qb12[-1] * np.ones(mesh_size)
        print(qb13)

        # qb21 -----------
        # range_21 is similar to range_12, the segment is split up to small segments of a defined mesh size.
        range_21 = np.linspace(0, l_sk, mesh_size)
        q21_z = lambd_array[0] * (t_spar * (radius - z_bar) * range_21)
        q21_y = lambd_array[1] * (t_spar * -range_12 * range_21)
        # prepare matrix for the stringers' contribution
        print(radius)
        print(l_sk, d_ip, d_i, d_st)
        print(range_21)
        q21_stif = np.where(range_21 < d_ip, 0, ayaz_mat[n_st - 1, 0])
        range_21_stif = (range_12 - d_ip) / d_st

        # print(q21_stif)
        print(range_21_stif)
        return

    def get_shortest_normal():
        '''
        Solves for the shortest distance from the trailing edge skin to the centre of spar.
        :return:
        '''
        theta = m.atan(aircraft_class.h / 2 / z_tr)
        return z_tr * m.sin(theta)

    def get_shear_center(aircraft_class):
        """
        Solves for shear centre, eta. It's measured from the leading edge of the aileron cross-section.
        :param aircraft_class:
        :return:
        """
        # Since aileron is symmetric on the z axis, only a shear of S_y = 1 is applied.

        lambdas = get_constants([0, -1])
        # epses = get_idealised_shear_flow(aircraft_class)
        qb_list = get_shear_flow_base_v1(lambdas)

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
            b_1 = -1 * ((qb_list[0, 0] + qb_list[0, 2]) / t_skin * (0.5 * m.pi * radius) +
                        qb_list[0, 1] * h_spar / t_spar)
            b_2 = -1 * ((qb_list[1, 0] + qb_list[1, 2]) * l_sk / t_skin +
                        qb_list[1, 1] * h_spar / t_spar)
            b_matrix = np.array([[b_1], [b_2]])

            # Matrix c is a 2x1 matrix containing the base shear flows
            return np.linalg.solve(a_matrix, b_matrix)

        c_matrix = get_sc_twist().flatten()
        print('c_matrix: ', c_matrix)
        d = get_shortest_normal()

        a_i = aircraft_class.A1
        a_ii = aircraft_class.A2
        # print(a_i, a_ii)
        # solve for shear centre
        eta = - ((qb_list[0, 0] + qb_list[0, 2]) * (1 / 2 * m.pi * radius * radius) + (
                qb_list[1, 0] + qb_list[1, 2]) * l_sk * d + 2 * a_i * c_matrix[0] + 2 * a_ii * c_matrix[1])
        # print('eta is: ', eta)

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
        lambds = get_constants(szy_magnit)
        q_base = get_shear_flow_base_v1(lambds)

        # get shortest distance perpendicular to the line
        d = get_shortest_normal()

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
        twist_a1 = 1 / (2 * aircraft_class.A1)
        twist_a2 = 1 / (2 * aircraft_class.A2)
        twist_11 = ((m.pi * radius) / t_skin + h_spar / t_spar) * twist_a1
        twist_12 = - h_spar / t_spar * twist_a1
        twist_21 = - h_spar / t_spar * twist_a2
        twist_22 = (2 * l_sk / t_skin + h_spar / t_spar) * twist_a2

        twist_rhs_1 = -1 * ((q_base[0, 0] + q_base[0, 2]) / t_skin * (1 / 2 * m.pi * radius) +
                            q_base[0, 1] * h_spar / t_spar)
        twist_rhs_2 = -1 * ((q_base[1, 0] + q_base[1, 2]) * l_sk / t_skin +
                            q_base[1, 1] * h_spar / t_spar)

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

    def torsional_stiffness(aircraft_class):
        # aileron dimensions and thicknesses imported through self.
        r = aircraft_class.h / 2
        z_tr = aircraft_class.C_a - r
        ds_halfcircle = m.pi * r
        ds_spar = aircraft_class.h
        ds_skin = m.sqrt((r * r) + (z_tr * z_tr))
        t_sp = aircraft_class.t_sp
        t_sk = aircraft_class.t_sk

    a1 = aircraft_class.A1
    a2 = aircraft_class.A2
    # print(r, a1,a2)
    #
    # print(t_sk, t_sp)

        # terms from the twist compatibility equation
        twist_a1 = 1 / (2 * a1)
        twist_a2 = 1 / (2 * a2)
        twist_11 = ((m.pi * r) / t_sk + ds_spar / t_sp) * twist_a1
        twist_12 = - ds_spar / t_sp * twist_a1
        twist_21 = - ds_spar / t_sp * twist_a2
        twist_22 = (2 * ds_skin / t_sk + ds_spar / t_sp) * twist_a2

        a_mat = np.array([[2 * a1, 2 * a2, 0],
                          [twist_11, twist_12, -1],
                          [twist_21, twist_22, -1]])
        b_mat = np.array([[1], [0], [0]])

        # solve for the G_d/dz_theta
        g_dtdz = np.linalg.solve(a_mat, b_mat)
        # print('g_dtdz', g_dtdz)
        # divide the following formula in part with A1 and A2
        # Both parts should consist of a spar part and the rest
        # Gdthetadz = 1 / (2 * a1) * (
        #         (shearflows[0, 0] + shearflows[0, 2]) * ds_quartercircle / t_sk +
        #         shearflows[0, 1] * ds_spar / t_sp) + 1 / (2 * a2) * (
        #                     (shearflows[1, 0] + shearflows[1, 2]) * ds_skin / t_sk + shearflows[1, 1] * ds_spar / t_sp)

        # T = 1  # given assumption
        J = 1 / g_dtdz[-1]

        return J

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
        j = torsional_stiffness(aircraft_class)
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

    def grapher(x, y, z):
        plt.title("Shear in y axis")
        plt.xticks(x, rotation=40)
        plt.ylabel('Shear force [N]')
        plt.xlabel('x-location [m]')
        plt.plot(x, y)
        plt.show()

        plt.title("Shear in z axis")
        plt.xticks(x, rotation=40)
        plt.ylabel('Shear force [N]')
        plt.xlabel('x-location[m]')
        plt.plot(x, z)
        plt.show()


# szy_list = list(szy_gen)
# print('szy_list: \n', szy_list)
# get_shear_distr(A320, [0, 1], [0, 0])

    # get_idealised_shear_contribution(A320)
    # mesh_points, szy_outputs, qtot = get_shear_analysis(A320)
    # mesh_size = 200
    x = np.linspace(0, A320.l_a, num=mesh_size)  # Subsequent functions accept numpy-arrays
    # x = [0, 0.5, 0.6, 2.7]
    # y = veri.aileron.Sy(x)
    # z = veri.aileron.Sz(x)
    y = []
    z = []
    # z = [-142927.3208521, 50482.49383354, 47337.01848392, 28483.05827346,
    #      81931.44007969, -47482.17520149, -40350.65987382, -50588.59840216,
    #      -52326.56806738, 109083.76089247]
    # y = [-40780.67606372, -24640.70289515, -23940.76485954, -27257.63765062,
    #      21286.85433852, 14219.16544322, 13295.9750672, 8606.43415073,
    #      9175.02022091, 90179.46565151]

    # mesh_points, szy_outputs, qtot = get_shear_analysis(A320, x, y, z)
    # j = torsional_stiffness(A320)
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


shear_calc_env(A320, [0, 0])
