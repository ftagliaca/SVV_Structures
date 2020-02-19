import numpy as np
from typing import Tuple
from geometricalProperties import Aileron

# http://fourier.eng.hmc.edu/e176/lectures/ch7/node7.html
# https://en.wikipedia.org/wiki/Bicubic_interpolation

class AerodynamicLoad:
    
    def __init__(self, aileron: Aileron, filename: str):
        self.x_i, self.z_i = aileron.x_i, aileron.z_i # [] -> [m]

        self.data = np.genfromtxt(filename, delimiter=',') # [kN/m^2]
        self.n_z, self.n_x = self.data.shape

        # Coordinate of every data point [z, x]
        self.z_coordinates = np.zeros(self.n_z)
        self.x_coordinates = np.zeros(self.n_x)

        for i in range(self.n_z):
            self.z_coordinates = self.z_i(i, N_z=self.n_z)
        
        for j in range(self.n_x):
            self.x_coordinates = self.x_i(j, N_x=self.n_x)

        # np.ndarray that contains the 16 a_ij coefficients for every square
        self.grid_rectangles = np.zeros((self.n_z - 1, self.n_x - 1, 16))

        A_inv = self.__generate_interpolation_matrix__()

        for i in range(self.n_z):
            for j in range(self.n_x):
                f = self.__get_f_array__(i, j)
                self.grid_rectangles[i, j] = A_inv.dot(f)
    
    def __generate_interpolation_matrix__(self) -> np.ndarray:
        """Generates the array used to calculate the a_ij coefficients.

        Args:

        Returns:
            np.ndarray: A_inv. 
        """
        A = np.zeros((16, 16))
        
        # Helper for readability of array indices
        a = np.zeros((4, 4, 16))
        for i in range(4):
            for j in range(4):
                a[i, j] = np.zeros(16)
                a[i, j, i + j * 4] = True

        # The rows in the array are: (where the columns are a_00 til a_33)
        # f(0, 0) = a_00
        A[0] = a[0, 0]
        # f(1, 0) = a_00 + a_10 + a_20 + a_30
        A[1] = a[0, 0] + a[1, 0] + a[2, 0] + a[3, 0]
        # f(0, 1) = a_00 + a_01 + a_02 + a_03
        A[2] = a[0, 0] + a[0, 1] + a[0, 2] + a[0, 3]
        # f(1, 1) = sum a_ij
        A[3] = np.ones(16)
        

        # f_z(0, 0) = a_10
        A[4] = a[1, 0]
        # f_z(1, 0) = a_10 + 2 * a_20 + 3 * a_30
        A[5] = a[1, 0] + 2 * a[2, 0] + 3 * a[3, 0]
        # f_z(0, 1) = a_10 + a_11 + a_12 + a_13
        A[6] = a[1, 0] + a[1, 1] + a[1, 2] + a[1, 3]
        # f_z(1, 1) = sum i * a_ij
        A[7] = np.array([0, 1, 2, 3] * 4)
         
        # f_x(0, 0) = a_01
        A[8] = a[0, 1]
        # f_x(1, 0) = a_01 + a_11 + a_21 + a_31
        A[9] = a[0, 1] + a[1, 1] + a[2, 1] + a[3, 1]
        # f_x(0, 1) = a_01 + 2 * a_02 + 3 * a_03
        A[10] = a[0, 1] + 2 * a[0, 2] + 3 * a[0, 3]
        # f_x(1, 1) = sum j * a_ij
        A[11] = np.array([np.floor(i / 4) for i in range(4 * 4)])

        # f_zx(0, 0) = a_11
        A[12] = a[1, 1]
        # f_zx(1, 0) = a_11 + 2 * a_21 + 3 * a_31
        A[13] = a[1, 1] + 2 * a[2, 1] + 3 * a[3, 1]
        # f_zx(0, 1) = a_11 + 2 * a_12 + 3 * a_13
        A[14] = a[1, 1] + 2 * a[1, 2] + 3 * a[1, 3]
        # f_zx(1, 1) = sum i * j * a_ij
        for i in range(4):
            for j in range(4):
                A[15, i * 4 + j] = i * j

        return  np.linalg.inv(A)

    def __get_f_array__(self, z_i: int, x_i: int) -> np.ndarray:
        """Generates the f array that, by multiplying it with A_inv, results in the coefficients needed for interpolation.

        Args:
            z_i (int): z index of point (left upper corner of interpolation rectangle)
            x_i (int): x index of point (left upper corner of interpolation rectangle)

        Returns:
            np.ndarray: Array with 16 elements representing the f values
        """

        f = np.zeros(16)

        # 00 01 10 11

        for index, delta_z_i, delta_x_i in enumerate([(0, 0), (0, 1), (1, 0), (1, 1)]):
            local_z_i, local_x_i = z_i + delta_z_i, x_i + delta_x_i

            f[index] = self.data[local_z_i, local_x_i]

            # df / d_z

            h1, h2 = 1, 1

            if local_z_i == 0:
                h1 = 0
            elif local_z_i == self.n_z - 1:
                h2 = 0
            
            f[4 + index] = (self.data[local_z_i - h1, local_x_i] - self.data[local_z_i + h2, local_x_i]) / (self.z_coordinates[local_z_i - h1] - self.z_coordinates[local_z_i + h2])


            # df / d_x
            k1, k2 = 1, 1

            if local_x_i == 0:
                k1 = 0
            elif local_x_i == self.n_x - 1:
                k2 = 0
            
            f[8 + index] = (self.data[local_z_i, local_x_i - k1] - self.data[local_z_i, local_x_i + k2]) / (self.x_coordinates[local_x_i - k1] - self.x_coordinates[local_x_i + k2])

            # https://math.stackexchange.com/q/3296431
            # df / d_z d_x
            f[12 + index] = self.data[local_z_i + h2, local_x_i + k2] \
                          - self.data[local_z_i + h2, local_x_i - k1] \
                          - self.data[local_z_i - h1, local_x_i + k2] \
                          + self.data[local_z_i - h1, local_x_i - k1] \
                          / \
                          ( \
                            (self.x_coordinates[local_x_i - k1] - self.x_coordinates[local_x_i + k2]) \
                            * \
                            (self.z_coordinates[local_z_i - h1] - self.z_coordinates[local_z_i + h2])  \
                          )

        return f

    def __find_grid_square__(self, z: float, x: float) -> Tuple[int, int]:
        """Finds row and column of square that contains point with coordinates z, x.

        Args:
            z (float): z coordinate of point
            x (float): x coordinate of point

        Returns:
            Tuple[int, int]: Row and column of the square.
        """

        pass

    def get_value_at(self, z: float, x: float) -> float:
        """Finds the value of the aerodynamic load at the given position (z, x).

        Args:
            z (float): z coordinate of point of interest
            x (float): x coordinate of point of interest

        Returns:
            float: Aerodynamic load in N/m^2 at point (z, x).
        """

        z_i, x_i = self.__find_grid_square__(z, x)

        result = 0

        a_ijs = self.grid_rectangles[z_i, x_i]
        for i in range(4):
            for j in range(4):
                result += a_ijs[i + j * 4] * z**i * x**j

        return result * 1e3





        