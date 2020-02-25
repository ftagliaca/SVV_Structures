from typing import Tuple, Optional, List, Union, IO
from aileronProperties import Aileron
import numpy as np

# http://fourier.eng.hmc.edu/e176/lectures/ch7/node7.html
# https://en.wikipedia.org/wiki/Bicubic_interpolation

class AerodynamicLoad:

    def __init__(self, aileron: Aileron, filename: Union[str, IO], correction_factor: float = 1e3):
        """Initialises a new aerodynamic load object

        Args:
            aileron (Aileron): The aileron object
            filename (str or file-like object): The file containing the aerodynamic load data
            correct_factor (float): A factor by which the data is multiplied to account e.g. for different units (Default: 1000 N/kN)

        Returns:
            AerodynamicLoad: The load object.
        """

        self.z_i, self.x_i = aileron.z_i, aileron.x_i # [] -> [m]

        self.data = np.genfromtxt(filename, delimiter=',') # [kN/m^2]
        self.data *= correction_factor

        self.n_z, self.n_x = self.data.shape

        # Coordinate of every data point [z, x]
        self.grid_z_coordinates = self.z_i(np.arange(self.n_z) + 1, N_z=self.n_z)
        self.grid_x_coordinates = self.x_i(np.arange(self.n_x) + 1, N_x=self.n_x)

        # Invert order of data if coordinates are in decreasing order
        if self.grid_z_coordinates[0] > self.grid_z_coordinates[-1]:
            self.grid_z_coordinates = self.grid_z_coordinates[::-1]
            self.data = self.data[::-1, :]

        if self.grid_x_coordinates[0] > self.grid_x_coordinates[-1]:
            self.grid_x_coordinates = self.grid_x_coordinates[::-1]
            self.data = self.data[:, ::-1]

        # 2D array that contains all the tiles composing the interpolation.
        self.tiles: np.ndarray = np.empty((self.n_z, self.n_x), dtype='object')

        # Matrix necessary to calculate the coefficients of interpolation.
        self.A_inv = self.__generate_interpolation_matrix__()

        for i in range(self.n_z - 1):
            for j in range(self.n_x - 1):
                self.tiles[i, j] = self.init_tile(i, j)


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

        return np.linalg.inv(A)

    def init_tile(self, tile_z_i: int, tile_x_i: int) -> 'Tile':
        """Initialises a new tile object, which then can be used to interpolate that part of the graph.

        Args:
            tile_z_i (int): z index of tile
            tile_x_i (int): x index of tile

        Returns:
            Tile: The tile instance.
        """

        outer_self = self

        z_exponents: np.ndarray = np.tile(np.arange(0, 4), (4, 1)).flatten()
        x_exponents: np.ndarray = np.arange(0, 4).repeat(4)

        class Tile:

            # z_i and x_i are also the index of the corner with smallest z and x, respectively
            def __init__(self, z_i: int, x_i: int):
                self.z_i = z_i
                self.x_i = x_i

                # tile edge points
                self.z0 = outer_self.grid_z_coordinates[z_i]
                self.z1 = outer_self.grid_z_coordinates[z_i + 1]

                if self.z0 > self.z1:
                    tmp = self.z0
                    self.z0 = self.z1
                    self.z1 = tmp

                self.x0 = outer_self.grid_x_coordinates[x_i]
                self.x1 = outer_self.grid_x_coordinates[x_i + 1]

                if self.x0 > self.x1:
                    tmp = self.x0
                    self.x0 = self.x1
                    self.x1 = tmp

                # tile dimensions
                self.dz = self.z1 - self.z0
                self.dx = self.x1 - self.x0

                # f array
                self.f = np.zeros(16)
                self.__populate_f_array__()

                self.a = outer_self.A_inv.dot(self.f)


            def __get_value_and_derivatives__(self, delta_z_i: int, delta_x_i: int) -> Tuple[float, float, float, float]:
                """Internal function to get the values and approximate derivatives at the points relative to the smallest edge (0, 0), (1, 0), (0, 1), and (1, 1)
                The derivatives are approximated using finite differences, that, in term, use the two adjacent points. In case only
                one neighbouring point is available, the point itself is used for the finite difference calculation.

                Args:
                    delta_z_i (int): relative position in z direction
                    delta_x_i (int): relative position in x direction
                Returns:
                    f   : The value of the point.
                    fz  : A finite difference approximation of the partial derivative in z direction of that point.
                    fx  : A finite difference approximation of the partial derivative in x direction of that point.
                    fzx : A finite difference approximation of the partial derivative in x and z direction of that point.
                """

                def f(z_i, x_i):
                    return outer_self.data[z_i, x_i]

                z_i = self.z_i + delta_z_i
                x_i = self.x_i + delta_x_i

                minus_hz = 1 if not z_i == 0 else 0
                plus_hz = 1 if not z_i >= outer_self.n_z - 1 else 0
                delta_z = (minus_hz + plus_hz)

                minus_hx = 1 if not x_i == 0 else 0
                plus_hx = 1 if not x_i >= outer_self.n_x - 1 else 0
                delta_x = (minus_hx + plus_hx)

                f_ = f(z_i, x_i)
                fz = (f(z_i + plus_hz, x_i) - f(z_i - minus_hz, x_i)) / delta_z
                fx = (f(z_i, x_i + plus_hx) - f(z_i, x_i - minus_hx)) / delta_x

                fzx = (  f(z_i + plus_hz, x_i + plus_hx)
                       - f(z_i - minus_hz, x_i + plus_hx)
                       - f(z_i + plus_hz, x_i - minus_hx)
                       + f(z_i - minus_hz, x_i - minus_hx)
                      ) / (delta_x * delta_z)

                return f_, fz, fx, fzx


            def __populate_f_array__(self):
                """Internal function to populate this tile's 'f' array.

                The 'f' array is the array, that, in combination with the previously calculated inverse A array,
                can be used to calculate the coefficients necessary for interpolation.

                Args:

                Returns:

                """

                for index, (delta_z_i, delta_x_i) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)]):
                    f, fz, fx, fzx = self.__get_value_and_derivatives__(delta_z_i, delta_x_i)
                    self.f[index + 0 ] = f
                    self.f[index + 4 ] = fz * self.dz
                    self.f[index + 8 ] = fx * self.dx
                    self.f[index + 12] = fzx * self.dx * self.dz


            def get_relative_dir(self, z: float, x: float) -> int:
                """This function returns the direction (in term of indices) of the tile, in which the point (z, x) can be found.

                Args:
                    z (float): z-coordinate of the point
                    x (float): x-coordinate of the point

                Returns:
                    (z, x) Tuple[int, int]: The index direction in z and x.
                """

                dz = 0
                if z > self.z1 and self.z_i != outer_self.n_z - 2:
                    dz = 1
                elif z < self.z0 and self.z_i != 0:
                    dz = -1

                dx = 0
                if x > self.x1 and self.x_i != outer_self.n_x - 2:
                    dx = 1
                elif x < self.x0 and self.x_i != 0:
                    dx = -1

                return dz, dx


            def get_value_at(self, z: float, x: float) -> Optional[float]:
                """This function returns the interpolated value at the point (z, x).

                Args:
                    z (float): z-coordinate of the point
                    x (float): x-coordinate of the point

                Returns:
                    float: The value at the point (z, x).
                """


                z_bar = (z - self.z0) / (self.z1 - self.z0)
                x_bar = (x - self.x0) / (self.x1 - self.x0)

                return (self.a * z_bar ** z_exponents * x_bar ** x_exponents).sum()


            def __repr__(self):
                return "T{" + f"{self.z0, self.x0} -> {self.z1, self.x1}" + "}"

        return Tile(tile_z_i, tile_x_i)


    def __find_grid_square__(self, z: float, x: float) -> 'Tile':
        """Finds the tile that contains point with coordinates z, x.

        Args:
            z (float): z coordinate of point
            x (float): x coordinate of point

        Returns:
            Tile: Tile containing the point.
        """

        # exempt last value as its index does not correspond to a tile
        id_z = (np.abs(self.grid_z_coordinates[:-1] - z)).argmin()
        id_x = (np.abs(self.grid_x_coordinates[:-1] - x)).argmin()

        tile = self.tiles[id_z][id_x]

        dz, dx = tile.get_relative_dir(z, x)
        tile = self.tiles[id_z + dz][id_x + dx]

        return tile

    def __find_grid_squares__(self, z: Union[float, np.ndarray], x: Union[float, np.ndarray]) -> 'Tile':
        """Finds the tile that contains point with coordinates z, x.

        Args:
            z (float): z coordinate of point
            x (float): x coordinate of point

        Returns:
            Tile: Tile containing the point.
        """

        z: np.ndarray = np.array(z, ndmin=1)
        x: np.ndarray = np.array(x, ndmin=1)

        def find_closest(reference_coords: np.ndarray, coords: np.ndarray) -> np.ndarray:
            ref_coords_array = np.tile(reference_coords[:, np.newaxis], (1, len(coords)))

            diff = coords - ref_coords_array
            diff[diff < 0] = diff.max()
            indices = diff.argmin(axis=0)

            return indices

        # exempt last value as its index does not correspond to a tile
        id_z = find_closest(self.grid_z_coordinates[:-1], z)
        id_x = find_closest(self.grid_x_coordinates[:-1], x)

        tiles = self.tiles[id_z, id_x]

        return tiles

    def get_values_grid(self, z_coordinates: np.ndarray, x_coordinates: np.ndarray) -> np.ndarray:
        """Finds the values of the aerodynamic load in the given (rectilinear) grid with the given z and x grid coordinates.

        Args:
            z (float or np.ndarray): z coordinate(s) of point of interest
            x (float or np.ndarray): x coordinate(s) of point of interest

        Returns:
            float or np.ndarray: Aerodynamic load in N/m^2 at point(s) (z, x).
        """
        Z, X = np.meshgrid(z_coordinates, x_coordinates, indexing='ij')
        z_coordinates: np.ndarray = np.array(Z, ndmin=1)
        x_coordinates: np.ndarray = np.array(X, ndmin=1)
        return self.get_value_at(Z.flatten(), X.flatten()).reshape((z_coordinates.shape[0], x_coordinates.shape[0]))

    def get_value_at(self, z: Union[float, np.ndarray], x: Union[np.ndarray, float]) -> Union[float, np.ndarray]:
        """Finds the value of the aerodynamic load at the given position (z, x).

        Args:
            z (float or np.ndarray): z coordinate(s) of point of interest
            x (float or np.ndarray): x coordinate(s) of point of interest

        Returns:
            float or np.ndarray: Aerodynamic load in N/m^2 at point(s) (z, x).
        """

        z: np.ndarray = np.array(z, ndmin=1)
        x: np.ndarray = np.array(x, ndmin=1)

        tiles = self.__find_grid_squares__(z, x)
        print("-"*10)
        print(np.size(tiles),np.size(x),np.size(z))
        result = [tiles[i].get_value_at(z[i], x[i]) for i, _ in enumerate(z)]

        return np.array(result)
