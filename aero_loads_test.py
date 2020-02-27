from mpl_toolkits.mplot3d import Axes3D
from aero_loads import AerodynamicLoad
from aileronProperties import A320
from tempfile import TemporaryFile
import matplotlib.pyplot as plt
from unittest import TestCase
from typing import Tuple
import numpy as np
import unittest
import random


class AeroLoadsTest(TestCase):

    def __init__(self, methodName='runTest'):
        super().__init__(methodName=methodName)
        self.f1 = lambda x1, x2: np.sin(x1 * x2)
        self.f2 = lambda x1, x2: np.sin((x1 - 1) * 3.2) * np.sin((x2 - 7) * 0.3)
        self.f3 = lambda x1, x2: np.sin((x1 - 0.5) * 3)**2 + np.cos((x2 - 0.5) * 0.14)**2

        self.fs = [self.f1, self.f2, self.f3]


    def generic_test_aero_loads(self, f: int, label: str = 'Aerodynamic Loading', grid: Tuple = (np.arange(-0.5, 0.5, 0.011), np.arange(-7, 7, 0.11), np.arange(-0.5 + 0.02, 0.5, 0.01), np.arange(-7 + 0.1, 7, 0.04))):
        def f_(x: np.ndarray, y: np.ndarray):
            return f(x[:, np.newaxis], y[np.newaxis, :])

        # Source data for interpolation
        
        x, y, interp_x, interp_y = grid

        data_x, data_y = np.meshgrid(x, y, indexing='ij')
        data = f_(x, y)

        scaling_factor = data.max() - data.min()
        data /= scaling_factor

        def z_i(i, N_z = -1):
            return x[i - 1]

        def x_i(i, N_x = -1):
            return y[i - 1]

        A320.z_i, A320.x_i = z_i, x_i

        # Target and verification data
        data_verification = f_(interp_x, interp_y) / scaling_factor
        data_interpolation = data_verification * 0
        mesh_x, mesh_y = np.meshgrid(interp_x, interp_y, indexing='ij')

        # interpolation
        with TemporaryFile() as file:
            np.savetxt(file, data, delimiter=',')
            file.seek(0, 0)

            aero_loads = AerodynamicLoad(A320, file, correction_factor=1)
        
        data_interpolation = aero_loads.get_values_grid(interp_x, interp_y)

        error = np.average(np.abs(data_verification - data_interpolation))
        print(f"average error: {error}")

        fig = plt.figure(label)
        ax = plt.axes(projection='3d')
        

        ax.plot_wireframe(mesh_x, mesh_y, data_verification, color='black', label="Verification data")
        ax.plot_wireframe(mesh_x, mesh_y, data_interpolation, color='blue', label=f"Interpolated data (Avg Error: {error})")
#        ax.plot_wireframe(data_x, data_y, data, color='green', label="Source data")

        ax.set_xlabel("z")
        ax.set_ylabel("x")
        ax.set_zlabel("y")

        plt.legend()

        #fig.show()

        self.assertAlmostEqual(np.average(data_verification - data_interpolation), 0, 2)


    def run_tests(self, grid: Tuple, label = "Generic Grid"):
        for i, f in enumerate(self.fs):
            self.generic_test_aero_loads(f, label=f'Test {i} ({label})', grid=grid)


    def test_grid_1(self):
        print("Test grid 1:")
        x = np.arange(-0.5, 0.5, 0.011)
        y = np.arange(-7, 7, 0.11)
        interp_x = np.arange(-0.5 + 0.02, 0.5, 0.01)
        interp_y = np.arange(-7 + 0.1, 7, 0.04)

        grid = (x, y, interp_x, interp_y)

        self.run_tests(grid, label="Grid 1")


    def test_grid_2(self):
        print("Test grid 2:")
        x = np.arange(-0.5, 0.5, 0.011)
        y = np.arange(-7, 7, 0.11)
        interp_x = np.arange(-0.5 + 0.02, 0.5, 0.011 / 2.01)
        interp_y = np.arange(-7 + 0.1, 7, 0.11 / 2.01)

        grid = (x, y, interp_x, interp_y)

        self.run_tests(grid, label="Grid 2")
        

    def test_grid_3(self):
        print("Test grid 3:")
        x = np.arange(-0.5, 0.5, 0.011)
        y = np.arange(-7, 7, 0.11)
        interp_x = np.arange(-0.5 + 0.05, 0.5, 0.11)
        interp_y = np.arange(-7 + 0.01, 7, 0.011)

        grid = (x, y, interp_x, interp_y)

        self.run_tests(grid, label="Grid 3")
        
        plt.show()

if __name__ == "__main__":
    unittest.main()
    



                    
