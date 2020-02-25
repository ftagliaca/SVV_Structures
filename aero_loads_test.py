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

# Using 0.011 and 0.11
# Test 1:
# average diff: 9.320488189361016e-08
# average diff: 6.095636171884944e-07
# average diff: 8.289327858409918e-07
# average diff: 1.4716753667975696e-06
# average diff: -5.028727649139186e-07
# .Test 2:
# average diff: -1.4372062435501226e-06
# average diff: 6.712648547564665e-08
# average diff: -1.4372062435501226e-06
# average diff: 1.3260976223282184e-06
# average diff: 2.0125387562366698e-07

class AeroLoadsTest(TestCase):

    def generic_test_aero_loads(self, f: int):
        def f_(x: np.ndarray, y: np.ndarray):
            return f(x[:, np.newaxis], y[np.newaxis, :])

        # Source data for interpolation
        
        x = np.arange(-5, 5 + 0.5, 0.011)
        y = np.arange(-7, 7 + 0.5, 0.11)
        data = f_(x, y)

        def z_i(i, N_z = -1):
            return x[i - 1]

        def x_i(i, N_x = -1):
            return y[i - 1]

        A320.z_i, A320.x_i = z_i, x_i

        # Target and verification data
        interp_x = np.arange(-5 + 0.25, 5 + 0.25, 0.03)
        interp_y = np.arange(-7 + 0.25, 7 + 0.25, 0.03)
        data_verification = f_(interp_x, interp_y)
        data_interpolation = data_verification * 0
        mesh_x, mesh_y = np.meshgrid(interp_x, interp_y, indexing='ij')


        # interpolation
        with TemporaryFile() as file:
            np.savetxt(file, data, delimiter=',')
            file.seek(0, 0)

            aero_loads = AerodynamicLoad(A320, file)
            
            for x_i, x in enumerate(interp_x):
                for y_i, y in enumerate(interp_y):
                    data_interpolation[x_i, y_i] = aero_loads.get_value_at(x, y)
        
        print(f"average diff: {np.average(data_verification - data_interpolation)}")
        self.assertAlmostEqual(np.average(data_verification - data_interpolation), 0, 3)

#        plt.figure('Aerodynamic Loading')
#        ax = plt.axes(projection='3d')
#        ax.plot_wireframe(mesh_x, mesh_y, data_verification, color='black')
#        ax.plot_wireframe(mesh_x, mesh_y, data_interpolation, color='green')
#
#
#        ax.set_xlabel("z")
#        ax.set_ylabel("x")
#        ax.set_zlabel("y")
#
#        plt.show()


    def test_1(self):
        print("Test 1:")
        for _ in range(5):
            r = random.random()
            f = lambda x1, x2: np.sin(x1 * x2 * 0.25 * r)
            self.generic_test_aero_loads(f)

    def test_2(self):
        print("Test 2:")
        for _ in range(5):
            r = random.random()
            r1 = 1 if random.random() >= 0.5 else r
            r2 = 1 if random.random() < 0.5 else r

            f = lambda x1, x2: np.sin(x1 * r2) * np.sin(x2 * r1)
            self.generic_test_aero_loads(f)

if __name__ == "__main__":
    unittest.main()



                    
