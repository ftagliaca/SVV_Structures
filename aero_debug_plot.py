import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from aero_loads import AerodynamicLoad
from aileronProperties import A320

aero_load = AerodynamicLoad(A320, 'aerodynamicloada320.dat')

## Constants
C_a = 0.547  # [m]
l_a = 2.771  # [m]

## Methods
z_i, x_i = A320.z_i, A320.x_i

data = np.genfromtxt('aerodynamicloada320.dat', delimiter=',')
print("data:", data.shape)
n_z, n_x = data.shape

dim_z = aero_load.grid_z_coordinates[-1] - aero_load.grid_z_coordinates[0]
dim_x = aero_load.grid_x_coordinates[-1] - aero_load.grid_x_coordinates[0]
dz = dim_z / n_z
dx = dim_x / n_x
regular_gridded_data = data * 0

Z_reg = np.arange(aero_load.grid_z_coordinates[0], aero_load.grid_z_coordinates[-1], dz)
X_reg = np.arange(aero_load.grid_x_coordinates[0], aero_load.grid_x_coordinates[-1], dx)

print("Z", Z_reg.shape)
print("X", X_reg.shape)

for i, z in enumerate(Z_reg):
    for j, x in enumerate(X_reg):
        regular_gridded_data[i, j] = aero_load.get_value_at(z, x)

Z_reg, X_reg = np.meshgrid(Z_reg, X_reg)

X, Z = np.meshgrid(x_i(np.arange(n_x) + 1), z_i(np.arange(n_z) + 1))

fig = plt.figure('Aerodynamic Loading')
ax = plt.axes(projection='3d')


print(Z_reg.shape, X_reg.shape, regular_gridded_data.shape)
ax.plot_wireframe(Z_reg, X_reg, regular_gridded_data.T, color='green')

print(Z.shape, X.shape, data.shape)
ax.plot_wireframe(Z, X, data, color='black')






plt.show()