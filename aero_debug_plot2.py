import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from aero_loads import AerodynamicLoad
from aileronProperties import A320
import time

np.set_printoptions(linewidth=None)

# Graph setup
fig = plt.figure('Aerodynamic Loading')
ax = plt.axes(projection='3d')
ax.set_xlabel("z")
ax.set_ylabel("x")
ax.set_zlabel("y")

# interpolation stuff
# Overwriting the grid functions because of simpler grid

aero_load = AerodynamicLoad(A320, 'aerodynamicloada320.dat')

z_sizes = -aero_load.grid_z_coordinates[:-1] + aero_load.grid_z_coordinates[1:]
x_sizes = -aero_load.grid_x_coordinates[:-1] + aero_load.grid_x_coordinates[1:]
print(f"Tile sizes in z-dir avg: {np.average(z_sizes)}, min: {z_sizes.min()}, max: {z_sizes.max()}")
print(f"Tile sizes in x-dir avg: {np.average(x_sizes)}, min: {x_sizes.min()}, max: {x_sizes.max()}")

n_z, n_x = aero_load.n_z, aero_load.n_x
print(n_z, n_x)


z = np.linspace(aero_load.grid_z_coordinates[0], aero_load.grid_z_coordinates[-1], 60)
x = np.linspace(aero_load.grid_x_coordinates[0], aero_load.grid_x_coordinates[-1], 60)



start_time = time.time()
data = aero_load.get_values_grid(z, x)
print(f"Time taken: {time.time() - start_time}")

Z, X = np.meshgrid(z, x, indexing='ij')
ax.plot_wireframe(Z, X, data, color='green')


# plotting
data = np.genfromtxt('aerodynamicloada320.dat', delimiter=',') * 1e3
print("data:", data.shape)
n_z, n_x = data.shape

print("avergae:", np.average(data))

Z, X = np.meshgrid(A320.z_i(np.arange(n_z) + 1), A320.x_i(np.arange(n_x) + 1), indexing='ij')
print("n_z", n_z, "n_x", n_x)

ax.plot_wireframe(Z, X, data, color='black')




plt.show()