from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Circle, PathPatch
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

inp_file = 'data/B737.inp'

###   Loading the data   ###
# Basis points
points = np.loadtxt(inp_file, skiprows=9, max_rows=6597-9, delimiter=',')
points[:, 1:] /= 1000.  # Scale from mm to m

# Rectangular elements
s4r_elements = np.loadtxt(inp_file, skiprows=6598, dtype=int, max_rows=13232-6598, delimiter=',')

# Assembly elements, where loads are applied
assembly_nodes = np.zeros((16, 4))

for i in range(16):
    assembly_nodes[i, :] = np.loadtxt(inp_file, skiprows=14146 + 2 * i, max_rows=1, delimiter=',')

assembly_nodes[:, 1:] /= 1000.  # Scale from mm to m

###   End loading data   ###

###   Setting up the plot   ###
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

ax.set_xlim(0, 1.25 * 2)
ax.set_ylim(0 - 0.25, 0.25)
ax.set_zlim(0 - 0.25, 0.25)
###   End setting up the plot ###



elements = []

for index, element_points in enumerate(s4r_elements[:, 1:]):
    if index % 6 != 0: continue

    elements.append(tuple(map(tuple, points[element_points - 1, 1:])))

ax.add_collection3d(Poly3DCollection(elements))


# assembly_nodes indices [4, 15)
# 4 and 14: load = -0.737
# [3, 13]: load = -1.474

ax.quiver(*np.split(assembly_nodes[4:15, 1:].T, 3), [0] * 11, [-0.737] + [-1.474]*9 + [-0.737], [0] * 11, length=0.05, normalize=False, pivot='tip', color='red')
ax.quiver(*np.split(assembly_nodes[3, 1:].T, 3), [0], [0], [-97.4], length=0.005, normalize=False, pivot='tip', color='red')

#ax.quiver(*np.split(assembly_nodes[4:15, 1:], 3, axis=1), [0] * 11, [-0.737] + [-1.474]*9 + [-0.737], [0] * 11, length=0.05, normalize=False, pivot='tip', color='red')


plt.show()