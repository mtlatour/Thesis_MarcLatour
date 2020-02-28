# File to create camber normal values for openFOAM-SU2 comparison test case
import numpy as np
import cmath as m
import pandas as pd
import matplotlib.pyplot as plt

# Latex template for plots
plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 18,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

num_points = 100
start_bf = 0.0
end_bf = 1.0
coord_list = np.linspace(start_bf, end_bf, num_points)
nx_list = []
ny_list = []
ycoord_list = [0]
beta_list = []

beta_1 = 35.0
beta_2 = 25.0
delta_beta = beta_2 - beta_1
delta_coord = end_bf - start_bf
last_coord = coord_list[len(coord_list)-1]
first_coord = coord_list[0]

for i in range(0, len(coord_list)):
    beta_deg = beta_1 + (delta_beta/delta_coord) * (coord_list[i]-start_bf)
    beta = np.deg2rad(beta_deg)
    beta_list.append(beta)
    # nx = m.sin(m.pi + beta).real
    # ny = m.cos(beta).real
    # nx_list.append(nx)
    # ny_list.append(ny)

print(beta_list)

for i in range(1, len(coord_list)):
    dy = np.tan(beta_list[i-1]) * (coord_list[i] - coord_list[i-1])
    ycoord_list.append(ycoord_list[i-1]+dy)

plt.figure(figsize=(10,8))
plt.plot(coord_list, ycoord_list, 'k', linewidth=3)
plt.xlabel("x")
plt.ylabel("y")
# plt.axis('off')
plt.show()


# filename_nx = "Nx.dat"
# filename_ny = "Ny.dat"
#
# nx_file = open(filename_nx, "w")
# ny_file = open(filename_ny, "w")
#
# for item in range(0, len(coord_list)):
#     nx_file.write("{}, {}\n".format(coord_list[item], nx_list[item]))
#     ny_file.write("{}, {}\n".format(coord_list[item], ny_list[item]))
#
# # Create list of camber normals for use in c++
# filename = "Whittle_2D_cambernormals.txt"
# whittle = pd.read_csv(filename, sep="\s+", header=None)
# whittle.columns = ["x", "Nx", "Ny"]
# whittle_x = []
# whittle_nx = []
# whittle_ny = []
# for i in range(0,len(whittle)):
#     whittle_x.append(whittle["x"].iloc[i])
#     whittle_nx.append(whittle["Nx"].iloc[i])
#     whittle_ny.append(whittle["Ny"].iloc[i])
# # print(whittle_x)
# # print(whittle_nx)
# # print(whittle_ny)