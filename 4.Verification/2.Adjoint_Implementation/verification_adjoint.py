# Marc Latour - Thesis - 29/01/2020
# Visualize verification adjoint
import pandas as pd
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib import cm

# Latex template for plots
plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 18,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

# Columns for slices
columns = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]
columns_history_direct = ["Iteration","CL","CD","CSF","CMx","CMy","CMz","CFx","CFy","CFz","CL/CD","AoA","Custom_ObjFunc","Res_Flow[0]","Res_Flow[1]","Res_Flow[2]","Res_Flow[3]","Res_Flow[4]","Avg_MassFlow","Avg_Mach","Avg_Temp","Avg_Press","Avg_Density","Avg_Enthalpy","Avg_NormalVel","Uniformity","Secondary_Strength","Momentum_Distortion","Secondary_Over_Uniformity","Avg_TotalTemp","Avg_TotalPress","Pressure_Drop","Linear_Solver_Iterations","CFL_Number","Time(min)"]
columns_history_adjoint = ["Iteration","Sens_Geo","Sens_Mach","Sens_AoA","Sens_Press","Sens_Temp","Sens_AoS","Res_AdjFlow[0]","Res_AdjFlow[1]","Res_AdjFlow[2]","Res_AdjFlow[3]","Res_AdjFlow[4]","Linear_Solver_Iterations","CFL_Number","Time(min)"]

def create_dataframe_history(aoa, type, zone):
    filename = "aoa" + str(aoa) + "_" + str(type) + "_history" + str(zone) + ".dat"
    df = pd.read_csv(filename, sep=",", header=None, skiprows=3)
    if type == "direct":
        df.columns = columns_history_direct
    else:
        df.columns = columns_history_adjoint
    return df


def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors



# Analyze convergence history of simulations
plt.figure(figsize=(10, 8))
colors = set_viridis(0.2, 0.8, 2)
typerange = ["direct", "adjoint"]
for i in range(0, 2):
    data = create_dataframe_history(20, typerange[i], 1)
    if typerange[i] == "direct":
        plt.plot(data["Iteration"], data["Res_Flow[0]"], color=colors[i], label="Direct")
    else:
        plt.plot(data["Iteration"], data["Res_AdjFlow[0]"], color=colors[i], label="Adjoint")
plt.legend(loc="best")
plt.grid()
plt.xlim([0, 3000])
# plt.ylim([-11, -3.5])
plt.xlabel("Iteration")
plt.ylabel("Res[0]")

plt.show()
