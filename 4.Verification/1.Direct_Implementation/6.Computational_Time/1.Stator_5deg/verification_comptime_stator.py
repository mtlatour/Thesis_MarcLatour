# Marc Latour - Thesis - 06/09/2019
# Visualize verification part 1
import pandas as pd
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib import cm

# Latex template for plots
plt.rcParams['text.latex.preamble'] = [r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 20,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

# Columns for slices
columns = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]
columns_history = ["Iteration", "CL", "CD", "CSF", "CMx", "CMy", "CMz", "CFx","CFy","CFz","CL/CD","AoA","Custom_ObjFunc","Res_Flow[0]","Res_Flow[1]","Res_Flow[2]","Res_Flow[3]","Res_Flow[4]","Linear_Solver_Iterations","CFL_Number","Time(min)"]


# Functions for creating and plotting dataframes
def calculate_flowangle(df):
    vel_x = df["Momentum_x"] / df["Density"]
    vel_y = df["Momentum_y"] / df["Density"]
    df["Velocity_x"] = vel_x
    df["Velocity_y"] = vel_y
    df["Flow_Angle"] = np.arctan2(vel_y, vel_x) * 180 / m.pi


def calculate_totalpressure(df):
    vel_x = df["Momentum_x"] / df["Density"]
    vel_y = df["Momentum_y"] / df["Density"]
    vel = (vel_x * vel_x + vel_y * vel_y)**0.5
    df["Total_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * vel * vel
    df["p/pt"] = df["Pressure"] / df["Total_Pressure"]


def create_dataframe_bf(cells):
    filename = "bf_cells" + str(cells) + "_meridional.dat"
    df = pd.read_csv(filename, sep="\s+", header=None, skiprows=13)
    df.columns = columns
    calculate_flowangle(df)
    calculate_totalpressure(df)
    return df


def create_dataframe_history(cells, blade):
    if blade:
        filename = "history_blade_cells" + str(cells) + ".dat"
    else:
        filename = "history_cells" + str(cells) + ".dat"
    df = pd.read_csv(filename, sep="\s+", header=None, skiprows=3)
    df.columns = columns_history
    return df


def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors


def create_dataframe_blade():
    whittle_blade_flow_filename = "blade.csv"
    whittle_blade_flow = pd.read_csv(whittle_blade_flow_filename, skiprows=1, header=None)
    whittle_blade_flow = whittle_blade_flow.drop(whittle_blade_flow.columns[11], axis=1)
    whittle_blade_flow.columns = columns
    xlocations = whittle_blade_flow.x.unique()
    whittle_blade_averages = pd.DataFrame(columns=["x"])
    whittle_blade_averages.x = xlocations

    for column in columns[2:]:
        column_array = np.zeros(len(xlocations), dtype='float64')
        for i in range(0, len(xlocations)):
            x = whittle_blade_averages["x"].iloc[i]
            flow_at_x = whittle_blade_flow[whittle_blade_flow["x"] == x]
            col_average = flow_at_x[column].mean()
            column_array[i] = col_average
        column_df = pd.DataFrame(column_array)
        whittle_blade_averages[column] = column_df

    calculate_flowangle(whittle_blade_averages)
    calculate_totalpressure(whittle_blade_averages)
    return whittle_blade_averages


# Comparison of blade flow and bf flow with decreasing mesh size
cells_range = [189, 362, 553, 796, 959, 1335, 1980, 2462, 3220, 4356, 8824, 13377]
cells_range_blades = [3159, 5037, 9176, 12903]
blade_flow = create_dataframe_blade()
blade_history = create_dataframe_history(12903, True)
blade_cellnum = 12903

bf_results = np.zeros((len(cells_range), 2))
blade_results = np.zeros((len(cells_range_blades), 2))
colors = set_viridis(0.2, 1, len(cells_range))

# Compare history files for run time
for case in range(0, len(cells_range)):
    bf_history = create_dataframe_history(cells_range[case], False)
    bf_results[case, 0] = cells_range[case]
    bf_results[case, 1] = bf_history["Time(min)"].iloc[-1] * 60
for case in range(0, len(cells_range_blades)):
    blade_history = create_dataframe_history(cells_range_blades[case], True)
    blade_results[case, 0] = cells_range_blades[case]
    blade_results[case, 1] = blade_history["Time(min)"].iloc[-1] * 60

# Flow Angle comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow.Flow_Angle, label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(cells_range[case])
    plt.plot(bf_flow.x, bf_flow.Flow_Angle, label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.5, 1.5])
plt.xlabel("x [-]")
plt.ylabel(r"Flow Angle [$^\circ$]")

# Total pressure comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow.Total_Pressure, label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(cells_range[case])
    plt.plot(bf_flow.x, bf_flow.Total_Pressure, label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.5, 1.5])
plt.xlabel("x [-]")
plt.ylabel(r"$p_{t}$ [Pa]")

# X-momentum comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow.Momentum_x, label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(cells_range[case])
    plt.plot(bf_flow.x, bf_flow.Momentum_x, label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.5, 1.5])
plt.xlabel("x [-]")
plt.ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")

# pressure ratio comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow["p/pt"], label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(cells_range[case])
    plt.plot(bf_flow.x, bf_flow["p/pt"], label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.5, 1.5])
plt.xlabel("x [-]")
plt.ylabel(r"$p/p_{t}$ [-]")

plt.figure(figsize=(12, 9.6))
plt.loglog(bf_results[:, 0], bf_results[:, 1], linewidth=2.5, color="r", marker="*", label="Body Force")
plt.loglog(blade_results[:, 0], blade_results[:, 1], linewidth=2.5, color="k", marker="*", label="Blade")
# plt.plot(blade_cellnum, blade_history["Time(min)"].iloc[-1] * 60, marker="*", color="b", label="Blade")
plt.grid()
plt.legend(loc="best")
plt.xlabel("Number of Cells [-]")
plt.ylabel("Simulation Time [s]")


plt.show()
