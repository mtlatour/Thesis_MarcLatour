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
          'font.size': 18,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

# Columns for slices
columns = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]
columns_blade = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Grid_Velocity_x", "Grid_Velocity_y", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]
columns_history = ["Iteration", "CL", "CD", "CSF", "CMx", "CMy", "CMz", "CFx","CFy","CFz","CL/CD","AoA","Custom_ObjFunc","Res_Flow[0]","Res_Flow[1]","Res_Flow[2]","Res_Flow[3]","Res_Flow[4]","Linear_Solver_Iterations","CFL_Number","Time(min)"]
columns_history_blade = ["Iteration","TotalPressureLoss_1","KineticEnergyLoss_1","EntropyGen_1","EulerianWork_1","PressureRatio_1","FlowAngleIn_1","FlowAngleOut_1","AbsFlowAngleIn_1","AbsFlowAngleOut_1","MassFlowIn_1","MassFlowOut_1","MachIn_1","MachOut_1","TotalEfficiency_1","TotalStaticEfficiency_1","TotalPressureLoss_2","KineticEnergyLoss_2","EntropyGen_2","EulerianWork_2","PressureRatio_2","FlowAngleIn_2","FlowAngleOut_2","AbsFlowAngleIn_2","AbsFlowAngleOut_2","MassFlowIn_2","MassFlowOut_2","MachIn_2","MachOut_2","TotalEfficiency_2","TotalStaticEfficiency_2","TotalPressureLoss_3","KineticEnergyLoss_3","EntropyGen_3","EulerianWork_3","PressureRatio_3","FlowAngleIn_3","FlowAngleOut_3","AbsFlowAngleIn_3","AbsFlowAngleOut_3","MassFlowIn_3","MassFlowOut_3","MachIn_3","MachOut_3","TotalEfficiency_3","TotalStaticEfficiency_3","TotalPressureLoss_4","KineticEnergyLoss_4","EntropyGen_4","EulerianWork_4","PressureRatio_4","FlowAngleIn_4","FlowAngleOut_4","AbsFlowAngleIn_4","AbsFlowAngleOut_4","MassFlowIn_4","MassFlowOut_4","MachIn_4","MachOut_4","TotalEfficiency_4","TotalStaticEfficiency_4","Res_Flow[0]","Res_Flow[1]","Res_Flow[2]","Res_Flow[3]","Res_Flow[4]","Linear_Solver_Iterations","CFL_Number","Time(min)"]

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


def calculate_totalrelpressure(df, N):
    U = (-N / 60.0) * 2 * m.pi * 0.1625
    df["Tangential_Vel"] = U
    vel_x = df["Velocity_x"]
    vel_y = df["Velocity_y"] - U
    rel_vel = (vel_x * vel_x + vel_y * vel_y) ** 0.5
    df["Relative_Vel"] = rel_vel
    df["Total_Rel_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * rel_vel * rel_vel


def create_dataframe_bf(n, cells):
    filename = "bf_n" + str(n) + "_cells" + str(cells) + ".dat"
    df = pd.read_csv(filename, sep="\s+", header=None, skiprows=13)
    df.columns = columns
    calculate_flowangle(df)
    calculate_totalpressure(df)
    calculate_totalrelpressure(df, n)
    return df


def create_dataframe_history(n, cells, blade):
    if blade:
        filename = "history_blade_n" + str(n) + ".dat"
        df = pd.read_csv(filename, sep=",", header=None, skiprows=3)
        df.columns = columns_history_blade
    else:
        filename = "history_bf_n" + str(n) + "_cells" + str(cells) + ".dat"
        df = pd.read_csv(filename, sep=",", header=None, skiprows=3)
        df.columns = columns_history
    return df


def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors


def create_dataframe_blade(n):
    whittle_blade_flow_filename = "blade_n" + str(n) + ".csv"
    whittle_blade_flow = pd.read_csv(whittle_blade_flow_filename, skiprows=1, header=None)
    whittle_blade_flow = whittle_blade_flow.drop(whittle_blade_flow.columns[11], axis=1)
    whittle_blade_flow.columns = columns_blade
    xlocations = whittle_blade_flow.x.unique()
    whittle_blade_averages = pd.DataFrame(columns=["x"])
    whittle_blade_averages.x = xlocations

    for column in columns_blade[2:]:
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
    calculate_totalrelpressure(whittle_blade_averages, n)
    return whittle_blade_averages


# Comparison of blade flow and bf flow with decreasing mesh size
blade_cellnum = 39441
rot_speed = 2362
cells_range = [2714, 4802, 6930, 10807, 19122, 38922]
blade_flow = create_dataframe_blade(rot_speed)
bf_results = np.zeros((len(cells_range), 2))
colors = set_viridis(0.2, 1, len(cells_range))

# Determine total run time blade simulation
blade_history = create_dataframe_history(rot_speed, 12903, True)
total_iter = blade_history["Iteration"].iloc[-1]
num_thousands = total_iter/1000
blade_runtime = 0
for i in range(1,num_thousands+1):
    iter = 1000 * i
    runtime = blade_history["Time(min)"].iloc[iter]
    blade_runtime = blade_runtime + runtime

# Compare history files for run time
for case in range(0, len(cells_range)):
    bf_history = create_dataframe_history(rot_speed, cells_range[case], False)
    bf_results[case, 0] = cells_range[case]
    bf_results[case, 1] = bf_history["Time(min)"].iloc[-1] * 60

# Flow Angle comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow.Flow_Angle, label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(rot_speed, cells_range[case])
    plt.plot(bf_flow.x, bf_flow.Flow_Angle, label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.04, 0.08])
plt.xlabel("x [-]")
plt.ylabel(r"Flow Angle [$^\circ$]")

# Total pressure comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow.Total_Pressure, label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(rot_speed, cells_range[case])
    plt.plot(bf_flow.x, bf_flow.Total_Pressure, label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.04, 0.08])
plt.xlabel("x [-]")
plt.ylabel(r"$p_{t}$ [Pa]")

# X-momentum comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow.Momentum_x, label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(rot_speed, cells_range[case])
    plt.plot(bf_flow.x, bf_flow.Momentum_x, label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.04, 0.08])
plt.xlabel("x [-]")
plt.ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")

# pressure ratio comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow["p/pt"], label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(rot_speed, cells_range[case])
    plt.plot(bf_flow.x, bf_flow["p/pt"], label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.04, 0.08])
plt.xlabel("x [-]")
plt.ylabel(r"$p/p_{t}$ [-]")

# total relative pressure comparison between BF simulations
plt.figure(figsize=(12, 9.6))
plt.plot(blade_flow.x, blade_flow["Total_Rel_Pressure"], label="blade", color='k', linewidth=2.5)
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(rot_speed, cells_range[case])
    plt.plot(bf_flow.x, bf_flow["Total_Rel_Pressure"], label="bf,"+str(cells_range[case]), color=colors[case], linewidth=2.5)
plt.grid()
plt.legend(loc="best")
plt.xlim([-0.04, 0.08])
plt.xlabel("x [-]")
plt.ylabel(r"$p_{t_{rel}}$ [Pa]")

# Total enthalpy difference between BF simulations
dht_blade = blade_flow["Tangential_Vel"].iloc[0] * (blade_flow["Velocity_y"].iloc[-1] - blade_flow["Velocity_y"].iloc[0])
print("dht blade = " + str(dht_blade))
for case in range(0, len(cells_range)):
    bf_flow = create_dataframe_bf(rot_speed, cells_range[case])
    dht_bf = bf_flow["Tangential_Vel"].iloc[0] * (bf_flow["Velocity_y"].iloc[-1] - bf_flow["Velocity_y"].iloc[0])
    print("dht (cells=" + str(cells_range[case]) + ") =" + str(dht_bf))


# Plot of computational time blade vs. body force
plt.figure(figsize=(12, 9.6))
plt.loglog(bf_results[:, 0], bf_results[:, 1], linewidth=2.5, color="r", marker="*", label="Body Force")
# plt.loglog(blade_results[:, 0], blade_results[:, 1], linewidth=2.5, color="k", marker="*", label="Blade")
plt.loglog(blade_cellnum, blade_runtime * 60, marker="*", color="k", label="Blade")
plt.grid()
plt.legend(loc="best")
plt.xlabel("Number of Cells [-]")
plt.ylabel("Simulation Time [s]")


plt.show()
