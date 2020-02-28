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
columns_history = ["Iteration","CL","CD","CSF","CMx","CMy","CMz","CFx","CFy","CFz","CL/CD","AoA","Custom_ObjFunc","Res_Flow[0]","Res_Flow[1]","Res_Flow[2]","Res_Flow[3]","Res_Flow[4]","Linear_Solver_Iterations","CFL_Number","Time(min)"]


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
    vel = (vel_x * vel_x + vel_y * vel_y)**(0.5)
    df["Total_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * vel * vel
    df["p/pt"] = df["Pressure"] / df["Total_Pressure"]


def create_dataframe(aoa, blades, location):
    filename = "aoa" + str(aoa) + "_blade" + str(blades) + "_" + location + ".dat"
    df = pd.read_csv(filename, sep="\s+", header=None, skiprows=13)
    df.columns = columns
    calculate_flowangle(df)
    calculate_totalpressure(df)
    return df


def create_dataframe_history(aoa, blades, zone):
    filename = "aoa" + str(aoa) + "_blade" + str(blades) + "_history" + str(zone) + ".dat"
    df = pd.read_csv(filename, sep=",", header=None, skiprows=3)
    df.columns = columns_history
    return df


def create_plot(df, x, y, color, label, xlabel, ylabel):
    plt.plot(df[x], df[y], label=label, color=color, linewidth=2.5)
    plt.legend(loc='best')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    axes = plt.gca()
    axes.set_xlim([-0.5, 1.5])


def create_multi_plot(df, x, y1, y2, color1, color2, label1, label2, xlabel, y1label, y2label):
    fig, ax1 = plt.subplots(figsize=(10,8))
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(y1label)
    ax1.plot(df[x], df[y1], color=color1, label=label1, linewidth = 2.5)
    plt.grid()
    plt.legend(loc="center left")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel(y2label)  # we already handled the x-label with ax1
    ax2.plot(df[x], df[y2], color=color2, label=label2, linewidth = 2.5)
    # fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.legend(loc="center right")
    axes = plt.gca()
    axes.set_xlim([-0.5, 1.5])


def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors

# # Plots per angle of attack for range of blades
aoa = [1, 2, 5, 10, 20]
aoarange = [1, 2, 5, 10, 20]
blades = [1, 5, 10, 20]

# ### Plot used in thesis ###
# # Flow angle for range of aoa and fixed number of blades
# plt.figure(figsize=(10,8))
# for i in range(0, len(aoa)):
#     colors = set_viridis(0.2, 0.8, len(aoa))
#     meridional = create_dataframe(aoa[i], 20, "meridional")
#     create_plot(meridional, "x", "Flow_Angle", colors[i], "AoA: " + str(aoa[i]), "x", "Flow Angle [deg]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.grid()
#
# ### Plot used in thesis ###
# # Flow angle for fixed aoa and range of blades
# plt.figure(figsize=(10,8))
# for i in range(0, len(blades)):
#     colors = set_viridis(0.2, 0.8, len(blades))
#     meridional = create_dataframe(5, blades[i], "meridional")
#     create_plot(meridional, "x", "Flow_Angle", colors[i], "Blades: " + str(blades[i]), "x", "Flow Angle [deg]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.grid()
#
# # AoA 5, blades 20: 1) pt and p/pt, 2) flow angle and mom_x
# meridional = create_dataframe(5, 20, "meridional")
# colors = set_viridis(0.2, 0.8, 2)
# #1) pt and p/pt,
# fig1, ax1 = plt.subplots(figsize=(12, 9.6))
# ln1 = ax1.plot(meridional["x"], meridional["Total_Pressure"], color=colors[0], label=r"$p_{t}$", linewidth=2.5)
# plt.grid()
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# ln2 = ax2.plot(meridional["x"], meridional["p/pt"], color=colors[1], label=r"$p/p_{t}$", linewidth=2.5)
# # make one legend for both lines
# lns = ln1+ln2
# labs = [l.get_label() for l in lns]
# ax1.legend(lns, labs, loc=0)
# ax1.set_xlabel("x")
# ax1.set_ylabel(r"$p$ [Pa]")
# ax1.set_ylim([102024,102044])
# ax2.set_ylabel(r"$p/p_{t}$ [-]")  # we already handled the x-label with ax1
# axes = plt.gca()
# axes.set_xlim([-0.5, 1.5])
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
#
# #2) flow angle and mom_x
# fig2, ax1 = plt.subplots(figsize=(12, 9.6))
# ln1 = ax1.plot(meridional["x"], meridional["Flow_Angle"], color=colors[0], label="Flow Angle", linewidth=2.5)
# plt.grid()
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# ln2 = ax2.plot(meridional["x"], meridional["Momentum_x"], color=colors[1], label=r"$\rho u$", linewidth=2.5)
# # make one legend for both lines
# lns = ln1+ln2
# labs = [l.get_label() for l in lns]
# ax1.legend(lns, labs, loc=0)
# ax1.set_xlabel("x")
# ax1.set_ylabel(r"Flow Angle [$^\circ$]")
# ax2.set_ylim([41.6606,41.6806])
# ax2.set_ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")  # we already handled the x-label with ax1
# axes = plt.gca()
# axes.set_xlim([-0.5, 1.5])
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
#
# # Plot of total pressure for range of AoA
# plt.figure(figsize=(10,8))
# for i in range(0, len(aoarange)):
#     colors = set_viridis(0.2, 0.8, len(aoarange))
#     meridional = create_dataframe(aoarange[i], 20, "meridional")
#     create_plot(meridional, "x", "Total_Pressure", colors[i], "AoA: " + str(aoarange[i]), "x", r"$p_{t}$ [Pa]")
#     plt.ylim([101000,103000])
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.grid()
#
# # Plot of p/pt for range of AoA
# plt.figure(figsize=(10,8))
# for i in range(0, len(aoarange)):
#     colors = set_viridis(0, 1, len(aoarange))
#     meridional = create_dataframe(aoarange[i], 20, "meridional")
#     create_plot(meridional, "x", "p/pt", colors[i], "AoA: " + str(aoarange[i]), "x", r"$p/p_{t}$ [-]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.ylim([0.992, 0.9935])
# plt.grid()
#
# # Plot of p/pt for range of blade numbers
# plt.figure(figsize=(10,8))
# for i in range(0, len(blades)):
#     colors = set_viridis(0, 1, len(blades))
#     meridional = create_dataframe(5, blades[i], "meridional")
#     create_plot(meridional, "x", "p/pt", colors[i], "Blades: " + str(blades[i]), "x", r"$p/p_{t}$ [-]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
# # plt.ylim([0.992, 0.9935])
# plt.grid()


# # Plot of velocity for range of AoA
# plt.figure()
# for i in range(0, len(aoarange)):
#     colors = set_viridis(0, 1.0, len(aoarange))
#     meridional = create_dataframe(aoarange[i], 20, "meridional")
#     create_plot(meridional, "x", "Velocity_x", colors[i], "AoA: " + str(aoarange[i]), "Vx over Body Force Domain", "x", "Vx [m/s]")


# ### Plot used in thesis ###
# # Plot of simulated vs analytical values of pressure rise
# aoarange = [1, 2, 5, 10, 20]
# analytical_dp = [-0.216363133617, -0.865980146513, -5.43555494521, -22.0789221222, -94.0745310168]
# simulation_dp = [0, 0, 0, 0, 0]
# for i in range(0, len(aoarange)):
#     meridional = create_dataframe(aoarange[i], 20, "meridional")
#     total_points = len(meridional)
#     delta_p = meridional.iloc[total_points-1]["Pressure"] - meridional.iloc[0]["Pressure"]
#     simulation_dp[i] = delta_p
#
# print(simulation_dp)
#
# x = np.arange(len(aoarange))  # the label locations
# width = 0.35  # the width of the bars
# colors = set_viridis(0.0, 1.0, 2)
#
# fig, ax = plt.subplots(figsize=(10,8))
# rects1 = ax.bar(x - width/2, analytical_dp, width, label='Analytical', color=colors[0])
# rects2 = ax.bar(x + width/2, simulation_dp, width, label='Simulation', color=colors[1])
# # Add some text for labels, title and custom x-axis tick labels, etc.
# ax.set_ylabel('Pressure Difference [Pa]')
# ax.set_xlabel('Angle of Attack')
# ax.set_xticks(x)
# ax.set_xticklabels(aoarange)
# ax.legend(loc="lower left")
# ax.grid()
# ax.set_axisbelow(True)
#
# fig.tight_layout()


#--------------------------- Compare NACA data with body force results ---------------------------#
aoa5_blade_flow_filename = "aoa5_blade20_naca.csv"
aoa5_blade_flow = pd.read_csv(aoa5_blade_flow_filename, skiprows=1, header=None)
aoa5_blade_flow = aoa5_blade_flow.drop(aoa5_blade_flow.columns[11], axis=1)
aoa5_blade_flow.columns = columns
xlocations = aoa5_blade_flow.x.unique()
aoa5_blade_averages = pd.DataFrame(columns=["x"])
aoa5_blade_averages.x = xlocations

for column in columns[2:]:
    column_array = np.zeros(len(xlocations), dtype='float64')
    for i in range(0, len(xlocations)):
        x = aoa5_blade_averages["x"].iloc[i]
        flow_at_x = aoa5_blade_flow[aoa5_blade_flow["x"] == x]
        col_average = flow_at_x[column].mean()
        column_array[i] = col_average
    column_df = pd.DataFrame(column_array)
    aoa5_blade_averages[column] = column_df

calculate_flowangle(aoa5_blade_averages)
calculate_totalpressure(aoa5_blade_averages)

aoa5_bf_flow = create_dataframe(5, 20, "meridional")

# Plots comparing blade and body force results
colors = set_viridis(0.2,0.8,2)
plt.figure(figsize=(10,8))
plt.plot(aoa5_blade_averages.x, aoa5_blade_averages.Flow_Angle,label="Physical Blade",color=colors[0],linewidth=2.5)
plt.plot(aoa5_bf_flow.x, aoa5_bf_flow.Flow_Angle, label="Body Force",color=colors[1],linewidth=2.5)
plt.legend(loc='best')
plt.xlabel("x [-]")
plt.ylabel(r"Flow Angle [$^\circ$]")
plt.grid()

colors = set_viridis(0.2,0.8,2)
plt.figure(figsize=(10,8))
plt.plot(aoa5_blade_averages.x, aoa5_blade_averages["p/pt"],label="Physical Blade",color=colors[0],linewidth=2.5)
plt.plot(aoa5_bf_flow.x, aoa5_bf_flow["p/pt"], label="Body Force",color=colors[1],linewidth=2.5)
plt.legend(loc='best')
plt.xlabel("x [-]")
plt.ylabel(r"$\frac{p}{p_t}$ [$-$]")
plt.grid()

colors = set_viridis(0.2,0.8,2)
plt.figure(figsize=(10,8))
plt.plot(aoa5_blade_averages.x, aoa5_blade_averages.Total_Pressure,label="Physical Blade",color=colors[0],linewidth=2.5)
plt.plot(aoa5_bf_flow.x, aoa5_bf_flow.Total_Pressure, label="Body Force",color=colors[1],linewidth=2.5)
plt.legend(loc='best')
plt.xlabel("x [-]")
plt.ylabel(r"Total Pressure [$Pa$]")
plt.grid()

colors = set_viridis(0.2,0.8,2)
plt.figure(figsize=(10,8))
plt.plot(aoa5_blade_averages.x, aoa5_blade_averages.Momentum_x,label="Physical Blade",color=colors[0],linewidth=2.5)
plt.plot(aoa5_bf_flow.x, aoa5_bf_flow.Momentum_x, label="Body Force",color=colors[1],linewidth=2.5)
plt.legend(loc='best')
plt.xlabel("x [-]")
plt.ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")
plt.grid()

# # Compare structured with unstructured grid results (aoa10, blade20)
# colors = set_viridis(0.2, 0.8, 2)
# structured = create_dataframe(10, 20, "meridional")
# unstructured = create_dataframe(10, 20, "meridional_unstructured")
#
# plt.figure(figsize=(10, 8))
# plt.plot(structured.x, structured.Flow_Angle, label="Structured", color=colors[0], linewidth=2.5)
# plt.plot(unstructured.x, unstructured.Flow_Angle, label="Unstructured", color=colors[1], linewidth=2.5)
# # plt.ylabel(r"Flow Angle [$^\circ$]")
# plt.xlabel("x [-]")
# plt.grid()
# plt.legend(loc="best")

# # Analyze convergence history of simulations
# plt.figure(figsize=(10, 8))
# colors = set_viridis(0.2, 0.8, 3)
# for i in range(0, 3):
#     data = create_dataframe_history(10, 20, i)
#     plt.plot(data["Iteration"], data["Res_Flow[0]"], color=colors[i], label="Zone "+str(i+1))
# plt.legend(loc="best")
# plt.grid()
# plt.xlim([10, 13500])
# plt.ylim([-11, -3.5])
# plt.xlabel("Iteration")
# plt.ylabel("Res[0]")

plt.show()
