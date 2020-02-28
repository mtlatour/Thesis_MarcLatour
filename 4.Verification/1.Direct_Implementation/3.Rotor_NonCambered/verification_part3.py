# Marc Latour - Thesis - 16/09/2019
# Visualize verification part 3 - Rotating flat plate
import pandas as pd
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib import cm
import os.path

#Latex template for plots
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 18,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

# Columns for slices
columns = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]


# Functions for creating and plotting dataframes
def calculate_flowangle(df):
    vel_x = df["Momentum_x"] / df["Density"]
    vel_y = df["Momentum_y"] / df["Density"]
    df["Velocity_x"] = vel_x
    df["Velocity_y"] = vel_y
    df["Flow_Angle"] = np.arctan2(vel_y, vel_x) * 180 / m.pi


def calculate_totalpressure(df):
    vel_x = df["Velocity_x"]
    vel_y = df["Velocity_y"]
    vel = (vel_x * vel_x + vel_y * vel_y)**0.5
    df["Total_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * vel * vel
    df["p/pt"] = df["Pressure"] / df["Total_Pressure"]


def calculate_totalrelpressure(df, N, R):
    vel_tangential = (-N/60.0) * 2 * m.pi * R
    df["Tangential_Vel"] = vel_tangential
    vel_x = df["Velocity_x"]
    vel_y = df["Velocity_y"] - vel_tangential
    rel_vel = (vel_x * vel_x + vel_y * vel_y) ** 0.5
    df["Relative_Vel"] = rel_vel
    df["Total_Rel_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * rel_vel * rel_vel


def create_dataframe(aoa, n, location):
    filename = "aoa" + str(aoa) + "_n" + str(n) + "_" + location + ".dat"
    df = pd.read_csv(filename, sep="\s+", header=None, skiprows=13)
    df.columns = columns
    calculate_flowangle(df)
    calculate_totalpressure(df)
    calculate_totalrelpressure(df, n, 1)
    return df


def create_plot(df, x, y, color, label, title, xlabel,  ylabel):
    plt.plot(df[x], df[y], label=label, color=color, linewidth=2.5)
    plt.legend(loc='best')
    if title:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    axes = plt.gca()
    axes.set_xlim([-0.5, 1.5])


def create_multi_plot(df, x, y1, y2, color1, color2, label1, label2, title, xlabel, y1label, y2label):
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(y1label)
    ax1.plot(df[x], df[y1], color=color1, label=label1)
    plt.legend(loc="center left")
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel(y2label)  # we already handled the x-label with ax1
    ax2.plot(df[x], df[y2], color=color2, label=label2)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.title(title)
    plt.legend(loc="center right")
    axes = plt.gca()
    axes.set_xlim([-0.5, 1.5])


def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors


def check(mdot1, mdot2, specific_work, p1, p2, V1, Vt1, rho1, U, aoa):
    # percentage error mass flow
    mdot_delta = mdot2 - mdot1
    percentage_error_mdot = mdot_delta/mdot1 * 100
    print("Mass Flow Error: " + str(percentage_error_mdot) + "%")

    # percentage error specific work
    Vt2 = U - m.tan(aoa) * V1
    specific_work_theoretical = U * (Vt2-Vt1)
    specific_work_difference = abs(abs(specific_work) - specific_work_theoretical)
    percentage_error_sp_work = specific_work_difference/specific_work_theoretical * 100
    print("Specific Work Error: " + str(percentage_error_sp_work) + "%")

    # percentage error static pressure rise
    V2 = m.sqrt(V1 ** 2 + Vt2 ** 2)
    del_p_theoretical = rho1 * U * (Vt2 - Vt1) - 0.5 * rho1 * (V2 ** 2 - V1 ** 2)
    del_p = p2 - p1
    del_p_difference = abs(del_p) - abs(del_p_theoretical)
    del_p_error = del_p_difference / del_p_theoretical * 100
    print("Static Pressure Rise Error: " + str(del_p_error) + "%")


def delta_ht_sim(df, rot, R):
    vt1 = df["Velocity_y"].iloc[0]
    vt2 = df["Velocity_y"].iloc[-1]
    U = (rot/60) * 2 * m.pi * R
    deltaht = U*(vt2-vt1)
    return deltaht


def delta_ht_theory(rot, R, v_meridional, alpha):
    u = (rot/60) * 2 * m.pi * R
    vt2 = -(u - (v_meridional*m.sin(m.radians(alpha))))
    dht = u * vt2
    return dht


# Plot settings for varying aoa and N
n_range = [100, 200, 500, 1000]
aoa_range = [10, 20, 30, 40, 50]

# # Flow angle for range of aoa and N
# plt.figure(figsize=(10, 8))
# for i in range(0, len(n_range)):
#     flow = create_dataframe(10, n_range[i], "meridional")
#     colors = set_viridis(0.2, 0.8, len(n_range))
#     create_plot(flow, "x", "Flow_Angle", colors[i], r"$\Omega$=" + str(n_range[i]), False, "x [-]", r"Flow Angle [$^\circ$]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()
#
# plt.figure(figsize=(10, 8))
# for i in range(0, len(aoa_range)):
#     flow = create_dataframe(aoa_range[i], 500, "meridional")
#     colors = set_viridis(0.2, 0.8, len(aoa_range))
#     create_plot(flow, "x", "Flow_Angle", colors[i], r"$\alpha$=" + str(aoa_range[i]), False, "x [-]", r"Flow Angle [$^\circ$]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()
#
# # Ptrel for range of aoa and N
# plt.figure(figsize=(10, 8))
# for i in range(0, len(n_range)):
#     flow = create_dataframe(10, n_range[i], "meridional")
#     colors = set_viridis(0.2, 0.8, len(n_range))
#     create_plot(flow, "x", "Total_Rel_Pressure", colors[i], r"$\Omega$=" + str(n_range[i]), False, "x [-]", r"$p_{t_{rel}}$ [Pa]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()
#
# plt.figure(figsize=(10, 8))
# for i in range(0, len(aoa_range)):
#     flow = create_dataframe(aoa_range[i], 500, "meridional")
#     colors = set_viridis(0.2, 0.8, len(aoa_range))
#     create_plot(flow, "x", "Total_Rel_Pressure", colors[i], r"$\alpha$=" + str(aoa_range[i]), False, "x [-]", r"$p_{t_{rel}}$ [Pa]")
# plt.ylim([102300, 105000])
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()
#
# # Pt for range of aoa and N
# plt.figure(figsize=(10, 8))
# for i in range(0, len(n_range)):
#     flow = create_dataframe(10, n_range[i], "meridional")
#     colors = set_viridis(0.2, 0.8, len(n_range))
#     create_plot(flow, "x", "Total_Pressure", colors[i], r"$\Omega$=" + str(n_range[i]), False, "x [-]", r"$p_{t}$ [Pa]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()
#
# plt.figure(figsize=(10, 8))
# for i in range(0, len(aoa_range)):
#     flow = create_dataframe(aoa_range[i], 500, "meridional")
#     colors = set_viridis(0.2, 0.8, len(aoa_range))
#     create_plot(flow, "x", "Total_Pressure", colors[i], r"$\alpha$=" + str(aoa_range[i]), False, "x [-]", r"$p_{t}$ [Pa]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()
#
# # p/pt for range of aoa and N
# plt.figure(figsize=(10, 8))
# for i in range(0, len(n_range)):
#     flow = create_dataframe(10, n_range[i], "meridional")
#     colors = set_viridis(0.2, 0.8, len(n_range))
#     create_plot(flow, "x", "p/pt", colors[i], r"$\Omega$=" + str(n_range[i]), False, "x [-]", r"$p/p_{t}$ [-]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()
#
# plt.figure(figsize=(10, 8))
# for i in range(0, len(aoa_range)):
#     flow = create_dataframe(aoa_range[i], 500, "meridional")
#     colors = set_viridis(0.2, 0.8, len(aoa_range))
#     create_plot(flow, "x", "p/pt", colors[i], r"$\alpha$=" + str(aoa_range[i]), False, "x [-]", r"$p/p_{t}$ [-]")
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()

# # Single test case analysis of all flow characteristics aoa=30, N=500, B=20
# angle = 30
# rot = 500
# colors = set_viridis(0.2, 0.8, 3)
# meridional = create_dataframe(angle, rot, "meridional")
# fig1, ax1 = plt.subplots(figsize=(12, 9.6))
# ln1 = ax1.plot(meridional["x"], meridional["Total_Pressure"], color=colors[0], label=r"$p_{t}$", linewidth=2.5)
# ln2 = ax1.plot(meridional["x"], meridional["Total_Rel_Pressure"], color=colors[1], label=r"$p_{t_{rel}}$", linewidth=2.5)
# plt.grid()
# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# ln3 = ax2.plot(meridional["x"], meridional["p/pt"], color=colors[2], label=r"$p/p_{t}$", linewidth=2.5)
# # make one legend for both lines
# lns = ln1+ln2+ln3
# labs = [l.get_label() for l in lns]
# ax1.legend(lns, labs, loc="best")
# ax1.set_xlabel("x")
# ax1.set_ylabel(r"$p$ [Pa]")
# # ax1.set_ylim([101500,104500])
# ax2.set_ylabel(r"$p/p_{t}$ [-]")  # we already handled the x-label with ax1
# axes = plt.gca()
# axes.set_xlim([-0.5, 1.5])
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
#
# colors = set_viridis(0.2, 0.8, 2)
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
# ax2.set_ylim([42.65, 46.65])
# ax2.set_ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")  # we already handled the x-label with ax1
# axes = plt.gca()
# axes.set_xlim([-0.5, 1.5])
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

# Calculate total enthalpy rise
# R = 1.0
# flow = create_dataframe(angle, rot, "meridional")
# vt1 = flow["Velocity_y"].iloc[0]
# vt2 = flow["Velocity_y"].iloc[-1]
# U = (rot/60) * 2 * m.pi * R
# delta_ht = U * (vt2-vt1)
# delta_ht_theoretical = delta_ht_theory(U, flow["Velocity_x"].iloc[0], angle)
# print("U: " + str(U))
# print("V1: " + str(flow["Velocity_x"].iloc[0]))
# print("p1: " + str(flow["Pressure"].iloc[0]))
# print("rho1: " + str(flow["Density"].iloc[0]))
# print("Vt1: " + str(vt1))
# print("Vt2: " + str(vt2))
# print("dht: " + str(delta_ht))
# print("Cm: " + str(flow["Velocity_x"].iloc[0]))
# print("dht_theory: " + str(delta_ht_theoretical))


# Compare total enthalpy rise for all test cases
total_enthalpy = pd.DataFrame(columns=["AoA", "N", "Simulation", "Theoretical", "Delta", "Percentage"])
total_relative_enthalpy = pd.DataFrame(columns=["AoA", "N", "Simulation", "Theoretical", "Delta", "Percentage"])
for i in range(0, len(aoa_range)):
    row_input_ht = np.zeros(6)
    row_input_htrel = np.zeros(6)
    for j in range(0, len(n_range)):
        filename = "aoa" + str(aoa_range[i]) + "_n" + str(n_range[j]) + "_meridional.dat"
        exists = os.path.exists(filename)
        if exists:
            flow = create_dataframe(aoa_range[i], n_range[j], "meridional")
            dht_sim = delta_ht_sim(flow, n_range[j], 1.0)
            dht_th = delta_ht_theory(n_range[j], 1.0, flow["Velocity_x"].iloc[0], aoa_range[i])
            delta_sim_th = abs(abs(dht_sim) - abs(dht_th))
            delta_sim_th_percentage = delta_sim_th / abs(dht_sim) * 100
            row_input_ht[0] = aoa_range[i]
            row_input_ht[1] = n_range[j]
            row_input_ht[2] = dht_sim
            row_input_ht[3] = dht_th
            row_input_ht[4] = delta_sim_th
            row_input_ht[5] = delta_sim_th_percentage
            new_row = pd.DataFrame(data=[row_input_ht], columns=["AoA", "N", "Simulation", "Theoretical", "Delta", "Percentage"])
            total_enthalpy = total_enthalpy.append(new_row, ignore_index=True)

print(total_enthalpy)

# Plot of simulated vs analytical values of pressure rise
x = np.arange(len(total_enthalpy["AoA"]))  # the label locations
width = 0.75  # the width of the bars
colors = set_viridis(0.0, 1.0, 100)

fig, ax = plt.subplots(figsize=(10, 8))
rects1 = ax.bar(x, total_enthalpy["Percentage"], width, color="k")
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel(r"Difference Simulation and Theory [$\%$]")
ax.set_xticks(x)
ax.set_xticklabels((r"$\alpha10, \Omega100$", r"$\alpha10, \Omega200$", r"$\alpha10, \Omega500$", r"$\alpha10, \Omega1000$", r"$\alpha20, \Omega200$", r"$\alpha20, \Omega500$", r"$\alpha20, \Omega1000$", r"$\alpha30, \Omega200$", r"$\alpha30, \Omega500$", r"$\alpha30, \Omega1000$", r"$\alpha40, \Omega500$", r"$\alpha40, \Omega1000$", r"$\alpha50, \Omega500$", r"$\alpha50, \Omega1000$"), rotation=70, fontsize=10)
ax.grid()
ax.set_axisbelow(True)

fig.tight_layout()

# Compare simulation and theoretical total enthalpy rise
colors = set_viridis(0.2, 0.8, 2)
plt.figure(figsize=(10,8))
plt.plot(-total_enthalpy["Theoretical"], linewidth=2.5, color=colors[0], marker="^", label="Theoretical")
plt.plot(-total_enthalpy["Simulation"], linewidth=2.5, color=colors[1], marker="^", label="Simulation")
plt.ylabel(r"$h_{t}$ [$J$]")
locs, labels = plt.xticks()
plt.xticks(np.arange(0, 13, step=1))
# plt.xticks(np.arange(14), (r"$\alpha=10, N=100$", r"$\alpha=10, N=200$", r"$\alpha=10, N=500$", r"$\alpha=10, N=1000$", r"$\alpha=20, N=200$", r"$\alpha=20, N=500$", r"$\alpha=20, N=1000$", r"$\alpha=30, N=200$", r"$\alpha=30, N=500$", r"$\alpha=30, N=1000$", r"$\alpha=40, N=500$", r"$\alpha=40, N=1000$", r"$\alpha=50, N=500$", r"$\alpha=50, N=1000$"), rotation=70, fontsize=8)
plt.xticks(np.arange(14), (r"$\alpha10, \Omega100$", r"$\alpha10, \Omega200$", r"$\alpha10, \Omega500$", r"$\alpha10, \Omega1000$", r"$\alpha20, \Omega200$", r"$\alpha20, \Omega500$", r"$\alpha20, \Omega1000$", r"$\alpha30, \Omega200$", r"$\alpha30, \Omega500$", r"$\alpha30, \Omega1000$", r"$\alpha40, \Omega500$", r"$\alpha40, \Omega1000$", r"$\alpha50, \Omega500$", r"$\alpha50, \Omega1000$"), rotation=70, fontsize=10)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()



# ### Plot used in thesis ###
# # Compare relative total pressure values for different mesh sizes for numerical error
# colors = set_viridis(0.2, 0.8, 3)
# coarse = create_dataframe(40, 500, "meridional_coarse")
# refined = create_dataframe(40, 500, "meridional_refined")
# locally_refined = create_dataframe(40, 500, "meridional_locallyrefined")
# calculate_totalrelpressure(coarse, 500, 1.0)
# calculate_totalrelpressure(refined, 500, 1.0)
# calculate_totalrelpressure(locally_refined, 500, 1.0)
# plt.figure(figsize=(10, 8))
# plt.plot(coarse.x, coarse.Total_Rel_Pressure, label="9560 Cells", linewidth=2.5, color=colors[0])
# plt.plot(refined.x, refined.Total_Rel_Pressure, label="11220 Cells", linewidth=2.5, color=colors[1])
# plt.plot(locally_refined.x, locally_refined.Total_Rel_Pressure, label="16200 Cells", linewidth=2.5, color=colors[2])
# # plt.plot(coarse.x, coarse.Total_Pressure, label="pt-coarse")
# # plt.plot(refined.x, refined.Total_Pressure, label="pt-refined")
# # plt.plot(locally_refined.x, locally_refined.Total_Pressure, label="pt-loc-refined")
# # plt.title("Pressure over BF Domain (Coarse vs. Refined)")
# plt.legend(loc="best")
# plt.xlabel("x [-]")
# plt.ylabel(r"$p_{t_{rel}}$ [Pa]")
# plt.xlim([-0.1, 0.25])
# plt.ylim([103650, 103750])
# plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
# plt.grid()


plt.show()
