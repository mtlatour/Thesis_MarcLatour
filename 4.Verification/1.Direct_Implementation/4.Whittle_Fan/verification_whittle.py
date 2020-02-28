# Comparison of Whittle fan 2D section physical blade and body force simulations
# Marc Latour; TU Delft; 03/10/2019

import pandas as pd
import numpy as np
import cmath as m
import matplotlib.pyplot as plt
from matplotlib import cm

#Latex template for plots
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 18,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

# Columns for reading data
columns_bf = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]
columns_blade = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Grid_Velocity_x", "Grid_Velocity_y", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]


# Functions for dataframe calculations
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


def calculate_totalrelpressure(df, N):
    U = (-N / 60.0) * 2 * m.pi * 0.1625
    df["Tangential_Vel"] = U
    vel_x = df["Velocity_x"]
    vel_y = df["Velocity_y"] - U
    rel_vel = (vel_x * vel_x + vel_y * vel_y) ** 0.5
    df["Relative_Vel"] = rel_vel
    df["Total_Rel_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * rel_vel * rel_vel


def calculate_dht(df, N):
    vt1 = df["Velocity_y"].iloc[0]
    vt2 = df["Velocity_y"].iloc[-1]
    U = (-N / 60.0) * 2 * m.pi * 0.1625
    dht = U * (vt2-vt1)
    return dht


def create_dataframe(x, section, N):
    filename = "whittle_" + str(x) + "_" + str(section) + ".dat"
    if x == "bf":
        df = pd.read_csv(filename, sep="\s+", header=None, skiprows=13)
        df.columns = columns_bf
    else:
        df = pd.read_csv(filename, sep="\s+", header=None, skiprows=15)
        df.columns = columns_blade
    calculate_flowangle(df)
    calculate_totalpressure(df)
    calculate_totalrelpressure(df, N)
    return df


def create_dataframe_blade(N):
    whittle_blade_flow_filename = "whittle_blade_n" + str(N) + ".csv"
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
    calculate_totalrelpressure(whittle_blade_averages, rot_speed)
    return whittle_blade_averages


def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors


# Set case rotational speed
rot_speed = 2898
rot_speed_range = [1803, 2362, 2898]

# Load body force and blade flow
whittle_bf_flow = create_dataframe("bf", "n" + str(rot_speed) + "_cells38922", rot_speed)
dht_bf = calculate_dht(whittle_bf_flow, rot_speed)
whittle_blade_flow = create_dataframe_blade(rot_speed)
dht_blade = calculate_dht(whittle_blade_flow, rot_speed)

# Set number of colors
colors = set_viridis(0.2, 0.8, 2)

# Plots comparing bf and blade results
plt.figure(figsize=(10, 8))
plt.plot(whittle_blade_flow.x, whittle_blade_flow.Flow_Angle, label="Physical Blade", color=colors[0], linewidth=2.5)
plt.plot(whittle_bf_flow.x, whittle_bf_flow.Flow_Angle, label="Body Force", color=colors[1], linewidth=2.5)
plt.xlim([-0.025, 0.075])
plt.xlabel("x [-]")
plt.ylabel(r"Flow Angle [$^\circ$]")
plt.legend()
plt.grid()
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.04569681, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

plt.figure(figsize=(10, 8))
plt.plot(whittle_blade_flow.x, whittle_blade_flow.Total_Pressure, label="Physical Blade", color=colors[0], linewidth=2.5)
plt.plot(whittle_bf_flow.x, whittle_bf_flow.Total_Pressure, label="Body Force", color=colors[1], linewidth=2.5)
plt.xlim([-0.025, 0.075])
plt.xlabel("x [-]")
plt.ylabel(r"$p_{t}$ [Pa]")
plt.legend()
plt.grid()
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.04569681, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

plt.figure(figsize=(10, 8))
plt.plot(whittle_blade_flow.x, whittle_blade_flow.Total_Rel_Pressure, label="Physical Blade", color=colors[0], linewidth=2.5)
plt.plot(whittle_bf_flow.x, whittle_bf_flow.Total_Rel_Pressure, label="Body Force", color=colors[1], linewidth=2.5)
plt.xlim([-0.025, 0.075])
plt.ylim([103800, 104200])
plt.xlabel("x [-]")
plt.ylabel(r"$p_{t_{rel}}$ [Pa]")
plt.legend()
plt.grid()
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.04569681, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

plt.figure(figsize=(10, 8))
plt.plot(whittle_blade_flow.x, whittle_blade_flow.Momentum_x, label="Physical Blade", color=colors[0], linewidth=2.5)
plt.plot(whittle_bf_flow.x, whittle_bf_flow.Momentum_x, label="Body Force", color=colors[1], linewidth=2.5)
plt.xlim([-0.025, 0.075])
plt.xlabel("x [-]")
plt.ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")
plt.legend()
plt.grid()
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.04569681, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

# Plot difference in results for all three test cases
colors = set_viridis(0.2, 0.8, len(rot_speed_range))
flowangle_comp = np.zeros((3, 3))
pt_comp = np.zeros((3, 3))
plt.figure(figsize=(12, 9.6))
for case in range(0, len(rot_speed_range)):
    bf_flow = create_dataframe("bf", "n" + str(rot_speed_range[case]) + "_cells38922", rot_speed_range[case])
    blade_flow = create_dataframe_blade(rot_speed_range[case])
    plt.plot(bf_flow.x, bf_flow.Flow_Angle, label="BF-" + str(rot_speed_range[case]), color=colors[case], linewidth=2.5)
    plt.plot(blade_flow.x, blade_flow.Flow_Angle, label="Blade-" + str(rot_speed_range[case]), color=colors[case], linewidth=2.5, linestyle="--")
    flowangle_comp[case, 0] = bf_flow["Flow_Angle"].iloc[-1]
    flowangle_comp[case, 1] = blade_flow["Flow_Angle"].iloc[-1]
    flowangle_comp[case, 2] = flowangle_comp[case, 1] - flowangle_comp[case, 0]
    pt_comp[case, 0] = bf_flow["Total_Pressure"].iloc[-1]
    pt_comp[case, 1] = blade_flow["Total_Pressure"].iloc[-1]
    pt_comp[case, 2] = pt_comp[case, 1] - pt_comp[case, 0]

plt.xlabel("x [-]")
plt.ylabel(r"Flow Angle [$^\circ$]")

plt.xlim([-0.025, 0.075])
plt.grid()
plt.legend(loc="best")
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.04569681, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

print(flowangle_comp)
print(pt_comp)

plt.show()
