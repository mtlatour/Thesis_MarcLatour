# Marc Latour - Thesis - 12/09/2019
# Visualize verification part 2
import pandas as pd
import numpy as np
import math as m
import matplotlib.pyplot as plt
from matplotlib import cm

#Latex template for plots
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
params = {'text.usetex': True,
          'font.size': 16,
          'font.family': 'lmodern',
          'text.latex.unicode': True,
          }
plt.rcParams.update(params)

# Columns for slices
columns = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]

rho_inf = 1.225
p_inf = 101325
v_inf = 34.05

# Functions for creating and plotting dataframes
def calculate_flowangle(df):
    vel_x = df["Momentum_x"] / df["Density"]
    vel_y = df["Momentum_y"] / df["Density"]
    df["Flow_Angle"] = np.arctan2(vel_y, vel_x) * 180 / m.pi
    return df
def calculate_totalpressure(df):
    vel_x = df["Momentum_x"] / df["Density"]
    vel_y = df["Momentum_y"] / df["Density"]
    vel = (vel_x * vel_x + vel_y * vel_y)**(0.5)
    df["Total_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * vel * vel
    return df
def calculate_cp(df):
    denom = 0.5 * rho_inf * v_inf * v_inf
    df["Cp_calc"] = (df["Pressure"] - p_inf) / denom
    return df
def create_dataframe(blades, location):
    filename = "blades" + str(blades) + "_" + location + ".dat"
    df = pd.read_csv(filename, sep="\s+", header=None, skiprows=13)
    df.columns = columns
    calculate_flowangle(df)
    calculate_totalpressure(df)
    calculate_cp(df)
    return df
def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors


# Creating dataframes and plotting data
bladesrange = [1, 5, 10, 20]
colors = set_viridis(0.2,0.8,len(bladesrange))
for i in range(0, len(bladesrange)):
    meridional = create_dataframe(bladesrange[i], "meridional")
    plt.plot(meridional["x"], meridional["Flow_Angle"], color=colors[i], label="Blades: " + str(bladesrange[i]))

plt.legend(loc="best")
plt.xlim([-0.5, 1.5])
plt.grid(linewidth=0.3)
plt.title("Flow Angle Variation over Body Force Domain")
plt.xlabel("x")
plt.ylabel("Flow Angle [deg]")
plt.show()


