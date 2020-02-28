# Marc Latour - Thesis - 13/09/2019
# Visualize verification part 2 - compare with part 1
import pandas as pd
import numpy as np
import math as m
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

# Columns for slices
columns = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]

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
    df["p/pt"] = df["Pressure"] / df["Total_Pressure"]
    return df
def create_dataframe(aoa, blades, location, directory):
    if directory == 1:
        dir = "1.No_Rot_Flat"
    elif directory == 2:
        dir = "2.No_Rot_Cambered"
    elif directory == 3:
        dir = "3.Rot_Flat"
    elif directory == 4:
        dir = "4.Rot_Cambered"
    filename = "/Users/marclatour/Documents/Thesis/Thesis/Analysis/Verification/" + dir + "/aoa" + str(aoa) + "_blade" + str(blades) + "_" + location + ".dat"
    df = pd.read_csv(filename, sep="\s+", header=None, skiprows=13)
    df.columns = columns
    calculate_flowangle(df)
    calculate_totalpressure(df)
    return df
def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors

### Plot used in thesis ###
# Creating dataframes and plotting data
bladesrange = [1, 5, 10, 20]
aoalist = [10, "camb"]
aoalabel = ["flat", "camb"]
linestyles = ['-', '--']
colors = set_viridis(0.2,0.8,len(bladesrange))
plt.figure(figsize=(10,8))
for j in range(0, len(aoalist)):
    for i in range(0, len(bladesrange)):
        meridional = create_dataframe(aoalist[j], bladesrange[i], "meridional", j+1)
        plt.plot(meridional["x"], meridional["Flow_Angle"], color=colors[i], label= "B=" + str(bladesrange[i]) + ", " + aoalabel[j]  , linestyle=linestyles[j], linewidth=2.5)
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
plt.legend(loc="best", prop = {"size": 12})
plt.xlim([-0.5, 1.5])
plt.grid(linewidth=0.3)
plt.xlabel("x")
plt.ylabel("Flow Angle [deg]")


# Aoa 10, blades 5 comparison with cambered airfoil
bladenum = 5
aoa = 10
colors = set_viridis(0.2, 0.8, 2)
meridional_cambered = create_dataframe("camb", bladenum, "meridional", 2)
meridional_flat = create_dataframe(aoa, bladenum, "meridional", 1)
fig, ax1 = plt.subplots(figsize=(12, 9.6))
ln1 = ax1.plot(meridional_cambered["x"], meridional_cambered["Total_Pressure"], color=colors[0], label=r"$p_{t,camb}$", linewidth=2.5, linestyle="--")
ln2 = ax1.plot(meridional_flat["x"], meridional_flat["Total_Pressure"], color=colors[0], label=r"$p_{t,flat}$", linewidth=2.5)
plt.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ln3 = ax2.plot(meridional_cambered["x"], meridional_cambered["p/pt"], color=colors[1], label=r"$(p/p_{t})_{camb}$", linewidth=2.5, linestyle="--")
ln4 = ax2.plot(meridional_flat["x"], meridional_flat["p/pt"], color=colors[1], label=r"$(p/p_{t})_{flat}$", linewidth=2.5)
# make one legend for both lines
lns = ln1+ln2+ln3+ln4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)
ax1.set_xlabel("x")
ax1.set_ylabel(r"$p$ [Pa]")
ax1.set_ylim([102024, 102044])
ax2.set_ylabel(r"$p/p_{t}$ [-]")  # we already handled the x-label with ax1
axes = plt.gca()
axes.set_xlim([-0.5, 1.5])
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))
plt.axvline(x=1.0, color="black", linestyle="--", linewidth=1, dashes=(10,15))






plt.show()


