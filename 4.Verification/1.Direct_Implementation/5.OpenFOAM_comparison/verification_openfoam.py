# Comparison of openFOAM and SU2 results for stator and rotor test case
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

def calculations_su2(df):
    vel_x = df["Momentum_x"] / df["Density"]
    vel_y = df["Momentum_y"] / df["Density"]
    df["Velocity_x"] = vel_x
    df["Velocity_y"] = vel_y
    df["Flow_Angle"] = np.arctan2(vel_y, vel_x) * 180 / m.pi
    df["p_norm"] = df["Pressure"] - df["Pressure"].iloc[0]
    vel = (vel_x * vel_x + vel_y * vel_y)**0.5
    df["Total_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * vel * vel
    df["pt_norm"] = df["Total_Pressure"] - df["Total_Pressure"].iloc[0]
    df["p/pt"] = df["Pressure"] / df["Total_Pressure"]
def calculations_openfoam(df, init):
    dp = init - df["Pressure"].iloc[0] * 1.225
    df["Pressure_new"] = df["Pressure"] * 1.225 + dp
    vel_x = df["Velocity_x"]
    vel_y = df["Velocity_y"]
    vel = (vel_x * vel_x + vel_y * vel_y) ** 0.5
    df["Total_Pressure"] = df["Pressure_new"] + 0.5 * 1.225 * vel * vel
    # df["pt_norm"] = df["Total_Pressure"] - df["Total_Pressure"].iloc[0]
    df["Flow_Angle"] = np.arctan2(df["Velocity_y"], df["Velocity_x"]) * 180 / m.pi
    # df["p_norm"] = ( df["Pressure"] - df["Pressure"].iloc[0] ) * 1.225
    df["p/pt"] = df["Pressure_new"] / df["Total_Pressure"]
def calculate_totalrelpressure(df, N, R, type):
    vel_tangential = (-N/60.0) * 2 * m.pi * R
    df["Tangential_Vel"] = vel_tangential
    vel_x = df["Velocity_x"]
    vel_y = df["Velocity_y"] - vel_tangential
    df["Rel_Velocity_y"] = vel_y
    rel_vel = (vel_x * vel_x + vel_y * vel_y) ** 0.5
    df["Relative_Vel"] = rel_vel
    if type == "OF":
        df["Total_Rel_Pressure"] = df["Pressure_new"] + 0.5 * 1.225 * rel_vel * rel_vel
    else:
        df["Total_Rel_Pressure"] = df["Pressure"] + 0.5 * df["Density"] * rel_vel * rel_vel
def set_viridis(start, stop, num_lines):
    cm_subsection = np.linspace(start, stop, num_lines)
    colors = [cm.viridis(x) for x in cm_subsection]
    return colors

# Columns for slices
columns_su2 = ["x", "y", "Density", "Momentum_x", "Momentum_y", "Energy", "Pressure", "Temperature", "Mach",
           "Cp", "Roe_Dissipation"]
columns_openfoam = ["Velocity_x", "Velocity_y", "Velocity_z", "Pressure", "vtkvalidpointmark", "arc_length", "x", "y", "z", "Point_ID"]

# Read files
rotor_su2_filename = "rotor_meridional.dat"
rotor_openfoam_filename = "rotor.csv"
stator_su2_filename = "stator_meridional.dat"
stator_openfoam_filename = "stator.csv"

rotor_su2 = pd.read_csv(rotor_su2_filename, sep="\s+", header=None, skiprows=13)
rotor_su2.columns = columns_su2
stator_su2 = pd.read_csv(stator_su2_filename, sep="\s+", header=None, skiprows=13)
stator_su2.columns = columns_su2
rotor_openfoam = pd.read_csv(rotor_openfoam_filename)
rotor_openfoam.columns = columns_openfoam
stator_openfoam = pd.read_csv(stator_openfoam_filename, skiprows=1)
stator_openfoam.columns = columns_openfoam

# Calculate flow angle and total pressure
calculations_su2(stator_su2)
calculations_su2(rotor_su2)
calculations_openfoam(stator_openfoam, stator_su2["Pressure"].iloc[0])
calculations_openfoam(rotor_openfoam, rotor_su2["Pressure"].iloc[0])
calculate_totalrelpressure(rotor_su2, 150, 1.5, "SU2")
calculate_totalrelpressure(rotor_openfoam, 150, 1.5, "OF")


## Plots comparing stator openFOAM and SU2
colors = set_viridis(0.2, 0.8, 2)
#1) pt and p/pt
fig1, ax1 = plt.subplots(figsize=(12, 9.6))
ln1 = ax1.plot(stator_su2["x"], stator_su2["Total_Pressure"], color=colors[0], label=r"$p_{t,SU2}$", linewidth=2.5)
ln2 = ax1.plot(stator_openfoam["x"], stator_openfoam["Total_Pressure"], color=colors[0], label=r"$p_{t,OF}$", linewidth=2.5, linestyle="--")
plt.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ln3 = ax2.plot(stator_su2["x"], stator_su2["p/pt"], color=colors[1], label=r"$(p/p_{t})_{SU2}$", linewidth=2.5)
ln4 = ax2.plot(stator_openfoam["x"], stator_openfoam["p/pt"], color=colors[1], label=r"$(p/p_{t})_{OF}$", linewidth=2.5, linestyle="--")
# make one legend for both lines
lns = ln1+ln2+ln3+ln4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)
ax1.set_xlabel("x")
ax1.set_ylabel(r"$p$ [Pa]")
ax1.set_ylim([101400, 101600])
ax2.set_ylabel(r"$p/p_{t}$ [-]")  # we already handled the x-label with ax1
axes = plt.gca()
axes.set_xlim([-0.5, 1.0])
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.5, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

#2) flow angle and mom_x
fig2, ax1 = plt.subplots(figsize=(12, 9.6))
ln1 = ax1.plot(stator_su2["x"], stator_su2["Flow_Angle"], color=colors[0], label="Flow Angle SU2", linewidth=2.5)
ln2 = ax1.plot(stator_openfoam["x"], stator_openfoam["Flow_Angle"], color=colors[0], label="Flow Angle OF", linewidth=2.5, linestyle="--")
plt.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ln3 = ax2.plot(stator_su2["x"], stator_su2["Momentum_x"], color=colors[1], label=r"$(\rho u)_{SU2}$", linewidth=2.5)
ln4 = ax2.plot(stator_openfoam["x"], stator_openfoam["Velocity_x"]*1.225, color=colors[1], label=r"($\rho u)_{OF}$", linewidth=2.5, linestyle="--")
# make one legend for both lines
lns = ln1+ln2+ln3+ln4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)
ax1.set_xlabel("x")
ax1.set_ylabel(r"Flow Angle [$^\circ$]")
ax2.set_ylim([20.36, 21.36])
ax2.set_ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")  # we already handled the x-label with ax1
axes = plt.gca()
axes.set_xlim([-0.5, 1.0])
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.5, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

## Plots comparing rotor openFOAM and SU2
#1) pt and p/pt
fig3, ax1 = plt.subplots(figsize=(12, 9.6))
ln1 = ax1.plot(rotor_su2["x"], rotor_su2["Total_Pressure"], color=colors[0], label=r"$p_{t,SU2}$", linewidth=2.5)
ln2 = ax1.plot(rotor_openfoam["x"], rotor_openfoam["Total_Pressure"], color=colors[0], label=r"$p_{t,OF}$", linewidth=2.5, linestyle="--")
plt.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ln3 = ax2.plot(rotor_su2["x"], rotor_su2["p/pt"], color=colors[1], label=r"$(p/p_{t})_{SU2}$", linewidth=2.5)
ln4 = ax2.plot(rotor_openfoam["x"], rotor_openfoam["p/pt"], color=colors[1], label=r"$(p/p_{t})_{OF}$", linewidth=2.5, linestyle="--")
# make one legend for both lines
lns = ln1+ln2+ln3+ln4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)
ax1.set_xlabel("x")
ax1.set_ylabel(r"$p$ [Pa]")
# ax1.set_ylim([102024,102044])
ax2.set_ylabel(r"$p/p_{t}$ [-]")  # we already handled the x-label with ax1
ax2.set_ylim([0.9962, 0.9974])
axes = plt.gca()
axes.set_xlim([-0.5, 1.0])
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.5, color="black", linestyle="--", linewidth=1, dashes=(10, 15))

#2) flow angle and mom_x
fig4, ax1 = plt.subplots(figsize=(12, 9.6))
ln1 = ax1.plot(rotor_su2["x"], rotor_su2["Flow_Angle"], color=colors[0], label="Flow Angle SU2", linewidth=2.5)
ln2 = ax1.plot(rotor_openfoam["x"], rotor_openfoam["Flow_Angle"], color=colors[0], label="Flow Angle OF", linewidth=2.5, linestyle="--")
plt.grid()
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ln3 = ax2.plot(rotor_su2["x"], rotor_su2["Momentum_x"], color=colors[1], label=r"$(\rho u)_{SU2}$", linewidth=2.5)
ln4 = ax2.plot(rotor_openfoam["x"], rotor_openfoam["Velocity_x"]*1.225, color=colors[1], label=r"($\rho u)_{OF}$", linewidth=2.5, linestyle="--")
# make one legend for both lines
lns = ln1+ln2+ln3+ln4
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)
ax1.set_xlabel("x")
ax1.set_ylabel(r"Flow Angle [$^\circ$]")
ax2.set_ylim([20, 32])
ax2.set_ylabel(r"$\rho u$ [$\frac{kg}{m^2 \cdot s}$]")  # we already handled the x-label with ax1
axes = plt.gca()
axes.set_xlim([-0.5, 1.0])
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.5, color="black", linestyle="--", linewidth=1, dashes=(10, 15))


plt.figure(figsize=(10,8))
plt.plot(rotor_su2["x"], rotor_su2["Total_Rel_Pressure"], color=colors[0], label='SU2', linewidth=2.5)
plt.plot(rotor_openfoam["x"], rotor_openfoam["Total_Rel_Pressure"], color=colors[1], label='openFOAM', linewidth=2.5)
plt.ylim([101800, 101850])
plt.xlim([-0.5, 1.0])
plt.xlabel("x [-]")
plt.ylabel(r"$p_{t_{rel}}$ [Pa]")
plt.legend()
plt.grid()
plt.axvline(x=0.0, color="black", linestyle="--", linewidth=1, dashes=(10, 15))
plt.axvline(x=0.5, color="black", linestyle="--", linewidth=1, dashes=(10, 15))


# Calculate total enthalpy rise for rotor
delta_ht_rotor_su2 = 23.56 * ( rotor_su2["Velocity_y"].iloc[-1] - rotor_su2["Velocity_y"].iloc[0] )
delta_ht_rotor_openfoam = 23.56 * ( rotor_openfoam["Velocity_y"].iloc[-1] - rotor_openfoam["Velocity_y"].iloc[0] )

print("Total enthalpy rise rotor SU2: " + str(delta_ht_rotor_su2))
print("Total enthalpy rise rotor openfoam: " + str(delta_ht_rotor_openfoam))

plt.show()
