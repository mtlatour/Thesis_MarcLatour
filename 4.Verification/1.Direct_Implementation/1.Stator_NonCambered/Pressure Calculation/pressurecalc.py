# File calculating pressure change analytically
import cmath as m

# Flat Plae
V1 = 34.3 #m/s
Vt1 = 0 #m/s
plate_angle = 5 * m.pi / 180 #deg
p1 = 101325 #Pa
rho = 1.225 #k/m3
N = 0 #RPM
R = 1.0 #m
U = (N/60.0) * m.pi * 2 * R #m/s

Vt2_beta = m.tan(plate_angle) * V1
Vt2 = U-Vt2_beta
V2 = m.sqrt(V1 ** 2 + Vt2 ** 2)
del_p_st = rho * U * (Vt2-Vt1) - 0.5 * rho * (V2 ** 2 - V1 ** 2)
W2 = V2 - U
specific_work = U * (Vt2-Vt1)
print("Vt1: " + str(Vt1))
print("Vt2: " + str(Vt2))
print("Specific Work: " + str(specific_work.real))

if U > 0:
    specific_work = U * Vt2
    flow_coeff = V1 / U
    work_coeff = specific_work / ((U*U)/2.0)
    work_coeff = work_coeff.real
    print("Flow Coefficient: " + str(flow_coeff))
    print("Work Coefficient: " + str(work_coeff))

print("Change in static pressure: " + str(del_p_st) + "Pa")
print("Set outlet back pressure to: " + str(p1+del_p_st) + "Pa")