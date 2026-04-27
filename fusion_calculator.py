import numpy as np
import matplotlib.pyplot as plt
import math
print("=== Fusion Reaction Calculator ===")

# Pyhysical constants
K_B = 1.380649e-23      # Boltzmann constant (J/K)
EV_TO_J =1.602e-19      # eV - Joules
KEV_TO_J = 1.602e-16    # keV - Joules

# Get user input
user_temp = float(input("Enter temperature(million C): "))

# Convert million C to keV
T_keV = user_temp / 11.6

# Temperature range for graph
max_temp = max(200, user_temp + 50)
temperatures = np.linspace(1, max_temp, 200)

density = float(input("Enter particle density(particles per m3, e.g. 1e20). "))
confinement_time = float(input("Enter confinement time(seconds). "))

ntau = density  * confinement_time

LAWSON_DT = 1e20 * (10 / T_keV) 

LAWSON_THRESHOLD = 1e21 # simplified

# S-factor models

def S_dt(T_keV):
    # D-T is relatively flat
    return 1e-28

def S_dd(T_keV):
    return 1e-30

def S_pb(T_keV):
    return 1e-33

# Reactivity

def reactivity_saddle(T_keV, Z1, Z2, mu, S_func):
   
    T = max(T_keV, 1e-6)
   
    S_val = S_func(T)

    B = 20 * (Z1 * Z2)**(2/3) * (mu**(1/3)) 

    return S_val * (T**2.5) * np.exp(-B / np.sqrt(T + 1e-6)) 

# Approx Values

# D-T
dt_params = (1, 1, 3.34e-27, S_dt)

# D-D
dd_params = (1, 1, 6.68e-27, S_dd)

# p_B11
pb_params = (1, 5, 1.83e-26, S_pb) 

# Calculate values at user temperature 
dt_value = reactivity_saddle(T_keV, *dt_params)
dd_value = reactivity_saddle(T_keV, *dd_params)
pb_value = reactivity_saddle(T_keV, *pb_params)

# Energy per reaction (MeV)
DT_ENERGY = 17.6
DD_ENERGY = 4.0
PB_ENERGY = 8.7

# Convert MeV to Joules
MEV_TO_J = 1.602e-13

# Energy output (Joules, relative)
dt_energy = dt_value * DT_ENERGY * MEV_TO_J
dd_energy = dd_value * DD_ENERGY * MEV_TO_J
pb_energy = pb_value * PB_ENERGY * MEV_TO_J

# Confinement efficiency 
efficiency = float (input("Enter confinement efficiency (0-1): "))

dt_energy *= efficiency
dd_energy *= efficiency 
pb_energy *= efficiency

# Density
n = 1e20  # paricles per m3

# Energy per reaction (Joules)
E_DT = 17.6e6 * 1.602e-19
E_DD = 3.6e6 * 1.602e-19
E_PB = 8.7e6 * 1.602e-19

# Temperature in keV for full range
T_keV_curve = temperatures / 11.6

# Radiation losses 
radiation_loss = 1e-37 * n**2 * np.sqrt(T_keV_curve)

# Plasma composition 
fuel_mix = float(input("Enter Fuel mix (0-1). "))

total_energy = (fuel_mix * dt_energy * (1 - fuel_mix) * dd_energy)

 

# Print results
print (f"\nAt {user_temp} million C:") 
print (f"D-T fusion rate: {dt_value:.3e} sigma-v")
print (f"D-D fusion rate: {dd_value:.3e} sigma-v")
print (f"p-B11 fusion rate: {pb_value:.3e} sigma-v")
print ("\n=== Energy Output (relative Joules) ===")
print (f"D-T energy: {dt_energy:.6e} J")
print (f"D-D energy: {dd_energy:.6e} J") 
print (f"p-B11 energy: {pb_energy:.6e} J")
print ("\n=== Lawson Criterion ===")
print (f"Temperature: {T_keV:.2f} keV")
print (f"nt value: {ntau:.3e}")
if ntau >= LAWSON_DT and 10 <= T_keV <= 20:
    print ("Fusion conditions ideal for D-T")
elif ntau >= LAWSON_DT:
    print ("Fusion possible, but temperature not optimal")
else: 
    print ("Fusion not achieved") 
print (f"Radiation loss min:{radiation_loss.min():.3e}")
print (f"Radiation loss max:{radiation_loss.max():.3e}")
print ("\n===Power Balance===")
print (f"Fusion power: {total_energy:.3e}")
radiation_loss_point = 1e-37 * n**2 * np.sqrt(T_keV)
print (f"Losses: {radiation_loss_point:.3e}")
net_power_point = n**2 * dt_value * E_DT - radiation_loss_point
if net_power_point > 0:
    print ("Net positive fusion acieved")
else: 
    print ("Net negative (losses dominate)")

# Calculate curves
dt_rates = np.array([reactivity_saddle(T, *dt_params) for T in T_keV_curve])
dd_rates = np.array([reactivity_saddle(T, *dd_params) for T in T_keV_curve])
pb_rates = np.array([reactivity_saddle(T, *pb_params) for T in T_keV_curve]) 
print ("dt_rates min:", dt_rates.min())
print ("dt_rates max:", dt_rates.max())

dt_power = n**2 * dt_rates * E_DT
dd_power = n**2 * dd_rates * E_DD
pb_power = n**2 * pb_rates * E_PB

total_power = dt_power + dd_power + pb_power
net_power = total_power - radiation_loss 

# Energy curves
dt_energy_curve = dt_rates * density**2 * confinement_time * DT_ENERGY * MEV_TO_J * efficiency
dd_energy_curve = dd_rates * density**2 * confinement_time * DD_ENERGY * MEV_TO_J * efficiency
pb_energy_curve = pb_rates * density**2 * confinement_time * PB_ENERGY * MEV_TO_J * efficiency

# Temperature in keV for full range
T_keV_curve = temperatures / 11.6

# Convert keV to million C
T_million_C = T_keV_curve * 11.6

# Radiation loss curve
radiation_curve = 1e-38 * density**2 * np.sqrt(T_keV_curve)

# Total fusion (based on mix)
fusion_total_curve = (dt_energy_curve + dd_energy_curve + pb_energy_curve) 

# Net power
net_power_curve = fusion_total_curve - radiation_curve

print ("n\=== DEBUG GET POWER ===")
print ("Min net:", np.min(net_power_curve))
print ("Max net:", np.max(net_power_curve))

# Plotting section
plt.figure(figsize=(12, 7))

#===== TOP GRAPH: FUSION RATES =====
plt.subplot(2, 1, 1)

# Fusion Rates
plt.plot(T_keV_curve, dt_rates, label="D-T Rate", linewidth=2)
plt.plot(T_keV_curve, dd_rates, label="D-D Rate", linewidth=2)
plt.plot(T_keV_curve, pb_rates, label="p-B11 rate", linewidth=2)

# Scatter
plt.scatter(T_keV, dt_value, s=120)
plt.scatter(T_keV, dd_value, s=120)
plt.scatter(T_keV, pb_value, s=120)

# Text Sizes
plt.text(T_keV, dt_value, " D-T", fontsize=9)
plt.text(T_keV, dd_value, " D-D", fontsize=9)
plt.text(T_keV, pb_value, " p-B11", fontsize=9)

# Log Scale
plt.yscale("log")

# Axis
plt.xlim(0, max(T_keV_curve))

# Labels
plt.xlabel("Temperature (keV))")
plt.ylabel("Reactivity")
plt.title("Fusion Reaction Rates")

plt.legend()
plt.grid()

#===== Bottom Graph: Energy + Losses =====
plt.subplot(2, 1, 2)

# Fusion Energy Curves
plt.plot(T_million_C, fusion_total_curve, label="Fusion Power", linewidth=2)

# Radiation Loss
plt.plot(T_million_C, radiation_curve, "--", label="Radiation Loss", linewidth=2)

# Net Power
plt.plot(T_million_C, net_power_curve, linewidth=3, label="Net Power")

# Break even line
zero_index = np.where(np.diff(np.sign(net_power_curve)) !=0)[0]
if len(zero_index) > 0:
    i = zero_index[0]
   
    
    x1, x2 = T_million_C[i], T_million_C[i+1]
    y1, y2 = net_power_curve[i], net_power_curve[i+1]
    break_even_temp = x1 - y1 * (x2 - x1) / (y2 - y1)
 
        
    plt.axvline(break_even_temp, color='red', linestyle='--', linewidth=2)
    y_position = net_power_curve[i]   
    
    plt.text(break_even_temp, y1, " Break-even", color='red', fontsize=10)
    print (f"Break-even at ~{break_even_temp:.2f} million C")
else: 
    print ("No Break-even point found")

# Axis + styling
plt.xlim(0, max_temp)
plt.axhline(0, linestyle=":", linewidth=1)
plt.axvline(user_temp, linestyle=":", linewidth=1)
plt.yscale("log")
plt.xlabel("Temperature (million C)")
plt.ylabel("Power")
plt.title("Fusion Power vs Losses vs Net Output")

plt.legend()
plt.grid()

# Layout Fix
plt.tight_layout()

# Save
plt.savefig("fusion_dashboard.png", dpi=300)

print ("Dashboard saved as fusion_dashboard.png")
 

 
