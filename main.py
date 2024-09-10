import math
from tabulate import tabulate

# Constants
g_0 = 9.8  # m/s^2, gravitational acceleration
payload_weight = 362.87  # Payload weight in kg
crew = 1
crew_weight = crew * 90
SPEED_SOUND_1 = 296  # m/s, speed of sound at altitude 1
SPEED_SOUND_2 = 296  # m/s, speed of sound at altitude 2
H2_A1_ratio = 142/43
cryogenic_h2_density = 71

# Assumptions
MAC = 1.6  # meter, mean aerodynamic chord
WINGSPAN = 13  # meter, wingspan
K_LD = 15.5  # empirical constant for L/D max
F_TAXI = 1 - (1 - 0.99) / H2_A1_ratio  # fuel fraction for taxi and take-off
F_DECENT = 1 - (1 - 0.952) / H2_A1_ratio  # fuel fraction for descent
F_LANDING = 1 - (1 - 0.995) / H2_A1_ratio  # fuel fraction for landing
SFC_cruise = 19 / H2_A1_ratio / 1000000  # Specific fuel consumption, kg/s
SFC_loiter = SFC_cruise * (1 - 0.128)
S = 6.26  # Wetted area, m^2

# Flight data
cruise_distance_1 = 600 * 1852  # 600 NM converted to meters (Cruise 1)
cruise_distance_2 = 200 * 1852  # 200 NM converted to meters (Cruise 2, diversion)
cruise_speed = 360 * 0.514444  # 360 knots converted to m/s
loiter_time_1 = 30  # 30 minutes loiter (first phase)
loiter_time_2 = 45  # 45 minutes loiter (second phase)

# Empty weight constants
a_ewf = 1.3  # constant for empty weight calculation
b_ewf = -0.096  #  exponent for empty weight calculation


# Tank-penalty constants
a_tank = 16.53
b_tank = -4.818
c_tank = 0.4216
d_tank = -0.1132

# Functions

# Tank penalty function
def tank_penalty(A_t):
    """Calculate the fuel-to-total mass ratio G_i for the tank."""
    G_i = a_tank * math.exp(b_tank * A_t) + c_tank * math.exp(d_tank * A_t)
    return G_i

# Calculate wing area
def calc_wing_area(MAC, wingspan):
    return MAC * wingspan

# Calculate aspect ratio
def cal_AP_wet(Wingspan, A_wing):
    return Wingspan**2 / A_wing

# Calculate L/D max
def L_D_MAX(AP_wet, S, K_LD):
    return K_LD * math.sqrt(AP_wet / S)

# Weight fraction for climb
def weight_fraction_climb(speed, a):
    M = speed / a  # Mach number
    return 1.0065 - 0.0325 * M  # weight fraction during climb

# Weight fraction for cruise
def weight_fraction_cruise(SFC, L_D, Distance, Speed):
    R = Distance  # Cruise range, m
    g = g_0
    V = Speed  # Cruise speed, m/s
    return math.exp((-R * SFC * g) / (V * (L_D * 0.886)))

# Weight fraction for loiter
def weight_fraction_loiter(SFC, L_D, Time, Speed):
    E = Time * 60  # Time converted from minutes to seconds
    g = g_0
    return math.exp((-E * SFC * g) / (L_D))

def solve_tank_mass(fuel_mass, G_i):
    """Calculate the tank mass given the fuel mass and the fuel-to-total mass ratio G_i."""
    if G_i <= 0 or G_i >= 1:
        raise ValueError("G_i must be between 0 and 1.")
    tank_mass = (1 / G_i - 1) * fuel_mass
    return tank_mass

# Empty weight fraction formula
def empty_weight_fraction(w0, a, b, GI, total_fuel_weight_fraction):
    """Calculate the empty weight fraction using the formula: EWF = a * W0^b + Wtank / W0."""
    if w0 <= 0:
        raise ValueError("W0 must be positive")
    w_tank = solve_tank_mass((w0 * total_fuel_weight_fraction), GI)
    ewf = a * (w0 ** b) +  (w_tank / w0)
    return ewf, w_tank

    
    tank_mass = (fuel_mass / G_i) - fuel_mass
    return tank_mass

# Calculating total weight fraction (Wf/W0)
def total_weight_fraction(MAC, WINGSPAN, cruise_distance_1, cruise_distance_2, cruise_speed, loiter_time_1, loiter_time_2):
    # Wing and AP_wet calculations
    A_wing = calc_wing_area(MAC, WINGSPAN)
    AR_wet = cal_AP_wet(WINGSPAN, A_wing)
    L_D_max_value = L_D_MAX(AR_wet, S, K_LD)

    # Weight fractions for mission phases
    Wf_W0_taxi = F_TAXI
    Wf_W0_climb_1 = weight_fraction_climb(cruise_speed, SPEED_SOUND_1)
    Wf_W0_cruise_1 = weight_fraction_cruise(SFC_cruise, L_D_max_value, cruise_distance_1, cruise_speed)
    Wf_W0_descent_1 = F_DECENT
    Wf_W0_loiter_1 = weight_fraction_loiter(SFC_loiter, L_D_max_value, loiter_time_1, cruise_speed)
    Wf_W0_landing_1 = F_LANDING

    Wf_W0_climb_2 = weight_fraction_climb(cruise_speed, SPEED_SOUND_2)
    Wf_W0_cruise_2 = weight_fraction_cruise(SFC_cruise, L_D_max_value, cruise_distance_2, cruise_speed)
    Wf_W0_descent_2 = F_DECENT
    Wf_W0_loiter_2 = weight_fraction_loiter(SFC_loiter, L_D_max_value, loiter_time_2, cruise_speed)
    Wf_W0_landing_2 = F_LANDING

    # Total weight fraction 
    Wf_W0_total = (Wf_W0_taxi * Wf_W0_climb_1 * Wf_W0_cruise_1 * Wf_W0_descent_1 * Wf_W0_loiter_1 * Wf_W0_landing_1 * Wf_W0_climb_2 * Wf_W0_cruise_2 * Wf_W0_descent_2 *Wf_W0_loiter_2 * Wf_W0_landing_2)

    return 1.06*(1 - Wf_W0_total), L_D_max_value, WINGSPAN, AR_wet, MAC

# Iteration function to find correct W0
def find_correct_w0(W0_initial, payload_weight, a, b, tolerance=1e-6, max_iterations=1000, A_t=3):
    W0 = W0_initial
    G_i = tank_penalty(A_t)  # Calculate tank penalty
    
    # Calculate total weight fraction
    total_fuel_weight_fraction, L_D_max_value, wingspan, AR_wet, MAC_value = total_weight_fraction(
            MAC, WINGSPAN, cruise_distance_1, cruise_distance_2, cruise_speed, loiter_time_1, loiter_time_2
        )

    for iteration in range(max_iterations):
        # Calculate empty weight fraction
        EWF, tank_mass = empty_weight_fraction(W0, a, b, G_i, total_fuel_weight_fraction)

        # Update W0 to account for tank penalty
        W0_new = (crew_weight + payload_weight) / (1 - EWF - total_fuel_weight_fraction)

        # Check for convergence
        if abs(W0_new - W0) < tolerance:
            # Return converged W0 and other results
            return W0_new, L_D_max_value, wingspan, AR_wet, MAC_value, G_i, tank_mass, total_fuel_weight_fraction, EWF

        # Update W0 for next iteration
        W0 = W0_new

    raise ValueError("Did not converge within the maximum number of iterations")

# Initial guess for W0
W0_initial = 10000  # Example initial guess in kg

# Call the iteration function
A_t = 2.32  
correct_W0, L_D_max_value, wingspan, AR_wet, MAC_value, G_i, tank_mass, total_fuel_weight_fraction, empty_weight_fraction = find_correct_w0(
    W0_initial, payload_weight, a_ewf, b_ewf, A_t=A_t
)

fuel_mass = correct_W0 * total_fuel_weight_fraction

tank_volume = fuel_mass / cryogenic_h2_density

empty_weight = empty_weight_fraction * correct_W0

# Prepare data for tabulation
data = [
    ["Takeoff Weight (W0)", f"{correct_W0:.2f} kg"],
    ["Empty Weight (WE)", f"{empty_weight:.2f} kg"],
     ["Fuel Weight", f"{fuel_mass:.2f} kg"],
    ["L/D_max", f"{L_D_max_value:.2f}"],
    ["Wingspan", f"{wingspan:.2f} meters"],
    ["Aspect Ratio (AR_wet)", f"{AR_wet:.2f}"],
    ["Mean Aerodynamic Chord (MAC)", f"{MAC_value:.2f} meters"],
    ["Tank Penalty (G_i)", f"{G_i:.4f}"],
    ["Tank Mass", f"{tank_mass:.2f} kg"],
    ["Tank Volume", f"{tank_volume:.2f} m^3"]
]

# Print the table
table = tabulate(data, headers=["Parameter", "Value"], tablefmt="grid")
print(table)
