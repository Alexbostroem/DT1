import math
from config import *  # Import constants

class Tank:
    def __init__(self, AR_tank, shell_thickness):
        self.A_t = AR_tank
        self.t_shell = shell_thickness
    
    def tank_penalty(self):
        """Calculate the fuel-to-total mass ratio G_i for the tank."""
        G_i = a_tank * math.exp(b_tank * self.A_t) + c_tank * math.exp(d_tank * self.A_t)
        return G_i
    
    def calc_tank_mass(self, fuel_mass, G_i):
        """Calculate the tank mass given the fuel mass and the fuel-to-total mass ratio G_i."""
        if G_i <= 0 or G_i >= 1:
            raise ValueError("G_i must be between 0 and 1.")
        return (1 / G_i - 1) * fuel_mass
    
    def calc_tank_volume(self, fuel_mass, cryo_density):
        """Calculate tank volume based on fuel mass and cryogenic density."""
        return fuel_mass / cryo_density

    def tank_dimensions(self, fuel_mass):
        """Calculate the inner and outer dimensions of the tank based on fuel mass."""
        tank_volume_per_tank = self.calc_tank_volume(fuel_mass / 2, cryogenic_h2_density)  # 2 tanks
        D_inner = ( (4 * tank_volume_per_tank) / (math.pi * self.A_t) ) ** (1 / 3)  # Inner diameter
        L_inner = self.A_t * D_inner  # Inner length
        D_outer = D_inner + 2 * self.t_shell  # Outer diameter
        L_outer = L_inner + 2 * self.t_shell  # Outer length
        V_outer_tank = math.pi * (D_outer / 2) ** 2 * L_outer  # Outer volume
        return D_inner, L_inner, D_outer, L_outer, V_outer_tank

class Aircraft:
    def __init__(self, payload_weight, crew_weight, wingspan, MAC):
        self.payload_weight = payload_weight
        self.crew_weight = crew_weight
        self.wingspan = wingspan
        self.MAC = MAC
        self.total_weight_fraction = None
    
    def calc_wing_area(self):
        return self.MAC * self.wingspan
    
    def calc_aspect_ratio(self, wing_area):
        return self.wingspan ** 2 / wing_area
    
    def calc_ld_max(self, AR_wet, S):
        return K_LD * math.sqrt(AR_wet / S)
    
    def weight_fraction_climb(self, speed, sound_speed):
        M = speed / sound_speed
        return 1.0065 - 0.0325 * M

class MissionProfile:
    def __init__(self, cruise_distance_1, cruise_distance_2, cruise_speed, loiter_time_1, loiter_time_2):
        self.cruise_distance_1 = cruise_distance_1
        self.cruise_distance_2 = cruise_distance_2
        self.cruise_speed = cruise_speed
        self.loiter_time_1 = loiter_time_1
        self.loiter_time_2 = loiter_time_2
    
    def weight_fraction_cruise(self, SFC, L_D, distance, speed):
        R = distance  # Cruise range, m
        g = g_0
        V = speed  # Cruise speed, m/s
        return math.exp((-R * SFC * g) / (V * (L_D * 0.886)))

    def weight_fraction_loiter(self, SFC, L_D, time):
        E = time * 60  # Convert minutes to seconds
        g = g_0
        return math.exp((-E * SFC * g) / (L_D))

# Function to calculate total weight fraction for a complete mission
def total_weight_fraction(aircraft, mission, tank, SFC_cruise, SFC_loiter):
    A_wing = aircraft.calc_wing_area()
    AR_wet = aircraft.calc_aspect_ratio(A_wing)
    L_D_max_value = 0.9 * aircraft.calc_ld_max(AR_wet, S)

    Wf_W0_taxi = F_TAXI
    Wf_W0_climb_1 = aircraft.weight_fraction_climb(cruise_speed, SPEED_SOUND_1)
    Wf_W0_cruise_1 = mission.weight_fraction_cruise(SFC_cruise, L_D_max_value, mission.cruise_distance_1, mission.cruise_speed)
    Wf_W0_loiter_1 = mission.weight_fraction_loiter(SFC_loiter, L_D_max_value, mission.loiter_time_1)
    Wf_W0_landing_1 = F_LANDING

    Wf_W0_climb_2 = aircraft.weight_fraction_climb(cruise_speed, SPEED_SOUND_2)
    Wf_W0_cruise_2 = mission.weight_fraction_cruise(SFC_cruise, L_D_max_value, mission.cruise_distance_2, mission.cruise_speed)
    Wf_W0_loiter_2 = mission.weight_fraction_loiter(SFC_loiter, L_D_max_value, mission.loiter_time_2)
    Wf_W0_landing_2 = F_LANDING

    Wf_W0_total = (Wf_W0_taxi * Wf_W0_climb_1 * Wf_W0_cruise_1 * Wf_W0_loiter_1 * Wf_W0_landing_1 *
                   Wf_W0_climb_2 * Wf_W0_cruise_2 * Wf_W0_loiter_2 * Wf_W0_landing_2)

    return 1.06 * (1 - Wf_W0_total), L_D_max_value, AR_wet

def empty_weight_fraction(w0, a, b, GI, total_fuel_weight_fraction):
    """Calculate the empty weight fraction using the formula: EWF = a * W0^b + Wtank / W0."""
    if w0 <= 0:
        raise ValueError("W0 must be positive")
    w_tank = tank.calc_tank_mass(w0 * total_fuel_weight_fraction, GI)
    ewf = a * (w0 ** b) +  (w_tank / w0)
    return ewf, w_tank

# Iteration function to find correct W0
def find_correct_w0(W0_initial, aircraft, mission, tank, tolerance=1e-6, max_iterations=1000):
    W0 = W0_initial
    G_i = tank.tank_penalty()  # Calculate tank penalty
    for iteration in range(max_iterations):
        total_fuel_weight_fraction, L_D_max_value, AR_wet = total_weight_fraction(
            aircraft, mission, tank, SFC_cruise, SFC_loiter
        )
        EWF, tank_mass = empty_weight_fraction(W0, a_ewf, b_ewf, G_i, total_fuel_weight_fraction)
        W0_new = (aircraft.crew_weight + aircraft.payload_weight) / (1 - EWF - total_fuel_weight_fraction)
        if abs(W0_new - W0) < tolerance:
            return W0_new, L_D_max_value, AR_wet, G_i, tank_mass, total_fuel_weight_fraction, EWF
        W0 = W0_new
    raise ValueError("Did not converge within the maximum number of iterations")

# Initialize objects
tank = Tank(A_t, t_shell)
aircraft = Aircraft(payload_weight, crew_weight, WINGSPAN, MAC)
mission = MissionProfile(cruise_distance_1, cruise_distance_2, cruise_speed, loiter_time_1, loiter_time_2)

# Set initial guess for W0
W0_initial = 3000

# Call the iteration function
correct_W0, L_D_max_value, AR_wet, G_i, tank_mass, total_fuel_weight_fraction, empty_weight_fraction = find_correct_w0(
    W0_initial, aircraft, mission, tank
)

# Calculate fuel mass and tank volume
fuel_mass = correct_W0 * total_fuel_weight_fraction
tank_volume_total = tank.calc_tank_volume(fuel_mass, cryogenic_h2_density)

# Calculate tank dimensions
D_inner, L_inner, D_outer, L_outer, V_outer_tank = tank.tank_dimensions(fuel_mass)

# Print results
data = [
    ["Takeoff Weight (W0)", f"{correct_W0:.2f} kg"],
    ["Empty Weight Fraction (EWF)", f"{empty_weight_fraction:.2f}"],
    ["Fuel Mass", f"{fuel_mass:.2f} kg"],
    ["Fuel Fraction", f"{total_fuel_weight_fraction:.2f}"],
    ["L/D_max", f"{L_D_max_value:.2f}"],
    ["Aspect Ratio (AR_wet)", f"{AR_wet:.2f}"],
    ["Tank Penalty (G_i)", f"{G_i:.4f}"],
    ["Tank Mass", f"{tank_mass:.2f} kg"],
    ["Tank Volume", f"{tank_volume_total:.2f} m^3"],
    ["Inner Diameter", f"{D_inner:.2f} m"],
    ["Inner Length", f"{L_inner:.2f} m"],
    ["Outer Diameter", f"{D_outer:.2f} m"],
    ["Outer Length", f"{L_outer:.2f} m"],
    ["Outer Tank Volume", f"{V_outer_tank:.2f} m^3"],
]

from tabulate import tabulate
table = tabulate(data, headers=["Parameter", "Value"], tablefmt="grid")
print(table)
