# config.py

# Gravitational constants
g_0 = 9.8  # m/s^2, gravitational acceleration

# Crew and payload
crew_weight = 1 * 90  # 1 crew, 90 kg each
payload_weight = 362.87  # kg

# Aerodynamic data
SPEED_SOUND_1 = 296.5  # m/s
SPEED_SOUND_2 = 309.7  # m/s

# Fuel constants
H2_A1_ratio = 120 / 43
cryogenic_h2_density = 71  # kg/m^3

# Tank parameters
t_shell = 0.1  # m, shell thickness
A_t = 1.5  # Aspect ratio area of tank

# Mission data
SFC_cruise = 5.753 * 1 / 1000000  # Specific fuel consumption, kg/s
SFC_loiter = SFC_cruise * (1 - 0.128)
K_LD = 15.5  # Empirical constant for L/D max

# Aircraft design parameters
MAC = 2  # meter, mean aerodynamic chord
WINGSPAN = 16.5  # meter, wingspan
structual_factor = 0.95

# Wetted area
S = 7  # Wetted area, m^2

# Flight data
cruise_distance_1 = 1200 * 1852  # 600 NM converted to meters (Cruise 1)
cruise_distance_2 = 200 * 1852  # 200 NM converted to meters (Cruise 2, diversion)
cruise_speed = 360 * 0.514444  # 360 knots converted to m/s
loiter_time_1 = 30  # 30 minutes loiter (first phase)
loiter_time_2 = 30  # 45 minutes loiter (second phase)

# Empty weight constants
a_ewf = 1.5  # constant for empty weight calculation
b_ewf = -0.096  # exponent for empty weight calculation

# Tank-penalty constants
a_tank = 16.53
b_tank = -4.818
c_tank = 0.4216
d_tank = -0.1132

# Fuel fractions for different phases
F_TAXI = 1 - (1 - 0.97) / H2_A1_ratio  # fuel fraction for taxi and take-off
F_DECENT = 1 - (1 - 0.952) / H2_A1_ratio  # fuel fraction for descent
F_LANDING = 1 - (1 - 0.995) / H2_A1_ratio  # fuel fraction for landing

# Fuselage dimensions
D_fuselage = 2  # Fuselage diameter in meters
L_fuselage = 16  # Total fuselage length in meters
L_cabin = 4  # Assuming 80% of the fuselage length is used for cabin space
# Top-level requirements
required_cabin_volume_ft3 = 220  # Required cabin volume in cubic feet
required_cabin_volume_m3 = required_cabin_volume_ft3 / 35.3147  # Convert to cubic meters
