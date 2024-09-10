import math

'''Empty Weight Fraction'''
W_TANK = 800  # H2 tank weight (kg)
A = 1.5  # Curve fit parameter
B = -0.096  # Curve fit parameter

'''Segment Weight Fractions'''
F_TAXI = 0.985
F_DECENT = 0.952
F_LANDING = 0.99

SPEED_SOUND_1 = 296
SPEED_SOUND_2 = 296

S_WET_REF = 6.2
SFC = 6.5 / 1000000  # Kg/s
K_LD = 15.5
A_LD = 5.5

class Aircraft:
    def __init__(self):
        self.seats_cabin = 6
        self.seats_cockpit = 2
        self.weight_per_pax = 90  # kg
        self.cabin_volume_m3 = 6.23  # Minimum cabin volume in cubic meters
        self.payload_kg = 362.87  # 800 lbs converted to kg
        self.range_m = 2222400  # 1200 nautical miles converted to meters
        self.reserves_divert_m = 370400  # 200 nm converted to meters
        self.reserves_hold_min = 30  # Holding time in minutes
        self.reserve_fuel_percent = 0.05  # Reserve fuel 5%
        self.trapped_fuel_percent = 0.01  # Trapped fuel 1%
        self.balanced_field_length_m = 1219.2  # Balanced field length in meters
        self.lfl_m = 914.4  # Landing field length in meters
        self.max_cruise_speed_m_per_s = 216.07  # Max cruise speed in meters/second
        self.cruise_speed_m_per_s = 185.2  # Cruise speed in meters/second
        self.min_cruise_altitude_m = 10668  # Minimum cruise altitude in meters
        self.initial_roc_m_per_s = 20.32  # Initial rate of climb in meters/second
        self.roc_height_m = 121.92  # Rate of climb height in meters
        self.technology_freeze_year = 2035
        self.g_0 = 9.8  # Gravity in m/s²

        self.segment_weight_fractions = {
            "Taxi": F_TAXI,
            "Climb_1": None,
            "Cruise_1": None,
            "Decent_1": F_DECENT,
            "Loiter_1": None,
            "Landing_1": F_LANDING,
            "Climb_2": None,
            "Cruise_2": None,
            "Decent_2": F_DECENT,
            "Loiter_2": None,
            "Landing_2": F_LANDING,
            "Reserve": 1 - self.reserve_fuel_percent,  # Reserve fuel fraction
            "Trapped": 1 - self.trapped_fuel_percent, # Trapped fuel fraction
        }

    def empty_weight_fraction(self, w0, a, b, w_tank):
        """Calculate the empty weight fraction using the formula: EWF = a * W0^b + Wtank / W0."""
        ewf = a * (w0 ** b) + (w_tank / w0)
        return ewf
    
    def weight_fraction_climb(self, speed, a):
        M = speed / a 
        wfc = 1.0065 - 0.0325* M
        return wfc

    def weight_fraction_cruise(self, SFC, L_D, Distance, Speed):
        R = Distance
        g = self.g_0
        V = Speed
        wfc = math.exp((-R * SFC * g) / (V * (L_D * 0.886)))
        return wfc

    def weight_fraction_loiter(self, SFC, L_D , Time, Speed):
        E = Time * 60  # Minutes to seconds
        g = self.g_0
        V = Speed
        wfc = math.exp((-E * SFC * g) / (V * L_D))
        return wfc

    def calculate_L_D_max(self, S_WET_REF, A, K_LD):
        L_D_MAX = K_LD * math.sqrt(A / S_WET_REF)
        return L_D_MAX

    def mission_segment_weight_fractions(self):
        L_D = self.calculate_L_D_max(S_WET_REF, A_LD, K_LD)

        self.segment_weight_fractions["Climb_1"] = self.weight_fraction_climb(self.cruise_speed_m_per_s, SPEED_SOUND_1)
        self.segment_weight_fractions["Cruise_1"] = self.weight_fraction_cruise(SFC, L_D, self.range_m, self.cruise_speed_m_per_s)
        self.segment_weight_fractions["Loiter_1"] = self.weight_fraction_loiter(SFC, L_D, self.reserves_hold_min, self.cruise_speed_m_per_s*0.8)
        self.segment_weight_fractions["Climb_2"] = self.weight_fraction_climb(self.cruise_speed_m_per_s, SPEED_SOUND_2)
        self.segment_weight_fractions["Cruise_2"] = self.weight_fraction_cruise(SFC, L_D, self.reserves_divert_m, self.cruise_speed_m_per_s)
        self.segment_weight_fractions["Loiter_2"] = self.weight_fraction_loiter(SFC, L_D, self.reserves_hold_min, self.cruise_speed_m_per_s*0.8)

        mfw = 1
        # Iterate through all the values in the dictionary and multiply them
        for value in self.segment_weight_fractions.values():
            mfw *= value

        return mfw

def main():
    # Instantiate the Aircraft class
    aircraft = Aircraft()

    # Initial guess for W0
    W0 = 100000
    tolerance = 0.01  # Define a tolerance level for iteration convergence
    max_iterations = 1000  # Maximum number of iterations
    iteration = 0

    while iteration < max_iterations:
        iteration += 1
        
        # Calculate the empty weight fraction
        ewf = aircraft.empty_weight_fraction(W0, A, B, W_TANK)

        # Calculate the mission weight fraction (including all mission segments and reserves)
        mwf = aircraft.mission_segment_weight_fractions()

        # Update W0 based on the calculated fractions
        new_W0 = (aircraft.payload_kg + aircraft.payload_kg / (1 - mwf - ewf))

        # Check if the change in W0 is within the tolerance
        if abs(new_W0 - W0) < tolerance:
            break

        W0 = new_W0

    print(f"\nConverged W0: {W0:.2f}")
    print(f"Empty Weight Fraction: {ewf:.4f}")
    print(f"Mission Weight Fraction: {mwf:.4f}")
if __name__ == "__main__":
    main()
