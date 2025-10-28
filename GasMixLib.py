import numpy as np
import matplotlib.pyplot as plt

R = 8.314 # J/Mol*K

class GasPart():
    def __init__(self, name, Mi, A):
        self.R = 8.314 # J/Mol*K
        self.name = name
        self.Mi = Mi
        self.A = A
        self.Ri = self.R/self.Mi
        self.yi = None
    def set_yi(self, yi):
        self.yi = yi
    def get_xi(self, M):
        self.xi = self.yi*self.Mi/M
        return self.xi

    def get_Cpi_Cvi(self, T):
        Tz = T / 10000
        Cpi = sum(Aj * (Tz ** i) for i, Aj in enumerate(self.A))
        Cvi = Cpi - self.Ri/1000
        ki = Cpi / Cvi
        return Cpi, Cvi, ki

    def get_hi(self, T, Tref=298.15):

        A = self.A
        Tz = T / 10000
        Tz_ref = Tref / 10000

        hi = 0.0
        for i, Ai in enumerate(A):
            hi += Ai * 10000.0 / (i + 1) * (Tz ** (i + 1) - Tz_ref ** (i + 1))
        return hi

    def get_s0i(self, T, Tref=298.15):
        A = self.A
        Tz = T / 10000
        Tz_ref = Tref / 10000

        s0i = 0.0
        # First term: A0 * ln(T/Tref)
        s0i += A[0] * np.log(Tz / Tz_ref)

        # Remaining terms
        for j, Aj in enumerate(A[1:], start=1):
            s0i += Aj / j * (Tz ** j - Tz_ref ** j)

        return s0i

    def get_si(self, T, P, Tref=298.15, Pref=101.325e3):

        s0i = self.get_s0i(T, Tref)
        si = s0i - (self.Ri / 1000) * np.log(P / Pref)
        return si

    def get_ui(self, T, Tref=298.15):

        A = self.A
        Tz = T / 10000
        Tz_ref = Tref / 10000

        # Integrate Cp polynomial first
        ui = 0.0
        for i, Ai in enumerate(A):
            # Include *1000 because dT = 1000 * d(Tz)
            ui += Ai * 10000 / (i + 1) * (Tz ** (i + 1) - Tz_ref ** (i + 1))

        # Subtract the integral of Ri (constant) to get Cv integration
        ui -= (self.Ri / 1000) * (T - Tref)

        return ui

    def plot_Cp_vs_T(self, T_min=0, T_max=2000, n_points=100):

        T_range = np.linspace(T_min+273.15, T_max+273.15, n_points)
        Cp_values = []
        Cv_values = []
        k_values = []

        for T in T_range:
            Cpi, Cvi, k = self.get_Cpi_Cvi(T)
            Cp_values.append(Cpi)
            Cv_values.append(Cvi)
            k_values.append(k)

        # Plot Cp vs T
        plt.figure(figsize=(8, 5))
        plt.plot(T_range-273.15, Cp_values, label=f'{self.name}', lw=2, color='darkred')
        plt.xlabel("Temperature [C]", fontsize=12)
        plt.ylabel("Cp [kJ/kg·K]", fontsize=12)
        plt.title(f"Specific Heat Cp vs Temperature for {self.name}", fontsize=14)
        plt.grid(True)
        plt.legend()
        plt.show()


class Gas():
    def __init__(self, name, GasParts):
        self.name = name
        self.GasParts = GasParts
        M = 0
        for GasPart in GasParts:
            M += GasPart.Mi*GasPart.yi
        self.M = M
        for GasPart in GasParts:
            GasPart.get_xi(self.M)
    def get_properties(self, T, P, Yi):
        T = T+273.15
        P = P*1000
        u = 0
        h = 0
        Cp = 0
        Cv = 0
        s = 0
        s0=0
        for GasPart in self.GasParts:
            u += GasPart.get_ui(T)*GasPart.xi
            h += GasPart.get_hi(T)*GasPart.xi
            cpi, cvi, ki = GasPart.get_Cpi_Cvi(T)
            Cp += cpi*GasPart.xi
            Cv += cvi*GasPart.xi
            s += GasPart.get_si(T, P)*GasPart.xi
            s0 += GasPart.get_s0i(T)*GasPart.xi
        k = Cp/Cv
        return u, h, Cp, Cv, s, k, s0

    def plot_Cp_vs_T(self, T_min=0, T_max=2000, n_points=10):
        T_range = np.linspace(T_min, T_max, n_points)
        Cp_values = []

        # Compute Cp for each temperature
        for T in T_range:
            _, _, Cp, _, _, _, _ = self.get_properties(T, 101,
                                                   None)
            Cp_values.append(Cp)

        # Plot Cp vs T
        plt.figure(figsize=(8, 5))
        plt.plot(T_range, Cp_values, color='b', lw=2)
        plt.xlabel("Temperature [C]", fontsize=12)
        plt.ylabel("Cp [kJ/kg·K]", fontsize=12)
        plt.title(f"Specific Heat Cp vs Temperature for {self.name}", fontsize=14)
        plt.grid(True)
        plt.show()

# Provided coeffs
A_O2 = [
    1.006450,
    -1.047869,
    3.729558,
    -4.934172,
    3.284147,
    -1.095203,
    0.145737,
    0.0,
    0.0,
    -0.369790,
    0.000491
]

A_N2 = [
    1.006450,
    -1.047869,
    3.729558,
    -4.934172,
    3.284147,
    -1.095203,
    0.145737,
    0.0,
    0.0,
    -0.369790,
    0.000491
]

A_CO2 = [
    0.408089,
    2.027201,
    -2.405549,
    2.039166,
    -1.163088,
    0.381364,
    -0.052763,
    0.0,
    0.0,
    -0.366740,
    0.001736
]

A_H2O = [
    1.937043,
    -0.967916,
    3.338905,
    -3.652122,
    2.332470,
    -0.819451,
    0.118783,
    0.0,
    0.0,
    2.860773,
    -0.000219
]

A_DryAir = [
    0.992313,
    0.236688,
    -1.852148,
    6.083152,
    -8.893933,
    7.097112,
    -3.234725,
    0.794571,
    -0.081873,
    0.422178,
    0.001053
]

O2 = GasPart("O2", 31.999/1000, np.array(A_O2))
N2 = GasPart("N2", 28.0134/1000, np.array(A_N2))
CO2 = GasPart("CO2", 44.009/1000, np.array(A_CO2))
H2O = GasPart("H2O", 18.01528 /1000, np.array(A_H2O))
DryAir = GasPart("Dry Air", 29.96/1000, np.array(A_DryAir))

DryAir.set_yi(1)
DryAirMix = Gas("Dry Air", [DryAir])
DryAirMix.plot_Cp_vs_T(T_min=0, T_max=2000, n_points=10)
DryAir.plot_Cp_vs_T(T_min=0, T_max=2000, n_points=10)
N2.set_yi(0.79)
O2.set_yi(0.21)
DryAirApprox = Gas("Dry Approx", [N2, O2])
DryAirApprox.plot_Cp_vs_T(T_min=0, T_max=2000, n_points=10)

O2 = GasPart("O2", 31.999 / 1000, np.array(A_O2))
N2 = GasPart("N2", 28.0134 / 1000, np.array(A_N2))
CO2 = GasPart("CO2", 44.009 / 1000, np.array(A_CO2))
H2O = GasPart("H2O", 18.01528 / 1000, np.array(A_H2O))

# --- User input section ---
print("=== Gas Mixture Property Calculator ===")

T = float(input("Enter temperature [°C]: "))
P = float(input("Enter pressure [kPa]: "))

print("\nEnter molal composition (sum should be 1.0):")
y_O2 = float(input(" Mole fraction of O2: "))
y_N2 = float(input(" Mole fraction of N2: "))
y_CO2 = float(input(" Mole fraction of CO2: "))
y_H2O = float(input(" Mole fraction of H2O: "))

# Normalize to ensure total = 1 just in case
total_y = y_O2 + y_N2 + y_CO2 + y_H2O
y_O2 /= total_y
y_N2 /= total_y
y_CO2 /= total_y
y_H2O /= total_y

# --- Assign fractions to species ---
O2.set_yi(y_O2)
N2.set_yi(y_N2)
CO2.set_yi(y_CO2)
H2O.set_yi(y_H2O)

# --- Build the mixture ---
mixture = Gas("Custom Mixture", [O2, N2, CO2, H2O])

# --- Compute properties ---
u, h, Cp, Cv, s, k, s0 = mixture.get_properties(T, P, None)

# --- Display results ---
print("\n=== Thermodynamic Properties at Given State ===")
print(f"Temperature (°C): {T:.2f}")
print(f"Pressure (kPa): {P:.2f}")
print(f"--------------------------------------")
print(f"Cp     = {Cp:.3f} kJ/kg·K")
print(f"Cv     = {Cv:.3f} kJ/kg·K")
print(f"k      = {k:.4f}")
print(f"u      = {u:.3f} kJ/kg")
print(f"h      = {h:.3f} kJ/kg")
print(f"s⁰     = {s0:.3f} kJ/kg·K")
print(f"s      = {s:.3f} kJ/kg·K")
print("=======================================")

