import numpy as np
from fontTools.varLib.models import nonNone

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
        Tz = T / 1000
        Cpi = 0
        A = self.A
        for Aj, i in enumerate(A):
            Cpi += Aj * (Tz ** i)
        Cvi = Cpi - self.Ri
        return Cpi, Cvi
    def get_hi(self, T, Tref=0.5):
        A = self.A
        Mi = self.Mi  # [kg/mol]

        hi = 0.0
        for i, Ai in enumerate(A):
            term = Ai / (1000 ** i * (i + 1)) * (T ** (i + 1) - Tref ** (i + 1))
            hi += term

        # Convert from molar (J/mol) to specific (J/kg)
        hi /= Mi

        return hi
    def get_s0i(self, T, Tref=0.5):
        A = self.A
        Mi = self.Mi  # [kg/mol]

        s0i = 0.0

        # j = 0 term: integral(A0 / T) dT = A0 * ln(T/Tref)
        s0i += A[0] * np.log(T / Tref)

        # j >= 1 terms: ∫ (A_j * T^(j-1) / 1000^j) dT = A_j / (j * 1000^j) * (T^j - Tref^j)
        for j, Aj in enumerate(A[1:], start=1):
            s0i += Aj / (j * 1000 ** j) * (T ** j - Tref ** j)

        # Convert from molar (J/mol·K) to specific (J/kg·K)
        s0i /= Mi

        return s0i
    def get_si(self, T, P, Tref=0.5, Pref=101000):
        s0i = self.get_s0i(T, Tref)
        si = s0i-R*np.log(P/Pref)
        return si
    def get_ui(self, T, Tref=0.5):
        A = self.A
        Mi = self.Mi  # [kg/mol]
        Ri = self.Ri  # [J/kg·K]

        ui = 0.0
        # Integrate Cp(T) - Ri
        for i, Ai in enumerate(A):
            term = Ai / (1000 ** i * (i + 1)) * (T ** (i + 1) - Tref ** (i + 1))
            ui += term

        # Subtract the Ri * (T - Tref) term (integral of constant Ri)
        ui -= Ri * Mi * (T - Tref)  # temporarily in molar basis

        # Convert from molar (J/mol) to specific (J/kg)
        ui /= Mi

        return ui

class Gas():
    def __init__(self, name, GasParts):
        self.name = name
        self.GasParts = GasParts
        M = 0
        for GasPart in GasParts:
            M += GasPart.Mi()*GasPart.yi
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
        for GasPart in self.GasParts:
            u += GasPart.get_ui(T, GasPart.yi)*GasPart.xi
            h += GasPart.get_hi
            cpi, cvi = GasPart.get_Cpi_Cvi(T)
            Cp += cpi
            Cv += cvi
            s += GasPart.get_si(T, GasPart.yi)
        k = Cp/Cv
        return u, h, Cp, Cv, s, k

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
    0.0,
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
DryAir = GasPart("Air", 29.96/1000, np.array(A_DryAir))