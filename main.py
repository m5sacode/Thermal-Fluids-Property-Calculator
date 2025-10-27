import numpy as np

R = 8.314 # J/Mol*K

class GasPart():
    def __init__(self, name, Mi, A, yi):
        self.name = name
        self.Mi = Mi
        self.A = A
        self.Ri = R/self.Mi
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
        u = 0
        h = 0
        for GasPart in self.GasParts:
            u += GasPart.get_ui(T, GasPart.yi)*GasPart.xi
            h += GasPart.get_hi

        return U, h, Cp, Cv, S