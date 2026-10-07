"""
Infinite-medium k-infinity and critical boron — a runnable analytic cross-check
===============================================================================

A transparent, fast estimate that needs no nuclear-data library. It evaluates a
four-factor infinite medium with the pin-cell regions homogenised by volume:

        k_inf = eta * f * p * epsilon

        eta = nu * sigma_f(U235) / sigma_a(fuel)
        f   = Sigma_a(fuel) / [ Sigma_a(fuel) + Sigma_a(water) + Sigma_a(boron)
                                + Sigma_a(clad) ]            (volume-weighted)
        p, epsilon  : resonance-escape and fast-fission factors (typical values)

Sweeping the soluble-boron concentration shows reactivity falling with boron
and the **critical boron concentration** where k_inf crosses 1 — the same
behaviour the full OpenMC model reproduces. This is a teaching-level estimate
(constant p, epsilon; thermal one-group), not a substitute for the Monte Carlo
transport calculation in pin_cell.py.

Author: Wendlasida Rodrigue Kientega
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NA = 6.02214076e23
NAVY, GOLD, EMER = "#1c4571", "#b0822f", "#1f6b53"

# thermal (0.0253 eV) microscopic cross-sections, barns
sig = dict(U235_f=585.0, U235_a=681.0, U238_a=2.68,
           H_a=0.332, B_a=760.0, Zr_a=0.18)
NU = 2.43           # neutrons per thermal fission of U-235
P_RES = 0.78        # resonance-escape probability (typical PWR value)
EPS = 1.03          # fast-fission factor (typical)
BARN = 1e-24        # cm^2

# pin-cell geometry -> volume fractions (square pitch)
import numpy as _np
R_FUEL, R_CLAD, PITCH = 0.4096, 0.4750, 1.26
_A_cell = PITCH ** 2
VF_FUEL = _np.pi * R_FUEL ** 2 / _A_cell
VF_CLAD = _np.pi * (R_CLAD ** 2 - 0.4178 ** 2) / _A_cell
VF_WATER = 1.0 - _np.pi * R_CLAD ** 2 / _A_cell


def number_densities(enrich_wt=4.05, rho_uo2=10.5, rho_water=0.72, rho_zr=6.55):
    """Atomic densities (atoms/cm^3) within each material region."""
    M_UO2, M_Zr = 270.0, 91.22
    n_uo2 = rho_uo2 / M_UO2 * NA
    w = enrich_wt / 100.0
    a235 = (w / 235.0) / (w / 235.0 + (1 - w) / 238.0)   # atom fraction U-235
    N = dict(U235=a235 * n_uo2, U238=(1 - a235) * n_uo2)
    n_h2o = rho_water / 18.0 * NA
    N["H"] = 2 * n_h2o
    N["Zr"] = rho_zr / M_Zr * NA
    N["rho_water"] = rho_water
    return N


def k_inf(boron_ppm, N):
    """Four-factor k_inf of the volume-homogenised pin cell."""
    # macroscopic absorption/fission, volume-weighted (cm^-1)
    Sig_f = N["U235"] * sig["U235_f"] * BARN * VF_FUEL
    Sig_a_fuel = (N["U235"] * sig["U235_a"] + N["U238"] * sig["U238_a"]) * BARN * VF_FUEL
    Sig_a_water = N["H"] * sig["H_a"] * BARN * VF_WATER
    Sig_a_clad = N["Zr"] * sig["Zr_a"] * BARN * VF_CLAD
    n_B = (boron_ppm * 1e-6 * N["rho_water"]) / 10.81 * NA
    Sig_a_B = n_B * sig["B_a"] * BARN * VF_WATER

    eta = NU * Sig_f / Sig_a_fuel
    f = Sig_a_fuel / (Sig_a_fuel + Sig_a_water + Sig_a_clad + Sig_a_B)
    return eta * f * P_RES * EPS


def boron_worth(N, ppm_ref=1235.0):
    """Differential boron worth in pcm/ppm around a reference concentration."""
    def rho(b):
        k = k_inf(b, N)
        return (k - 1.0) / k
    dppm = 50.0
    dr = rho(ppm_ref + dppm) - rho(ppm_ref - dppm)
    return dr / (2 * dppm) * 1e5     # pcm per ppm


if __name__ == "__main__":
    N = number_densities()
    ppm = np.linspace(0, 2500, 400)
    k = np.array([k_inf(b, N) for b in ppm])

    print(f"Volume fractions:  fuel {VF_FUEL:.3f}  clad {VF_CLAD:.3f}  water {VF_WATER:.3f}")
    print(f"k_inf at    0 ppm boron = {k_inf(0, N):.4f}")
    print(f"k_inf at 1235 ppm boron = {k_inf(1235, N):.4f}")
    print(f"differential boron worth ~ {boron_worth(N):.1f} pcm/ppm "
          f"(textbook PWR value ~ -8 to -12 pcm/ppm)")

    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    ax.plot(ppm, k, color=NAVY, lw=1.9)
    ax.axvline(1235, color=GOLD, ls=":", lw=1.5, label="1235 ppm (model point)")
    ax.plot([1235], [k_inf(1235, N)], "o", color=EMER, ms=7)
    ax.set_xlabel("Soluble boron (ppm)"); ax.set_ylabel(r"k$_\infty$")
    ax.set_title("Infinite-medium reactivity vs soluble boron\n"
                 "(simplified four-factor model)", color=NAVY, fontweight="bold")
    ax.grid(alpha=.25); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("figures/kinf_vs_boron.png", dpi=130)
    print("saved figures/kinf_vs_boron.png")
