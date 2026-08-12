"""
Plot distributions from test_nhc_quartic NVT run.
Columns in nhc_quartic_T=300.dat:
    1: p       (momentum)
    2: x^2     (position squared)
    3: Ekin    (kinetic energy  = p^2/(2m))
    4: V       (potential energy = lambda*x^4)
Run from unit_tests/:  python plot_nhc_quartic.py
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.special import gamma as gamma_func

# ── physical parameters ──────────────────────────────────────────────────────
BoltzK = 3.16681520371153e-6   # Hartree / K
T      = 300.0                 # K
beta   = 1.0 / (BoltzK * T)
mass   = 1.0                   # au
lam    = 0.25                  # V = lam * x^4

# ── exact canonical targets ──────────────────────────────────────────────────
target_p2   = mass / beta      # <p^2>  = m/beta
target_Ekin = 0.5  / beta      # <Ekin> = kT/2
target_Epot = 0.25 / beta      # <V>    = kT/4   (virial theorem)

# ── load data ─────────────────────────────────────────────────────────────────
data = np.loadtxt("nhc_quartic_T=300.dat")
p    = data[:, 0]
x    = data[:, 1]
Ekin = data[:, 2]
V    = data[:, 3]

# ── figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(
    rf"NHC NVT — quartic oscillator $V(x)=\lambda x^4$,  $T={T:.0f}$ K"
    f"\n(mass={mass} au,  λ={lam},  β={beta:.3e} au⁻¹)",
    fontsize=12,
)

# ── 1. Momentum histogram vs Maxwell-Boltzmann ────────────────────────────────
ax = axes[0, 0]
sigma_p = np.sqrt(target_p2)
p_grid  = np.linspace(-4.5 * sigma_p, 4.5 * sigma_p, 400)
ax.hist(p, bins=80, density=True, alpha=0.7, color="steelblue", label="Sampled")
ax.plot(p_grid, norm.pdf(p_grid, 0.0, sigma_p), "r-", lw=2,
        label=r"$\mathcal{N}(0,\,\sqrt{m/\beta})$")
ax.set_xlabel("p  (au)", fontsize=11)
ax.set_ylabel("Density", fontsize=11)
ax.set_title(
    rf"Momentum:  $\langle p\rangle={np.mean(p):.2e}$"
    rf",  $\langle p^2\rangle={np.mean(p**2):.3e}$  (target {target_p2:.3e})"
)
ax.legend()

# ── 2. x histogram with exact theoretical P(x) ──────────────────────────────
#   P(x) = sqrt(2)*beta^(1/4)/Gamma(1/4) * exp(-beta*x^4/4)   [from image]
#   Column 2 stores signed x, so apply P(x) directly on symmetric range.
ax = axes[0, 1]
x_lo  = np.percentile(x, 0.5)
x_hi  = np.percentile(x, 99.5)
x_arr = np.linspace(x_lo, x_hi, 800)
px    = np.sqrt(2) * beta**0.25 / gamma_func(0.25) * np.exp(-beta * x_arr**4 / 4.0)
px   /= np.trapz(px, x_arr)   # renormalise to plotted range
ax.hist(x, bins=100, range=(x_lo, x_hi), density=True,
         color="steelblue", alpha=0.8, label="Sampled")
ax.plot(x_arr, px, "r-", lw=2,
        label=r"$\frac{\sqrt{2}\,\beta^{1/4}}{\Gamma(1/4)}\,e^{-\beta x^4/4}$")
ax.axvline(np.mean(x), color="k", lw=1.5, ls="--",
           label=rf"$\langle x\rangle={np.mean(x):.3e}$")
ax.set_xlabel(r"$x$  (au)", fontsize=12)
ax.set_ylabel("Density", fontsize=12)
ax.set_title(rf"$x$ distribution — NHC quartic oscillator, $T={T:.0f}$ K")
ax.legend(fontsize=9)


# ── 3. Running mean energies ──────────────────────────────────────────────────
ax = axes[1, 0]
N     = len(Ekin)
steps = np.arange(1, N + 1)
Ekin_run = np.cumsum(Ekin) / steps
Epot_run = np.cumsum(V)    / steps
ax.plot(steps, Ekin_run, color="steelblue", lw=1.0,
        label=r"$\langle E_{\rm kin}\rangle$ running avg")
ax.plot(steps, Epot_run, color="orange",    lw=1.0,
        label=r"$\langle V\rangle$ running avg")
ax.axhline(target_Ekin, color="steelblue", ls="--", lw=1.5,
           label=rf"$kT/2={target_Ekin:.3e}$")
ax.axhline(target_Epot, color="orange",    ls="--", lw=1.5,
           label=rf"$kT/4={target_Epot:.3e}$")
ax.set_xlabel("Sample index", fontsize=11)
ax.set_ylabel("Energy  (au)", fontsize=11)
ax.set_title("Running mean energies (convergence to equilibrium)")
ax.legend(fontsize=9)


# ── standalone x² histogram ───────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.hist(x, bins=100, range=(x_lo, x_hi), density=True,
         color="steelblue", alpha=0.8, label="Sampled")
ax2.plot(x_arr, px, "r-", lw=2,
         label=r"$\frac{\sqrt{2}\,\beta^{1/4}}{\Gamma(1/4)}\,e^{-\beta x^4/4}$")
ax2.axvline(np.mean(x), color="k", lw=2, ls="--",
            label=rf"$\langle x\rangle={np.mean(x):.4e}$")
ax2.set_xlabel(r"$x$  (au)", fontsize=12)
ax2.set_ylabel("Density", fontsize=12)
ax2.set_title(rf"$x$ distribution — NHC quartic oscillator, $T={T:.0f}$ K")
ax2.legend(fontsize=11)
plt.tight_layout()
out_x2 = "nhc_quartic_x_hist.png"
plt.savefig(out_x2, dpi=150, bbox_inches="tight")
print(f"Saved: {out_x2}")
