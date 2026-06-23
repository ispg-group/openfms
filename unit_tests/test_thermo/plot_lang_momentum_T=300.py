"""Plot Langevin OBABO canonical phase-space distribution for a 1-D harmonic oscillator.

Data file: lang_momentum_T=300.dat  (two columns: p  x)
Exact distributions:
    p ~ N(0, sqrt(m/beta))
    x ~ N(0, sqrt(1/(m*omega^2*beta)))
"""

import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── parameters ────────────────────────────────────────────────────────────────
mass   = 1.0
omega2 = 1.0
BoltzK = 3.16681520371153e-6   # Hartree / K
T      = 300.0
beta   = 1.0 / (BoltzK * T)

sigma_p = math.sqrt(mass / beta)
sigma_x = math.sqrt(1.0 / (mass * omega2 * beta))

# ── load data ─────────────────────────────────────────────────────────────────
data_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "lang_momentum_T=300.dat")
data = np.loadtxt(data_file)
p_data = data[:, 0]
x_data = data[:, 1]
N = len(p_data)

# ── print statistics ──────────────────────────────────────────────────────────
for label, d, sigma_exact in [("momentum p", p_data, sigma_p),
                               ("position x", x_data, sigma_x)]:
    mean = d.mean()
    var  = d.var()
    print(f"\n{label}:")
    print(f"  N       = {N}")
    print(f"  mean    = {mean:.6e}  (exact 0)")
    print(f"  var     = {var:.6e}  (exact {sigma_exact**2:.6e})")
    print(f"  std dev = {math.sqrt(var):.6e}  (exact {sigma_exact:.6e})")

# ── plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle(
    f"Langevin OBABO — 1-D HO canonical distribution  "
    f"($N={N:,}$,  $m=1$,  $\\omega=1$,  $T={T:.0f}$ K)",
    fontsize=12)

panels = [
    (axes[0], p_data, sigma_p,
     "momentum $p$",
     r"$\mathcal{N}(0,\,m/\beta)$  exact",
     "steelblue"),
    (axes[1], x_data, sigma_x,
     "position $x$",
     r"$\mathcal{N}(0,\,1/(m\omega^2\beta))$  exact",
     "seagreen"),
]

for ax, d, sigma, xlabel, exact_label, color in panels:
    lo, hi = -4.5 * sigma, 4.5 * sigma
    ax.hist(d, bins=80, density=True,
            color=color, alpha=0.55, label="samples")
    grid = np.linspace(lo, hi, 500)
    pdf  = np.exp(-0.5 * (grid / sigma)**2) / (sigma * math.sqrt(2 * math.pi))
    ax.plot(grid, pdf, color="tomato", linewidth=2, label=exact_label)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("density", fontsize=11)
    ax.legend(fontsize=9)

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "lang_momentum_T=300.png")
fig.savefig(out, dpi=150)
print(f"\nSaved {out}")
