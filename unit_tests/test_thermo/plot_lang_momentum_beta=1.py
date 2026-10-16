"""Plot histogram of Langevin momentum samples from lang_momentum.dat."""

import math
import os

data_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lang_momentum_beta=1.dat")
data = [float(line) for line in open(data_file)]
N = len(data)

mass = 1.0
beta = 1.0
sigma = math.sqrt(mass / beta)   # = 1.0

mean = sum(data) / N
var  = sum(x*x for x in data) / N - mean**2
print(f"N       = {N}")
print(f"mean    = {mean:.6f}  (exact 0)")
print(f"var     = {var:.6f}  (exact {mass/beta:.6f})")
print(f"std dev = {math.sqrt(var):.6f}  (exact {sigma:.6f})")


import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(data, bins=80, density=True, color="steelblue", alpha=0.6, label="Langevin samples")

p_grid = np.linspace(-4, 4, 400)
pdf = np.exp(-0.5 * (p_grid / sigma)**2) / (sigma * math.sqrt(2 * math.pi))
ax.plot(p_grid, pdf, color="tomato", linewidth=2,
        label=r"$\mathcal{N}(0,\,m/\beta)$  exact")

ax.set_xlabel("momentum $p$")
ax.set_ylabel("density")
ax.set_title(f"Langevin 1-D HO canonical momentum  ($N={N:,}$,  $m=\\beta=\\omega=1$)")
ax.legend()
fig.tight_layout()

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lang_momentum_beta=1.png")
fig.savefig(out, dpi=150)
