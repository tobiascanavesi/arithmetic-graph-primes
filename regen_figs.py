#!/usr/bin/env python3
"""Regenerate correlation-decay figure with the extended (exact-BC) data."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).parent
with open(here / "bc_results.json") as f:
    data = json.load(f)

Ns = [r["N"] for r in data]
rhos = [r["rho"] for r in data]

fig, ax = plt.subplots(figsize=(6.5, 4.2))
ax.axhline(0, color="grey", lw=0.7, ls="--", alpha=0.7)
ax.plot(Ns, rhos, "ko-", markersize=7, lw=1.5)
for n, r in zip(Ns, rhos):
    sign = "+" if r >= 0 else ""
    offset = (0, 12) if r >= 0 else (0, -16)
    ax.annotate(f"{sign}{r:.3f}", (n, r), textcoords="offset points",
                xytext=offset, ha="center", fontsize=8)

ax.set_xscale("log")
ax.set_xlabel(r"$N$")
ax.set_ylabel(r"$\rho(\mathrm{deg}, \mathrm{BC})$")
ax.set_title("Degree--Betweenness Correlation in $G_N$ (exact BC)")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(here / "fig_correlation_decay.pdf", dpi=300)
plt.close(fig)
print(f"Saved fig_correlation_decay.pdf ({len(data)} scales, rho range {min(rhos):.3f} to {max(rhos):.3f})")
