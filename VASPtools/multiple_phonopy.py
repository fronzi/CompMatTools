import matplotlib
matplotlib.use("Agg")   # HPC-safe
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------------
# Plot style (clean, Phonopy-like)
# -----------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.linewidth": 1.1,
    "figure.dpi": 120,
})

# -----------------------------
# Load Phonopy thermal data
# -----------------------------
def load_phonopy_thermal(filename):
    with open(filename, "r") as f:
        data = yaml.safe_load(f)
    tp = data["thermal_properties"]
    T  = np.array([x["temperature"] for x in tp])
    F  = np.array([x["free_energy"]  for x in tp])   # kJ/mol
    S  = np.array([x["entropy"]      for x in tp])   # J/K/mol
    Cv = np.array([x["heat_capacity"]for x in tp])   # J/K/mol
    E  = np.array([x["energy"]       for x in tp])   # kJ/mol
    return T, F, S, Cv, E

# -----------------------------
# Phases to compare
# -----------------------------
phase_dirs = [
    "bulk-alpha/phonons333/",
    "bulk-alpha-phenol/phonons333/",
    "bulk-alpha-phenol-2x1x1/phonons222/",   # ← missing comma fixed
    "bulk-delta/phonons222/",
    "bulk-delta-phenol/phonons222/",
]

labels = [
    "Alpha",
    "Alpha PEAI",
    "Alpha PEAI-211",
    "Delta",
    "Delta PEAI",
]

colors = [
    "#1f77b4",   # blue
    "#d62728",   # red
    "#2ca02c",   # green
    "#9467bd",   # purple
    "#ff7f0e",   # orange  ← 5th color added
]

linestyles = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]  # ← 5th linestyle added

# -----------------------------
# Figure (Phonopy layout)
# -----------------------------
fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
(axF, axS), (axCv, axE) = axes

for d, label, c, ls in zip(phase_dirs, labels, colors, linestyles):
    file = Path(d) / "thermal_properties.yaml"
    if not file.exists():
        print(f"WARNING: {file} not found — skipping {label}")
        continue
    try:
        T, F, S, Cv, E = load_phonopy_thermal(file)
    except Exception as e:
        print(f"WARNING: Could not load {file}: {e} — skipping {label}")
        continue
    axF.plot(T, F,  color=c, linestyle=ls, linewidth=2.0, label=label)
    axS.plot(T, S,  color=c, linestyle=ls, linewidth=2.0)
    axCv.plot(T, Cv, color=c, linestyle=ls, linewidth=2.0)
    axE.plot(T, E,  color=c, linestyle=ls, linewidth=2.0)

# -----------------------------
# Labels
# -----------------------------
axF.set_ylabel("Free energy (kJ/mol)")
axS.set_ylabel("Entropy (J/K/mol)")
axCv.set_ylabel("Heat capacity (J/K/mol)")
axE.set_ylabel("Energy (kJ/mol)")
axCv.set_xlabel("Temperature (K)")
axE.set_xlabel("Temperature (K)")

# -----------------------------
# Aesthetics
# -----------------------------
for ax in axes.flatten():
    ax.grid(True, linestyle=":", linewidth=0.6, alpha=0.7)
    ax.set_xlim(left=0)

axF.legend(loc="best", fontsize=10)
plt.tight_layout()
plt.subplots_adjust(hspace=0.18, wspace=0.25)

# -----------------------------
# Save
# -----------------------------
outname = "phonopy_style_thermal_properties_5phases"
plt.savefig(f"{outname}.pdf", dpi=300, bbox_inches="tight")
plt.savefig(f"{outname}.png", dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved {outname}.pdf and {outname}.png")
