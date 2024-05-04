#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 16:57:02 2024

@author: marco




"""
import matplotlib.pyplot as plt
import numpy as np

# Data
omp_levels = ['4-OMP', '8-OMP', '16-OMP']
scf_times = [32.129/4, 29.437/4, 29.643/4]
ref_time = 29.612/4  # SCF time for 1-OMP configuration

# Calculate normalized SCF times
normalized_times = np.array([time / ref_time for time in scf_times])

# Colors
colors = ['blue', 'darkred', 'green']

# Plot
plt.figure(figsize=(10, 8))

# Plot 1: Bar Plot
plt.subplot(1, 1, 1)
plt.bar(omp_levels, scf_times, color=colors)
plt.xlabel('OMP Parallelism', fontsize=12, fontweight='bold')
plt.ylabel('SCF Time (seconds)', fontsize=12, fontweight='bold')
plt.title('Impact of OMP Parallelism on Calculation Time', fontsize=14, fontweight='bold')
plt.xticks(fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')
plt.grid(axis='y', linestyle='--', alpha=0.7)


plt.axhline(y=1, color='gray', linestyle='--', linewidth=1)

# Add annotations to line plot
for i, txt in enumerate(normalized_times):
    plt.annotate(f'{txt:.3f}', (omp_levels[i], normalized_times[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=10, fontweight='bold')

# Adjust layout
plt.tight_layout()
plt.savefig("OMP.png", dpi=300)
# Show plot
plt.show()
