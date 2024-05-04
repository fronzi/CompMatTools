
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
markers = ['o', 's', '^']
colors = ['darkred', 'navy', 'green']
marker_size = 10
line_width = 2
line_styles = ['--', '-.', '-.']
label_fontsize = 22
tick_label_fontsize = 20

# Read data from the CSV file
df = pd.read_csv('SCF_222_111kp.csv')

# Ensure data types are correct, especially for numbers
df['SCF Time'] = pd.to_numeric(df['SCF Time'], errors='coerce')
df['Number of CPUs'] = pd.to_numeric(df['Number of CPUs'], errors='coerce')

# Sorting by the number of CPUs to make the plot easier to understand
df.sort_values('Number of CPUs', inplace=True)

# Assuming the base case is the smallest number of CPUs (using 1 CPU as the proxy for single-core)
base_time = df.loc[df['Number of CPUs'] == 1, 'SCF Time'].values[0]

# Calculate speedup and efficiency based on the provided definitions
df['Speedup'] = base_time / df['SCF Time']
max_speedup = df['Speedup'].max()
df['Normalized Speedup'] = df['Speedup']# / max_speedup
df['Efficiency'] = (df['Speedup'] / df['Number of CPUs']) * 100

# Plotting
fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(10, 12))  # creating two subplots

# Normalized Speedup and Efficiency plot
ax1.set_xlabel('Number of CPUs', fontsize=label_fontsize, fontweight='bold')
ax1.set_ylabel('Speedup', color=colors[0], fontsize=label_fontsize, fontweight='bold')
ax1.plot(df['Number of CPUs'], df['Normalized Speedup'], marker=markers[0], color=colors[0], markersize=marker_size, linewidth=line_width, linestyle=line_styles[0])
ax1.tick_params(axis='y', labelcolor=colors[0], labelsize=tick_label_fontsize)
ax1.set_xscale('log')  # Logarithmic scale for better visualization of wide range data
ax1.grid(True, which="both", linestyle='--', linewidth=0.5)

ax2 = ax1.twinx()
ax2.set_ylabel('Normalized Efficiency (%)', color=colors[1], fontsize=label_fontsize, fontweight='bold')
ax2.plot(df['Number of CPUs'], df['Efficiency'], marker=markers[1], color=colors[1], markersize=marker_size, linewidth=line_width, linestyle=line_styles[1])
ax2.tick_params(axis='y', labelcolor=colors[1], labelsize=tick_label_fontsize)
ax2.grid(True, which="both", linestyle=':', linewidth=0.5)

# Time vs Number of CPUs plot
ax3.set_xlabel('Number of CPUs', fontsize=label_fontsize, fontweight='bold')
ax3.set_ylabel('SCF Time (seconds)', color=colors[2], fontsize=label_fontsize, fontweight='bold')
ax3.plot(df['Number of CPUs'], df['SCF Time'], marker=markers[2], color=colors[2], markersize=marker_size, linewidth=line_width, linestyle=line_styles[2])
ax3.tick_params(axis='y', labelcolor=colors[2], labelsize=tick_label_fontsize)
ax3.set_xscale('log')  # Logarithmic scale for better visualization of wide range data
ax3.set_ylim(0, 200)
ax3.grid(True, which="both", linestyle=':', linewidth=0.5)

plt.title('Performance Metrics at 1 Node Across Various CPUs', fontsize=16, fontweight='bold', pad=20)
plt.tight_layout()  # Adjust layout to make room for plot

plt.savefig("metrics_1cpu_test.png", dpi=500)
plt.show()

