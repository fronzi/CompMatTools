import matplotlib.pyplot as plt
import pandas as pd

# Load data
df = pd.read_csv('SCF_222_111kp.csv')

# Convert columns to numeric, ensuring all data is correctly typed
df['SCF Time'] = pd.to_numeric(df['SCF Time'], errors='coerce')
df['Number of CPUs'] = pd.to_numeric(df['Number of CPUs'], errors='coerce')
df['Number of Nodes'] = pd.to_numeric(df['Number of Nodes'], errors='coerce')

# Sort the dataframe by the number of CPUs and Nodes to make the plot easier to understand
df.sort_values(['Number of Nodes', 'Number of CPUs'], inplace=True)

# Determine min_time for the baseline (SCF Time on a single node normalized by the number of CPUs)
base_time = df[df['Number of Nodes'] == 4]['SCF Time'].max()# / df[df['Number of Nodes'] == 1]['Number of CPUs'].iloc[0]


# Define colors for different node numbers to visually differentiate them
node_colors = {1: 'blue', 2: 'green', 4: 'red'}

# Marker setup for different metrics
metrics_markers = {'Speedup': 'o', 'Efficiency': 's', 'SCF Time': '^'}


marker_size=10



# Plotting setup
fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(10, 12))

# Group the data by Number of Nodes and plot
for node, group in df.groupby('Number of Nodes'):
    group = group.dropna(subset=['SCF Time'])
    if not group.empty:
        # min_time = group['SCF Time'].min()  # Tseq is the minimum SCF Time

        # Calculate Speedup and Efficiency
        group['Speedup'] = base_time / group['SCF Time']
        group['Efficiency'] = (base_time / (group['Number of CPUs'] * group['SCF Time'])) * 100


        # Normalize to maximum values
        max_speedup = group['Speedup'].max()
        max_efficiency = group['Efficiency'].max()
        group['Normalized Speedup'] = group['Speedup']# / max_speedup
        group['Normalized Efficiency'] = (group['Efficiency'] / max_efficiency) * 100

        # Plotting each metric with specific marker
        ax1.plot(group['Number of CPUs'], group['Normalized Speedup'], marker=metrics_markers['Speedup'], linestyle='--', linewidth=2, color=node_colors[node], label=f'{node} Nodes',markersize=marker_size)
        ax1.set_xlabel('Number of CPUs', fontsize=20, fontweight='bold')
        ax1.set_ylabel('Speedup', color='darkred', fontsize=22, fontweight='bold')
        ax1.tick_params(axis='y', labelcolor='darkred', labelsize=20)
        ax1.set_xscale('log')
        ax1.set_xlim(1, 140)
        # ax1.set_ylim(0, 1.05)
        ax1.grid(True, which="both", linestyle='--', linewidth=0.5)

        # Efficiency on twin axis
        ax2 = ax1.twinx()
        ax2.plot(group['Number of CPUs'], group['Normalized Efficiency'], marker=metrics_markers['Efficiency'], linestyle='-.', linewidth=2, color=node_colors[node], label=f'{node} Nodes',markersize=marker_size)
        ax2.set_ylabel('Normalized Efficiency (%)', color='navy', fontsize=22, fontweight='bold')
        ax2.set_xscale('log')
        ax2.tick_params(axis='y', labelcolor='navy', labelsize=20)
        ax2.set_xlim(1, 140)
        ax2.set_ylim(0, 105)
        ax2.grid(False)

        # SCF Time on a separate axis
        ax3.plot(group['Number of CPUs'], group['SCF Time'], marker=metrics_markers['SCF Time'], linestyle='-.', linewidth=2, color=node_colors[node], label=f'{node} Nodes',markersize=marker_size)
        ax3.set_xlabel('Number of CPUs', fontsize=22, fontweight='bold')
        ax3.set_ylabel('SCF Time (seconds)', color='green', fontsize=20, fontweight='bold')
        ax3.set_xscale('log')
        ax3.set_xlim(1, 140)
        ax3.set_ylim(0, max(group['SCF Time']) * 1.1)
        ax3.grid(True, which="both", linestyle='--', linewidth=0.5)
        ax3.set_xscale('log')
        ax3.set_xlim(0, 140)
        ax3.set_ylim(0, 70)
        ax3.tick_params(axis='y', labelcolor='green', labelsize=20)
        ax3.grid(True, which="both", linestyle='--', linewidth=0.5)

# Adding legends
ax1.legend(title='Speedup', loc='upper left')
ax2.legend(title='Norm. Efficiency', loc='lower left')
ax3.legend(title='SCF Time', loc='upper left')
plt.title('Performance Metrics Across Various Node Configurations', fontsize=16, fontweight='bold', pad=20)

plt.tight_layout()  # Adjust layout to make room for plot
# plt.suptitle('Performance Metrics Across Various Node Configurations', fontsize=16, fontweight='bold')
# fig.tight_layout()  # Adjust the subplots to fit the suptitle
# fig.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust the subplots to fit the suptitle

plt.savefig("metrics_Ncpu.png", dpi=300)

plt.show()
