import pandas as pd
import matplotlib.pyplot as plt

# --- 1. Load the Data ---
db_df = pd.read_csv('figures/database-size-timing.tsv', sep='\t')
th_df = pd.read_csv('figures/thread-count-timing.tsv', sep='\t')

# Convert Peak Memory from KB to GB
db_df['MEM_GB'] = db_df['MEM'] / 1048576.0
th_df['MEM_GB'] = th_df['MEM'] / 1048576.0

# --- 2. Setup the Figure ---
# 2 rows, 1 column, 6 inches wide by 6 inches tall
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6))

color_time = 'tab:blue'
color_mem = 'tab:red'

# --- 3. Panel 1 (Top): Database Size Scaling ---
# Left Axis: Wall Clock Time
ax1.plot(db_df['N'], db_df['TIME'], marker='o', color=color_time)
ax1.set_xlabel('Database Size (N)')
ax1.set_ylabel('Time (seconds)', color=color_time)
ax1.tick_params(axis='y', labelcolor=color_time)
ax1.set_title('Scaling versus database size, 64 threads')
ax1.set_ylim(0, db_df['TIME'].max() * 1.20)

# Time labels above data points
for x, y in zip(db_df['N'], db_df['TIME']):
    ax1.annotate(f"{y:.1f}s", (x, y), textcoords="offset points", 
                 xytext=(0, 7), ha='center', fontsize=6, color=color_time)

# Right Axis: Peak Memory
ax1b = ax1.twinx()
ax1b.plot(db_df['N'], db_df['MEM_GB'], marker='s', linestyle='--', color=color_mem)
ax1b.set_ylabel('Peak Memory (GB)', color=color_mem)
ax1b.tick_params(axis='y', labelcolor=color_mem)
ax1b.set_ylim(bottom=0, top=4)

# --- 4. Panel 2 (Bottom): Thread Count Scaling ---
# Left Axis: Wall Clock Time
ax2.plot(th_df['THREADS'], th_df['TIME'], marker='o', color=color_time)
ax2.set_xlabel('Thread Count')
ax2.set_ylabel('Time (seconds)', color=color_time)
ax2.tick_params(axis='y', labelcolor=color_time)
ax2.set_title('Scaling versus thread count, database size=10,000')
ax2.set_ylim(0, th_df['TIME'].max() * 1.20)

# Time labels above data points
for x, y in zip(th_df['THREADS'], th_df['TIME']):
    ax2.annotate(f"{y:.1f}s", (x, y), textcoords="offset points", 
                 xytext=(-3, 7), ha='center', fontsize=6, color=color_time)

# Right Axis: Peak Memory
ax2b = ax2.twinx()
ax2b.plot(th_df['THREADS'], th_df['MEM_GB'], marker='s', linestyle='--', color=color_mem)
ax2b.set_ylabel('Peak Memory (GB)', color=color_mem)
ax2b.tick_params(axis='y', labelcolor=color_mem)
ax2b.set_ylim(bottom=0, top=4)

# --- 5. Save Figure ---
plt.tight_layout()

plt.savefig('figures/scaling_performance.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/scaling_performance.pdf', bbox_inches='tight')

print("Figure successfully generated at figures/scaling_performance.png")
