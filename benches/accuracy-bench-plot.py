#/usr/bin/env 
try:
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ImportError:
    print("Missing required Python dependencies:")
    print("Please run: pip install matplotlib numpy pandas")
    sys.exit(1)



df = pd.read_csv('accuracy.tsv', sep='\t')

# Convert percentage columns to numeric
df['Mean Bias'] = df['Mean Bias'].str.replace('%', '').astype(float)
df['Mean Absolute Error'] = df['Mean Absolute Error'].str.replace('%', '').astype(float)

# Convert 'Mapped Proportion' from fraction string to float
df['Mapped Proportion'] = df['Mapped Proportion'].apply(lambda x: float(x.split('/')[0]) / float(x.split('/')[1]) if '/' in str(x) else float(x))

# Prepare data
plot_df = df.sort_values(['ANI', 'Length'])
anis = sorted(plot_df['ANI'].unique())
norm = plt.Normalize(min(anis), max(anis))

# You can easily try: 'viridis', 'plasma', 'inferno', 'magma', 'cividis'
cmap_name = 'magma'
colormap = plt.get_cmap(cmap_name)

# Create a 3-row figure - compact size
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(5, 7), sharex=True)

for i, ani in enumerate(anis):
    subset = plot_df[plot_df['ANI'] == ani]
    color = colormap(norm(ani))

    # Plot 1: Mapping Sensitivity
    ax1.plot(subset['Length'], subset['Mapped Proportion'],
             color=color, label=f'ANI {ani}', alpha=0.8, linewidth=2)

    # Plot 2: Top-hit coverage
    ax2.plot(subset['Length'], subset['Mean Length'] / subset['Length'],

             color=color, alpha=0.8, linewidth=2)

    # Plot 3: Mean ANI Error
    ax3.fill_between(subset['Length'],
                     subset['Mean Bias'] - subset['Mean Absolute Error'],
                     subset['Mean Bias'] + subset['Mean Absolute Error'],
                     color=color, alpha=0.1)
    ax3.plot(subset['Length'], subset['Mean Bias'], color=color, linewidth=1.5)

# Subplot 1 styling
ax1.set_ylabel('Mapping Sensitivity', fontsize=14)
#ax1.set_title(f'alamem performance metrics', fontsize=16, loc='left', alpha=0.6)
ax1.grid(True, linestyle=':', alpha=0.4)
ax1.legend(title='ANI', loc='lower right', fontsize='medium', ncol=2, framealpha=0.8, title_fontsize=12)

# Subplot 2 styling
ax2.set_ylabel('Top-hit coverage', fontsize=14)
ax2.grid(True, linestyle=':', alpha=0.4)

# Subplot 3 styling
ax3.axhline(0, color='black', linestyle='--', alpha=0.5)
ax3.set_ylabel('Mean ANI Error\n (SD ribbons)', fontsize=14)
ax3.set_xlabel('Element Length', fontsize=14)
ax3.set_ylim(-2, 2)

ax3.grid(True, linestyle=':', alpha=0.4)

plt.tight_layout()
#plt.show()
plt.savefig('accuracy.png')
