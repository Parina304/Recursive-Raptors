# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 20:13:08 2025

@author: Parina
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# === Step 1: Load and inspect data ===
# Replace 'your_data.csv' with the actual filename
filename = r"C:\Users\Parina\Desktop\Computing Concepts\Recursive-Raptors\Project_3\testfile_thickness_data_with_4_layers_.csv"

df = pd.read_csv(filename)

# Print a few rows to check structure
print(df.head())

# === Step 2: Optional: Filter for a single material ===
# For example, we look at only "Carbon" rows
material_filter = 'Glue'
df = df[df['material'] == material_filter]

# === Step 3: Pivot data into z (position) vs time vs thickness (x) ===
pivot_df = df.pivot_table(index='z', columns='time', values='x')

# Ensure it's sorted
pivot_df = pivot_df.sort_index()

# === Step 4: Prepare meshgrid ===
Z_vals = pivot_df.index.values
T_vals = pivot_df.columns.values
Z, T = np.meshgrid(Z_vals, T_vals, indexing='ij')
X = pivot_df.values

# === Step 5: Plot the surface ===
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(Z, T, X, cmap='viridis', edgecolor='k', linewidth=0.5)

ax.set_xlabel('Z Position (m)')
ax.set_ylabel('Time (s)')
ax.set_zlabel('Thickness (mm)')
ax.set_title(f'Thickness Evolution Over Time for {material_filter}')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10)
plt.show()
