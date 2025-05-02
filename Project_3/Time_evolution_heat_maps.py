# -*- coding: utf-8 -*-
"""Temperature Profile Heatmap Generator

This script reads a specially-formatted CSV file containing temperature profiles 
along a robot’s height measured over time. The file is divided into blocks 
separated by empty lines, where each block corresponds to a height `z` on the robot.

Each block:
- Contains time and multiple temperature sensor readings.
- Is processed to average temperature readings across sensors.
- Is stored in a dictionary mapping height `z` to (time, averaged temperatures).

Finally, the script constructs a temperature matrix and visualizes it as a 
heatmap showing the temperature evolution along the robot's height over a 
7-hour period.

Author: Parina
Created: May 1, 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Input CSV file containing temperature profile data
FILENAME = 'full_temp_profile_7hr.csv'

# Dictionary to store temperature profiles by height (z-value)
z_profiles = {}
current_block = []

# Read the CSV file block by block (blocks are separated by empty lines)
with open(FILENAME, 'r') as file:
    for line in file:
        line = line.strip()
        if not line:
            # Process the current block when an empty line is encountered
            if current_block:
                df_block = pd.DataFrame([row.split(',') for row in current_block])
                df_block = df_block.apply(pd.to_numeric, errors='coerce')

                if df_block.shape[1] >= 3:
                    z_val = df_block.iloc[0, 0]
                    time = df_block.iloc[:, 1]
                    temps = df_block.iloc[:, 2:].mean(axis=1)  # Average across sensors
                    z_profiles[z_val] = (time.values, temps.values)

                current_block = []
        else:
            current_block.append(line)

# Process the final block (in case the file doesn't end with an empty line)
if current_block:
    df_block = pd.DataFrame([row.split(',') for row in current_block])
    df_block = df_block.apply(pd.to_numeric, errors='coerce')

    if df_block.shape[1] >= 3:
        z_val = df_block.iloc[0, 0]
        time = df_block.iloc[:, 1]
        temps = df_block.iloc[:, 2:].mean(axis=1)
        z_profiles[z_val] = (time.values, temps.values)

# Sort the z-values to organize matrix rows from bottom to top of the robot
z_sorted = sorted(z_profiles.keys())

# Multiply time by 54 (unit correction or sampling interval factor)
time_array = z_profiles[z_sorted[0]][0] * 54

# Create a 2D temperature matrix: rows = z-values, columns = time samples
temp_matrix = np.array([z_profiles[z][1] for z in z_sorted])

# Plot the temperature profile as a heatmap
plt.figure(figsize=(12, 6))
plt.imshow(
    temp_matrix,
    aspect='auto',
    cmap='inferno',
    extent=[
        time_array[0],
        time_array[-1],
        z_sorted[0],
        z_sorted[-1]
    ],
    origin='lower'
)

plt.colorbar(label='Temperature (K)')
plt.xlabel('Time (s)')
plt.ylabel('Height along Robot (m)')
plt.title('Temperature Evolution Along Robot Height (7 hr)')
plt.tight_layout()
plt.show()
