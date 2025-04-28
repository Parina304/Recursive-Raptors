# -*- coding: utf-8 -*-
"""Visualize thickness evolution for multiple materials over time.

This script generates:
1. 3D surface plots for each material showing thickness vs. Z position and time.
2. 2D line plot of Z position vs. average thickness for each material.

Author: Layan
Created: 2025-04-24
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
FILENAME = r"C:\Users\layan\Downloads\testfile_thickness_data_with_4_layers_.csv"
MATERIALS = ['Glue', 'Carbon', 'Steel', 'Thermal Protection']


def load_data(filename: str) -> pd.DataFrame:
    """Loads the CSV data into a DataFrame.

    Args:
        filename: The full path to the CSV file.

    Returns:
        A pandas DataFrame containing the raw thickness data.
    """
    return pd.read_csv(filename)


def plot_3d_surfaces(df: pd.DataFrame, materials: list[str]) -> None:
    """Creates 3D surface plots for each material.

    Args:
        df: DataFrame containing all material data.
        materials: List of material names to plot.
    """
    fig = plt.figure(figsize=(16, 10))
    for i, material in enumerate(materials, start=1):
        material_df = df[df['material'] == material]

        pivot_df = material_df.pivot_table(index='z', columns='time', values='x')
        pivot_df = pivot_df.sort_index()

        z_vals = pivot_df.index.values
        t_vals = pivot_df.columns.values
        Z, T = np.meshgrid(z_vals, t_vals, indexing='ij')
        X = pivot_df.values

        ax = fig.add_subplot(2, 2, i, projection='3d')
        surf = ax.plot_surface(Z, T, X, cmap='viridis', edgecolor='k', linewidth=0.5)

        ax.set_xlabel('Z Position (m)')
        ax.set_ylabel('Time (s)')
        ax.set_zlabel('Thickness (mm)')
        ax.set_title(f'{material} Thickness')

    fig.suptitle('Thickness Evolution Over Time for All Materials', fontsize=18)
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.show()


def plot_avg_thickness_vs_z(df: pd.DataFrame, materials: list[str]) -> None:
    """Creates a 2D plot of average thickness vs. Z position for each material.

    Args:
        df: DataFrame containing all material data.
        materials: List of material names to plot.
    """
    plt.figure(figsize=(10, 6))

    for material in materials:
        material_df = df[df['material'] == material]
        grouped = material_df.groupby(['z', 'time'])['x'].mean().reset_index()
        pivot = grouped.pivot(index='z', columns='time', values='x')
        avg_thickness = pivot.mean(axis=1)

        plt.plot(avg_thickness.index, avg_thickness.values, label=material)

    plt.xlabel('Z Position (m)')
    plt.ylabel('Average Thickness (mm)')
    plt.title('Z Position vs Average Thickness for Each Material')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def main() -> None:
    """Main function to execute the visualization pipeline."""
    df = load_data(FILENAME)
    plot_3d_surfaces(df, MATERIALS)
    plot_avg_thickness_vs_z(df, MATERIALS)


if __name__ == '__main__':
    main()
