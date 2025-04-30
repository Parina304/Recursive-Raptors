# -*- coding: utf-8 -*-
"""Insulation Thickness Profile Plotter

This script visualizes how the thickness of an insulation layer changes
along the Z-position of a humanoid structure over time (1 hour, 3 hours, 
and 7 hours). It reads in three CSV files (with columns "z value" and 
"thickness (cm)"), each representing a snapshot of thickness distribution 
at a different time point.

The output is a line plot showing all three profiles for easy comparison.

Author: Layan
Modified by: Parina
Date: 2025-04-30
"""

import pandas as pd
import matplotlib.pyplot as plt

# File paths to CSV data for 1, 3, and 7 hours
FILENAME_1HR = (
    r"C:\Users\Parina\Desktop\Computing Concepts\Recursive-Raptors\Project_3\Thickness_1_hr.csv"
)
FILENAME_3HR = (
    r"C:\Users\Parina\Desktop\Computing Concepts\Recursive-Raptors\Project_3\Thickness_3_hr.csv"
)
FILENAME_7HR = (
    r"C:\Users\Parina\Desktop\Computing Concepts\Recursive-Raptors\Project_3\Thickness_7_hr.csv"
)


def load_data(filename: str) -> pd.DataFrame:
    """Loads a CSV file into a pandas DataFrame.

    Args:
        filename: Full path to the CSV file.

    Returns:
        A pandas DataFrame with the loaded data.
    """
    return pd.read_csv(filename)


def plot_multiple_thickness_profiles(dfs: dict[str, pd.DataFrame]) -> None:
    """Plots thickness vs. Z-position for multiple time points.

    Args:
        dfs: A dictionary mapping labels (e.g., "1 Hour") to their DataFrames.
    """
    plt.figure(figsize=(10, 6))

    # Plot each dataset with its label
    for label, df in dfs.items():
        plt.plot(df["z value"], df["thickness (cm)"], marker="o", label=label)

    # Configure plot aesthetics
    plt.xlabel("Z Position (m)")
    plt.ylabel("Insulation Thickness (cm)")
    plt.title("Insulation Thickness Distribution Along Z Over Time")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def main() -> None:
    """Main function to load data and generate the plot."""
    # Load each dataset into a DataFrame
    df_1hr = load_data(FILENAME_1HR)
    df_3hr = load_data(FILENAME_3HR)
    df_7hr = load_data(FILENAME_7HR)

    # Sanity check to ensure required columns are present
    for name, df in zip(["1 hr", "3 hr", "7 hr"], [df_1hr, df_3hr, df_7hr]):
        print(f"Columns in {name} dataset: {df.columns.tolist()}")
        if "z value" not in df.columns or "thickness (cm)" not in df.columns:
            raise ValueError(
                f"CSV for {name} must contain 'z value' and 'thickness (cm)' columns."
            )

    # Organize the DataFrames with labels for plotting
    dfs = {
        "1 Hour": df_1hr,
        "3 Hours": df_3hr,
        "7 Hours": df_7hr,
    }

    # Generate the plot
    plot_multiple_thickness_profiles(dfs)


if __name__ == "__main__":
    main()
