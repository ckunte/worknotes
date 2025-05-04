#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Hydrodynamic force on a tube
# C Kunte

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def create_depth_array(water_depth, num_points=200):
    # Creates an array of depths from seabed to surface
    return np.linspace(0, water_depth, num_points)


def interpolate_current_profile(current_data):
    # Interpolates the current profile based on depth and
    # velocity data
    depths = current_data[:, 0]
    velocities = current_data[:, 1]
    return interp1d(
        depths,
        velocities,
        kind="linear",
        fill_value="extrapolate",
    )


def calc_wave_number(
    wave_period,
    water_depth,
    tolerance=1e-6,
    max_iterations=100,
):
    # Calculates the wave number using the dispersion
    # relation (iterative)
    omega = 2 * np.pi / wave_period
    k = omega**2 / g  # Initial guess
    for _ in range(max_iterations):
        k_new = omega**2 / (g * np.tanh(k * water_depth))
        if np.abs(k_new - k) < tolerance:
            return k_new
        k = k_new
    return k  # Return the last estimate if no convergence


def wave_orbital_velo(
    depth,
    time,
    wave_height,
    wave_period,
    water_depth,
    wave_number,
):
    # Computes the wave orbital velocity at a given depth
    # and time
    omega = 2 * np.pi / wave_period
    amplitude = wave_height / 2
    return (
        amplitude
        * omega
        * np.cosh(wave_number * depth)
        / np.sinh(wave_number * water_depth)
        * np.cos(omega * time)
    )


def wave_orbital_accel(
    depth,
    time,
    wave_height,
    wave_period,
    water_depth,
    wave_number,
):
    # Computes the wave orbital acceleration at a given
    # depth and time
    omega = 2 * np.pi / wave_period
    amplitude = wave_height / 2
    return (
        -amplitude
        * omega**2
        * np.cosh(wave_number * depth)
        / np.sinh(wave_number * water_depth)
        * np.sin(omega * time)
    )


def hydrodynamic_force_at_depth(
    depth,
    time_points,
    wave_height,
    wave_period,
    water_depth,
    current_profile_func,
    wave_number,
    cd,
    cm,
    d_tube,
):
    # Calculates the hydrodynamic force per unit length at
    # a given depth over a wave cycle
    u_current = current_profile_func(depth)
    u_wave = wave_orbital_velo(
        depth,
        time_points,
        wave_height,
        wave_period,
        water_depth,
        wave_number,
    )
    du_dt_wave = wave_orbital_accel(
        depth,
        time_points,
        wave_height,
        wave_period,
        water_depth,
        wave_number,
    )
    u_total = u_wave + u_current
    f_drag = 0.5 * rho * cd * d_tube * u_total * np.abs(u_total)
    f_inertia = rho * cm * (np.pi * d_tube**2 / 4) * du_dt_wave
    f_total = f_drag + f_inertia
    return np.max(
        np.abs(f_total * 1e-3)
    )  # Peak absolute force in kN/m


def calc_force_profile(
    depth_array,
    time_points,
    wave_height,
    wave_period,
    water_depth,
    current_profile_func,
    wave_number,
    cd,
    cm,
    d_tube,
):
    # Calculates the hydrodynamic force profile over the
    # specified depth array.
    return [
        hydrodynamic_force_at_depth(
            z,
            time_points,
            wave_height,
            wave_period,
            water_depth,
            current_profile_func,
            wave_number,
            cd,
            cm,
            d_tube,
        )
        for z in depth_array
    ]


def plot_force_profile(
    force_profile_data,
    depth_array,
    wave_height,
    wave_period,
    water_depth,
    structure_data,
):
    # Plot hydrodynamic force profiles for multiple data sets
    # of Cd, Cm, and D_tube
    for i in range(len(force_profile_data)):
        force_profile = force_profile_data[i]
        cd = structure_data[i]["Cd"]
        cm = structure_data[i]["Cm"]
        d_tube = structure_data[i]["D_tube"]
        plt.plot(
            force_profile,
            depth_array,
            label=f"$C_d$: {cd}, $C_m$: {cm}, D: {d_tube} m",
        )

    plt.xlabel("Hydrodynamic force (kN/m)")
    plt.ylabel("Depth above seabed (m)")
    plt.title(
        f"Hydrodynamic load profile in {water_depth} m depth\nHmp: {wave_height} m; Tmp: {wave_period} s"
    )
    plt.grid(True)
    plt.legend()
    plt.savefig("hydforce.svg")
    plt.show()


if __name__ == "__main__":
    # --- INPUTS ---
    wave_height_input = 26.52  # wave height (m)
    wave_period_input = 13.58  # wave period (s)
    water_depth_input = 130.0  # total depth (m)
    current_data_input = np.array(
        [
            [0, 1.49],
            [5, 1.46],
            [10, 1.43],
            [20, 1.37],
            [30, 1.30],
            [40, 1.24],
            [50, 1.18],
            [60, 1.11],
            [70, 1.05],
            [80, 0.99],
            [90, 0.92],
            [100, 0.86],
            [110, 0.79],
            [120, 0.73],
            [130, 0.67],
        ]
    )
    # Define sets of Cd, Cm, and D_tube values
    struc_data_sets = [
        {"Cd": 1.05, "Cm": 1.20, "D_tube": 0.760},
        {"Cd": 2.0, "Cm": 1.2, "D_tube": 1.250},
    ]
    # --- END OF INPUTS ---

    # CONSTANTS
    rho = 1025  # seawater density (kg/m^3)
    g = 9.81  # gravitational acceleration (m/s^2)

    # Create depth array
    depths = create_depth_array(water_depth_input)

    # Interpolate current profile
    current_profile_interpolator = interpolate_current_profile(
        current_data_input
    )

    # Calculate wave number
    wave_number = calc_wave_number(
        wave_period_input, water_depth_input
    )

    # Time array for dynamic force
    time = np.linspace(0, wave_period_input, 200)  # one wave cycle

    # Compute force profiles for both sets of parameters
    all_force_profiles = []
    for data_set in struc_data_sets:
        force_profile = calc_force_profile(
            depths,
            time,
            wave_height_input,
            wave_period_input,
            water_depth_input,
            current_profile_interpolator,
            wave_number,
            data_set["Cd"],
            data_set["Cm"],
            data_set["D_tube"],
        )
        all_force_profiles.append(force_profile)

    # Plot the results
    plot_force_profile(
        all_force_profiles,
        depths,
        wave_height_input,
        wave_period_input,
        water_depth_input,
        struc_data_sets,
    )
