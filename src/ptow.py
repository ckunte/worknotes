#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Influence of vessel motions on transported piles in terms
of inertia forces, bending and shear stresses. 
2020 ckunte

Usage: ptow.py (--iner | --fb | --fv) [--tr=T1] [--tp=T2]
       ptow.py --help
       ptow.py --version

Options:
  --help      Show this help screen
  --iner      Plot inertia forces on pile during tow
  --fb        Plot bending stresses in pile during tow
  --fv        Plot shear stresses in pile during tow
  --tr=T1     Single amplitude roll period (s) [default: 10.]
  --tp=T2     Single amplitude pitch period (s) [default: 10.]
  --version   Show version.

"""
import numpy as np
import matplotlib.pyplot as plt
from docopt import docopt

args = docopt(
    __doc__,
    version="Influence of vessel motions on transported piles, v1.0.0",
)
# -------------------------------------------------------------------
#            Lp
#     |<-------------->|
#
#     +---+------------+-----------------+--------+----\
#     +---+---------------------------------------+----/
#     +---+------------+-----------------+--------+---\
#     +---+---------------------------------------+----\
#                (aft) +-----------------+-------------/
# ---------------------|-         Large barge       --/--------------
#                  --- +                              \  ---
#                   -   `-----------------------------/   -
#     .
#     .                    ELEVATION VIEW
#   ><>
# -------------------------------------------------------------------
# COMPUTE PILE SECTION PROPERTIES
def pipe_secprop(D, t):
    # Cross sectional ar, Lpea
    A = np.pi * (D - t) * t
    # Weight of overhung pile per unit length (in MN/m)
    # where 0.077 MN/m^3 => 7,850 kgf/m^3 (steel unit weight)
    w = A * 1.0 * 0.077
    # This is useful when MN/m^2 => MPa
    # Moment of inertia of the cross section
    I = (np.pi / 64.0) * (D ** 4 - (D - 2 * t) ** 4)
    # Polar moment of inertia of the cross section
    Ip = (np.pi / 32.0) * (D ** 4 - (D - 2 * t) ** 4)
    # Elastic section modulus
    Ze = I / (D / 2.0)
    # Plastic section modulus
    Zp = (1 / 6.0) * (D ** 3 - (D - 2 * t) ** 3)
    # Radius of gyration
    r = np.sqrt(I / A)
    return A, Ze, w


# COMPUTE INERTIA FORCES PER UNIT WEIGHT
def inertia(r, p):
    # r and p to be in radians
    r = r * (np.pi / 180.0)
    p = p * (np.pi / 180.0)
    # Angular acceleration: roll
    thta_r = r * (2 * np.pi / Tr) ** 2
    # Angular acceleration: pitch
    thta_p = p * (2 * np.pi / Tp) ** 2
    # Vertical force per unit mass
    Fvr = np.cos(r) + (L[1] / g) * thta_r + h * np.cos(r)
    Fvp = np.cos(p) + (L[0] / g) * thta_p + h * np.cos(p)
    # Horizontal force per unit mass
    Fhr = np.sin(r) + (L[2] / g) * thta_r + h * np.sin(r)
    Fhp = np.sin(p) + (L[2] / g) * thta_p + h * np.sin(p)
    # Resultant inertia force: roll
    Fr = np.sqrt(Fvr ** 2 + Fhr ** 2)
    # Resultant inertia force: pitch
    Fp = np.sqrt(Fvp ** 2 + Fhp ** 2)
    # Resultant inertia forces: quartering
    Fq = np.sqrt((0.6 * Fr) ** 2 + (0.6 * Fp) ** 2)
    # Return results
    return Fr, Fvr, Fhr, Fp, Fvp, Fhp, Fq


# COMPUTE BENDING STRESSES
def bending_stress(r, p):
    # Call results of inertia(r, p) function and multiply each value
    # of the tuple by an LRFD factor of 1.485 (= 1.1 * 1.35)
    # F = map(lambda x: x * 1.485, inertia(r, p)) # python2
    F = list(map(lambda x: x * 1.485, inertia(r, p)))  # for python3
    # Call results of pipe_secprop(D, t) function
    s = pipe_secprop(D, t)
    # For pile section overhung behind vessel stern (cantilever moment)
    # Bending stress (MPa) fb = M / Z, where M = (w * l) * l
    fbvr = ((s[2] * F[1]) * Lp ** 2) / s[1]
    fbhr = ((s[2] * F[2]) * Lp ** 2) / s[1]
    fbr = np.sqrt(fbvr ** 2 + fbhr ** 2)
    fbvp = ((s[2] * F[4]) * Lp ** 2) / s[1]
    fbhp = ((s[2] * F[5]) * Lp ** 2) / s[1]
    fbp = np.sqrt(fbvp ** 2 + fbhp ** 2)
    # where, F[1] => Fvr; F[2] => Fhr; F[4] => Fvp; F[5] => Fhp, and
    # where, s[1] => Ze; s[2] => w
    return fbr, fbvr, fbhr, fbp, fbvp, fbhp


# COMPUTE SHEAR STRESSES
def shear_stress(r, p):
    # F = map(lambda x: x * 1.485, inertia(r, p)) # python2
    F = list(map(lambda x: x * 1.485, inertia(r, p)))  # for python3
    s = pipe_secprop(D, t)
    # For pile section overhung behind vessel stern
    # Shear stress (MPa) fv = 2V / A, where V = (w * l)
    fvvr = (2.0 * ((s[2] * F[1]) * Lp)) / s[0]
    fvhr = (2.0 * ((s[2] * F[2]) * Lp)) / s[0]
    fvr = np.sqrt(fvvr ** 2 + fvhr ** 2)
    fvvp = (2.0 * ((s[2] * F[4]) * Lp)) / s[0]
    fvhp = (2.0 * ((s[2] * F[5]) * Lp)) / s[0]
    fvp = np.sqrt(fvvp ** 2 + fvhp ** 2)
    # where, s[0] => A
    return fvr, fvvr, fvhr, fvp, fvvp, fvhp


# PLOT FUNCTIONS
def misc():
    plt.xlabel("Motion angle (deg)")
    plt.legend(loc=0)
    plt.grid(True)
    pass


# PLOT MOTION v. INERTIA FORCE (ROLL)
def plot_roll_motion_inertia(r, p):
    F = inertia(r, p)
    plt.plot(
        r,
        F[0],
        label="$F_{r} = \\sqrt{F_{vr}^2 + F_{hr}^2}$",
        linewidth=2,
    )
    plt.plot(r, F[1], label="$F_{vr}$", linewidth=2)
    plt.plot(r, F[2], label="$F_{hr}$", linewidth=2)
    plt.ylabel("Inertia force in terms of unit weight of pile, W")
    misc()
    plt.savefig("pp-inertia-roll.png")
    plt.close()
    pass


# PLOT MOTION v. INERTIA FORCE (PITCH)
def plot_pitch_motion_inertia(r, p):
    F = inertia(r, p)
    plt.plot(
        p,
        F[3],
        label="$F_{p} = \\sqrt{F_{vp}^2 + F_{hp}^2}$",
        linewidth=2,
    )
    plt.plot(p, F[4], label="$F_{vp}$", linewidth=2)
    plt.plot(p, F[5], label="$F_{vp}$", linewidth=2)
    plt.ylabel("Inertia force in terms of unit weight of pile, W")
    misc()
    plt.savefig("pp-inertia-pitch.png")
    plt.close()
    pass


# PLOT MOTION v. BENDING STRESS (ROLL)
def plot_roll_motion_bendingstress(r, p):
    fb = bending_stress(r, p)
    plt.plot(
        r,
        fb[0],
        label="$f_{br} = \\sqrt{f_{bvr}^2 + f_{bhr}^2}$",
        linewidth=2,
    )
    plt.plot(r, fb[1], label="$f_{bvr}$", linewidth=2)
    plt.plot(r, fb[2], label="$f_{bhr}$", linewidth=2)
    plt.ylabel("Bending stress (MPa)")
    misc()
    plt.title(
        "Pile size: %0.0f$\\times$ %0.0f (D/t=%0.1f). Overhung span: %0.0fm aft"
        % (D * 1e3, t * 1e3, (D / t), Lp)
    )
    plt.axhspan(
        396.0, 500.0, linewidth=0, facecolor="r", alpha=0.18
    )  # where 396MPa is the bending strength of section
    plt.savefig("pp-fb-roll.png")
    plt.close()
    pass


# PLOT MOTION v. BENDING STRESS (PITCH)
def plot_pitch_motion_bendingstress(r, p):
    fb = bending_stress(r, p)
    plt.plot(
        p,
        fb[3],
        label="$f_{bp} = \\sqrt{f_{bvp}^2 + f_{bhp}^2}$",
        linewidth=2,
    )
    plt.plot(p, fb[4], label="$f_{bvp}$", linewidth=2)
    plt.plot(p, fb[5], label="$f_{bhp}$", linewidth=2)
    plt.ylabel("Bending stress (MPa)")
    misc()
    plt.title(
        "Pile size: %0.0f$\\times$ %0.0f (D/t=%0.1f). Overhung span: %0.0fm aft"
        % (D * 1e3, t * 1e3, (D / t), Lp)
    )
    plt.axhspan(396.0, 600.0, linewidth=0, facecolor="r", alpha=0.18)
    plt.savefig("pp-fb-pitch.png")
    plt.close()
    pass


# PLOT MOTION v. SHEAR STRESS (ROLL)
def plot_roll_motion_shearstress(r, p):
    fv = shear_stress(r, p)
    plt.plot(
        r,
        fv[0],
        label="$f_{vr} = \\sqrt{f_{vvr}^2 + f_{vhr}^2}$",
        linewidth=2,
    )
    plt.plot(r, fv[1], label="$f_{vvr}$", linewidth=2)
    plt.plot(r, fv[2], label="$f_{vhr}$", linewidth=2)
    plt.ylabel("Shear stress (MPa)")
    misc()
    plt.title(
        "Pile size: %0.0f$\\times$ %0.0f (D/t=%0.1f). Overhung span: %0.0fm aft"
        % (D * 1e3, t * 1e3, (D / t), Lp)
    )
    plt.savefig("pp-fv-roll.png")
    plt.close()
    pass


# PLOT MOTION v. SHEAR STRESS (PITCH)
def plot_pitch_motion_shearstress(r, p):
    fv = shear_stress(r, p)
    plt.plot(
        p,
        fv[3],
        label="$f_{vp} = \\sqrt{f_{vvp}^2 + f_{vhp}^2}$",
        linewidth=2,
    )
    plt.plot(p, fv[4], label="$f_{vvp}$", linewidth=2)
    plt.plot(p, fv[5], label="$f_{vhp}$", linewidth=2)
    plt.ylabel("Shear stress (MPa)")
    misc()
    plt.title(
        "Pile size: %0.0f$\\times$ %0.0f (D/t=%0.1f). Overhung span: %0.0fm aft"
        % (D * 1e3, t * 1e3, (D / t), Lp)
    )
    plt.savefig("pp-fv-pitch.png")
    plt.close()
    pass


def main():
    # -- BEGIN USER INPUTS --
    g = 9.81  # Acceleration due to gravity (m/s^2)
    # Cargo location w.r.t vessel:
    # Lever arm (x, y, z) between vessel C.O.R to overhung pile C.O.G (m):
    L = [95.00, 20.00, 15.25]
    # Steel pile properties:
    D = 3.000  # Pile diameter (m)
    t = 0.038  # Pile wall thickness (m)
    Lp = 35.0  # Overhung length of pile (m)
    # Motion characteristics:
    r = np.linspace(0, 35)  # roll angle range (to process)
    Tr = float(args["--tr"])  # Full cycle period: roll (default: 10s)
    p = np.linspace(0, 20)  # pitch angle range (to process)
    Tp = float(args["--tp"])  # Full cycle period: pitch (default: 10s)
    h = 0.20  # Heave amplitude (h) in terms of g
    # -- END USER INPUTS --
    # process option
    if args["--iner"]:
        plot_roll_motion_inertia(r, p)
        plot_pitch_motion_inertia(r, p)
    elif args["--fb"]:
        plot_roll_motion_bendingstress(r, p)
        plot_pitch_motion_bendingstress(r, p)
    elif args["--fv"]:
        plot_roll_motion_shearstress(r, p)
        plot_pitch_motion_shearstress(r, p)
    else:
        print("Please select plot option. Try: python ptow.py --help")
    pass


if __name__ == "__main__":
    main()