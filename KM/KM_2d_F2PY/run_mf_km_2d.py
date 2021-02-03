#!/usr/bin/env python
import os
import sys
import logging
import phasemap as pm
import matplotlib.pyplot as plt
import mf_km
from mf_km import *



plt.set_cmap("viridis")

POINT_SIZE = 1
num_step=6
init_mesh=4

def phase_fct(pos):
    t2,uloc = pos
    res = mf_km.mf_km_2d(t2,uloc)
    return res
    

def run(num_steps):
    os.makedirs("results", exist_ok=True)
    return pm.run(
        phase_fct,
        [(0, 6), (0, 6)],
        num_steps=num_steps,
        mesh=init_mesh,
        save_file="results/res_{}.json",
        save_interval=0.0,
        load=True
    )


def plot_boxes(res):
    pm.plot.boxes(res)
    plt.savefig("boxes.pdf", bbox_inches="tight")


def plot_points(res):
    pm.plot.points(res, s=POINT_SIZE, lw=0.0)
    plt.savefig("points.pdf", bbox_inches="tight")


def plot_combined(res):
    fig, ax = plt.subplots(figsize=[4.2, 4])
    ax.set_aspect(1.0)
    pm.plot.boxes(res, ax=ax, zorder=0, add_cbar=False, lw=0.1, edgecolor="k")
    pm.plot.points(res, ax=ax, edgecolors="k", lw=0.1, s=POINT_SIZE)
    plt.savefig("combined.pdf", bbox_inches="tight")


if __name__ == "__main__":
    res = run(num_step)
    plot_boxes(res)
    plot_points(res)
    plot_combined(res)
