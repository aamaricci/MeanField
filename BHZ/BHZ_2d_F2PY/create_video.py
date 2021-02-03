# © 2015-2018, ETH Zurich, Institut für Theoretische Physik
# Author: Dominik Gresch <greschd@gmx.ch>

import os
import functools
import contextlib

import numpy as np
import phasemap as pm

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize

RESULTS_DIR = "results"
POINT_SIZE = 2.0
VALS = (0,1, 2, 3)
CUTOFF = None


def get_idx(name):
    return int(name.split(".")[0].split("_")[1])


def get_filenames():
    return sorted(os.listdir(RESULTS_DIR), key=get_idx)[:CUTOFF]


def get_results():
    for name in get_filenames():
        yield pm.io.load(os.path.join(RESULTS_DIR, name))


def init():
    fig, (ax, cbar_ax) = plt.subplots(
        ncols=2, gridspec_kw=dict(width_ratios=(0.95, 0.05))
    )
    ax.set_aspect(1.0)
    fig.subplots_adjust(right=0.9)

    norm = Normalize()
    norm.autoscale(VALS)
    color_vals = [norm(c) for c in VALS]
    c_bar = ColorbarBase(
        cbar_ax,
        values=VALS,
        cmap=plt.get_cmap(),
        norm=norm,
        boundaries=range(5),
        ticklocation="right",
        ticks=[x + 0.5 for x in range(4)],
    )
    c_bar.solids.set_edgecolor("k")
    c_bar.set_ticklabels(VALS)
    return fig, ax


def plot(res, ax):
    ax.cla()
    pm.plot.boxes(
        res,
        ax=ax,
        zorder=0,
        add_cbar=False,
        lw=0.1,
        edgecolor="k",
        scale_val=VALS,
    )
    pm.plot.points(
        res, ax=ax, edgecolors="k", lw=0.1, add_cbar=False, scale_val=VALS, s=POINT_SIZE
    )


if __name__ == "__main__":
    with contextlib.suppress(FileNotFoundError):
        os.remove("video.mp4")
    fig, ax = init()
    plot_func = functools.partial(plot, ax=ax)
    ani = animation.FuncAnimation(
        fig, plot_func, get_results, interval=120, repeat=False
    )
    ani.save("video.mp4", dpi=300, bitrate=2000, writer="ffmpeg_file")
