"""Figures (RAST_Figuremaking.m) and the pictures saved by RAST_picture_segmentation.m.

picture_<i>_ST_<s>_linesselected.png: the smoothed picture with every line of score >= s drawn in the colour of
its best-matching family, and the theoretical trace directions of each family on the right.
Colours: {100} blue, {110} orange, {112} green, {123} yellow, special case magenta (bcc: matches only {112}
and/or {123}, i.e. unambiguously <111> slip; fcc: matches more than one {111} plane), no match red.
"""
import math
from pathlib import Path

import numpy as np

from . import io_mat, match

COLOURS = {100: (0, 0, 1), 110: (1, 128 / 255, 0), 112: (0, 1, 0), 123: (1, 1, 0), 111: (1, 0, 1), 0: (1, 0, 0)}
ARROW_COLOURS = [(0, 0, 1), (1, 128 / 255, 0), (0, 1, 0), (230 / 255, 230 / 255, 0)]
NAMES = {1: ["{100}", "{110}", "{112}", "{123}"], 2: ["{111}", "{-111}", "{1-11}", "{11-1}"]}


def _plt():
    import matplotlib
    matplotlib.use("Agg")           # files only, no window
    import matplotlib.pyplot as plt
    return plt


def save_picture(img, path):
    """Save a uint8 / uint16 / logical picture as PNG (MATLAB imwrite)."""
    from PIL import Image
    a = np.asarray(img)
    if a.dtype == bool:
        im = Image.fromarray((a * 255).astype(np.uint8))
    elif a.dtype == np.uint16:
        im = Image.fromarray(a, mode="I;16")
    else:
        im = Image.fromarray(a.astype(np.uint8))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    im.save(io_mat.long_path(path))


def lines_figure(I, picture_no, st, pixsize, ir, coefmat, code, numlines, euler, degrot, crystal, path):
    plt = _plt()
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(6, 9)
    ax = fig.add_subplot(gs[:, :6])
    ax.imshow(I, cmap="gray", interpolation="nearest")
    ax.set_axis_off()
    ax.set_title(f"Picture {picture_no}: lines with a score >= {st}\n"
                 f"euler = {euler} / N = {pixsize}*{pixsize} px / IR = {ir}° rotation correction = {degrot:g}°",
                 fontsize=11)
    cuth, cutl, _ = coefmat.shape
    for i in range(cuth):
        for j in range(cutl):
            for k in range(numlines):
                o = numlines + 6 * k
                if coefmat[i, j, o + 4] >= st:
                    c = COLOURS.get(int(code[i, j, k]))
                    if c is None:
                        continue
                    ax.plot([coefmat[i, j, o + 0], coefmat[i, j, o + 2]], [coefmat[i, j, o + 1], coefmat[i, j, o + 3]],
                            "-", lw=1, color=c)

    slips, _ = match.theoretical_traces(euler, crystal)
    for f in range(4):
        a2 = fig.add_subplot(gs[f + 1, 6])
        for a, b in slips[f]:
            if degrot == 0:
                x, y = a, b
            else:      # as in MATLAB: rotate the trace vectors by +degrot
                x = a * math.cos(math.radians(degrot)) - b * math.sin(math.radians(degrot))
                y = a * math.sin(math.radians(degrot)) + b * math.cos(math.radians(degrot))
            a2.annotate("", xy=(-x, y), xytext=(0, 0),
                        arrowprops=dict(arrowstyle="->", color=ARROW_COLOURS[f], lw=1))
        a2.add_patch(plt.Circle((0, 0), 1, fill=False, color=ARROW_COLOURS[f], lw=1.5))
        a2.set_title(NAMES[crystal][f], fontsize=10)
        a2.set_xlim(-1.05, 1.05)
        a2.set_ylim(-1.05, 1.05)
        a2.set_aspect("equal")
        a2.set_axis_off()
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(io_mat.long_path(path), dpi=200, bbox_inches="tight")
    plt.close(fig)
