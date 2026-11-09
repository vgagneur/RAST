"""`rast prepare`: make a RAST dataset from SEM pictures (port of RAST_Loading.m).

For each picture: convert to grey levels, crop (optional), mask the areas that must not be analysed (indent,
scratches, other grains, dust, strong straight features that are not slip traces), and give the Euler angles
of the grain (Bunge convention, degrees). The dataset is saved as a .mat that `rast run` and the MATLAB
RAST_Running.m both read.

Three ways to use it:
  interactive          rast prepare --out my_sample.mat [--pictures a.tif b.tif ...] [--eulers angles.csv]
                       a window opens for each picture (crop, then draw mask polygons); Euler angles are typed
                       in the terminal unless given in --eulers.
  from files (no GUI)  rast prepare --out my_sample.mat --pictures a.tif b.tif --masks mask_folder --eulers angles.csv
                       masks are PNGs named <picture stem>_mask.png (white = masked), same size as the picture.
  edit a dataset       rast prepare --edit my_sample.mat [--select 3 5] [--eulers angles.csv] [--out new.mat]
                       add mask areas to chosen pictures and/or change their Euler angles.

Euler CSV: a header line `picture,phi1,Phi,phi2`, then one line per picture; `picture` is the picture's file name
(or stem) or its 1-based number in the dataset.

Pictures of the same grain must be entered one after another (the results table numbers grains by
consecutive identical Euler angles). Masked pixels get the value 1, as in MATLAB.
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from . import io_mat

# MATLAB rgb2gray weights (the inverse of the NTSC YIQ matrix row, as rgb2gray.m uses them)
RGB2GRAY = (0.298936021293775, 0.587043074451121, 0.114020904255103)

INSTRUCTIONS_CROP = ("CROP: drag a rectangle around the part to keep (e.g. leave out the SEM legend).\n"
                     "Enter = accept   n = keep the whole picture")
INSTRUCTIONS_MASK = ("MASK: click the corners of a polygon around an area to hide (indent, scratch, other grain ...),\n"
                     "click the first corner again to close it, then press a to add it to the mask.\n"
                     "a = add polygon   u = undo last   c = clear all   Enter = done")


# ------------------------------------------------------------------ picture handling
def load_grey(path):
    """Read a picture file and return a 2-D uint8 / uint16 grey picture (MATLAB imread + im2gray)."""
    from PIL import Image
    im = Image.open(io_mat.long_path(path))
    a = np.asarray(im)
    if a.ndim == 3:
        a = a[:, :, :3]
        dt = a.dtype
        g = (RGB2GRAY[0] * a[:, :, 0].astype(np.float64) + RGB2GRAY[1] * a[:, :, 1] + RGB2GRAY[2] * a[:, :, 2])
        top = np.iinfo(dt).max if np.issubdtype(dt, np.integer) else 255
        a = np.clip(np.floor(g + 0.5), 0, top).astype(dt if dt in (np.uint8, np.uint16) else np.uint8)
    if a.dtype == bool:
        a = a.astype(np.uint8) * 255
    if a.dtype not in (np.uint8, np.uint16):
        if np.issubdtype(a.dtype, np.integer) and a.min() >= 0 and a.max() <= 65535:
            a = a.astype(np.uint16 if a.max() > 255 else np.uint8)
        else:
            raise ValueError(f"{path}: unsupported picture type {a.dtype}; save it as 8- or 16-bit grey or RGB")
    return np.ascontiguousarray(a)


def polygon_mask(shape, xs, ys):
    """Pixels whose centre is inside the polygon (x = column, y = row, 1-based as in MATLAB roipoly)."""
    from matplotlib.path import Path as MPath
    h, w = shape
    yy, xx = np.mgrid[1:h + 1, 1:w + 1]
    pts = np.column_stack([xx.ravel(), yy.ravel()])
    return MPath(np.column_stack([xs, ys])).contains_points(pts).reshape(h, w)


def apply_mask(I, mask):
    I = I.copy()
    I[mask] = 1            # MATLAB I(BWs)=0.5 stores 1 in an integer picture
    return I


def read_euler_csv(path):
    """-> dict key -> Euler text; keys are file names, stems and/or 1-based numbers (as text)."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        need = {"picture", "phi1", "Phi", "phi2"}
        if not r.fieldnames or not need <= set(r.fieldnames):
            raise SystemExit(f"{path}: the header must be picture,phi1,Phi,phi2 (got {r.fieldnames})")
        for row in r:
            key = row["picture"].strip()
            out[key] = io_mat.format_euler([row["phi1"], row["Phi"], row["phi2"]])
    return out


def euler_for(table, number, file=None):
    """Look a picture up in an Euler CSV table by file name, stem or number."""
    if table is None:
        return None
    for key in ([Path(file).name, Path(file).stem] if file else []) + [str(number)]:
        if key in table:
            return table[key]
    return None


def ask_euler(label):
    while True:
        s = input(f"Euler angles for {label} (Bunge, degrees) as [phi1,Phi,phi2]: ").strip()
        try:
            return io_mat.format_euler(io_mat.parse_euler(s))
        except ValueError as e:
            print(f"  {e}; try again, e.g. [90,50,131]")


# ------------------------------------------------------------------ interactive windows
def _gui():
    import matplotlib
    try:
        matplotlib.use("TkAgg")
    except Exception:
        pass
    import matplotlib.pyplot as plt
    if matplotlib.get_backend().lower() == "agg":
        raise SystemExit("No window can be opened here (no desktop session). Use --masks and --eulers instead.")
    return plt


def gui_crop(I, title):
    plt = _gui()
    from matplotlib.widgets import RectangleSelector
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.imshow(I, cmap="gray", interpolation="nearest")
    ax.set_title(f"{title}\n{INSTRUCTIONS_CROP}", fontsize=10, color="red")
    state = {"box": None, "done": False}

    def on_select(e1, e2):
        state["box"] = (e1.xdata, e1.ydata, e2.xdata, e2.ydata)

    sel = RectangleSelector(ax, on_select, useblit=True, interactive=True, button=[1])

    def on_key(ev):
        if ev.key == "enter":
            state["done"] = True
            plt.close(fig)
        elif ev.key == "n":
            state["box"] = None
            state["done"] = True
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()
    del sel
    if not state["box"]:
        return I, None
    x0, y0, x1, y1 = state["box"]
    c0, c1 = sorted((int(round(x0)), int(round(x1))))
    r0, r1 = sorted((int(round(y0)), int(round(y1))))
    h, w = I.shape
    c0, r0 = max(c0, 0), max(r0, 0)
    c1, r1 = min(c1, w - 1), min(r1, h - 1)
    return I[r0:r1 + 1, c0:c1 + 1].copy(), (r0, r1, c0, c1)


def gui_mask(I, title, start_mask=None):
    """Draw one or more polygons; returns (mask, xs, ys) with NaN between polygons in xs / ys."""
    plt = _gui()
    from matplotlib.widgets import PolygonSelector
    mask = np.zeros(I.shape, bool) if start_mask is None else start_mask.copy()
    polys = []
    fig, ax = plt.subplots(figsize=(11, 8))
    im = ax.imshow(I, cmap="gray", interpolation="nearest")
    overlay = ax.imshow(np.zeros((*I.shape, 4)), interpolation="nearest")
    ax.set_title(f"{title}\n{INSTRUCTIONS_MASK}", fontsize=10, color="red")
    state = {"verts": None}

    def redraw():
        m = mask.copy()
        for xs, ys in polys:
            m |= polygon_mask(I.shape, xs, ys)
        rgba = np.zeros((*I.shape, 4))
        rgba[m] = (1, 0, 0, 0.35)
        overlay.set_data(rgba)
        fig.canvas.draw_idle()

    def on_select(verts):
        state["verts"] = verts

    sel = PolygonSelector(ax, on_select, useblit=True)

    def on_key(ev):
        nonlocal sel
        if ev.key == "a" and state["verts"] and len(state["verts"]) >= 3:
            v = np.asarray(state["verts"])
            polys.append((v[:, 0] + 1, v[:, 1] + 1))     # matplotlib pixel centres are 0-based
            state["verts"] = None
            sel.clear()
            redraw()
        elif ev.key == "u" and polys:
            polys.pop()
            redraw()
        elif ev.key == "c":
            polys.clear()
            mask[:] = False if start_mask is None else start_mask
            redraw()
        elif ev.key == "enter":
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    redraw()
    plt.show()
    del sel, im
    m = mask.copy()
    xs, ys = [], []
    for px, py in polys:
        m |= polygon_mask(I.shape, px, py)
        xs += list(px) + [np.nan]
        ys += list(py) + [np.nan]
    return m, np.array(xs[:-1] if xs else []), np.array(ys[:-1] if ys else [])


def gui_confirm(I, mask, title):
    """Show the masked picture; returns 'keep', 'redo' or 'skip'."""
    plt = _gui()
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.imshow(apply_mask(I, mask), cmap="gray", interpolation="nearest")
    ax.set_title(f"{title}\nk = keep this picture   r = redo crop and mask   s = skip this picture",
                 fontsize=10, color="red")
    state = {"ans": "keep"}

    def on_key(ev):
        if ev.key in ("k", "r", "s", "enter"):
            state["ans"] = {"k": "keep", "r": "redo", "s": "skip", "enter": "keep"}[ev.key]
            plt.close(fig)

    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()
    return state["ans"]


# ------------------------------------------------------------------ the three modes
def prepare_interactive(files, euler_table, no_crop=False):
    pics, masks, eulers, xis, yis = [], [], [], [], []
    queue = list(files) if files else None
    while True:
        if queue is not None:
            if not queue:
                break
            f = queue.pop(0)
        else:
            from tkinter import Tk, filedialog
            root = Tk()
            root.withdraw()
            f = filedialog.askopenfilename(title="Select a picture (Cancel = finish)")
            root.destroy()
            if not f:
                break
        n = len(pics) + 1
        while True:
            I = load_grey(f)
            label = f"picture {n}: {Path(f).name}"
            if not no_crop:
                I, _ = gui_crop(I, label)
            m, xs, ys = gui_mask(I, label)
            ans = gui_confirm(I, m, label)
            if ans != "redo":
                break
        if ans == "skip":
            print(f"skipped {f}")
            continue
        e = euler_for(euler_table, n, f) or ask_euler(label)
        pics.append(apply_mask(I, m))
        masks.append(m)
        eulers.append(e)
        xis.append(xs)
        yis.append(ys)
        print(f"  {label}: {I.shape[1]} x {I.shape[0]} px, {int(m.sum())} masked px, Euler {e}")
    return pics, masks, eulers, xis, yis


def prepare_from_files(files, mask_dir, euler_table):
    from PIL import Image
    pics, masks, eulers, xis, yis = [], [], [], [], []
    for n, f in enumerate(files, start=1):
        I = load_grey(f)
        m = np.zeros(I.shape, bool)
        if mask_dir:
            mp = Path(mask_dir) / f"{Path(f).stem}_mask.png"
            if mp.exists():
                ma = np.asarray(Image.open(io_mat.long_path(mp)).convert("L")) > 127
                if ma.shape != I.shape:
                    raise SystemExit(f"{mp}: mask is {ma.shape[::-1]} px, picture is {I.shape[::-1]} px")
                m = ma
            else:
                print(f"  note: no mask {mp.name} for {Path(f).name}; nothing masked")
        e = euler_for(euler_table, n, f)
        if e is None:
            raise SystemExit(f"no Euler angles for {Path(f).name} (number {n}) in the --eulers file")
        pics.append(apply_mask(I, m))
        masks.append(m)
        eulers.append(e)
        xis.append(np.zeros(0))
        yis.append(np.zeros(0))
        print(f"  picture {n}: {Path(f).name}, {int(m.sum())} masked px, Euler {e}")
    return pics, masks, eulers, xis, yis


def edit_dataset(ds, select, euler_table, interactive):
    for i in select:
        if not 1 <= i <= ds["no"]:
            raise SystemExit(f"--select {i}: the dataset has pictures 1..{ds['no']}")
        k = i - 1
        e = euler_for(euler_table, i)
        if e is not None:
            print(f"  picture {i}: Euler {ds['eulers'][k]} -> {e}")
            ds["eulers"][k] = e
        if interactive:
            m, xs, ys = gui_mask(ds["pictures"][k], f"picture {i} (existing mask shown in red)", ds["masks"][k])
            added = int((m & ~ds["masks"][k]).sum())
            ds["masks"][k] = m
            ds["pictures"][k] = apply_mask(ds["pictures"][k], m)
            if xs.size:
                ds["xi"][k] = np.concatenate([np.asarray(ds["xi"][k], float).ravel(), [np.nan], xs])
                ds["yi"][k] = np.concatenate([np.asarray(ds["yi"][k], float).ravel(), [np.nan], ys])
            print(f"  picture {i}: {added} px added to the mask")
    return ds


def main(argv=None):
    ap = argparse.ArgumentParser(prog="rast prepare", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", help="dataset .mat to write")
    ap.add_argument("--pictures", nargs="+", help="picture files, in grain order (default: choose them one by one)")
    ap.add_argument("--eulers", help="CSV with columns picture,phi1,Phi,phi2")
    ap.add_argument("--masks", help="folder of <picture stem>_mask.png files (no window is opened)")
    ap.add_argument("--no-crop", action="store_true", help="skip the cropping step")
    ap.add_argument("--edit", help="existing dataset .mat to change")
    ap.add_argument("--select", type=int, nargs="+", help="with --edit: picture numbers to change (default all)")
    ap.add_argument("--no-window", action="store_true", help="with --edit: only change Euler angles (no mask window)")
    ap.add_argument("--overwrite", action="store_true", help="allow replacing an existing --out file")
    a = ap.parse_args(argv)

    table = read_euler_csv(a.eulers) if a.eulers else None
    if a.edit:
        out = Path(a.out or a.edit)
        if out.exists() and out.resolve() != Path(a.edit).resolve() and not a.overwrite:
            raise SystemExit(f"{out} exists; pass --overwrite to replace it")
        if out.resolve() == Path(a.edit).resolve() and not a.overwrite:
            raise SystemExit("writing back into the edited file needs --overwrite (or give a new --out)")
        ds = io_mat.read_dataset(a.edit)
        sel = a.select or list(range(1, ds["no"] + 1))
        ds = edit_dataset(ds, sel, table, interactive=not a.no_window)
        io_mat.write_dataset(out, ds["pictures"], ds["masks"], ds["eulers"], ds["xi"], ds["yi"])
        print(f"written {out} ({ds['no']} pictures)")
        return 0

    if not a.out:
        ap.error("--out is required")
    out = Path(a.out)
    if out.suffix.lower() != ".mat":
        out = out.with_suffix(".mat")
    if out.exists() and not a.overwrite:
        raise SystemExit(f"{out} exists; choose another --out or pass --overwrite")
    if a.masks:
        if not a.pictures or table is None:
            ap.error("--masks needs --pictures and --eulers")
        pics, masks, eulers, xis, yis = prepare_from_files(a.pictures, a.masks, table)
    else:
        pics, masks, eulers, xis, yis = prepare_interactive(a.pictures, table, a.no_crop)
    if not pics:
        raise SystemExit("no picture entered; nothing written")
    io_mat.write_dataset(out, pics, masks, eulers, xis, yis)
    print(f"written {out} ({len(pics)} pictures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
