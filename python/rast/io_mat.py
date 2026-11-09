"""Reading and writing RAST datasets (the .mat file made by RAST_Loading.m or `rast prepare`).

A dataset holds, for no pictures (MATLAB names):
    no        number of pictures
    Pictures  {1 x no} grey-level pictures (uint8, or uint16), cropped; masked pixels have the value 1
    BWss      {1 x no} logical masks (true = masked: indent, scratches, other grains ...)
    Eulers    {1 x no} each a {1 x 1} cell holding the Euler angles as text, e.g. '[345.057,160.014,265.877]'
              (Bunge convention, degrees)
    xi2, yi2  {1 x no} x and y coordinates of the mask polygon vertices (information only)

The same file works in MATLAB (RAST_Running.m) and in Python (`rast run`).
"""
import os
from pathlib import Path

import numpy as np
import scipy.io as sio


def long_path(p):
    """Path string that also works beyond Windows' 260-character limit."""
    s = str(Path(p).resolve())
    if os.name == "nt" and len(s) > 240 and not s.startswith("\\\\?\\"):
        s = "\\\\?\\" + s
    return s


def _cell_row(items):
    c = np.empty((1, len(items)), dtype=object)
    for i, x in enumerate(items):
        c[0, i] = x
    return c


def euler_text(x):
    """Euler entry of a dataset (nested MATLAB cell or text) -> text such as '[345.057,160.014,265.877]'."""
    while isinstance(x, np.ndarray) and x.dtype == object:
        x = x.ravel()[0]
    if isinstance(x, np.ndarray):
        x = x.ravel()[0] if x.size == 1 else "".join(x.ravel())
    return str(x).strip()


def parse_euler(text):
    """'[a,b,c]' or 'a b c' -> [a, b, c] (floats)."""
    vals = [float(v) for v in str(text).strip().strip("[]").replace(",", " ").split()]
    if len(vals) != 3:
        raise ValueError(f"Euler angles need 3 values (phi1, Phi, phi2), got {text!r}")
    return vals


def format_euler(vals):
    """[a, b, c] -> '[a,b,c]' with the shortest exact text of each value (90.0 -> '90')."""
    out = []
    for v in vals:
        s = repr(float(v))
        out.append(s[:-2] if s.endswith(".0") else s)
    return "[" + ",".join(out) + "]"


def read_dataset(path):
    """Read a RAST dataset. Returns dict(pictures, masks, eulers (text), no, xi, yi)."""
    try:
        d = sio.loadmat(long_path(path))
    except NotImplementedError:
        raise ValueError(f"{path} is a MATLAB v7.3 file, which RAST's Python reader can't open. In MATLAB, re-save it "
                         f"with: load('{Path(path).name}'); save('{Path(path).stem}_v7.mat', '-v7')") from None
    no = int(np.asarray(d["no"]).ravel()[0])
    pics = [np.asarray(d["Pictures"][0, i]) for i in range(no)]
    masks = [np.asarray(d["BWss"][0, i]).astype(bool) for i in range(no)]
    eulers = [euler_text(d["Eulers"][0, i]) for i in range(no)]
    xi = [np.asarray(d["xi2"][0, i], float) if "xi2" in d else np.zeros((0, 1)) for i in range(no)]
    yi = [np.asarray(d["yi2"][0, i], float) if "yi2" in d else np.zeros((0, 1)) for i in range(no)]
    for i, (p, m) in enumerate(zip(pics, masks)):
        if p.ndim != 2 or p.dtype not in (np.uint8, np.uint16):
            raise ValueError(f"picture {i + 1}: expected a 2-D uint8 or uint16 grey picture, got {p.dtype} {p.shape}")
        if m.shape != p.shape:
            raise ValueError(f"picture {i + 1}: mask shape {m.shape} differs from picture shape {p.shape}")
        parse_euler(eulers[i])
    return dict(pictures=pics, masks=masks, eulers=eulers, no=no, xi=xi, yi=yi)


def write_dataset(path, pictures, masks, eulers, xi=None, yi=None):
    """Write a RAST dataset in the layout of RAST_Loading.m (loads in MATLAB and Python)."""
    no = len(pictures)
    xi = xi if xi is not None else [np.zeros((0, 1)) for _ in range(no)]
    yi = yi if yi is not None else [np.zeros((0, 1)) for _ in range(no)]
    md = dict(
        no=np.array([[float(no)]]),
        Pictures=_cell_row([np.asarray(p) for p in pictures]),
        Eulers=_cell_row([_cell_row([np.array([str(e)])]) for e in eulers]),     # {1x1 cell of char}
        BWss=_cell_row([np.asarray(m, dtype=bool) for m in masks]),
        xi2=_cell_row([np.asarray(x, float).reshape(-1, 1) for x in xi]),
        yi2=_cell_row([np.asarray(y, float).reshape(-1, 1) for y in yi]),
    )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    sio.savemat(long_path(path), md, do_compression=True)
