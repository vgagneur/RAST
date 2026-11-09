"""RAST_picture_segmentation.m: segment grid, per-segment imadjust + Canny + Hough, coefmat packing.

MATLAB QUIRKS kept:
  * the grid has cuth = #{i : i*pixsize < h} rows of segments (cutl likewise); the right/bottom
    remainder is never analysed;
  * a segment is analysed only if it contains NO zero pixel (`ismember(0,J)==0`): masked areas,
    the zeroed picture border and natural black pixels all skip it;
  * an edge map is computed for EVERY segment (for the figure) from the raw segment; analysed
    segments overwrite it with the edge map of the imadjusted segment;
  * the figure edge map L keeps the original picture values on the one-pixel frame of each segment,
    so N_bw.png (= L(1:cf,1:df) > 0) shows the segment grid;
  * lines are accumulated in one flat list, the first two entries are dropped (`alllines(3:end)`),
    and the list is cut into consecutive slots of `numlines` per cell with the trailing remainder dropped.
"""
import math

import numpy as np

from . import canny


def grid_size(h, l, pixsize):
    cuth = 0
    seg = 0
    for _ in range(h):
        seg += pixsize
        if seg < h:
            cuth += 1
        else:
            break
    cutl = 0
    seg = 0
    for _ in range(l):
        seg += pixsize
        if seg < l:
            cutl += 1
        else:
            break
    return cuth, cutl


def segment_picture(I, pixsize, numlines, detect_lines=None, precision=canny.PRECISION):
    """Run the segmentation loop on the smoothed, border-zeroed picture I.

    detect_lines(Ji, numlines) -> list of dict(point1=(x,y), point2=(x,y)) in Ji's 1-based coords,
    or None to skip line detection (A4: figure maps only).
    Returns dict(cuth, cutl, Li (contrast picture), BW (edge figure), alllines (raw list)).
    """
    h, l = I.shape
    cuth, cutl = grid_size(h, l, pixsize)
    L = I.copy()
    Li = L.copy()
    alllines = []
    a, b, c, d = 1, 1, pixsize, pixsize
    e = f = pixsize
    cf = df = None
    for i in range(1, cuth + 1):
        for j in range(1, cutl + 1):
            J = L[a - 1:c, b - 1:d].copy()                           # MATLAB slices copy; numpy views alias
            Ji = canny.edge_canny(J[1:-1, 1:-1], precision)          # figure edge map (raw segment)
            if not (J == 0).any():
                J = canny.imadjust_uint8(J)
                Ji = canny.edge_canny(J[1:-1, 1:-1], precision)
                if detect_lines is not None:
                    for ln in detect_lines(Ji, numlines):
                        ln = dict(ln)
                        ln["xyco"] = (float(b), float(a))
                        ln["celx"], ln["cely"] = i, j
                        alllines.append(ln)
            L[a:c - 1, b:d - 1] = Ji.astype(L.dtype)                 # L(a+1:c-1,b+1:d-1) = Ji
            Li[a - 1:c, b - 1:d] = J
            b += f
            d += f
        b = 1
        df = d - f
        d = pixsize
        cf = c
        a += e
        c += e
    BW = L[:cf, :df] > 0 if cf and df else np.zeros((0, 0), bool)
    return dict(cuth=cuth, cutl=cutl, Li=Li, BW=BW, alllines=alllines)


def pack_coefmat(alllines, cuth, cutl, numlines, degrot):
    """RAST_picture_segmentation.m l.166-226: absolute endpoints, alllines(3:end), slot packing.

    coefmat[i, j, :] (0-based cell) holds, for slot l = 1..numlines (1-based):
        [l-1]                      orientation cur_deg + degrot
        [numlines + 6*(l-1) + 0..3] Ax, Ay, Bx, By (absolute, 1-based picture coords + segment offset)
        [numlines + 6*(l-1) + 4]   score (filled by score.score_lines)
        [numlines + 6*(l-1) + 5]   best plane (figure only)
    """
    lines = []
    for ln in alllines:
        bx, ay = ln["xyco"]
        p1 = (ln["point1"][0] + bx, ln["point1"][1] + ay)
        p2 = (ln["point2"][0] + bx, ln["point2"][1] + ay)
        lines.append((p1, p2, ln["celx"], ln["cely"]))
    lines = lines[2:]                                        # MATLAB QUIRK: alllines = alllines(3:end)
    coefmat = np.zeros((cuth, cutl, 7 * numlines))
    k = 0
    n = len(lines)
    with np.errstate(divide="ignore", invalid="ignore"):
        while k <= n - numlines:                             # MATLAB QUIRK: trailing remainder dropped
            for l in range(1, numlines + 1):
                k += 1
                (Ax, Ay), (Bx, By), X, Y = lines[k - 1]
                cur_deg = _atand(np.float64(Ay - By) / np.float64(Bx - Ax))
                o = numlines + 6 * (l - 1)
                coefmat[X - 1, Y - 1, l - 1] = cur_deg + degrot
                coefmat[X - 1, Y - 1, o + 0] = Ax
                coefmat[X - 1, Y - 1, o + 1] = Ay
                coefmat[X - 1, Y - 1, o + 2] = Bx
                coefmat[X - 1, Y - 1, o + 3] = By
                coefmat[X - 1, Y - 1, o + 4] = 0.0
                coefmat[X - 1, Y - 1, o + 5] = 0.0
    return coefmat


def _atand(v):
    return (180.0 / math.pi) * math.atan(float(v))
