"""MATLAB hough / houghpeaks / houghlines (Image Processing Toolbox), as RAST calls them:

    [H,theta,rho] = hough(Ji, 'Theta', -90:0.1:89);
    P = houghpeaks(H, numlines, 'threshold', 0.3*max(H(:)));
    lines = houghlines(Ji, theta, rho, P, 'FillGap', 100, 'MinLength', 4);

Coordinates: x = column - 1, y = row - 1 inside the transform; houghlines returns 1-based [x y] points.
"""
import math

import numpy as np


def matlab_colon(a, d, b):
    """MATLAB a:d:b (symmetric construction, so the last element is exact)."""
    n = int(math.floor((b - a) / d + 1e-10)) + 1
    k = np.arange(n, dtype=np.float64)
    hi = a + (n - 1) * d
    out = a + k * d
    half = n // 2
    out[n - half:] = hi - (n - 1 - np.arange(n - half, n)) * d
    return out


THETA = matlab_colon(-90.0, 0.1, 89.0)                    # 1791 values


def hough(BW, theta=THETA):
    M, N = BW.shape
    D = math.sqrt((M - 1) ** 2 + (N - 1) ** 2)
    q = math.ceil(D)                                       # RhoResolution 1
    nrho = 2 * q + 1
    rho = np.linspace(-q, q, nrho)
    th = theta * math.pi / 180.0
    cos_t, sin_t = np.cos(th), np.sin(th)
    rows, cols = np.nonzero(BW.T)                           # column-major pixel order
    x = rows.astype(np.float64)                             # (transposed) -> x = col-1
    y = cols.astype(np.float64)                             # y = row-1
    H = np.zeros((nrho, len(theta)), dtype=np.float64)
    if x.size:
        first = -float(q)
        slope = (nrho - 1) / (q - (-q))
        r = np.outer(x, cos_t) + np.outer(y, sin_t)
        idx = np.floor(slope * (r - first) + 0.5).astype(np.int64)     # houghmex: (int)(slope*(rho-first)+0.5)
        nt = len(theta)
        flat = (idx * nt + np.arange(nt)).ravel()          # integer vote counts: order-free, exact
        H = np.bincount(flat, minlength=nrho * nt).reshape(nrho, nt).astype(np.float64)
    return H, theta, rho


def houghpeaks(H, numpeaks, threshold):
    """Classic houghpeaks with the default NHoodSize and the theta = +/-90 antisymmetric wrap."""
    nh = np.array(H.shape, dtype=float) / 50.0
    nhood = np.maximum(2 * np.ceil(nh / 2.0) + 1, 1).astype(int)
    center = (nhood - 1) // 2
    nr, nc = H.shape
    hT = np.ascontiguousarray(H.T)                          # row-major over H.T = column-major over H
    peaks = []
    while True:
        mi = int(np.argmax(hT))                             # first max, column-major
        p, q = mi % nr + 1, mi // nr + 1                   # 1-based
        if hT[q - 1, p - 1] >= threshold:
            peaks.append((p, q))
            p1, p2 = p - center[0], p + center[0]
            q1, q2 = q - center[1], q + center[1]
            qq, pp = np.meshgrid(np.arange(q1, q2 + 1), np.arange(max(p1, 1), min(p2, nr) + 1))
            pp = pp.ravel(order="F").copy()
            qq = qq.ravel(order="F").copy()
            low = qq < 1
            qq[low] = nc + qq[low]
            pp[low] = nr - pp[low] + 1
            high = qq > nc
            qq[high] = qq[high] - nc
            pp[high] = nr - pp[high] + 1
            hT[qq - 1, pp - 1] = 0
            if len(peaks) == numpeaks:
                break
        else:
            break
    return peaks


def _round_half_away(v):
    return np.where(v >= 0, np.floor(v + 0.5), -np.floor(-v + 0.5))


def houghlines(BW, theta, rho, peaks, fillgap=100, minlength=4):
    """Returns a list of dict(point1=(x,y), point2=(x,y), theta, rho), 1-based points."""
    cols, rows = np.nonzero(BW.T)                           # [y, x] = find(BW) (column-major)
    x = cols.astype(np.float64)                             # careful: BW.T -> first index is column
    y = rows.astype(np.float64)
    # nonzeropix = [x-1, y-1]
    xs, ys = x, y                                           # already 0-based (col-1, row-1)
    nrho = len(rho)
    slope = (nrho - 1) / (rho[-1] - rho[0])
    lines = []
    fill_sq = fillgap ** 2
    min_sq = minlength ** 2
    for (pr, pt) in peaks:
        th = theta[pt - 1] * math.pi / 180.0
        rho_xy = xs * math.cos(th) + ys * math.sin(th)
        bins = _round_half_away(slope * (rho_xy - rho[0])) + 1
        sel = bins == pr
        if not sel.any():
            continue
        r = ys[sel] + 1
        c = xs[sel] + 1
        # reSortHoughPixels
        if (r.max() - r.min()) > (c.max() - c.min()):
            order = np.lexsort((c, r))                       # sort on r, then c
        else:
            order = np.lexsort((r, c))                       # sort on c, then r
        r, c = r[order], c[order]
        xy = np.stack([c, r], axis=1)
        dist_sq = np.sum(np.diff(xy, axis=0) ** 2, axis=1)
        idx = np.concatenate([[0], np.nonzero(dist_sq > fill_sq)[0] + 1, [xy.shape[0]]])
        for k in range(len(idx) - 1):
            p1 = xy[idx[k]]
            p2 = xy[idx[k + 1] - 1]
            if np.sum((p2 - p1) ** 2) >= min_sq:
                lines.append(dict(point1=(float(p1[0]), float(p1[1])), point2=(float(p2[0]), float(p2[1])),
                                  theta=float(theta[pt - 1]), rho=float(rho[pr - 1])))
    return lines


def detect_lines(Ji, numlines):
    """RAST's per-segment call sequence."""
    H, theta, rho = hough(Ji)
    P = houghpeaks(H, numlines, 0.3 * H.max())
    return houghlines(Ji, theta, rho, P, 100, 4)
