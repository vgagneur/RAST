"""MATLAB imadjust (default stretchlim) and edge(I, 'Canny') for uint8 images.

Follows the Image Processing Toolbox algorithms:
  imadjust(J): stretchlim with tol [0.01 0.99] on the 256-bin histogram, linear map (gamma 1) through a
      256-entry LUT built from linspace(0,1,256), converted with im2uint8 (round half up, saturate).
  edge(J,'Canny'): input converted with im2single (or im2double, see PRECISION), sigma sqrt(2),
      separable Gaussian + derivative-of-Gaussian filters (imfilter 'conv','replicate'), gradient
      magnitude hypot(dx,dy) normalised by its max, thresholds from imhist(mag, 64) with
      PercentOfPixelsNotEdges 0.7 and ThresholdRatio 0.4, directional non-maximum suppression with
      linear interpolation (cannyFindLocalMaxima), hysteresis by 8-connected selection of the weak
      map from the strong pixels (bwselect).
"""
import math

import numpy as np
from scipy import ndimage

PRECISION = np.float32      # MATLAB edge() runs on im2single for integer input


# ------------------------------------------------------------------ imadjust
def stretchlim_uint8(img, tol_low=0.01, tol_high=0.99):
    nbins = int(np.iinfo(img.dtype).max) + 1               # 256 (uint8) or 65536 (uint16)
    N = np.bincount(img.ravel(), minlength=nbins).astype(np.float64)
    cdf = np.cumsum(N) / np.sum(N)
    ilow = int(np.argmax(cdf > tol_low)) + 1        # find(cdf > tol_low, 1, 'first'), 1-based
    ihigh = int(np.argmax(cdf >= tol_high)) + 1
    if ilow == ihigh:
        lo, hi = 1, nbins
    else:
        lo, hi = ilow, ihigh
    return (lo - 1) / (nbins - 1), (hi - 1) / (nbins - 1)


def _linspace01_256(n=256):
    """MATLAB linspace(0,1,n)."""
    d1, d2 = 0.0, 1.0
    c = (d2 - d1) * (n - 1)
    k = np.arange(n, dtype=np.float64)
    # MATLAB: y = d1 + (0:n-1).*(d2-d1)/(n-1) ; endpoints forced
    y = d1 + k * (d2 - d1) / (n - 1)
    y[0], y[-1] = d1, d2
    return y


def imadjust_uint8(img):
    lo, hi = stretchlim_uint8(img)
    top = int(np.iinfo(img.dtype).max)
    lut = _linspace01_256(top + 1)
    lut = np.maximum(lo, np.minimum(hi, lut))
    lut = (lut - lo) / (hi - lo)                  # gamma 1
    lut = lut * (1.0 - 0.0) + 0.0
    lut8 = np.clip(np.floor(lut * top + 0.5), 0, top).astype(img.dtype)
    return lut8[img]


# ------------------------------------------------------------------ edge(..., 'Canny')
def _kernels(sigma=math.sqrt(2)):
    ext = math.ceil(4 * sigma)
    x = np.arange(-ext, ext + 1, dtype=np.float64)
    c = 1.0 / (math.sqrt(2 * math.pi) * sigma)
    g = c * np.exp(-(x ** 2) / (2 * sigma ** 2))
    g = g / np.sum(g)
    dg = np.gradient(g)                           # MATLAB gradient: central differences, one-sided at ends
    pos, neg = dg > 0, dg < 0
    dg[pos] = dg[pos] / np.sum(dg[pos])
    dg[neg] = dg[neg] / abs(np.sum(dg[neg]))
    return g, dg


def _conv1d_replicate(a, k, axis, dt):
    """imfilter(a, k, 'conv', 'replicate') along one axis (k odd length), accumulated in `dt`
    tap by tap in forward order. scipy.ndimage accumulates in double even for single input,
    which disagreed with MATLAB on 42 edge pixels of TiFe picture 1; this order leaves 6 (A4 notes)."""
    k = k.astype(dt)
    r = len(k) // 2
    pad = [(0, 0), (0, 0)]
    pad[axis] = (r, r)
    p = np.pad(a, pad, mode="edge")
    n = a.shape[axis]
    out = np.zeros(a.shape, dtype=dt)
    for t in range(len(k)):
        s = 2 * r - t                             # out[i] = sum_t k[t] * a[i + r - t]
        sl = [slice(None), slice(None)]
        sl[axis] = slice(s, s + n)
        out = (out + k[t] * p[tuple(sl)]).astype(dt)
    return out


def smooth_gradient(I, precision=PRECISION):
    g, dg = _kernels()
    # MATLAB passes the double kernels to imfilter on a single image -> filtered in single
    gx = _conv1d_replicate(I, g, 0, precision)    # gaussKernel' (column) -> along rows (axis 0)
    gx = _conv1d_replicate(gx, dg, 1, precision)  # derivGaussKernel (row) -> along columns
    gy = _conv1d_replicate(I, g, 1, precision)
    gy = _conv1d_replicate(gy, dg, 0, precision)
    return gx, gy


def _local_maxima(direction, ix, iy, mag):
    """cannyFindLocalMaxima on column-major linear indices (as MATLAB)."""
    m, n = mag.shape
    ixf = ix.ravel(order="F")
    iyf = iy.ravel(order="F")
    magf = mag.ravel(order="F")
    if direction == 1:
        sel = ((iyf <= 0) & (ixf > -iyf)) | ((iyf >= 0) & (ixf < -iyf))
    elif direction == 2:
        sel = ((ixf > 0) & (-iyf >= ixf)) | ((ixf < 0) & (-iyf <= ixf))
    elif direction == 3:
        sel = ((ixf <= 0) & (ixf > iyf)) | ((ixf >= 0) & (ixf < iyf))
    else:
        sel = ((iyf < 0) & (ixf <= iyf)) | ((iyf > 0) & (ixf >= iyf))
    idx = np.nonzero(sel)[0] + 1                   # 1-based linear indices
    if idx.size:
        v = np.mod(idx, m)
        ext = (v == 1) | (v == 0) | (idx <= m) | (idx > (n - 1) * m)
        idx = idx[~ext]
    if idx.size == 0:
        return idx
    ixv, iyv = ixf[idx - 1], iyf[idx - 1]
    gm = magf[idx - 1]
    M = lambda k: magf[k - 1]
    if direction == 1:
        d = np.abs(iyv / ixv)
        g1 = M(idx + m) * (1 - d) + M(idx + m - 1) * d
        g2 = M(idx - m) * (1 - d) + M(idx - m + 1) * d
    elif direction == 2:
        d = np.abs(ixv / iyv)
        g1 = M(idx - 1) * (1 - d) + M(idx + m - 1) * d
        g2 = M(idx + 1) * (1 - d) + M(idx - m + 1) * d
    elif direction == 3:
        d = np.abs(ixv / iyv)
        g1 = M(idx - 1) * (1 - d) + M(idx - m - 1) * d
        g2 = M(idx + 1) * (1 - d) + M(idx + m + 1) * d
    else:
        d = np.abs(iyv / ixv)
        g1 = M(idx - m) * (1 - d) + M(idx - m - 1) * d
        g2 = M(idx + m) * (1 - d) + M(idx + m + 1) * d
    return idx[(gm >= g1) & (gm >= g2)]


def edge_canny(img, precision=PRECISION):
    """MATLAB edge(img, 'Canny') for a uint8 image -> logical edge map."""
    a = img.astype(precision) / precision(np.iinfo(img.dtype).max)    # im2single / im2double
    dx, dy = smooth_gradient(a, precision)
    mag = np.hypot(dx, dy)
    magmax = mag.max()
    if magmax > 0:
        mag = mag / magmax
    # selectThresholds: imhist(mag, 64)
    m, n = mag.shape
    bins = np.floor(mag.astype(np.float64) * 63 + 0.5).astype(np.int64)
    counts = np.bincount(bins.ravel(), minlength=64)
    high = (int(np.argmax(np.cumsum(counts) > 0.7 * m * n)) + 1) / 64.0
    low = 0.4 * high
    # thinAndThreshold
    e = np.zeros(m * n, dtype=bool)
    magf = mag.ravel(order="F")
    strong = []
    for direction in (1, 2, 3, 4):
        lm = _local_maxima(direction, dx, dy, mag)
        weak = lm[magf[lm - 1] > low]
        e[weak - 1] = True
        strong.append(weak[magf[weak - 1] > high])
    E = e.reshape((m, n), order="F")
    strong = np.concatenate(strong) if strong else np.array([], int)
    if strong.size == 0:
        return np.zeros((m, n), dtype=bool)
    lab, _ = ndimage.label(E, structure=np.ones((3, 3), bool))
    keep = np.unique(lab.ravel(order="F")[strong - 1])
    keep = keep[keep > 0]
    return np.isin(lab, keep)
