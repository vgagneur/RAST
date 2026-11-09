"""Picture smoothing (RAST_function.m l.60-108) and the border zeroing at the top of
RAST_picture_segmentation.m (l.42-45).

MATLAB QUIRKS kept:
  * the "legend" search walks UP the LAST column from the bottom and stops one row above the
    5th non-zero pixel; only rows above that are smoothed (a leftover from SEM legends);
  * wiener2 works on im2double(uint8) and converts back to uint8 after EVERY pass, so the
    rounding happens nsmooth times;
  * the mask is applied after smoothing (masked pixels set to 0), then the picture's outer
    border is set to 0 by the segmentation step.
"""
import numpy as np


def _boxsum5(a):
    """filter2(ones(5), a) with zero padding ('same'), as a separable sum."""
    H, W = a.shape
    p = np.pad(a, 2)
    r = np.zeros((H + 4, W))
    for kj in range(5):
        r = r + p[:, kj:kj + W]
    out = np.zeros((H, W))
    for ki in range(5):
        out = out + r[ki:ki + H, :]
    return out


def wiener2_5x5_uint8(g8):
    """MATLAB wiener2(g, [5 5]) for a uint8 or uint16 image (noise = mean2 of the local variance);
    im2double divides by the class maximum and the result is converted back to the input class."""
    top = float(np.iinfo(g8.dtype).max)
    g = g8.astype(np.float64) / top                         # im2double
    local_mean = _boxsum5(g) / 25.0
    local_var = _boxsum5(g * g) / 25.0 - local_mean ** 2
    noise = np.mean(local_var)                             # mean2(localVar)
    f = g - local_mean
    g2 = np.maximum(local_var - noise, 0.0)
    local_var = np.maximum(local_var, noise)
    f = f / local_var
    f = f * g2
    f = f + local_mean
    return np.clip(np.floor(f * top + 0.5), 0, top).astype(g8.dtype)   # im2uint8/im2uint16 (round, saturate)


def legend_row(I):
    """1-based j: rows 1..j are smoothed (RAST_function.m l.79-89)."""
    Ix, Iy = I.shape
    l, j = 5, Ix
    while l != 0:
        if j < 1:     # MATLAB stops with an index error here
            raise ValueError("the last pixel column of the picture has fewer than 5 non-black pixels; "
                             "RAST needs them to find the region to smooth (crop the picture differently)")
        if I[j - 1, Iy - 1] != 0:
            l -= 1
        j -= 1
    return j


def smooth_picture(I, mask, nsmooth):
    """RAST_function.m step (4): returns the smoothed, masked uint8 picture."""
    I = I.copy()
    j = legend_row(I)
    L = I[:j, :]
    for _ in range(nsmooth):
        L = wiener2_5x5_uint8(L)
    I[:j, :] = L
    I[mask] = 0
    return I


def zero_border(I):
    """RAST_picture_segmentation.m l.42-45 (the picture written to N_picture.png)."""
    I = I.copy()
    I[:, 0] = 0
    I[0, :] = 0
    I[:, -1] = 0
    I[-1, :] = 0
    return I
