"""Line scoring, RAST_function.m step (6) (RAST v1.1).

For every cell (i, j) and slot lno, the score accumulates over neighbour orders o = 1..10, then k = -o..o,
l = -o..o (ring cells only: k or l = +/-o), then neighbour slots lno2 = 1..numlines:
    +1/o^2   if d < nei_deg                                    (similar orientation)
    -10/o    elif |co1| ~ 45 and d ~ 90                         (45-degree noise cross-hatch)
    -10/o    elif co1 ~ 0 or ~ +/-90, and |co2| ~ 45            (noise)
where ~ means within nei_deg, co = stored orientation - degrot, and d is the difference between co1 and
co2 on the 180 degree circle (0 to 90 deg) in the similarity test, and |co1 - co2| in the "~ 90" noise test.

v1.1 change (same as the MATLAB 1.1 files): the similarity test uses the difference on the 180 degree circle,
so lines at +89 and -89 deg count as similar. v1.0 used |co1 - co2| there too.

Details kept from the MATLAB code:
  * a neighbour offset (k, l) is used only if BOTH (i+k, j+l) and (i-k, j-l) are inside the grid;
  * a slot is "empty" when its Ax is 0, for the current line and for the neighbour;
  * the `(-45-nei_deg) < |co2| < (-45+nei_deg)` branch can never be true; kept as written;
  * scores are summed in exactly MATLAB's loop order (vectorised over cells only), because they are
    compared with integer score thresholds afterwards.
"""
import numpy as np


def score_lines(coefmat, numlines, nei_deg, degrot, max_order=10):
    cm = coefmat.copy()
    cuth, cutl, _ = cm.shape
    I, J = np.meshgrid(np.arange(1, cuth + 1), np.arange(1, cutl + 1), indexing="ij")
    co = cm[:, :, :numlines] - degrot                           # co_1b - degrot
    A = cm[:, :, [numlines + 6 * s for s in range(numlines)]]   # Ax per slot
    score = np.zeros((cuth, cutl, numlines))
    nd = nei_deg
    for o in range(1, max_order + 1):
        inc_sim = 1.0 / (o ** 2)
        inc_noise = -10.0 / o
        for k in range(-o, o + 1):
            for l in range(-o, o + 1):
                if not (k == o or l == o or k == -o or l == -o):
                    continue
                valid = ((I + k >= 1) & (J + l >= 1) & (I + k <= cuth) & (J + l <= cutl) &
                         (I - k >= 1) & (J - l >= 1) & (I - k <= cuth) & (J - l <= cutl))
                if not valid.any():
                    continue
                # neighbour arrays (values only used where valid)
                ii = np.clip(I + k, 1, cuth) - 1
                jj = np.clip(J + l, 1, cutl) - 1
                co2_all = co[ii, jj, :]
                A2_all = A[ii, jj, :]
                for lno in range(numlines):
                    co1 = co[:, :, lno]
                    A1 = A[:, :, lno]
                    for lno2 in range(numlines):
                        co2 = co2_all[:, :, lno2]
                        A2 = A2_all[:, :, lno2]
                        live = valid & (A1 != 0) & (A2 != 0)
                        d = np.abs(co1 - co2)
                        dw = np.where(d > 90.0, 180.0 - d, d)      # v1.1: similarity on the 180 deg circle
                        c_sim = dw < nd
                        c_n1 = (((45 - nd) < np.abs(co1)) & (np.abs(co1) < (45 + nd)) &
                                ((90 - nd) < d) & (d < (90 + nd)))
                        c_n2 = ((((0 - nd) < co1) & (co1 < (0 + nd))) | ((90 - nd) < co1) | (co1 < (-90 + nd))) & \
                               ((((45 - nd) < np.abs(co2)) & (np.abs(co2) < (45 + nd))) |
                                (((-45 - nd) < np.abs(co2)) & (np.abs(co2) < (-45 + nd))))
                        inc = np.where(c_sim, inc_sim, np.where(c_n1 | c_n2, inc_noise, 0.0))
                        inc = np.where(live, inc, 0.0)
                        hit = inc != 0
                        if hit.any():
                            s = score[:, :, lno]
                            s[hit] = s[hit] + inc[hit]          # sequential, in MATLAB's order
    for s in range(numlines):
        cm[:, :, numlines + 4 + 6 * s] = score[:, :, s]
    return cm
