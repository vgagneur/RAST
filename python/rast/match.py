"""Theoretical slip traces and line matching (RAST v1.1).

Ports of:
    RAST_euler_to_slip.m              -> euler_to_slip
    RAST_line_processing.m l.25-50    -> theoretical_traces (normalisation, slopes CO)
    RAST_slip_trace_analyser_single.m -> analyse_line

v1.1 change (same as the MATLAB 1.1 files): the angle between a line and a trace is measured on the
180 degree circle (0 to 90 deg). Which traces match a line (angle < IR) is decided by v1.0's own test, so
the match counts are those of v1.0; the angle used to rank the best / second match is now correct for
near-vertical lines (v1.0 gave ~170 deg instead of ~10 deg when the match was across +/-90 deg), which
only changes the best family (the line colour in the figures).

The rest follows the MATLAB code exactly, including its details:
  * the best / second logic finds a plane by FLOAT EQUALITY of its slope in the list of matched slopes,
    taking the LAST equal entry, and keeps `hasmatched` set once set within a family;
  * the index `x` survives across planes and families within one call;
  * the line angle goes through tand() then atand() before the difference;
  * theoretical traces are not rotated by degrot (only the detected lines carry degrot).
MATLAB's tand / atand / mod / round are reproduced so that angles agree to rounding level.
"""
import math

import numpy as np

# ------------------------------------------------------------------ MATLAB elementary functions
_R2D = 180.0 / math.pi
_D2R = math.pi / 180.0


def tand(x):
    """MATLAB tand (elfun/tand.m): exact at multiples of 45/90, reduction by 90-degree multiples."""
    x = float(x)
    if not math.isfinite(x):
        return math.nan
    n = round_half_away(x / 90.0)
    r = x - n * 90.0
    y = math.tan(_D2R * r)
    if int(n) % 2 != 0:                  # odd multiple of 90: tan -> -cot
        y = -math.copysign(math.inf, y) if y == 0.0 else -1.0 / y
    return y


def round_half_away(v):
    """MATLAB round: half away from zero."""
    return math.floor(v + 0.5) if v >= 0 else -math.floor(-v + 0.5)


def atand(x):
    """MATLAB atand: 180/pi * atan(x)."""
    return _R2D * math.atan(x)


def mmod(a, m):
    """MATLAB mod(a, m) for scalar doubles, m > 0: a - floor(a/m)*m, with MATLAB's
    round-off compensation (exact 0 when a/m is within eps of an integer)."""
    if m == 0:
        return a
    if not math.isfinite(a):
        return math.nan
    return a - math.floor(a / m) * m


# ------------------------------------------------------------------ RAST_euler_to_slip.m
P110 = [[1, 1, 0], [-1, 1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, -1, 1]]
P112 = [[1, 1, 2], [-1, 1, 2], [1, -1, 2], [1, 1, -2], [1, 2, 1], [-1, 2, 1], [1, -2, 1], [1, 2, -1],
        [2, 1, 1], [-2, 1, 1], [2, -1, 1], [2, 1, -1]]
P123 = [[1, 2, 3], [-1, 2, 3], [1, -2, 3], [1, 2, -3], [1, 3, 2], [-1, 3, 2], [1, -3, 2], [1, 3, -2],
        [2, 1, 3], [-2, 1, 3], [2, -1, 3], [2, 1, -3], [2, 3, 1], [-2, 3, 1], [2, -3, 1], [2, 3, -1],
        [3, 2, 1], [-3, 2, 1], [3, -2, 1], [3, 2, -1], [3, 1, 2], [-3, 1, 2], [3, -1, 2], [3, 1, -2]]
P100 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def planes_for(crystal):
    """Plane lists per family slot. crystal 1 = bcc/B2 ({100}, {110}, {112}, {123});
    crystal 2 = fcc, one {111} plane per slot (each entered twice, as in RAST_euler_to_slip.m)."""
    if crystal == 1:
        return [P100, P110, P112, P123]
    if crystal == 2:
        return [[[1, 1, 1]] * 2, [[-1, 1, 1]] * 2, [[1, -1, 1]] * 2, [[1, 1, -1]] * 2]
    raise ValueError(f"crystal must be 1 (bcc/B2) or 2 (fcc), not {crystal}")


def _rot_x(t):
    t = math.radians(t)  # MATLAB deg2rad: (pi/180)*t
    return np.array([[1, 0, 0], [0, math.cos(t), math.sin(t)], [0, -math.sin(t), math.cos(t)]])


def _rot_z(t):
    t = math.radians(t)
    return np.array([[math.cos(t), math.sin(t), 0], [-math.sin(t), math.cos(t), 0], [0, 0, 1]])


def parse_euler(s):
    """MATLAB str2num('[a,b,c]') -> 3 doubles."""
    return [float(v) for v in s.strip().strip("[]").replace(" ", ",").split(",") if v != ""]


def _mm(A, B):
    """Matrix product with plain sequential sums (closest to MATLAB's results in the A2 gate)."""
    A, B = np.asarray(A, float), np.asarray(B, float)
    out = np.zeros((A.shape[0], B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            s = 0.0
            for k in range(A.shape[1]):
                s += A[i, k] * B[k, j]
            out[i, j] = s
    return out


def euler_to_slip(euler_str, crystal):
    """Return [s100, s110, s112, s123], each (n_planes x 2) slip-trace direction (x, y)."""
    e = parse_euler(euler_str)
    om = _mm(_mm(_rot_z(e[2]), _rot_x(e[1])), _rot_z(e[0]))     # rotate_zxz
    out = []
    for fam in planes_for(crystal):
        rows = []
        for pl in fam:
            v = _mm(np.array([pl], dtype=float), om)[0]          # v * om (row vector)
            # cross(v, [0 0 1]) = [v2*1 - v3*0, v3*0 - v1*1, ...]; z column dropped
            rows.append([v[1] * 1.0 - v[2] * 0.0, v[2] * 0.0 - v[0] * 1.0])
        out.append(np.array(rows))
    return out


# ------------------------------------------------------------------ RAST_line_processing.m l.25-50
def theoretical_traces(euler_str, crystal):
    """Normalised trace directions `slips` and slopes `CO` (lists of 4 arrays, one per family slot)."""
    slips = euler_to_slip(euler_str, crystal)
    with np.errstate(divide="ignore", invalid="ignore"):
        # slipette(i,:)/norm(slipette(i,:)); hypot agrees with MATLAB's norm to rounding level
        slips = [s / np.array([[math.hypot(r[0], r[1])] for r in s]) for s in slips]
        CO = [s[:, 1] / (-s[:, 0]) for s in slips]          # slope = y / (-x)
    return slips, CO


# ------------------------------------------------------------------ RAST_slip_trace_analyser_single.m
def analyse_line(CO, line_deg, error):
    """Match one line (orientation in degrees) against the theoretical slopes CO with IR = error (deg).

    Returns best, second, matches: lists of 4 ints (one per plane family slot).
    """
    lineslope = tand(line_deg)
    linedeg = atand(lineslope)
    matched_slopes = []          # per family: slopes of the planes within IR (MATLAB sloped / slipXXXco)
    matched_degs = []            # per family: their angle to the line (MATLAB degoff / degoffXXX)
    for fam in CO:
        sl_f, de_f = [], []
        for sl in fam:
            sl = float(sl)
            refdeg = atand(sl)
            # matching test unchanged from v1.0: |line - trace| < IR, also across +/-90 deg
            if ((refdeg - error < linedeg < refdeg + error) or (refdeg - error + 180 < linedeg < refdeg + error + 180)
                    or (refdeg - error - 180 < linedeg < refdeg + error - 180)):
                deg = refdeg - linedeg                 # in (-180, 180)
                if deg >= 90.0:                        # v1.1: wrap across +/-90 deg; other angles are
                    deg -= 180.0                       # exactly those of v1.0
                elif deg < -90.0:
                    deg += 180.0
                sl_f.append(sl)
                de_f.append(abs(deg))
        matched_slopes.append(sl_f)
        matched_degs.append(de_f)

    co = [c for f in matched_slopes for c in f]
    de = [d for f in matched_degs for d in f]
    min_de = min(de) if de else None
    best = [0, 0, 0, 0]
    second = [0, 0, 0, 0]
    x = None
    for i, fam in enumerate(CO):
        min_sys = min(matched_degs[i]) if matched_degs[i] else None
        hasmatched = 0
        bestmatch = 0
        matchissime = 0
        for yj in fam:
            yj = float(yj)
            if any(yj == c for c in co):            # ismember(Y(j), co); NaN is never a member
                hasmatched = 1
            if hasmatched == 1:
                for z, c in enumerate(co):          # the LAST equal entry wins
                    if yj == c:
                        x = z
                if min_sys is not None and de[x] == min_sys:
                    bestmatch = 1
                    if de[x] == min_de:
                        matchissime = 1
        if matchissime == 1:
            best[i] = 1
        elif bestmatch == 1:
            second[i] = 1
    matches = [b + s for b, s in zip(best, second)]
    return best, second, matches
