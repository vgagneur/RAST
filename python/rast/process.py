"""RAST_line_processing.m: match every detected line to the theoretical traces and count per score threshold.

For one picture, returns:
    sol[st] for st = 0..20: the 82 RAST values for the lines with score >= st (same order as MATLAB's y:
        20 counts, 20 test-line counts, 20 percentages, 20 test-line percentages, nlines, tnlines);
    code[i, j, slot]: the colour code of each line for the figures (100 / 110 / 112 / 123 = best-matching
        family slot, 111 = special case, 0 = no match), as MATLAB writes it into coefmat.

Score-threshold rule (RAST v1.1, every output): at ST = st, a line is counted when its score is >= st.
Lines in a segment whose first slot is empty (masked area) are skipped, as in MATLAB.
"""
import math

import numpy as np

from . import match

N_ST = 21          # score thresholds 0..20
TEST_ANGLES = np.arange(1, 721) * 0.25     # generated test lines 0.25:0.25:180 deg

COUNT_NAMES = ["m100", "m110", "m112", "m123", "u100", "u110", "u112", "u123", "u111",
               "m100110", "m100112", "m100123", "m110112", "m110123", "m112123",
               "m110112123", "m100112123", "m100110123", "m100110112", "nomatch"]
SOL_NAMES = (COUNT_NAMES + ["t" + n for n in COUNT_NAMES] + ["p" + n for n in COUNT_NAMES] +
             ["pt" + n for n in COUNT_NAMES] + ["nlines", "tnlines"])
FAMILY_CODES = [100, 110, 112, 123]


def _counts(match_rows, crystal):
    """The 20 RAST counts (COUNT_NAMES order) for a list of 4-entry `matches` rows."""
    m = [0] * 4
    u = [0] * 4
    u111 = 0
    dbl = [0] * 6      # m100110, m100112, m100123, m110112, m110123, m112123
    tri = [0] * 4      # m110112123, m100112123, m100110123, m100110112
    nomatch = 0
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for mt in match_rows:
        s = sum(mt)
        for f in range(4):
            m[f] += mt[f] == 1
        if s == 0:
            nomatch += 1
        if crystal == 1 and s > 0 and not (mt[0] == 1 or mt[1] == 1):
            u111 += 1
        if s == 1:
            for f in range(4):
                u[f] += mt[f] == 1
        if s == 2:
            for q, (a, b) in enumerate(pairs):
                dbl[q] += mt[a] == 1 and mt[b] == 1
        if s == 3:
            for f in range(4):
                tri[f] += mt[f] == 0
    return [int(v) for v in m + u + [u111] + dbl + tri + [nomatch]]


def test_line_counts(CO, error):
    """Counts for the 720 generated test lines (identical for every score threshold).
    Test lines count u111 whatever the crystal mode, as in MATLAB."""
    rows = [match.analyse_line(CO, float(t), error)[2] for t in TEST_ANGLES]
    return _counts(rows, crystal=1), len(rows)


def line_code(best, matches, crystal):
    """Figure colour code of one line (RAST_line_processing.m l.84-97 and the 111 overrides)."""
    code = 0
    for f in range(4):
        if best[f] == 1:
            code = FAMILY_CODES[f]
            break
    s = sum(matches)
    if crystal == 1 and s > 0 and not (matches[0] == 1 or matches[1] == 1):
        code = 111          # bcc: unambiguously <111> slip (matches only {112} and/or {123})
    elif crystal == 2 and s > 1:
        code = 111          # fcc: matches more than one {111} plane (ambiguous)
    return code


def process_picture(euler_str, coefmat, numlines, error, crystal):
    _, CO = match.theoretical_traces(euler_str, crystal)
    cuth, cutl, _ = coefmat.shape
    code = np.zeros((cuth, cutl, numlines), dtype=int)
    scored = []                                        # (score, matches) per analysed line
    for i in range(cuth):
        for j in range(cutl):
            if coefmat[i, j, numlines] == 0:          # Ax of slot 1 == 0: masked / empty segment
                continue
            for ln in range(numlines):
                o = numlines + 6 * ln
                best, _, matches = match.analyse_line(CO, float(coefmat[i, j, ln]), error)
                sc = coefmat[i, j, o + 4]
                scored.append((sc, matches))
                # MATLAB applies the 111 overrides only to lines with score >= 0 (they are set inside the
                # ST loop, which starts at 0); figures are only drawn for ST >= 0, so this is equivalent.
                code[i, j, ln] = line_code(best, matches, crystal) if sc >= 0 else line_code(best, [0] * 4, 0)
    tcounts, tn = test_line_counts(CO, error)
    sol = []
    for st in range(N_ST):
        rows = [mt for sc, mt in scored if sc >= st]
        c = _counts(rows, crystal)
        n = len(rows)
        with np.errstate(divide="ignore", invalid="ignore"):
            pc = [np.float64(v) / n * 100 if n else math.nan for v in c]
        pt = [v / tn * 100 for v in tcounts]
        sol.append(np.array(c + tcounts + pc + pt + [n, tn], dtype=float))
    return dict(sol=sol, code=code)
