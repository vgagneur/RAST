"""The RAST results table (RAST_function.m step (8) and RAST_Running.m "Tablemaster").

One row per picture and score threshold ST = 0..20, in the column order of the MATLAB Excel sheet:
    Euler_1, Euler_2, Euler_3, grainID, pictureID, ST, IR, N,
    m100 ... nomatch (counts), tm100 ... tnomatch (test-line counts),
    pm100 ... pnomatch (% of nlines), ptm100 ... ptnomatch (% of tnlines), nlines, tnlines.
Column meanings: see data_analysis_tutorial.txt (m = matches the family, u = matches only that family,
u111 = bcc lines matching only {112} and/or {123}, mXXXYYY = matches exactly those families, t = test lines,
p = percentage). With crystal 2 (fcc) the four family columns are the four {111} planes, in the order
(111), (-111), (1-11), (11-1), under the names 100, 110, 112, 123.

grainID starts at 1 and goes up by one each time the Euler angles change from one picture to the next, so
pictures of the same grain must be consecutive in the dataset (as in MATLAB).
"""
import csv
import math
from pathlib import Path

from . import io_mat, process

COLUMNS = ["Euler_1", "Euler_2", "Euler_3", "grainID", "pictureID", "ST", "IR", "N"] + process.SOL_NAMES


def rows_for_run(eulers, sols, ir, n, picture_ids=None):
    """Rows for one (N, IR) run. sols[k] = list of 21 sol vectors for picture picture_ids[k] (1-based)."""
    picture_ids = picture_ids or list(range(1, len(sols) + 1))
    rows = []
    count = 1
    prev = None
    for e_text, sol, pid in zip(eulers, sols, picture_ids):
        e = io_mat.parse_euler(e_text)
        for st, y in enumerate(sol):
            if prev is not None and e != prev:
                count += 1
            prev = e
            rows.append(e + [count, pid, st, ir, n] + [_num(v) for v in y])
    return rows


def _num(v):
    v = float(v)
    if math.isnan(v):
        return math.nan
    return int(v) if v.is_integer() and abs(v) < 2 ** 53 else v


def write(rows, stem):
    """Write <stem>.xlsx and <stem>.csv. NaN (a percentage of 0 lines) is an empty cell in .xlsx, NaN in .csv."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    with open(io_mat.long_path(str(stem) + ".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(COLUMNS)
        for r in rows:
            w.writerow(["NaN" if isinstance(v, float) and math.isnan(v) else v for v in r])
    from openpyxl import Workbook       # imported here so that the rest of RAST works without it
    wb = Workbook(write_only=True)
    ws = wb.create_sheet("Sheet1")
    ws.append(COLUMNS)
    for r in rows:
        ws.append([None if isinstance(v, float) and math.isnan(v) else v for v in r])
    wb.save(io_mat.long_path(str(stem) + ".xlsx"))
