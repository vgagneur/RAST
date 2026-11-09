# Validation of RAST-python 1.1

RAST-python follows the MATLAB code step by step, including MATLAB's own image functions (wiener2, imadjust, Canny
edge detection, Hough transform), so that it gives the same lines, scores and match counts as MATLAB. This page
lists what has been checked so far, on test datasets of B2 Ti2AlMo, B2 TiFe, bcc Fe and fcc Ni (and others),
with N = 100 px, IR = 10°.

**Status:** steps 1 and 2 are done. Step 3 (direct comparison with MATLAB 1.1 output) is pending; until then this
version is a preview.

## 1. Image processing and line detection vs MATLAB

These parts are identical in RAST 1.0 and 1.1. They were checked against outputs of the MATLAB code on the paper
datasets.

| Step | MATLAB | Result |
|---|---|---|
| Dataset reading / writing | `load`, `save` | 8 datasets round-trip with no difference in values, classes or sizes |
| Smoothing | `wiener2` ×nsmooth | 286 of 286 pictures pixel-identical (8 datasets) |
| Contrast adjustment | `imadjust` | pixel-identical (TiFe, 9 pictures) |
| Edge detection | `edge(..., 'Canny')` | 38 of 3.49 million edge pixels differ (99.999 % identical, TiFe), from MATLAB's single-precision filter rounding |
| Hough transform | `hough`, `houghpeaks`, `houghlines` | identical on synthetic test lines |
| Theoretical traces | `RAST_euler_to_slip.m` | equal to 1.7 × 10⁻¹³ relative (rounding level) |
| Whole detection and scoring | `RAST_function.m` | 99.84 % (TiFe) and 99.89 % (Ti2AlMo) of detected lines identical; the rest follow from the few edge pixels above |

## 2. The 1.1 changes vs 1.0

On B2 Ti2AlMo (12 pictures), B2 TiFe (9) and bcc Fe (17), 34,485 detected lines and 27,360 test lines in total,
RAST-python 1.1 was compared with the same pipeline using the 1.0 rules:

| Check | Result |
|---|---|
| Smoothed pictures | identical to the 1.0 pipeline (38 of 38) |
| Which families each line matches (`m`, `u`, `nomatch` counts, test lines) | **identical to 1.0 for every line** |
| Best-matching family (line colour in the figures) | differs for 157 detected lines and 144 test lines, every one a line matching a trace across ±90° (the case 1.0 ranked wrongly) |
| Line scores | 1,116 of 34,485 change, every one with a near-vertical neighbour pair across ±90° (the case 1.0 missed); lines at ST = 10: 1,952 / 1,135 / 1,645, unchanged |
| Score-threshold rule | `nlines` at every ST equals the number of lines with score ≥ ST (all pictures) |
| `rast prepare` | rebuilding the Ti2AlMo dataset from its pictures, mask PNGs and an Euler CSV gives identical pictures, masks and Euler angles; MATLAB-made datasets run unchanged |

## 3. Pending: RAST-python 1.1 vs RAST-matlab 1.1

The MATLAB 1.1 files will be run on Ti2AlMo, TiFe, Fe and fcc Ni, and the Python results compared line by line
(scores, match counts per picture and ST). The acceptance criterion is the one met in section 1: angles equal to
rounding level, ≥ 99.8 % identical lines, the rest traced to the edge-pixel rounding. This page will be updated
with the results.
