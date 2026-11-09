# Parameters and conventions (agent note)

| Parameter | Meaning | Default | Notes |
|---|---|---|---|
| `N` | segment size (px): the picture is cut into N × N segments and lines are found per segment | 100 | keep above ~90 px; a list = one run per value |
| `IR` | identification range (deg): a line matches a trace if their angle is < IR | 10 | a list = one run per value |
| `nsmooth` | Wiener smoothing passes before edge detection | 10 | 15–20 for dirty or scratched pictures |
| `degrot` | anticlockwise rotation of the SEM picture relative to the EBSD frame (deg) | 0 | added to every detected line angle |
| `nei_deg` | neighbouring lines within this angle count as "similar" when scoring | 3.6 | |
| `numlines` | lines detected per segment (Hough peaks) | 3 | |
| `crystal` | 1 = bcc/B2 ({100}, {110}, {112}, {123}); 2 = fcc ({111}, one plane per column) | 1 | |
| `scorefigs` | ST values for which figures are drawn | [0, 10] | tables always contain ST = 0..20 |

## Conventions
- **Euler angles**: Bunge (phi1, Phi, phi2), in degrees, in the frame of the SEM picture. RAST was developed with
  MTEX exports. A wrong convention or frame gives plausible-looking but wrong matches.
- **Score**: each line gains 1/o² for every similar line among its neighbours of order o (o = 1..10), and loses
  10/o for noise patterns (a 45° line with a neighbour at 90° from it; a 0° or 90° line with a 45° neighbour).
  Orientation differences are taken on the 180° circle (RAST 1.1).
- **ST**: at ST = s, a line is counted if its score is >= s, in tables and figures alike.
- **Matching**: the angle between a line and a trace is taken on the 180° circle (0–90°). A family matches if one of
  its traces is within IR. The family with the smallest angle is the "best" match (its colour in the figure); the
  other matching families are "second" matches; `m*` counts both.
- **Grain numbering**: consecutive pictures with identical Euler text belong to one grain.
- **Masks**: a segment containing any masked (or black) pixel is skipped entirely. The right and bottom strips of
  the picture that don't fill a whole segment are not analysed, and the picture border is set to zero.
