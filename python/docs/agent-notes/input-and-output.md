# Input and output formats (agent note)

## Dataset (.mat), written by `rast prepare` or MATLAB `RAST_Loading.m`
| Variable | Content |
|---|---|
| `no` | number of pictures |
| `Pictures` | 1 × no cell; each a 2-D uint8 (or uint16) grey picture, already cropped; masked pixels set to 1 |
| `BWss` | 1 × no cell of logical masks, same size as the picture (true = masked) |
| `Eulers` | 1 × no cell; each a 1 × 1 cell holding text like `'[345.057,160.014,265.877]'` (Bunge, degrees) |
| `xi2`, `yi2` | 1 × no cell of mask-polygon vertex coordinates (information only; NaN separates polygons) |

In Python: `rast.io_mat.read_dataset(path)` returns dict(pictures, masks, eulers, no, xi, yi), and
`rast.io_mat.write_dataset(path, pictures, masks, eulers, xi, yi)` writes one. MATLAB v7.3 files (over 2 GB) can't
be read: re-save them with `-v7` in MATLAB.

Other inputs of `rast prepare`:
- picture files: any format Pillow reads; colour is converted with MATLAB's rgb2gray weights;
- mask PNGs named `<picture stem>_mask.png` (white = masked);
- an Euler CSV with the header `picture,phi1,Phi,phi2`, where `picture` is the file name, its stem, or the 1-based
  number.

## Results
```
NAME/NAME_Tablemaster.xlsx|.csv           all rows of all (N, IR) runs
NAME/NAME_<N>_N_<IR>_IR/
    NAME_<N>_N_<IR>_IR_Table.xlsx|.csv    one row per picture × ST (0..20)
    <i>_picture.png                       smoothed, masked picture with zeroed border (what was segmented)
    <i>_contrast.png                      picture after per-segment contrast adjustment
    <i>_bw.png                            Canny edges per segment (shows the segment grid)
    picture_<i>_ST_<s>_linesselected.png  lines with score >= s, coloured by best family (s in scorefigs only)
    params.json                           rast_version, input, pictures, N, IR and all other parameters
```
`<i>` is the 1-based picture number in the dataset.

## Table columns (same as the MATLAB Excel sheet)
`Euler_1..3`; `grainID` (goes up by one when the Euler angles change between consecutive pictures); `pictureID`;
`ST`; `IR`; `N`. Then come 20 counts, the same 20 for the test lines (prefix `t`), 20 percentages of `nlines`
(prefix `p`), 20 percentages of `tnlines` (prefix `pt`), and finally `nlines`, `tnlines`.

The 20 counts:
- `m100 m110 m112 m123`: lines matching the family (within IR, as best or second-best family);
- `u100 u110 u112 u123`: lines matching only that family;
- `u111`: bcc only; lines matching {112} and/or {123} but neither {100} nor {110} (unambiguous <111> slip);
- `m100110 … m112123`: lines matching exactly those two families; `m110112123 … m100110112`: exactly three;
- `nomatch`: lines matching no family.

Test lines are 720 lines at 0.25° steps (0.25° to 180°) matched against the same grain, which gives the chance
level for that orientation. For fcc (crystal 2) the four family columns are the four {111} planes. `pXXX` is NaN
(an empty cell in .xlsx) when `nlines` = 0.
