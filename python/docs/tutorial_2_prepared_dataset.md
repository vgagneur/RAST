# Tutorial 2: analysing a dataset that is already prepared

Use this tutorial when the masks are already drawn and the Euler angles already entered: a `.mat` dataset made with
`rast prepare` ([Tutorial 1](tutorial_1_new_dataset.md)) or with the MATLAB `RAST_Loading.m`. Both work the same.
Installation: see [Tutorial 1, section 1](tutorial_1_new_dataset.md#1-install-and-check).

## 1. Run with the defaults

```
rast run --input my_sample.mat
```

This writes the folder `my_sample/` next to where you run the command (or choose with `--out-dir D:\results`), with
the results table `my_sample_Tablemaster.xlsx` (and `.csv`) and one sub-folder per (N, IR) with figures. RAST
won't write into an existing results folder. Use `--name` to give each run its own name:

```
rast run --input my_sample.mat --name my_sample_nsmooth20 --nsmooth 20
```

## 2. Choose the parameters

On the command line:

```
rast run --input my_sample.mat --crystal 2 --nsmooth 15 --scorefigs 10 20 --workers 4
```

or in a parameter file, convenient to keep the settings of a study together (copy
[`examples/params_example.toml`](../examples/params_example.toml)):

```
rast run --input my_sample.mat --params my_params.toml
```

| Parameter | Meaning | Default |
|---|---|---|
| `N` | segment size, px (several values = one analysis each) | 100 |
| `IR` | identification range: max angle between a line and a trace for a match, deg | 10 |
| `nsmooth` | smoothing passes before line detection | 10 |
| `degrot` | anticlockwise rotation of the SEM picture relative to the EBSD frame, deg | 0 |
| `nei_deg` | max orientation difference for neighbouring lines to count as similar, deg | 3.6 |
| `numlines` | lines detected per segment | 3 |
| `crystal` | 1 = bcc / B2 ({100}, {110}, {112}, {123}), 2 = fcc ({111}) | 1 |
| `scorefigs` | score thresholds for which figures are drawn | 0 10 |

Command-line values win over the file. `params.json` in each output folder records what was used.

## 3. Ranges of N and IR

```
rast run --input my_sample.mat --N 100 120 140 --IR 6 8 10
```

runs nine analyses. Each (N, IR) gets its own folder, and the Tablemaster holds them all with `N` and `IR`
columns: in an Excel pivot chart, put `N` and `IR` in the filter box and select one value of each (see
`data_analysis_tutorial.txt`). Line detection is done once per N and reused for every IR.

## 4. Fix a mask or an Euler angle

Pictures can be re-masked or given new angles without starting again:

```
rast prepare --edit my_sample.mat --select 3 5 --out my_sample_v2.mat
```

opens pictures 3 and 5 with their current mask in red; draw extra polygons (**a** to add, **Enter** when done). To
change only Euler angles, without any window:

```
rast prepare --edit my_sample.mat --select 3 --eulers new_angles.csv --no-window --out my_sample_v2.mat
```

where `new_angles.csv` has the line `3,phi1,Phi,phi2` (the 1-based picture number) under the header
`picture,phi1,Phi,phi2`. A mask can be extended this way but not reduced, because masked pixels are overwritten in
the dataset. To un-mask an area, prepare that picture again from the original file.

Then run the new dataset under a new name and compare.

## 5. Run only some pictures

```
rast run --input my_sample.mat --pictures 1 2 7 --name test_three
```

Picture numbers in the outputs stay those of the dataset. Useful to try parameters quickly before a full run.

## 6. Tables only, or faster runs

- `--no-figures`: no PNG output (much faster to write, and the folders stay small).
- `--workers N`: N pictures at a time. Each worker holds one picture and its working arrays in memory; use fewer
  workers if the computer runs out of memory.

## 7. Compare two runs

The `.csv` tables can be read by any software. In Python, for example:

```python
import pandas as pd
a = pd.read_csv("run1/run1_Tablemaster.csv")
b = pd.read_csv("run2/run2_Tablemaster.csv")
st10 = lambda t: t[t.ST == 10].groupby("grainID")[["m110", "nlines"]].sum()
print(st10(a).assign(f110=lambda d: d.m110 / d.nlines), st10(b).assign(f110=lambda d: d.m110 / d.nlines))
```

## Notes

- Datasets saved by MATLAB in the v7.3 format (files over 2 GB) can't be read by Python. RAST tells you if that is
  the case: re-save in MATLAB with `save('name_v7.mat', '-v7')` after loading it.
- 16-bit pictures are supported.
- On Windows, very long folder paths (over 260 characters) are handled.
