# RAST 1.1: Python version

RAST (Rapid Analysis of Slip Traces) detects slip traces on SEM pictures of indents and matches them to the slip-plane
families allowed by each grain's EBSD orientation. This is the Python version of RAST 1.1. It runs the same default
pipeline as the MATLAB 1.1 files in [`../matlab`](../matlab), and it doesn't need a MATLAB licence.

> **Status: preview.** RAST-python 1.1 is a port of RAST-matlab 1.1 made with the AI model Claude Opus 5.5
> (Anthropic) and reviewed by V. Gagneur. It is being checked against MATLAB 1.1 outputs; see
> [VALIDATION.md](VALIDATION.md) for what has been checked so far.

## Install

You need Python 3.11 or newer.

```
git clone https://github.com/vgagneur/RAST
cd RAST/python
pip install .
```

Check the installation with the demo (downloads a 46 MB B2 Ti2AlMo dataset, then runs for a few minutes):

```
rast demo
```

In `rast_demo/demo_TiAlMo/demo_TiAlMo_Tablemaster.xlsx`, at ST = 10, `pm110` (lines matching {110}) should be close
to 100 % in every picture, well above the test-line value `ptm110`.

## Use

| Step | Command | Tutorial |
|---|---|---|
| Make a dataset from your pictures: crop, mask indents and scratches, give the Euler angles | `rast prepare --out my_sample.mat` | [Tutorial 1](docs/tutorial_1_new_dataset.md) |
| Analyse a dataset (yours, or one made with MATLAB `RAST_Loading.m`) | `rast run --input my_sample.mat` | [Tutorial 2](docs/tutorial_2_prepared_dataset.md) |
| Work with an AI coding agent | open the repository in the agent | [Tutorial 3](docs/tutorial_3_ai_agent.md) |

`rast prepare -h` and `rast run -h` list every option. If the `rast` command isn't found after installing, use
`python -m rast` instead.

## Outputs

As in MATLAB, a run named `NAME` writes:

```
NAME/NAME_Tablemaster.xlsx (.csv)            all (N, IR) runs together
NAME/NAME_100_N_10_IR/                       one folder per segmentation size N and identification range IR
    NAME_100_N_10_IR_Table.xlsx (.csv)       one row per picture and score threshold ST = 0..20
    1_picture.png, 1_contrast.png, 1_bw.png  smoothed picture, contrast-adjusted segments, detected edges
    picture_1_ST_10_linesselected.png        lines with score >= 10, coloured by best-matching family
    params.json                              parameters of this run
```

The table columns are the same as in the MATLAB Excel sheet, so
[`data_analysis_tutorial.txt`](../matlab-1.0-legacy/data_analysis_tutorial.txt) (pivot charts, box plots, checks)
applies unchanged.

## How it relates to the MATLAB code

| MATLAB file | Python |
|---|---|
| `RAST_Loading.m` | `rast prepare` (`rast/prepare.py`) |
| `RAST_Running.m`, `RAST_function.m` | `rast run` (`rast/run.py`, `rast/pipeline.py`, `rast/smooth.py`, `rast/score.py`, `rast/table.py`) |
| `RAST_picture_segmentation.m` | `rast/segment.py`, `rast/canny.py` (imadjust, Canny), `rast/hough.py` (Hough transform) |
| `RAST_line_processing.m` | `rast/process.py` |
| `RAST_euler_to_slip.m`, `RAST_slip_trace_analyser_single.m` | `rast/match.py` |
| `RAST_Figuremaking.m` | `rast/figures.py` |

The code follows the MATLAB steps one to one, including MATLAB's own image functions (wiener2, imadjust, Canny,
Hough), so that results agree with MATLAB to rounding level. Module docstrings point to the MATLAB lines they port.

## Citation

V. Gagneur, J. P. Magnussen, K. Ma, M. J. Lloyd, N. Martin, A. J. Knowles, *Automated indentation-based slip trace
analysis for bcc and B2 intermetallic plasticity*, Materials Characterization 229 (2025) 115524.
https://doi.org/10.1016/j.matchar.2025.115524
