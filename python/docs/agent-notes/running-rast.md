# Running RAST (agent note)

## Install
```
cd python
pip install .                 # or: pip install -e .  (editable, for development)
python -m rast --version      # 1.1.x
```
Dependencies: numpy, scipy, matplotlib, pillow, openpyxl. No MATLAB needed. Behind a proxy, set `HTTPS_PROXY`
(and `HTTP_PROXY`) first.

## Check the install
`rast demo [DIR]` downloads `demo_B2_tialmo.mat` (12 pictures of 5 B2 Ti2AlMo grains) from the v1.0 GitHub release
into `DIR` (default `rast_demo/`) and runs it with the defaults. Expected: at ST = 10, `pm110` close to 100 % in
every picture, with `ptm110` (test lines) much lower. Extra arguments go to `rast run`, e.g.
`rast demo rast_demo --workers 4`.

## Typical runs
```
rast run --input sample.mat --name sample_default
rast run --input sample.mat --name sample_fcc --crystal 2 --nsmooth 15 --scorefigs 10 20
rast run --input sample.mat --name sample_Nscan --N 80 100 120 140 --no-figures --workers 4
rast run --input sample.mat --name quick --pictures 1 2 --no-figures
rast run --input sample.mat --params study.toml --name study_run1
```
- One output folder per run name. RAST stops if it exists: use a new `--name`, and pass `--overwrite` only when the
  user asks.
- Line detection runs once per N and is reused for every IR value. Figures are written in each (N, IR) folder.
- `--workers k` processes k pictures in parallel (separate processes). Use fewer if memory runs out.
- Runtime: about 30–60 s per picture per N value with one worker, for ~3000 × 2000 px pictures.

## From Python code
```python
from rast import io_mat, pipeline, table
ds = io_mat.read_dataset("sample.mat")
p = dict(nsmooth=10, pixsize=100, numlines=3, nei_deg=3.6, degrot=0.0, crystal=1)
res = pipeline.run_pictures(ds["pictures"], ds["masks"], ds["eulers"], p, irs=[10.0], indices=range(ds["no"]))
rows = table.rows_for_run(ds["eulers"], [res[i]["per_ir"][10.0]["sol"] for i in range(ds["no"])], ir=10.0, n=100)
```
`res[i]["per_ir"][ir]["sol"][st]` is the 82-value vector of picture i at ST = st (names in `rast.process.SOL_NAMES`).
`res[i]["det"]["coefmat"]` holds every detected line (orientation, end points, score) per segment. Guard scripts that
use `workers > 1` with `if __name__ == "__main__":` (Windows).

## MATLAB instead
Run `matlab/RAST_Loading.m` to make a dataset. Then set the dataset name and the parameters in section (1) of
`matlab/RAST_Running.m`, and run it. It needs MATLAB with the Image Processing Toolbox. The file formats are the same
as in Python.
