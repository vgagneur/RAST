"""`rast run`: the RAST analysis of a prepared dataset (RAST_Running.m + RAST_function.m).

Output layout (as in MATLAB), for a run named NAME:
    NAME/NAME_Tablemaster.xlsx / .csv                  all (N, IR) runs together
    NAME/NAME_<N>_N_<IR>_IR/                            one folder per (N, IR)
        NAME_<N>_N_<IR>_IR_Table.xlsx / .csv            results table (one row per picture and ST)
        <i>_picture.png, <i>_contrast.png, <i>_bw.png   smoothed picture, contrast-adjusted segments, edges
        picture_<i>_ST_<s>_linesselected.png            figures for the ST values in --scorefigs
        params.json                                     parameters and RAST version of this run
"""
import argparse
import json
import sys
import time
from pathlib import Path

from . import __version__, figures, io_mat, pipeline, table

DEFAULTS = dict(N=[100], IR=[10.0], nsmooth=10, degrot=0.0, nei_deg=3.6, numlines=3, crystal=1, scorefigs=[0, 10])

PARAM_HELP = """\
Parameters (defaults as in RAST_Running.m):
  N         segmentation size in px; several values = one run each (MATLAB loop over pixsize)   [100]
  IR        identification range in deg: max angle between a line and a trace to call it a match [10]
  nsmooth   number of smoothing passes; raise (e.g. 20) for dirty or poorly contrasted pictures [10]
  degrot    anticlockwise rotation (deg) between SEM and EBSD frames, e.g. after changing microscope [0]
  nei_deg   max orientation difference (deg) for two neighbouring lines to count as similar       [3.6]
  numlines  lines detected per segment                                                            [3]
  crystal   1 = bcc / B2 ({100}, {110}, {112}, {123}); 2 = fcc (the four {111} planes)          [1]
  scorefigs score thresholds ST for which figures are drawn (tables always hold ST = 0..20)    [0 10]
"""


def build_parser():
    ap = argparse.ArgumentParser(prog="rast run", description="Run RAST on a prepared dataset (.mat).",
                                 epilog=PARAM_HELP, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="dataset .mat made by `rast prepare` or RAST_Loading.m")
    ap.add_argument("--name", help="name of the run (output folder); default: the dataset file name")
    ap.add_argument("--out-dir", default=".", help="folder in which the NAME folder is created [current folder]")
    ap.add_argument("--params", help="TOML file with any of the parameters below (command-line values win)")
    ap.add_argument("--N", type=int, nargs="+")
    ap.add_argument("--IR", type=float, nargs="+")
    ap.add_argument("--nsmooth", type=int)
    ap.add_argument("--degrot", type=float)
    ap.add_argument("--nei-deg", dest="nei_deg", type=float)
    ap.add_argument("--numlines", type=int)
    ap.add_argument("--crystal", type=int, choices=[1, 2])
    ap.add_argument("--scorefigs", type=int, nargs="*")
    ap.add_argument("--pictures", type=int, nargs="+", help="run only these pictures (1-based numbers)")
    ap.add_argument("--no-figures", action="store_true", help="skip all PNG output (tables only, faster)")
    ap.add_argument("--workers", type=int, default=1, help="pictures processed in parallel [1]")
    ap.add_argument("--overwrite", action="store_true", help="allow writing into an existing NAME folder")
    return ap


def resolve_params(a):
    p = dict(DEFAULTS)
    if a.params:
        import tomllib
        with open(a.params, "rb") as f:
            user = tomllib.load(f)
        unknown = set(user) - set(DEFAULTS)
        if unknown:
            raise SystemExit(f"unknown parameter(s) in {a.params}: {sorted(unknown)}; allowed: {sorted(DEFAULTS)}")
        p.update(user)
    for k in DEFAULTS:
        v = getattr(a, k, None)
        if v is not None:
            p[k] = v
    for k in ("N", "IR", "scorefigs"):
        if not isinstance(p[k], (list, tuple)):
            p[k] = [p[k]]
    p["N"] = [int(v) for v in p["N"]]
    p["IR"] = [float(v) for v in p["IR"]]
    p["scorefigs"] = [int(v) for v in p["scorefigs"]]
    if p["crystal"] not in (1, 2):
        raise SystemExit("crystal must be 1 (bcc/B2) or 2 (fcc)")
    return p


def _fmt(v):
    """Number as MATLAB num2str writes it in folder names (10.0 -> '10')."""
    return str(int(v)) if float(v).is_integer() else f"{v:g}"


def main(argv=None):
    a = build_parser().parse_args(argv)
    p = resolve_params(a)
    t0 = time.time()
    ds = io_mat.read_dataset(a.input)
    name = a.name or Path(a.input).stem
    main_dir = Path(a.out_dir) / name
    if main_dir.exists() and not a.overwrite:
        raise SystemExit(f"{main_dir} already exists; choose another --name, or pass --overwrite")
    main_dir.mkdir(parents=True, exist_ok=True)
    idx = list(range(ds["no"])) if not a.pictures else [i - 1 for i in a.pictures]
    bad = [i + 1 for i in idx if not 0 <= i < ds["no"]]
    if bad:
        raise SystemExit(f"--pictures {bad}: the dataset has pictures 1..{ds['no']}")
    print(f"RAST {__version__}: {ds['no']} pictures in {a.input}; running {len(idx)}; parameters {p}", flush=True)

    master = []
    for n in p["N"]:
        folders = {ir: main_dir / f"{name}_{n}_N_{_fmt(ir)}_IR" for ir in p["IR"]}
        for f in folders.values():
            f.mkdir(parents=True, exist_ok=True)
        pp = dict(nsmooth=p["nsmooth"], pixsize=n, numlines=p["numlines"], nei_deg=p["nei_deg"],
                  degrot=p["degrot"], crystal=p["crystal"])

        def on_result(i, r, folders=folders, n=n):
            if a.no_figures:
                return
            ii = i + 1
            for ir, f in folders.items():
                figures.save_picture(r["det"]["picture"], f / f"{ii}_picture.png")
                figures.save_picture(r["det"]["contrast"], f / f"{ii}_contrast.png")
                figures.save_picture(r["det"]["bw"], f / f"{ii}_bw.png")
                sol = r["per_ir"][ir]["sol"]
                for st in p["scorefigs"]:
                    if 0 <= st <= 20 and sum(sol[st][0:5]) != 0:     # MATLAB: bnlines = sum(y(1:5)) ~= 0
                        figures.lines_figure(r["smoothed"], ii, st, n, _fmt(ir), r["det"]["coefmat"],
                                             r["per_ir"][ir]["code"], p["numlines"], ds["eulers"][i],
                                             float(p["degrot"]), p["crystal"],
                                             f / f"picture_{ii}_ST_{st}_linesselected.png")

        print(f"N = {n} px, IR = {p['IR']}", flush=True)
        res = pipeline.run_pictures(ds["pictures"], ds["masks"], ds["eulers"], pp, p["IR"], idx,
                                    workers=a.workers, on_result=on_result, log=lambda s: print(s, flush=True))
        for ir, f in folders.items():
            rows = table.rows_for_run([ds["eulers"][i] for i in idx], [res[i]["per_ir"][ir]["sol"] for i in idx],
                                      ir, n, picture_ids=[i + 1 for i in idx])
            table.write(rows, f / f"{f.name}_Table")
            master += rows
            json.dump(dict(rast_version=__version__, input=str(a.input), pictures=[i + 1 for i in idx], N=n, IR=ir,
                           **{k: p[k] for k in ("nsmooth", "degrot", "nei_deg", "numlines", "crystal", "scorefigs")}),
                      open(f / "params.json", "w"), indent=1)
            print(f"  written {f}", flush=True)
    table.write(master, main_dir / f"{name}_Tablemaster")
    print(f"done in {(time.time() - t0) / 60:.1f} min; results in {main_dir}")


if __name__ == "__main__":
    sys.exit(main())
