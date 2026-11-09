# Notes for AI coding agents working with RAST

You are in the RAST repository (Rapid Analysis of Slip Traces; Gagneur et al., Materials Characterization 229 (2025)
115524). RAST detects straight slip traces on SEM pictures of indents and matches each line to the slip-plane families
allowed by the grain's EBSD orientation. These notes are for any coding agent (Claude Code, Codex, Cursor, Copilot
...) a RAST user opens this folder with. They are the agent's memory for this project. Read this file fully, then
the note in `python/docs/agent-notes/` that fits the task.

## Layout

| Folder | What |
|---|---|
| `python/` | **RAST 1.1, Python.** Use this unless the user asks for MATLAB. Package `rast`, command `rast` |
| `matlab/` | RAST 1.1, MATLAB (needs the Image Processing Toolbox) |
| `matlab-1.0-legacy/` | v1.0, the exact code of the 2025 paper. Don't edit; use it only to reproduce the paper |
| `python/docs/` | user tutorials 1-3 and `agent-notes/` |

## Quick start (Python)

```
cd python && pip install .            # Python >= 3.11
rast demo                             # install check: downloads the Ti2AlMo demo (46 MB) and runs it
rast prepare --out sample.mat         # interactive: crop, mask, Euler angles (needs a desktop session)
rast run --input sample.mat           # analysis with the paper defaults
rast run -h / rast prepare -h         # all options
```

If `rast` is not on PATH, use `python -m rast ...`. Behind a proxy, set `HTTPS_PROXY` before `pip` and `rast demo`.

## Facts you can't guess

- **Input dataset** = a `.mat` with `no`, `Pictures` (uint8/uint16, cropped, masked pixels = 1), `BWss` (logical
  masks), `Eulers` (each a 1x1 cell with text `'[phi1,Phi,phi2]'`), `xi2`, `yi2`. Made by `rast prepare` or MATLAB
  `RAST_Loading.m`; both formats are the same. Details: `python/docs/agent-notes/input-and-output.md`.
- **Euler angles are Bunge, in degrees**, in the same frame as the SEM picture. A wrong convention gives wrong
  results with no error. `degrot` = anticlockwise rotation (deg) between SEM and EBSD frames.
- **Pictures of one grain must be consecutive**: `grainID` increments whenever the Euler angles change.
- **crystal 1** = bcc/B2, families {100}, {110}, {112}, {123}; **crystal 2** = fcc: the four columns 100/110/112/123
  are the planes (111), (-111), (1-11), (11-1).
- **Score threshold ST**: at ST = s, a line counts if its score >= s (every output). Tables hold ST = 0..20.
- **Parameters / defaults**: N 100 px, IR 10 deg, nsmooth 10, degrot 0, nei_deg 3.6, numlines 3, crystal 1,
  scorefigs [0, 10]. Meanings: `python/docs/agent-notes/parameters-and-conventions.md`.
- **Output**: `NAME/NAME_Tablemaster.xlsx|.csv` and `NAME/NAME_<N>_N_<IR>_IR/` with the per-run table, figures and
  `params.json`. Table columns: `python/docs/agent-notes/input-and-output.md`.
- Runtime: ~30-60 s per picture per N value (single worker); `--workers` runs pictures in parallel.

## Rules

1. **Never draw, change or remove masks, and never choose or change Euler angles yourself.** They are the user's
   experimental decisions. You may open the `rast prepare` windows for the user, read angles from a file the user
   gives you, or convert a file format, and you must say exactly what you did.
2. **Never overwrite** a dataset or a results folder. Use a new `--name` / `--out` for every run (RAST refuses
   existing folders unless `--overwrite` is given; only pass it when the user asks).
3. **Report the parameters with every result** (or point to `params.json`), and the ST used.
4. **Flag suspicious results instead of interpreting them**: many "no match" lines (`pnomatch`), real-line
   percentages close to the test-line ones (`pm*` vs `ptm*`), or an unexpected answer on a known sample. Likely
   causes: Euler convention or frame, `degrot`, scratches not masked. See `python/docs/agent-notes/checking-results.md`.
5. Don't change the algorithm (`python/rast/*.py`, `matlab/*.m`) unless the user asks. Results are meant to match
   the MATLAB version; any change must be applied to both and noted in `CHANGELOG.md`.
6. `rast prepare` without `--masks` opens windows. In a remote shell or CI there is no display: use the file mode
   (`--pictures`, `--masks`, `--eulers`), or ask the user to run that step on their desktop.

## Notes for specific tasks (`python/docs/agent-notes/`)

| Note | Read when |
|---|---|
| `running-rast.md` | installing, running, batch runs, parameter studies |
| `input-and-output.md` | reading/writing datasets, table columns, figure files |
| `parameters-and-conventions.md` | choosing parameters, Euler / crystal / ST conventions |
| `checking-results.md` | judging whether an analysis worked, summarising results for the user |
| `known-traps.md` | errors and platform issues |

User tutorials: `python/docs/tutorial_1_new_dataset.md` (new dataset), `tutorial_2_prepared_dataset.md` (prepared
dataset), `tutorial_3_ai_agent.md` (working with you).
