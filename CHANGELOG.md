# Changelog

## 1.1 (in preparation)

### Fixed
- **Best-match choice for near-vertical lines.** In 1.0, a line close to ±90° could match a trace on the other side
  of the vertical (e.g. a line at −86° and a trace at +84°), but the angle between them was then recorded as ~170°
  instead of ~10°, so that trace was ranked last when choosing the best-matching family. 1.1 measures that angle on
  the 180° circle (always between 0° and 90°). Which families a line matches is unchanged, so **all match counts in
  the results table are the same as in 1.0**. Only the best family, i.e. the colour of the line in the figures,
  can change (checked on the Ti2AlMo, TiFe and Fe datasets of the paper: 157 of ~34,500 detected lines, all
  near-vertical).
- **Scoring of near-vertical lines.** In 1.0, two neighbouring lines at +89° and −89° were not counted as similarly
  oriented when scoring. 1.1 measures that difference on the 180° circle. Near-vertical slip traces can therefore
  get slightly higher scores than in 1.0 (on the same three datasets: 1,116 of ~34,500 line scores change; the number
  of lines at ST = 10 is unchanged in each dataset).

### Changed
- Default parameters in `RAST_Running.m`: one segmentation size (N = 100 px), IR = 10°, figures at ST = 0 and 10.
- The score-threshold rule is now stated in the code: at threshold ST, a line is counted when its score is ≥ ST,
  in the tables and in the figures (as in 1.0).

### Added
- Python version (`python/`) of the full default pipeline: dataset preparation (cropping, masking, Euler angles,
  typed or from a CSV file), analysis, figures, and the results table as `.xlsx` and `.csv`. Datasets prepared in
  MATLAB run in Python and the reverse.
- Tutorials for new and already-prepared datasets, and for running RAST with an AI coding agent (`AGENTS.md`).

## 1.0 (2025)

Version used in the paper (Materials Characterization 229 (2025) 115524). Files in `matlab-1.0-legacy/`, and the
`v1.0` tag.
