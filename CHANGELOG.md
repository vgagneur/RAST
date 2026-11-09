# Changelog

## 1.1 (in preparation)

### Fixed
- **Matching of near-vertical lines.** In 1.0, a line close to ±90° could match a trace on the other side of the
  vertical (e.g. a line at −86° and a trace at +84°), but the angle between them was then recorded as ~170° instead of
  ~10°, so that trace was ranked as the worst candidate when choosing the best and second-best match. 1.1 measures
  every line-to-trace angle on the 180° circle (always between 0° and 90°).
- **Scoring of near-vertical lines.** In 1.0, two neighbouring lines at +89° and −89° were not counted as similarly
  oriented when scoring. 1.1 measures that difference on the 180° circle too. Near-vertical slip traces can therefore
  get higher scores than in 1.0.

### Changed
- Every output uses the same score-threshold rule: at threshold ST, a line is counted when its score is ≥ ST.
- Default parameters in `RAST_Running.m`: one segmentation size (N = 100 px), IR = 10°, figures at ST = 0 and 10.

### Added
- Python version (`python/`) of the full default pipeline: dataset preparation (cropping, masking, Euler angles,
  typed or from a CSV file), analysis, figures, and the results table as `.xlsx` and `.csv`. Datasets prepared in
  MATLAB run in Python and the reverse.
- Tutorials for new and already-prepared datasets, and for running RAST with an AI coding agent (`AGENTS.md`).

## 1.0 (2025)

Version used in the paper (Materials Characterization 229 (2025) 115524). Files in `matlab-1.0-legacy/`, and the
`v1.0` tag.
