# Tutorial 1: a new dataset, from SEM pictures to results

This tutorial takes you from raw SEM pictures of indents to the RAST results table. It covers installing RAST,
masking indents and scratches, entering the EBSD orientations, running the analysis and checking that it worked.
For a dataset that is already prepared (masks drawn, angles entered), go to
[Tutorial 2](tutorial_2_prepared_dataset.md).

## 1. Install and check

You need Python 3.11 or newer ([python.org](https://www.python.org/downloads/); on Windows, tick "Add python.exe to
PATH" during installation).

```
git clone https://github.com/vgagneur/RAST
cd RAST/python
pip install .
rast demo
```

`rast demo` downloads a B2 Ti2AlMo dataset and analyses it. It should end with `done in ... min`. If it does, RAST
works on your computer. (Behind a proxy, set `HTTPS_PROXY` first, e.g. `set HTTPS_PROXY=http://proxy:port`.)

## 2. What you need

1. **SEM pictures of the indents** (TIFF, PNG, JPEG ...; 8- or 16-bit, grey or colour). One or more pictures per
   grain, with the slip traces around each indent clearly visible. A single indent per picture is simplest.
2. **The EBSD orientation of each grain**, as Euler angles in the **Bunge convention, in degrees**: (phi1, Phi, phi2).
   RAST was developed with orientations exported from MTEX. In MTEX, for a grain `id`:

   ```matlab
   o = grains(id).meanOrientation;
   fprintf('[%.3f,%.3f,%.3f]\n', o.phi1/degree, o.Phi/degree, o.phi2/degree)
   ```

   The SEM picture and the EBSD map must be in the same reference frame (same sample orientation, same axes). If the
   sample was rotated between the two (e.g. moved to another microscope), measure that rotation: it is the
   `degrot` parameter (anticlockwise, degrees) in step 5.
3. **The order of your pictures.** Pictures of the same grain must be entered one after the other: the results table
   numbers the grains (`grainID`) by consecutive identical Euler angles.

> **Check the convention first.** If your EBSD software uses another Euler convention or reference frame than
> RAST, every match will be wrong without any error message. Before analysing a new material, run one sample whose
> answer you know (section 7): B2 Ti2AlMo matches {110} close to 100 %; pure fcc Ni matches {111} close to 100 %.

## 3. Make the dataset: `rast prepare`

```
rast prepare --out my_sample.mat
```

For each picture, RAST:

1. **Asks for the picture file** (a file dialog). Choose the pictures in grain order. Cancel when you have entered
   them all. Or give them all at once: `rast prepare --out my_sample.mat --pictures p1.tif p2.tif p3.tif`.
2. **Crop.** Drag a rectangle around the part to keep, then press **Enter**. Leave out the SEM information bar
   (magnification, voltage ...) and anything outside the grain. Press **n** to keep the whole picture.
3. **Mask.** Click the corners of a polygon around an area that RAST must ignore, click the first corner again to
   close it, then press **a** to add it to the mask (shown in red). Draw as many polygons as needed. **u** undoes the
   last polygon, **c** clears them all, **Enter** finishes.
4. **Confirm.** The masked picture is shown: **k** keeps it, **r** redoes the crop and mask, **s** skips the picture.
5. **Euler angles.** Type them in the terminal as `[phi1,Phi,phi2]`, e.g. `[345.057,160.014,265.877]`.

The dataset `my_sample.mat` is then saved. It works in Python and in MATLAB (`RAST_Running.m`).

### What to mask, and why

RAST finds straight lines in every 100 × 100 px segment of the picture and keeps those that line up with similar
lines in neighbouring segments. Anything straight that is not a slip trace becomes a false line and dilutes the
statistics. **A segment that touches the mask is not analysed at all**, so mask generously but not more than needed.

| Mask | Don't mask |
|---|---|
| the indent itself, with its edges and pile-up rim | the slip traces around the indent |
| scratches from polishing (long straight lines crossing the picture) | faint traces: RAST is designed to pick them up |
| other grains, and grain boundaries | small dust specks far from any trace (they rarely give straight lines) |
| cracks, dust, contamination, image labels or scale bars left after cropping | |
| strong straight features that are not slip (e.g. second-phase lamellae) | |

### Euler angles from a file

Instead of typing the angles, give a CSV file (see [`examples/eulers_example.csv`](../examples/eulers_example.csv)):

```
picture,phi1,Phi,phi2
indent_01.tif,345.057,160.014,265.877
indent_02.tif,345.057,160.014,265.877
indent_03.tif,12.5,40.2,301.7
```

`picture` is the file name (or its 1-based number). Then: `rast prepare --out my_sample.mat --eulers angles.csv`.

### Masks drawn in another program

If you prefer to draw masks in ImageJ, GIMP or Photoshop, save each one as a black-and-white PNG of the same size as
the picture (white = masked), named `<picture name>_mask.png`, all in one folder. Then no window is needed:

```
rast prepare --out my_sample.mat --pictures p1.tif p2.tif --masks mask_folder --eulers angles.csv
```

(Crop the pictures before drawing masks in this case.)

## 4. Check the dataset

Open the dataset again to look at the masks, without changing anything:

```
rast prepare --edit my_sample.mat --out my_sample_checked.mat
```

Each picture opens with its mask in red; press **Enter** to go to the next. To fix a picture later, see
[Tutorial 2, section 4](tutorial_2_prepared_dataset.md#4-fix-a-mask-or-an-euler-angle).

## 5. Run the analysis: `rast run`

```
rast run --input my_sample.mat
```

The defaults are those of the paper: N = 100 px segments, IR = 10°, 10 smoothing passes, bcc/B2 (crystal 1),
figures at ST = 0 and 10. Common options:

| Option | When |
|---|---|
| `--crystal 2` | fcc samples ({111} slip) |
| `--degrot 157` | the SEM picture is rotated by 157° (anticlockwise) relative to the EBSD frame |
| `--nsmooth 20` | dirty or poorly contrasted pictures, where microstructure gets detected as lines |
| `--N 100 120 140` | test the effect of the segment size (one analysis per value) |
| `--IR 6 8 10` | test the effect of the identification range |
| `--workers 4` | analyse 4 pictures at a time (faster on multi-core computers) |
| `--name run1` | name of the output folder (default: the dataset name) |

Each picture takes about 30 s to a minute. All parameters: `rast run -h`.

## 6. Read the results

Open the output folder (`my_sample/`):

- **Figures** `picture_<i>_ST_10_linesselected.png`: the lines with a score ≥ 10, coloured by the family they match
  best: {100} blue, {110} orange, {112} green, {123} yellow, magenta for lines that match only {112} and/or {123}
  (unambiguous <111> slip), red for no match. On the right, the theoretical trace directions of each family.
  For fcc, the four colours are the four {111} planes, and magenta marks lines that match more than one plane.
- **`<i>_bw.png`**: the detected edges. If they follow scratches or microstructure, mask those or raise `--nsmooth`.
- **The table** `my_sample_Tablemaster.xlsx`: one row per picture and score threshold ST = 0 to 20.
  `pm110` is the percentage of lines matching {110}, `ptm110` the same for 720 evenly spread test lines, i.e. what
  a random line would give for that grain orientation. `nlines` is the number of lines.
  [`data_analysis_tutorial.txt`](../../matlab-1.0-legacy/data_analysis_tutorial.txt) shows how to make the pivot
  charts and box plots of the paper in Excel.

**The score threshold ST.** A line's score grows with the number of similar lines around it and drops for patterns
typical of noise. Higher ST keeps fewer, more reliable lines. ST = 10 is a good starting point for bcc/B2 pictures
(it is the example used in `data_analysis_tutorial.txt`); for scratched samples a higher ST may be needed.

## 7. Did it work?

From the v1.0 tutorial, three checks:

1. **"No match" lines** (`pnomatch`) should be rare at the ST you use. Many no-match lines mean the Euler angles,
   their convention, or `degrot` are probably wrong, or that IR is too small.
2. **Real lines vs test lines.** Compare `pm100`, `pm110` ... with `ptm100`, `ptm110` .... If the real lines look
   like the test lines, RAST found no preference beyond what the grain orientation gives by chance: check the
   detection figures and the angles.
3. **A sample with a known answer.** B2 Ti2AlMo (50Ti-25Al-25Mo at.%, homogenised at 1200 °C) matches {110} close to
   100 %; clean pure fcc Ni matches {111} close to 100 %. Get that result with your microscope and EBSD settings
   before analysing materials with more complex slip, such as bcc Fe.

## 8. When to change parameters

| Symptom | Try |
|---|---|
| lines follow scratches, polishing marks or microstructure | mask them; `--nsmooth 20` |
| very few lines on clearly visible traces | check the edges in `<i>_bw.png`; lower ST when reading the table |
| traces are short and dense | `--N 80` (smaller segments; below ~90 px results get noisier) |
| results change a lot between N or IR values | the data is borderline: report the spread, use more indents |
| systematic offset between lines and traces in the figures | check `degrot` and the Euler convention |
