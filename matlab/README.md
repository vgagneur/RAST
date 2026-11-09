# RAST 1.1: MATLAB version

Requires MATLAB with the Image Processing Toolbox. Keep all `.m` files in one folder, together with your dataset.

> **Status: preview.** These files are v1.0 with the fixes listed in [CHANGELOG.md](../CHANGELOG.md). They are being
> run on the reference datasets before release. For the exact code of the 2025 paper, use
> [`../matlab-1.0-legacy`](../matlab-1.0-legacy).

1. **Make a dataset:** run `RAST_Loading.m`. For each picture, choose the file, crop it (drag the area to keep, then
   right-click and choose *Crop Image*), mask the indent (left-click the corners of a polygon, right-click to
   complete it, then right-click and choose *Create Mask*), and type the Euler angles (Bunge convention, degrees) as `[phi1,Phi,phi2]`.
   Enter pictures of the same grain one after another. The dataset is saved as a `.mat` file. To mask several areas
   of one picture, use the Python `rast prepare` instead (it reads and writes the same files).
2. **Run:** open `RAST_Running.m`, set `matfile` to your dataset and the parameters in section (1) (comments explain
   each one), and run it. The defaults are N = 100 px, IR = 10°, 10 smoothing passes, bcc/B2.
3. **Results:** a folder with the figures and the Excel table `<name>_Tablemaster.xlsx`.
   [`data_analysis_tutorial.txt`](../matlab-1.0-legacy/data_analysis_tutorial.txt) explains how to analyse it.

The Python version ([`../python`](../python)) runs the same steps without MATLAB, and its tutorials
([`../python/docs`](../python/docs)) also explain what to mask and how to check that an analysis worked.
