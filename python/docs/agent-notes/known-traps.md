# Known traps (agent note)

- **No display.** `rast prepare` (without `--masks`) and `rast prepare --edit` (without `--no-window`) open
  matplotlib / Tk windows. Over SSH, in WSL without a display, in CI or in a container, they stop with "No window can
  be opened here". Use the file mode, or ask the user to run that step on their desktop.
- **Proxy.** `pip install` and `rast demo` need internet access; behind a proxy set `HTTPS_PROXY` / `HTTP_PROXY`.
- **MATLAB v7.3 datasets** (over 2 GB, HDF5 inside) can't be read; RAST prints the MATLAB command that re-saves them
  as v7.
- **Long Windows paths** (over 260 characters) are handled inside RAST, but other tools (Excel, zip) may fail on
  them. Keep output folders short.
- **Picture whose last pixel column is almost all black.** Smoothing looks for the SEM legend in the last column and
  needs at least 5 non-black pixels there; otherwise RAST stops with a message. Crop differently.
- **16-bit pictures** work (smoothing, contrast and edges use the 16-bit range).
- **Existing output folder**: RAST refuses to write into it; choose a new `--name`.
- **Same grain, different Euler text** (e.g. `[90,50,131]` and `[90.0,50,131]`) counts as two grains in `grainID`.
  Use identical text, or the same CSV values, for all pictures of a grain.
- **Windows and `--workers`**: works from the command line, and from scripts guarded by
  `if __name__ == "__main__":`.
- Runtime grows with picture size and the number of N values; use `--pictures` and `--no-figures` for quick tests.
