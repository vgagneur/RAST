"""`rast demo [DIR]`: download the B2 Ti2AlMo demo dataset (if needed) and run RAST on it with the defaults.

Ti2AlMo slips almost only on {110}: in the results table, pm110 should be close to 100 % at ST = 10, well above
the test-line value ptm110. Use it to check an installation, or as a first example.
"""
import sys
import urllib.request
from pathlib import Path

DEMO_URL = "https://github.com/vgagneur/RAST/releases/download/v1.0/demo_B2_tialmo.mat"
DEMO_FILE = "demo_B2_tialmo.mat"


def fetch(folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / DEMO_FILE
    if target.exists():
        print(f"using {target}")
        return target
    print(f"downloading {DEMO_URL} (46 MB) ...", flush=True)
    tmp = target.with_suffix(".part")
    urllib.request.urlretrieve(DEMO_URL, tmp)       # honours the HTTPS_PROXY environment variable
    tmp.replace(target)
    return target


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    folder = Path(argv[0]) if argv else Path("rast_demo")
    ds = fetch(folder)
    from . import run
    name = "demo_TiAlMo"
    if (folder / name).exists():
        print(f"{folder / name} already exists: delete it or look at the results there")
        return 1
    return run.main(["--input", str(ds), "--name", name, "--out-dir", str(folder)] + argv[1:])


if __name__ == "__main__":
    sys.exit(main())
