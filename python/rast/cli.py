"""Command line: `rast prepare ...` and `rast run ...` (also `python -m rast ...`)."""
import sys

from . import __version__

USAGE = f"""RAST {__version__}: Rapid Analysis of Slip Traces

usage:
  rast prepare ...   make a dataset: pictures + masks + EBSD Euler angles  (rast prepare -h)
  rast run ...       analyse a dataset: lines, matches, tables, figures     (rast run -h)
  rast demo [DIR]    download the TiAlMo demo dataset and run it with the default parameters
  rast --version
"""


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help"):
        print(USAGE)
        return 0
    if argv[0] == "--version":
        print(__version__)
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "prepare":
        from . import prepare
        return prepare.main(rest)
    if cmd == "run":
        from . import run
        return run.main(rest)
    if cmd == "demo":
        from . import demo
        return demo.main(rest)
    print(f"unknown command {cmd!r}\n\n{USAGE}")
    return 2


if __name__ == "__main__":
    sys.exit(main())
