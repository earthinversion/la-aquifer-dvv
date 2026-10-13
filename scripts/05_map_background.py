"""Shaded-relief background for the site map, drawn with GMT.

Cuts GMT's 15 arc-second Earth relief to the map region (downloaded once by GMT and
cached), illuminates it from the northwest, colors it with the muted palette in
data/relief_muted.cpt, and adds the coastline. The image has no text or symbols; step 6 draws those as vector graphics.

Requires the GMT command-line tool (`gmt`) on PATH, or GMT=/path/to/gmt.

Usage:
    python scripts/05_map_background.py
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config as C


def main():
    gmt = os.environ.get("GMT") or shutil.which("gmt")
    if not gmt:
        sys.exit("GMT not found: install it (conda install -c conda-forge gmt) or set GMT=")
    region = "-R{}/{}/{}/{}".format(*C.MAP_REGION)
    with tempfile.TemporaryDirectory(prefix="la-map-") as tmp:
        env = dict(os.environ, GMT_SESSIONDIR=tmp, GMT_TMPDIR=tmp,
                   GMT_SESSION_NAME=f"la-map-{os.getpid()}")

        def run(*args):
            subprocess.run([gmt, *args], env=env, check=True)

        if not C.RELIEF_GRID.exists():
            run("grdcut", "@earth_relief_15s", region, f"-G{C.RELIEF_GRID}")
        shade = Path(tmp) / "shade.nc"
        run("grdgradient", str(C.RELIEF_GRID), "-A315", "-Nt0.35", f"-G{shade}")
        run("begin", str(C.MAP_FILE.with_suffix("")), "png", "E1200")
        run("grdimage", str(C.RELIEF_GRID), region, f"-JM{C.MAP_WIDTH_IN}i",
            f"-C{C.RELIEF_CPT}", f"-I{shade}", "-B+n")
        run("coast", "-W0.25p,gray45", "-Df", "-A10")
        run("end")
    print(f"wrote {C.MAP_FILE}")


if __name__ == "__main__":
    main()
