"""RAST_function.m steps (4)-(7) for one picture: smooth -> segment + detect -> score -> match and count.

Segmentation and scoring do not depend on the crystal mode or on IR; `detect` can therefore be run once and
`process.process_picture` called for every IR value.
"""
import time
from concurrent.futures import ProcessPoolExecutor

from . import hough, process, score, segment, smooth


def smooth_picture(I, mask, nsmooth):
    """Step (4): the smoothed picture as RAST keeps it (Pictures{i} after smoothing, mask set to 0)."""
    return smooth.smooth_picture(I, mask, nsmooth)


def detect(S, pixsize, numlines, nei_deg, degrot):
    """Steps (5)-(6) on a smoothed picture: segmentation, line detection, coefmat packing, scoring.

    Returns dict(coefmat (scored), cuth, cutl, picture (border-zeroed), contrast, bw) where the last three are
    the N_picture / N_contrast / N_bw pictures that RAST_picture_segmentation.m saves.
    """
    Sz = smooth.zero_border(S)
    seg = segment.segment_picture(Sz, pixsize, numlines, detect_lines=hough.detect_lines)
    cm = segment.pack_coefmat(seg["alllines"], seg["cuth"], seg["cutl"], numlines, degrot)
    cm = score.score_lines(cm, numlines, nei_deg, degrot)
    return dict(coefmat=cm, cuth=seg["cuth"], cutl=seg["cutl"], picture=Sz, contrast=seg["Li"], bw=seg["BW"])


def run_picture(I, mask, euler, p, irs):
    """All steps for one picture, for every IR in `irs` (one detection pass).

    p: dict(nsmooth, pixsize, numlines, nei_deg, degrot, crystal).
    Returns dict(smoothed, det (see detect), per_ir={IR: dict(sol, code)}).
    """
    S = smooth_picture(I, mask, p["nsmooth"])
    det = detect(S, p["pixsize"], p["numlines"], p["nei_deg"], p["degrot"])
    per_ir = {ir: process.process_picture(euler, det["coefmat"], p["numlines"], ir, p["crystal"]) for ir in irs}
    return dict(smoothed=S, det=det, per_ir=per_ir)


def _job(args):
    i, I, mask, euler, p, irs = args
    t0 = time.time()
    r = run_picture(I, mask, euler, p, irs)
    return i, r, time.time() - t0


def run_pictures(pictures, masks, eulers, p, irs, indices, workers=1, on_result=None, log=print):
    """Run the listed 0-based pictures; calls on_result(i, result) as each one finishes (in any order)."""
    jobs = [(i, pictures[i], masks[i], eulers[i], p, list(irs)) for i in indices]
    results = {}

    def _done(i, r, dt):
        n = int((r["det"]["coefmat"][:, :, p["numlines"]] != 0).sum()) * p["numlines"]
        log(f"  picture {i + 1}: {n} lines detected ({dt:.1f} s)")
        if on_result is not None:
            on_result(i, r)
        # keep only the small parts; the pictures are written by on_result and dropped to save memory
        results[i] = dict(det=dict(coefmat=r["det"]["coefmat"], cuth=r["det"]["cuth"], cutl=r["det"]["cutl"]),
                          per_ir=r["per_ir"])

    if workers > 1 and len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for i, r, dt in ex.map(_job, jobs):
                _done(i, r, dt)
    else:
        for j in jobs:
            _done(*_job(j))
    return results
