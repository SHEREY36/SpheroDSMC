#!/usr/bin/env python3
"""
Compare SpheroDSMC with the original Python implementation of the method.

    python3 tests/compare_reference.py --fortran 'runs/validate/f*.ts' \
        --python '/path/to/py_s*.txt'

The Python runs (DSMC_V2, DSMC_0D_v2.simulation.run_simulation with the
encounter-unit model) write rows: t, tau_enc, T_tr, T_rot, T.  Both ensembles
are interpolated on a common encounter clock and compared with their standard
errors.
"""
import argparse
import glob

import numpy as np

GRID = np.array([1.0, 2.0, 4.0, 6.0, 8.0, 10.0])


def fortran(path):
    a = np.loadtxt(path)
    t, tau, ttr, trot = a[:, 0], a[:, 1], a[:, 3], a[:, 4]
    return np.array([np.interp(GRID, tau, x) for x in (t, ttr, trot, ttr / trot)])


def python(path):
    a = np.loadtxt(path)
    t, tau, ttr, trot = a[:, 0], a[:, 1], a[:, 2], a[:, 3]
    return np.array([np.interp(GRID, tau, x) for x in (t, ttr, trot, ttr / trot)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fortran", required=True)
    ap.add_argument("--python", required=True)
    a = ap.parse_args()
    F = np.array([fortran(p) for p in sorted(glob.glob(a.fortran))])
    P = np.array([python(p) for p in sorted(glob.glob(a.python))])
    worst = 0.0
    for i, name in enumerate(("t", "T_tr", "T_rot", "theta")):
        print(name)
        for j, g in enumerate(GRID):
            fm, fs = F[:, i, j].mean(), F[:, i, j].std(ddof=1) / np.sqrt(len(F))
            pm, ps = P[:, i, j].mean(), P[:, i, j].std(ddof=1) / np.sqrt(len(P))
            z = (fm - pm) / np.hypot(fs, ps)
            worst = max(worst, abs(z))
            print(f"  tau_enc = {g:5.1f}   fortran {fm:9.4f} +- {fs:.4f}   "
                  f"python {pm:9.4f} +- {ps:.4f}   z = {z:+.2f}")
    print(f"largest |z| = {worst:.2f}")


if __name__ == "__main__":
    main()
