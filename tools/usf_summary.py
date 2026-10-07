#!/usr/bin/env python3
"""Steady uniform-shear-flow results of a SpheroDSMC campaign, one row per (AR, alpha).

    python3 tools/usf_summary.py runs/C1/usf [--dem DEM_TABLE] [--out summary.csv]

Reads the logs <dir>/AR<AR>_a<alpha>_<cold|hot>_s<seed>.log (flow = 'usf').  The cold and
hot runs of one seed are averaged first, so the standard error (sample std / sqrt(n)) has
one degree of freedom per independent seed; start_gap is the largest |cold - hot| over the
seeds in units of the block standard errors of the runs (a check that both starts reached
the same steady state).  Every run must be contact model C1 unless --any-contact-model.
With --dem (e.g. DSMC_V2/DSMC_0D_v2/reference/usf_dem_C1.csv) the relative difference
from DEM is added for each quantity.
"""
import argparse
import csv
import glob
import math
import os
import re
from collections import defaultdict

NAMES = ["Tstar", "theta", "Pk_xx", "Pk_yy", "Pk_zz", "Pk_xy", "N1", "N2", "a2_tr", "a2_rot"]
DEM = {"Tstar": "Tstar", "theta": "theta", "Pk_xx": "Pk_xx", "Pk_yy": "Pk_yy", "Pk_zz": "Pk_zz",
       "Pk_xy": "Pk_xy", "contacts_per_strain": "nu_over_gdot", "a2_tr": "a2_tr"}


def read_log(path):
    txt = open(path).read()
    if "wall seconds total" not in txt:
        return None
    get = lambda label: [float(x) for x in re.search(label + r"\s*(.*)", txt).group(1).split()]
    return {"mean": dict(zip(NAMES, get(r"steady mean:"))),
            "se": dict(zip(NAMES, get(r"steady se \(10 blocks\):"))),
            "contacts_per_strain": get(r"contacts per particle per unit strain:")[0],
            "energy_balance": get(r"\(shear \+ collisional\)/shear:")[-1],
            "contact_model": re.search(r"contact model:\s*(\S+)", txt).group(1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs")
    ap.add_argument("--dem", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--any-contact-model", action="store_true")
    a = ap.parse_args()
    runs = defaultdict(dict)
    for path in sorted(glob.glob(os.path.join(a.runs, "AR*_a*_*_s*.log"))):
        m = re.search(r"AR([\d.]+)_a([\d.]+)_(cold|hot)_s(\d+)\.log$", path)
        r = read_log(path)
        if not m or r is None:
            continue
        if not a.any_contact_model and r["contact_model"] != "C1":
            raise SystemExit(f"{path}: contact model {r['contact_model']}")
        runs[(float(m.group(1)), float(m.group(2)))][(m.group(3), int(m.group(4)))] = r
    dem = {}
    if a.dem:
        dem = {(float(r["AR"]), round(float(r["alpha"]), 4)): r for r in csv.DictReader(open(a.dem))}
    rows = []
    for (ar, alpha), by in sorted(runs.items()):
        seeds = sorted({s for _, s in by if ("cold", s) in by and ("hot", s) in by})
        if len(seeds) < 2:
            continue
        row = {"AR": ar, "alpha": alpha, "seeds": len(seeds)}
        gap = 0.0
        for q in NAMES + ["contacts_per_strain", "energy_balance"]:
            per_seed = []
            for s in seeds:
                c, h = by[("cold", s)], by[("hot", s)]
                vc = c["mean"][q] if q in NAMES else c[q]
                vh = h["mean"][q] if q in NAMES else h[q]
                per_seed.append(0.5 * (vc + vh))
                if q in NAMES:
                    se = math.hypot(c["se"][q], h["se"][q])
                    gap = max(gap, abs(vc - vh) / se if se > 0 else 0.0)
            mean = sum(per_seed) / len(per_seed)
            sd = math.sqrt(sum((x - mean) ** 2 for x in per_seed) / (len(per_seed) - 1))
            row[q], row[q + "_se"] = mean, sd / math.sqrt(len(per_seed))
        row["start_gap_in_se"] = gap
        ref = dem.get((ar, round(alpha, 4)))
        if ref:
            for q, col in DEM.items():
                row[f"{q}_vs_DEM"] = row[q] / float(ref[col]) - 1.0
        rows.append(row)
    if not rows:
        raise SystemExit(f"no finished USF runs in {a.runs}")
    out = a.out or os.path.join(a.runs, "summary.csv")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"{out}: {len(rows)} cases")
    if dem:
        for q in DEM:
            d = [abs(r[f"{q}_vs_DEM"]) for r in rows if f"{q}_vs_DEM" in r]
            if d:
                print(f"  {q:20s} |DSMC/DEM - 1|: median {sorted(d)[len(d) // 2]:.2%}  max {max(d):.2%}")


if __name__ == "__main__":
    main()
