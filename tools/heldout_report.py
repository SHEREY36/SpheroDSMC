#!/usr/bin/env python3
"""Closure predictions against CTC encounters the closure was never fitted to.

    python3 tools/heldout_report.py --dsmc-v2 ../DSMC_V2 --planes models_C1 --out runs/C1/heldout

For every held-out CTC shard (DSMC_V2 results/ctc/nodes/holdout/*/*: nodes of independent
seed at fitted theta, and nodes at withheld theta between the fitted ones) the incoming
encounters of the shard are passed through the closure of the exported plane (alpha, AR)
exactly as SpheroDSMC samples it (src/collisions.f90, src/closure_model.f90: linear
stencil in theta, loss eps ~ Beta(kappa mu, kappa (1 - mu)) with mu = c_t z + c_r (1 - z),
exchange by quantile interpolation of the node laws, scattering cosine from the node
mixture), and the outcome moments are compared with those of the CTC:

  eps_mean   <eps>                              loss per encounter
  c_t, c_r   least squares eps = c_t z + c_r (1 - z)
  exchange   <E ((1 - eps) z' - z)> / <E>       translational energy gained per encounter
  stress     <E (1 - eps) z' P2(cos chi)> / <E> anisotropy carried out of the encounter
  cos_chi    <E' cos chi> / <E'>, E' = (1 - eps) E   deflection, weighted by the outgoing energy:
             the angular law is fitted with this weight (fit_angular_memory_table.py), so it is
             its target; cos_chi_plain (unweighted) is reported as a diagnostic only
  k_mean     <k>                                contacts per encounter (closure: table at theta)

CTC standard errors come from the 128 attempt blocks of the shard; the closure draw is one
sample per CTC encounter, so its Monte Carlo error is of the same size and is combined.
ks_eps, ks_zout, ks_cos are two-sample Kolmogorov-Smirnov distances of the outcome laws
(ks_cos between the E'-weighted laws).
Writes <out>/heldout.csv and prints a summary per set and quantity.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

QUANTITIES = ["eps_mean", "c_t", "c_r", "exchange", "stress", "cos_chi", "k_mean"]


class Plane:
    """One exported (alpha, AR) plane with the sampling of SpheroDSMC."""

    def __init__(self, path):
        z = np.load(path)
        self.theta = z["theta_nodes"]
        self.loss_ct, self.loss_cr, self.kappa = z["loss_ct"], z["loss_cr"], z["loss_kappa"]
        self.mean_loss = z["mean_loss"]
        self.lambda1, self.lambda4 = z["lambda1"], z["lambda4"]
        self.memory_form, self.memory_coef = z["memory_form"], z["memory_coef"]
        self.memory_center, self.memory_scale = z["memory_center"], z["memory_scale"]
        self.nu = len(z["quantile_probability"])
        off = z["energy_a_offsets"]
        self.rows = [(int(off[k]), int(off[k + 1]) - 1) for k in range(len(self.theta))]
        self.a_grid, self.q = z["energy_a_grid"], z["energy_quantiles"]     # q[row, j]
        self.angular = z["angular"]                                          # angular[node, p]
        self.k_theta, self.k_mean = z["k_theta"], z["k_mean"]
        self.contact_model = "C1" if "contact_model_c1" in z.files and z["contact_model_c1"][0] > 0.5 else "R1-legacy"

    def stencil(self, theta):
        t = self.theta
        if theta <= t[0]:
            return [0], [1.0]
        if theta >= t[-1]:
            return [len(t) - 1], [1.0]
        k = int(np.searchsorted(t, theta, side="left")) - 1
        f = (theta - t[k]) / (t[k + 1] - t[k])
        return [k, k + 1], [1.0 - f, f]

    def memory(self, nd, z):
        if self.memory_form[nd] == 1:
            zc = np.clip(z, 1e-12, 1 - 1e-12)
            x = np.tanh(((np.log(zc / (1 - zc)) - self.memory_center[nd])
                         / max(self.memory_scale[nd], 1e-8)) / 2)
            c = self.memory_coef[nd]
            return c[0] * x + c[1] * x ** 2 + c[2] * x ** 3
        return self.memory_coef[nd][0] * z

    def sample(self, theta, z, rng):
        nodes, w = self.stencil(theta)
        ct = sum(wk * self.loss_ct[n] for n, wk in zip(nodes, w))
        cr = sum(wk * self.loss_cr[n] for n, wk in zip(nodes, w))
        kap = sum(wk * self.kappa[n] for n, wk in zip(nodes, w))
        sml = sum(wk * self.mean_loss[n] for n, wk in zip(nodes, w))
        # loss
        mu = np.minimum(ct * z + cr * (1 - z), 1 - 1e-8)
        eps = np.zeros_like(z)
        ok = (mu > 1e-12) & (kap > 0)
        eps[ok] = rng.beta(kap * mu[ok], kap * (1 - mu[ok]))
        # exchange: quantile interpolation at a shared U
        x = rng.random(len(z)) * (self.nu - 1)
        j = np.minimum(x.astype(int), self.nu - 2)
        fu = x - j
        zp = np.zeros_like(z)
        clamps = 0
        for n, wk in zip(nodes, w):
            cov = eps * self.mean_loss[n] / sml if (sml > 0 and self.mean_loss[n] > 0) else eps
            a = self.lambda1[n] + self.memory(n, z) + self.lambda4[n] * cov
            lo, hi = self.rows[n]
            g = self.a_grid[lo:hi + 1]
            clamps += int(np.sum((a < g[0]) | (a > g[-1])))
            a = np.clip(a, g[0], g[-1])
            r = np.clip(np.searchsorted(g, a, side="right") - 1, 0, len(g) - 2)
            blend = (a - g[r]) / (g[r + 1] - g[r])
            rr = r + lo
            qlo = (1 - fu) * self.q[rr, j] + fu * self.q[rr, j + 1]
            qhi = (1 - fu) * self.q[rr + 1, j] + fu * self.q[rr + 1, j + 1]
            zp += wk * ((1 - blend) * qlo + blend * qhi)
        # scattering cosine from the node mixture
        nd = np.full(len(z), nodes[0])
        if len(nodes) == 2:
            nd[rng.random(len(z)) >= w[0]] = nodes[1]
        p = self.angular[nd]
        a = np.clip(p[:, 0] + p[:, 1] * zp + p[:, 2] * z, p[:, 6], p[:, 7])
        b = np.clip(p[:, 3] + p[:, 4] * zp + p[:, 5] * z, p[:, 8], p[:, 9])
        fmax = np.maximum(a + b, -a + b)
        v = np.where(b < 0, -a / np.where(b < 0, 3 * b, 1.0), 2.0)
        inside = (b < 0) & (np.abs(v) < 1)
        fmax[inside] = np.maximum(fmax[inside], a[inside] * v[inside] + b[inside] * 0.5 * (3 * v[inside] ** 2 - 1))
        c = np.empty_like(z)
        todo = np.arange(len(z))
        while len(todo):
            cc = 2 * rng.random(len(todo)) - 1
            r = a[todo] * cc + b[todo] * 0.5 * (3 * cc * cc - 1) - fmax[todo]
            acc = rng.random(len(todo)) <= np.exp(r)
            c[todo[acc]] = cc[acc]
            todo = todo[~acc]
        k = float(np.interp(theta, self.k_theta, self.k_mean))
        return eps, zp, c, k, clamps


def moments(z, eps, zp, c, E, k):
    p2 = 0.5 * (3 * c * c - 1)
    X = np.column_stack([z, 1 - z])
    ct, cr = np.linalg.lstsq(X, eps, rcond=None)[0]
    Ef = (1 - eps) * E
    return {"eps_mean": eps.mean(), "c_t": ct, "c_r": cr,
            "exchange": np.mean(E * ((1 - eps) * zp - z)) / E.mean(),
            "stress": np.mean(Ef * zp * p2) / E.mean(),
            "cos_chi": np.sum(Ef * c) / np.sum(Ef), "cos_chi_plain": c.mean(), "k_mean": np.mean(k)}


def block_se(z, eps, zp, c, E, k, block, nb=128):
    vals = []
    for b in range(nb):
        m = block == b
        if m.sum() > 20:
            vals.append(moments(z[m], eps[m], zp[m], c[m], E[m], k[m] if np.ndim(k) else k))
    return {q: np.std([v[q] for v in vals], ddof=1) / np.sqrt(len(vals)) for q in QUANTITIES}


def ks(x, y, wx=None, wy=None):
    """Two-sample KS distance of (optionally weighted) empirical laws."""
    wx = np.ones(len(x)) if wx is None else np.asarray(wx, float)
    wy = np.ones(len(y)) if wy is None else np.asarray(wy, float)
    ix, iy = np.argsort(x), np.argsort(y)
    x, y, cx, cy = x[ix], y[iy], np.cumsum(wx[ix]) / wx.sum(), np.cumsum(wy[iy]) / wy.sum()
    grid = np.concatenate([x, y])
    fx = np.concatenate([[0.0], cx])[np.searchsorted(x, grid, side="right")]
    fy = np.concatenate([[0.0], cy])[np.searchsorted(y, grid, side="right")]
    return float(np.max(np.abs(fx - fy)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dsmc-v2", required=True, help="DSMC_V2 checkout (contracts and results/ctc)")
    ap.add_argument("--shards", default="results/ctc/nodes/holdout/*/*")
    ap.add_argument("--planes", default="models_C1")
    ap.add_argument("--out", default="runs/C1/heldout")
    ap.add_argument("--seed", type=int, default=20261007)
    a = ap.parse_args()
    v2 = os.path.abspath(a.dsmc_v2)
    sys.path.insert(0, os.path.join(v2, "contracts", "python"))
    from dsmc_v2_contracts import load_run
    from dsmc_v2_contracts.io import AI, OI
    rng = np.random.default_rng(a.seed)
    rows = []
    for shard in sorted(glob.glob(os.path.join(v2, a.shards))):
        if not os.path.isfile(os.path.join(shard, "_SUCCESS")):
            continue
        run = load_run(shard)
        md = run.metadata
        alpha, theta, ar = float(md["alpha"]), float(md["theta"]), float(md["aspect_ratio"])
        tag = f"{ar:.2f}"
        tag = tag[:-1] if tag.endswith("0") else tag                  # 2.0, 1.25, 1.1 (export_planes.sh)
        plane_path = os.path.join(a.planes, f"model_AR{tag}_a{alpha:.2f}.npz")
        plane = Plane(plane_path)
        if plane.contact_model != md.get("contact_model_id", "R1-legacy"):
            raise SystemExit(f"{shard}: shard contact model {md.get('contact_model_id')} "
                             f"but plane {plane.contact_model}")
        o = np.asarray(run.outcomes)
        v = o["values"]
        att = np.asarray(run.attempts)
        lookup = {(int(e), int(i)): n for n, (e, i) in enumerate(zip(att["event_id"], att["attempt_index"]))}
        idx = np.array([lookup[(int(e), int(i))] for e, i in zip(o["event_id"], o["attempt_index"])])
        av = att["values"]
        c1 = av[idx][:, [AI[f"c1_{x}"] for x in "xyz"]]
        c2 = av[idx][:, [AI[f"c2_{x}"] for x in "xyz"]]
        vcm = 0.5 * (c1 + c2)
        E = v[:, OI["e_initial"]]
        z = float(md["mass"]) * (np.sum((c1 - vcm) ** 2, 1) + np.sum((c2 - vcm) ** 2, 1)) / E
        eps = v[:, OI["delta_total"]] / E
        zp = v[:, OI["et_inelastic"]] / (E - v[:, OI["delta_total"]])
        bad = (zp <= 0) | (zp >= 1)
        et = v[bad, OI["et_inelastic"]]
        zp[bad] = et / (et + v[bad, OI["er1_inelastic"]] + v[bad, OI["er2_inelastic"]])
        gpre = v[:, [OI[f"ghat_pre_{x}"] for x in "xyz"]]
        gpost = v[:, [OI[f"ghat_post_{x}"] for x in "xyz"]]
        cos = np.clip(np.einsum("ni,ni->n", gpre, gpost), -1, 1)
        k = o["n_contact"].astype(float)
        block = o["block_id"].astype(int) % 128
        # closure draw on the same incoming encounters
        eps_c, zp_c, cos_c, k_c, clamps = plane.sample(theta, z, rng)
        ctc, cl = moments(z, eps, zp, cos, E, k), moments(z, eps_c, zp_c, cos_c, E, np.full(len(z), k_c))
        se_ctc = block_se(z, eps, zp, cos, E, k, block)
        se_cl = block_se(z, eps_c, zp_c, cos_c, E, np.full(len(z), k_c), block)
        se_cl["k_mean"] = 0.0
        row = {"set": os.path.basename(os.path.dirname(shard)), "alpha": alpha, "theta": theta, "AR": ar,
               "events": len(z), "exchange_clamps": clamps,
               "ks_eps": ks(eps, eps_c) if alpha < 1 else 0.0, "ks_zout": ks(zp, zp_c),
               "ks_cos": ks(cos, cos_c, (1 - eps) * E, (1 - eps_c) * E),
               "cos_chi_plain_ctc": ctc["cos_chi_plain"], "cos_chi_plain_closure": cl["cos_chi_plain"]}
        for q in QUANTITIES:
            se = float(np.hypot(se_ctc[q], se_cl[q]))
            row[f"{q}_ctc"], row[f"{q}_closure"], row[f"{q}_se"] = ctc[q], cl[q], se
            row[f"{q}_z"] = (cl[q] - ctc[q]) / se if se > 0 else 0.0
        rows.append(row)
        print(f"{row['set']:22s} a={alpha:.2f} th={theta:<6g} AR={ar:<4g} "
              + " ".join(f"{q}:{row[q + '_z']:+.1f}" for q in QUANTITIES)
              + f"  KS eps {row['ks_eps']:.3f} z' {row['ks_zout']:.3f} cos {row['ks_cos']:.3f}", flush=True)
    os.makedirs(a.out, exist_ok=True)
    import csv
    path = os.path.join(a.out, "heldout.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.8g}" if isinstance(v, float) else v) for k, v in r.items()})
    summary = {}
    for s in sorted({r["set"] for r in rows}):
        g = [r for r in rows if r["set"] == s]
        summary[s] = {"nodes": len(g)}
        for q in QUANTITIES:
            zz = np.abs([r[f"{q}_z"] for r in g])
            rel = [abs(r[f"{q}_closure"] / r[f"{q}_ctc"] - 1) for r in g if abs(r[f"{q}_ctc"]) > 1e-9]
            summary[s][q] = {"median_abs_z": float(np.median(zz)), "max_abs_z": float(np.max(zz)),
                             "median_abs_rel": float(np.median(rel)) if rel else None}
        for q in ("ks_eps", "ks_zout", "ks_cos"):
            summary[s][q] = {"median": float(np.median([r[q] for r in g])), "max": float(np.max([r[q] for r in g]))}
    json.dump(summary, open(os.path.join(a.out, "heldout_summary.json"), "w"), indent=1)
    print(f"{path}: {len(rows)} nodes")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
