#!/usr/bin/env python3
"""
Export the collision closures of one (alpha, AR) plane to a single .npz model
file read by the Fortran DSMC (src/closure_model.f90).

    python3 tools/export_model.py --artifact closure_v2.npz \
        --angular angular_memory.json --cross-section encounter_cross_section.json \
        --loss-memory loss_memory.json --alpha 0.8 --AR 2.0 --out models/model_AR2.0_a0.80.npz

At fixed restitution alpha and aspect ratio AR every closure depends on the gas
only through theta = T_tr/T_rot, so the file holds, for each theta node:

  theta_nodes (K)              node temperature ratios (increasing)
  mean_loss (K)                fitted mean loss per encounter (rescales the loss
                               entering the exchange law at each node)
  loss_ct, loss_cr (K)         loss law E[eps|z] = c_t z + c_r (1 - z) (with --loss-memory)
  loss_kappa (K)               Beta concentration of the bounded loss law (loss-memory-v2)
  contact_model_c1             1 for contact model C1 data, 0 for the archived model R1
  loss_eps_mean (K)            plain per-encounter mean loss of the loss table
  xi_grid (NX), xi_enhancement (K, NX)   spin factor G(Xi)
  lambda1, lambda4 (K)         dissipative tilt of the exchange kernel
  memory_form (K)              0: memory lambda3 z;  1: bounded-logit cubic
  memory_coef (K, 3), memory_center (K), memory_scale (K)
  quantile_probability (NU)    uniform grid on [0, 1]
  energy_a_offsets (K+1), energy_a_grid (NA), energy_quantiles (NA, NU)
                               quantile function Q(U; a) of the exchange kernel
  angular (K, 10)              eta10 eta11 eta12 eta20 eta21 eta22 a_lo a_hi b_lo b_hi
  cs_theta, cs_dynamic_factor  sigma_c(theta) = A_bar D(theta)
  k_theta, k_mean              contacts per encounter (diagnostic)
  alpha, aspect_ratio, loss_beta_a, loss_beta_b, anchor (2)

All arrays are float64 (offsets and forms int64) and the archive is written
uncompressed (np.savez), which is what the Fortran reader expects.
"""
import argparse
import json

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--artifact", required=True)
    ap.add_argument("--angular", required=True, help="angular-memory JSON table")
    ap.add_argument("--cross-section", required=True, help="encounter cross-section JSON")
    ap.add_argument("--loss-memory", default=None,
                    help="loss-memory table: E[eps|z] = c_t z + c_r (1 - z) per node (v1), plus the "
                         "Beta concentration kappa of the bounded law (v2)")
    ap.add_argument("--contact-model", required=True, choices=("C1", "R1-legacy"),
                    help="contact model of the CTC data behind the artifact; checked against the "
                         "model_card.json next to the artifact when there is one")
    ap.add_argument("--elastic-plane", action="store_true",
                    help="export the learned alpha = 1 plane (elastic-equilibrium test only)")
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--AR", type=float, required=True)
    ap.add_argument("--beta-a", type=float, default=1.21,
                    help="loss shape: B ~ Beta(zeta_n/2, zeta_s/2), zeta_n = 2.42")
    ap.add_argument("--beta-b", type=float, default=3.67, help="zeta_s = 7.34")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    art = np.load(a.artifact, allow_pickle=False)
    coords = np.asarray(art["surface_coordinates"], float)
    on_plane = np.flatnonzero(np.isclose(coords[:, 0], a.alpha, atol=1e-9)
                              & np.isclose(coords[:, 2], a.AR, atol=1e-9))
    if on_plane.size < 2:
        raise SystemExit(f"no fitted (alpha={a.alpha}, AR={a.AR}) plane in {a.artifact}; "
                         "the Fortran code samples fitted planes only")
    nodes = on_plane[np.argsort(coords[on_plane, 1])]
    theta = coords[nodes, 1]
    if a.alpha >= 1.0 and not a.elastic_plane:
        raise SystemExit("alpha = 1 uses the elastic Borgnakke-Larsen block, not a fitted plane "
                         "(--elastic-plane exports the learned alpha = 1 plane for the equilibrium test)")
    import os
    card = os.path.join(os.path.dirname(os.path.abspath(a.artifact)), "model_card.json")
    if os.path.exists(card):
        declared = json.load(open(card)).get("contact_model_id", "R1-legacy")
        if declared != a.contact_model:
            raise SystemExit(f"{card} declares contact model {declared}, not {a.contact_model}")

    forms = np.asarray(art["energy_kernel_forms"]).astype(str)[nodes]
    memory_form = np.array([1 if f == "conditional_logit_cubic_v3" else 0 for f in forms], np.int64)
    params = np.asarray(art["energy_parameters"], float)[nodes]          # lambda1..4
    coef = np.asarray(art["energy_memory_coefficients"], float)[nodes]
    center = np.asarray(art["energy_memory_center"], float)[nodes]
    scale = np.asarray(art["energy_memory_scale"], float)[nodes]

    prob = np.asarray(art["quantile_probability"], float)
    if not np.allclose(np.diff(prob), prob[1] - prob[0]):
        raise SystemExit("quantile probability grid must be uniform")
    off = np.asarray(art["energy_a_offsets"], np.int64)
    a_grid_all = np.asarray(art["energy_a_grid"], float)
    q_all = art["energy_quantiles"]
    grids, tables, offsets = [], [], [0]
    for i in nodes:
        grids.append(a_grid_all[off[i]:off[i + 1]])
        tables.append(np.asarray(q_all[off[i]:off[i + 1]], float))
        offsets.append(offsets[-1] + off[i + 1] - off[i])

    # angular law with incoming-partition memory, per node
    ang_json = json.load(open(a.angular))
    if ang_json.get("schema") != "angular-memory-v1":
        raise SystemExit("not an angular-memory-v1 table")
    key = lambda al, th, ar: (round(al, 4), round(th, 5), round(ar, 4))
    by_node = {key(n["alpha"], n["theta"], n["AR"]): n for n in ang_json["nodes"]}
    angular = []
    for th in theta:
        n = by_node[key(a.alpha, th, a.AR)]
        angular.append([n["eta1"], n["xi"], n["rho1"], n["eta2"], n["zeta"], n["rho2"],
                        *n["a_support"], *n["b_support"]])

    loss_tables = {}
    if a.loss_memory:
        lm = json.load(open(a.loss_memory))
        if lm.get("schema") not in ("loss-memory-v1", "loss-memory-v2"):
            raise SystemExit("not a loss-memory table")
        by = {key(n["alpha"], n["theta"], n["AR"]): n for n in lm["nodes"]}
        rows = [by[key(a.alpha, th, a.AR)] for th in theta]
        loss_tables = dict(loss_ct=np.array([r["c_t"] for r in rows], float),
                           loss_cr=np.array([r["c_r"] for r in rows], float),
                           loss_eps_mean=np.array([r["eps_mean"] for r in rows], float))
        if lm["schema"] == "loss-memory-v2":
            loss_tables["loss_kappa"] = np.array([r["kappa"] or 0.0 for r in rows], float)
        elif a.contact_model == "C1":
            raise SystemExit("contact model C1 needs a loss-memory-v2 table (bounded loss law)")

    cs = json.load(open(a.cross_section))
    if cs.get("schema") != "encounter-cross-section-v2":
        raise SystemExit("need an encounter-cross-section-v2 table")
    row = cs["cross_section"][f"{a.AR:.4f}"]
    order = np.argsort(row["theta"])
    cs_theta = np.asarray(row["theta"], float)[order]
    cs_dyn = np.asarray(row["dynamic_factor"], float)[order]
    kk = sorted((e["theta"], e["k_mean"]) for e in cs["contacts_per_encounter"]
                if abs(e["alpha"] - a.alpha) < 1e-9 and abs(e["AR"] - a.AR) < 1e-9)
    if not kk:
        raise SystemExit("no contacts-per-encounter entries for this plane")

    np.savez(a.out,
             alpha=np.array([a.alpha]), aspect_ratio=np.array([a.AR]),
             loss_beta_a=np.array([a.beta_a]), loss_beta_b=np.array([a.beta_b]),
             theta_nodes=theta,
             mean_loss=np.asarray(art["energy_mean_loss"], float)[nodes],
             xi_grid=np.asarray(art["xi_grid"], float),
             xi_enhancement=np.asarray(art["xi_enhancement"], float)[nodes],
             lambda1=params[:, 0].copy(), lambda4=params[:, 3].copy(),
             memory_form=memory_form, memory_coef=coef, memory_center=center,
             memory_scale=scale, quantile_probability=prob,
             energy_a_offsets=np.asarray(offsets, np.int64),
             energy_a_grid=np.concatenate(grids), energy_quantiles=np.vstack(tables),
             angular=np.asarray(angular, float),
             cs_theta=cs_theta, cs_dynamic_factor=cs_dyn,
             k_theta=np.array([t for t, _ in kk]), k_mean=np.array([k for _, k in kk]),
             anchor=np.asarray(art["energy_anchor"], float)[nodes[0]],
             contact_model_c1=np.array([1.0 if a.contact_model == "C1" else 0.0]),
             **loss_tables)
    print(f"{a.out}: alpha={a.alpha} AR={a.AR} theta nodes {list(theta)}  "
          f"energy rows {offsets[-1]} x {prob.size}  forms {list(forms)}")


if __name__ == "__main__":
    main()
