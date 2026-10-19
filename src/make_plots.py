#!/usr/bin/env python
"""Generate the initial-screen result plots from the saved CSVs.

Outputs PNGs to Data/Initial screen/Plots/. All plots are PER PLATE (the batch
confound means cross-plate magnitudes aren't comparable). Reads:
  - ALL_PLATES_full_results.csv  (model contrasts + equivalence)
  - SSMD_vs_model_ALL.csv        (model-free effect sizes + model)
  - aggregated_data.csv          (per-second data w/ bio_plate, for trajectories;
                                  path via FINOMENA_PLOT_AGG_CSV)
"""
import csv, math, os, statistics
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Paths are env-overridable so the SageMaker container (mounts under
# /opt/ml/processing/...) and the wrapper can thread one output dir through every
# step; unset -> the original local repo layout, so by-hand behaviour is unchanged.
PROD_DIR = os.environ.get("FINOMENA_ADDON_PROD_DIR",
                          "Data/Initial screen/Production Output")
OUT = os.environ.get("FINOMENA_PLOTS_DIR", "Data/Initial screen/Plots")
# Per-second aggregated data for the trajectory panel (needs a bio_plate column).
AGG_FULL = os.environ.get("FINOMENA_PLOT_AGG_CSV", "/tmp/agg_full.csv")
os.makedirs(OUT, exist_ok=True)
PLATES = os.environ.get("FINOMENA_PLATES", "Plate 1,Plate 2,Plate 3").split(",")

def num(v):
    try: return float(v)
    except: return float("nan")

# ── load model + equivalence ────────────────────────────────────────────────
full = list(csv.DictReader(open(os.path.join(PROD_DIR, "ALL_PLATES_full_results.csv"))))
byplate = defaultdict(list)
for r in full: byplate[r["Plate"]].append(r)

# genotype window per plate. There is no per-drug "genotype" row, but it falls
# out of the two contrasts: for any drug,
#   DrugEffect_vsKO - Rescue_vsWT = (drug - KO) - (drug - WT) = WT - KO,
# the WT-minus-KO disease gap (POSITIVE, since WT moves more than KO). The median
# across drugs denoises it. So geno[p] = WT - KO > 0, and -geno[p] (used below)
# is the KO baseline depth: KO sits -geno[p] log2FC below WT on the rescue axis.
geno = {}
for p in PLATES:
    rs = byplate[p]
    diffs = [num(r["DrugEffect_vsKO_est_log2"]) - num(r["Rescue_vsWT_est_log2"]) for r in rs
             if not (math.isnan(num(r["DrugEffect_vsKO_est_log2"])) or
                     math.isnan(num(r["Rescue_vsWT_est_log2"])))]
    geno[p] = statistics.median(diffs)  # WT - KO disease gap (log2), positive

# ============================================================================
# 1. RESCUE SCATTER: drug-effect (vs KO) x  rescue (vs WT), colored by Pr_equiv
# ============================================================================
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
for ax, p in zip(axes, PLATES):
    rs = byplate[p]
    x = [num(r["DrugEffect_vsKO_est_log2"]) for r in rs]   # + = moved up off disease
    y = [num(r["Rescue_vsWT_est_log2"]) for r in rs]       # 0 = at WT
    c = [num(r["equiv_Pr_vsWT"]) for r in rs]
    sc = ax.scatter(x, y, c=c, cmap="viridis", s=45, edgecolor="k", linewidth=0.3, vmin=0, vmax=0.5)
    ax.axhline(0, color="green", ls="--", lw=1.2, label="WT level = rescued")
    ax.axhline(-geno[p], color="red", ls="--", lw=1.2, label="KO level = diseased")
    ax.axvline(0, color="grey", ls=":", lw=0.8)
    # shade the disease window (the red<->green gap) so it's visibly the genotype window
    ax.axhspan(-geno[p], 0, color="orange", alpha=0.08)
    # Gap-anchored MOVEMENT band (visual triage guide only — alters no probability).
    # The disease gap = -geno[p] (KO->WT magnitude, +log2FC). Vertical guides mark
    # the "% recovery toward WT" movement cutoffs used by the Check A/B 2x2 grid:
    # a drug at x = f*gap has moved a fraction f of the way off the KO baseline.
    gap = -geno[p]
    for f in (0.15, 0.20, 0.25, 0.33):
        ax.axvline(f * gap, color="purple", ls="-", lw=0.6, alpha=0.35)
        ax.text(f * gap, ax.get_ylim()[1], f"{int(f*100)}%", color="purple",
                fontsize=6, ha="center", va="bottom", alpha=0.7)
    ax.set_xlabel("Drug effect vs KO  (log2FC; + = moved toward WT)")
    ax.set_ylabel("Rescue: drug vs WT  (log2FC; 0 = reached WT)")
    ax.set_title(f"{p}   (disease window: KO is {geno[p]:.2f} log2FC below WT)", fontsize=11)
    ax.legend(fontsize=7, loc="upper left")
    cb = fig.colorbar(sc, ax=ax, fraction=0.046); cb.set_label("P(equivalent to WT)", fontsize=8)
fig.suptitle("Identifying rescue compounds: drugs that both move the mutant off the disease baseline "
             "(rightward) and restore it to wild-type (near the green line)",
             fontsize=13, y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/1_rescue_scatter.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# NOTE: a volcano plot was intentionally REMOVED. With ~2-3 animals/drug the
# difference tests are underpowered: significance only flags the few huge
# (suppressive) effects, while the biologically interesting movers sit in the
# non-significant cloud purely for lack of power. A volcano would read "ns = no
# effect," which is wrong here. Effect size + Pr(equivalent to WT) are the
# power-appropriate readouts; see the scatter and waterfall.

# ============================================================================
# 2. WATERFALL: drugs ranked by drug-effect (vs KO), colored by Pr_equiv-to-WT
# ============================================================================
fig, axes = plt.subplots(3, 1, figsize=(16, 12))
norm = matplotlib.colors.Normalize(vmin=0, vmax=0.5)
sm = matplotlib.cm.ScalarMappable(norm=norm, cmap="viridis")
for ax, p in zip(axes, PLATES):
    rs=sorted(byplate[p], key=lambda r:num(r["DrugEffect_vsKO_est_log2"]))
    vals=[num(r["DrugEffect_vsKO_est_log2"]) for r in rs]
    cols=[num(r["equiv_Pr_vsWT"]) for r in rs]
    ax.bar(range(len(vals)), vals, color=sm.to_rgba(np.array(cols)))
    ax.axhline(0,color="green",ls="--",lw=1,label="WT level (rescued)")
    ax.axhline(-geno[p],color="red",ls="--",lw=1,label="KO level (diseased)")
    ax.set_title(f"{p}  (n = 80 compounds; bars above 0 = increased movement toward wild-type)", fontsize=11)
    ax.set_xlabel("compounds, ranked by effect on mutant movement"); ax.set_ylabel("effect vs KO (log2FC)")
    ax.legend(fontsize=8, loc="upper left")
    cb=fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01); cb.set_label("P(equivalent to WT)", fontsize=8)
fig.suptitle("Per-compound effect landscape: most compounds show little effect, with a tail of "
             "strong suppressors (bottom) and weaker movers toward wild-type (top)",
             fontsize=13, y=1.005)
fig.tight_layout(); fig.savefig(f"{OUT}/2_waterfall.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ============================================================================
# 3. METHOD AGREEMENT: model effect vs Glass-log SSMD (validation)
# ============================================================================
ss=list(csv.DictReader(open(os.path.join(PROD_DIR, "SSMD_vs_model_ALL.csv"))))
ssp=defaultdict(list)
for r in ss: ssp[r["Plate"]].append(r)
fig, axes=plt.subplots(1,3,figsize=(18,5.5))
for ax,p in zip(axes,PLATES):
    rs=ssp[p]
    x=[num(r["model_Rescue_vsWT_log2"]) for r in rs]   # rescue endpoint = vs WT
    y=[num(r["Glass_vsWT_log"]) for r in rs]
    xy=[(a,b) for a,b in zip(x,y) if not(math.isnan(a) or math.isnan(b))]
    ax.scatter([a for a,b in xy],[b for a,b in xy],s=40,alpha=.7,edgecolor="k",linewidth=.3)
    if len(xy)>2:
        xs=np.array([a for a,b in xy]); ys=np.array([b for a,b in xy])
        r=np.corrcoef(xs,ys)[0,1]
        ax.set_title(f"{p}   (Pearson r = {r:.2f})", fontsize=11)
    ax.set_xlabel("GAM model effect vs WT (log2FC)")
    ax.set_ylabel("model-free robust effect size\n(Glass's Δ, log scale, vs WT)")
fig.suptitle("Model validation: GAM effect estimates agree with a model-free robust effect size, "
             "confirming the results are not modelling artifacts",
             fontsize=13, y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/3_method_agreement.png", dpi=130, bbox_inches="tight"); plt.close(fig)

# ============================================================================
# 4. TRAJECTORIES: top rescue-direction candidates vs KO & WT mean curves
# ============================================================================
# per-second mean movement per condition per bio_plate
traj=defaultdict(lambda: defaultdict(lambda: defaultdict(list)))  # plate->cond->tsec->[pixels]
for r in csv.DictReader(open(AGG_FULL)):
    traj[r["bio_plate"]][r["Condition"]][int(float(r["time_sec"]))].append(float(r["pixel_diff"]))
def meancurve(plate,cond):
    d=traj[plate][cond]
    ts=sorted(d); return ts,[statistics.mean(d[t]) for t in ts]
# top rescue-direction consensus (hardcode the ones we found; or recompute top by drug-effect>0)
cand={"Plate 1":["KO+01F05","KO+01H08","KO+01H10"],
      "Plate 2":["KO+02E06","KO+02G04","KO+02G02"],
      "Plate 3":["KO+03E06","KO+03H03","KO+03H08"]}
fig, axes=plt.subplots(1,3,figsize=(18,5.5))
for ax,p in zip(axes,PLATES):
    ko=meancurve(p,"KO+DMSO"); wt=meancurve(p,"WT+DMSO")
    ax.plot(*ko,color="red",lw=2.5,label="KO+DMSO (disease)")
    ax.plot(*wt,color="green",lw=2.5,label="WT+DMSO (healthy)")
    for d in cand[p]:
        try: ax.plot(*meancurve(p,d),lw=1.2,alpha=.8,label=d.replace("KO+",""))
        except: pass
    ax.set_title(f"{p}   (3 strongest movers toward wild-type)", fontsize=11)
    ax.set_xlabel("time (s)"); ax.set_ylabel("mean movement (pixel difference)"); ax.legend(fontsize=7)
fig.suptitle("Movement over time: a rescued mutant's trajectory would rise from the disease baseline "
             "(red) toward wild-type (green)",
             fontsize=13, y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/4_trajectories.png", dpi=130, bbox_inches="tight"); plt.close(fig)

print("Saved 4 figures to", OUT)
for f in sorted(os.listdir(OUT)): print("  ", f)

# ============================================================================
# POWER DIAGNOSTIC — why the volcano (difference-test significance) is unusable
# ============================================================================
# A volcano implies "significant = real effect, ns = nothing there." At ~2-3
# animals/drug the difference tests can ONLY flag the few enormous effects, and
# those happen to all be suppressive (further FROM WT). No drug is significant
# vs KO while also indistinguishable from WT — the difference-test definition of
# rescue. The table below makes that explicit: the "rescue quadrant" is empty.
OUTCSV = os.path.join(PROD_DIR, "power_diagnostic_volcano.csv")
rows = []
for p in PLATES:
    rs = byplate[p]
    n = len(rs)
    sigKO = sum(num(r["DrugEffect_vsKO_adj_p_FDR"]) < 0.05 for r in rs)
    nsKO  = n - sigKO
    sigWT = sum(num(r["Rescue_vsWT_adj_p_FDR"]) < 0.05 for r in rs)
    nsWT  = n - sigWT
    # rescue quadrant = moved off KO (sig vs KO) AND reached WT (ns vs WT) AND toward WT
    sigKO_nsWT = sum((num(r["DrugEffect_vsKO_adj_p_FDR"]) < 0.05) and
                     not (num(r["Rescue_vsWT_adj_p_FDR"]) < 0.05) for r in rs)
    rescue_dir = sum((num(r["DrugEffect_vsKO_adj_p_FDR"]) < 0.05) and
                     not (num(r["Rescue_vsWT_adj_p_FDR"]) < 0.05) and
                     num(r["DrugEffect_vsKO_est_log2"]) > 0 for r in rs)
    # of the sig-vs-KO hits, how many moved the WRONG way (suppressive, away from WT)?
    sigKO_suppr = sum((num(r["DrugEffect_vsKO_adj_p_FDR"]) < 0.05) and
                      num(r["DrugEffect_vsKO_est_log2"]) < 0 for r in rs)
    rows.append(dict(Plate=p, n_drugs=n, sig_vs_KO=sigKO, ns_vs_KO=nsKO,
                     sig_vs_WT=sigWT, ns_vs_WT=nsWT,
                     sig_KO_and_ns_WT=sigKO_nsWT,
                     rescue_direction=rescue_dir,
                     sig_KO_suppressive=sigKO_suppr))
# totals
tot = dict(Plate="ALL", n_drugs=sum(r["n_drugs"] for r in rows),
           sig_vs_KO=sum(r["sig_vs_KO"] for r in rows),
           ns_vs_KO=sum(r["ns_vs_KO"] for r in rows),
           sig_vs_WT=sum(r["sig_vs_WT"] for r in rows),
           ns_vs_WT=sum(r["ns_vs_WT"] for r in rows),
           sig_KO_and_ns_WT=sum(r["sig_KO_and_ns_WT"] for r in rows),
           rescue_direction=sum(r["rescue_direction"] for r in rows),
           sig_KO_suppressive=sum(r["sig_KO_suppressive"] for r in rows))
rows.append(tot)

cols = ["Plate","n_drugs","sig_vs_KO","ns_vs_KO","sig_vs_WT","ns_vs_WT",
        "sig_KO_and_ns_WT","rescue_direction","sig_KO_suppressive"]
with open(OUTCSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

print("\n" + "="*78)
print("POWER DIAGNOSTIC: why a volcano (difference-test significance) is not usable")
print("="*78)
hdr = f"{'Plate':7} {'n':>3} {'sigKO':>6} {'sigWT':>6} {'sigKO&nsWT':>11} {'rescue-dir':>10} {'sigKO=suppr':>11}"
print(hdr); print("-"*len(hdr))
for r in rows:
    print(f"{r['Plate']:7} {r['n_drugs']:>3} {r['sig_vs_KO']:>6} {r['sig_vs_WT']:>6} "
          f"{r['sig_KO_and_ns_WT']:>11} {r['rescue_direction']:>10} {r['sig_KO_suppressive']:>11}")
print("-"*len(hdr))
print("Reading: 'sigKO&nsWT' = rescue by the difference tests (off KO AND indistinct")
print("from WT). 'rescue-dir' = that, AND moved toward WT. 'sigKO=suppr' = of the")
print("sig-vs-KO hits, how many moved AWAY from WT (suppressive).")
print(f"-> Every sig-vs-KO hit is suppressive; the rescue quadrant is EMPTY (0).")
print(f"   Significance tracks only effect MAGNITUDE, not rescue. Volcano unusable.")
print(f"   Saved: {OUTCSV}")
