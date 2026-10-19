#!/usr/bin/env python
"""Model-free robust effect sizes (SSMD + robust Glass's delta) per drug per plate.

Reconstructed from the original inline step (run ad-hoc on 2026-06-11, never
previously committed). Verified to reproduce the committed SSMD_vs_model_ALL.csv
to algorithmic precision (differences only from aggregated-data drift; control
medians match to 4 decimals, formula is exact).

Per (bio_plate, condition), collapses the per-frame movement to a per-animal
mean, then takes robust summaries (median + MAD, MAD scaled by 1.4826) over the
animals. For each drug vs each control (KO+DMSO, WT+DMSO), on both raw and
log2(mean+1) scales:

  SSMD (two-MAD, robust):  (med_drug - med_ctrl) / sqrt(MAD_drug^2 + MAD_ctrl^2)
  robust Glass's delta:    (med_drug - med_ctrl) / MAD_ctrl

Glass uses the CONTROL spread only (the screen's reference scale); SSMD uses the
pooled drug+control spread. The script also copies the model's Drug_effect /
Rescue log2FC + adjusted p (read-only) alongside, for the method-agreement plot.

NOTE: the "mean_*"/"logmean_*" output columns actually hold the per-condition
MEDIAN of per-animal means (robust center), not the mean — column names kept for
backward compatibility with the committed file's schema. This is read DOWNSTREAM
as-is per the rescue-screen spec: Glass's delta is never recomputed or shrunk by
the addons; this script merely regenerates the same frozen quantity in-pipeline.

Inputs (env-overridable):
  FINOMENA_ADDON_AGG_CSV   aggregated_data.csv (per-frame long; needs bio_plate,
                           Condition, animal_id, pixel_diff). bio_plate is
                           derived from `plate` if the column is absent.
  per-plate master_results.csv (for the copied model_* columns) via OUTPUT_DIR.
"""
import csv
import math
import os
import re
import statistics
from collections import defaultdict

# Single source of truth: aggregated_data.csv lives wherever the BAM writes it
# (FINOMENA_OUTPUT_DIR), so its default is DERIVED from that dir rather than a
# divergent standalone path. Each is still independently env-overridable.
OUTPUT_DIR = os.environ.get(
    "FINOMENA_OUTPUT_DIR",
    "Data/Initial screen/Production Output")
PROD_DIR = os.environ.get("FINOMENA_ADDON_PROD_DIR", OUTPUT_DIR)
AGG_CSV = os.environ.get(
    "FINOMENA_ADDON_AGG_CSV",
    os.path.join(OUTPUT_DIR, "aggregated_data.csv"))  # derived from output dir
PLATES = os.environ.get("FINOMENA_PLATES", "Plate 1,Plate 2,Plate 3").split(",")
POS = os.environ.get("FINOMENA_POS_CONTROL", "KO+DMSO")  # KO+DMSO
NEG = os.environ.get("FINOMENA_NEG_CONTROL", "WT+DMSO")   # WT+DMSO


def med(x):
    return statistics.median(x) if x else float("nan")


def mad(x):
    if not x:
        return float("nan")
    m = med(x)
    return 1.4826 * med([abs(v - m) for v in x])


def main():
    os.makedirs(PROD_DIR, exist_ok=True)

    # bio_plate -> condition -> animal_id -> [per-frame pixel_diff]
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    with open(AGG_CSV) as fh:
        rdr = csv.DictReader(fh)
        has_bio = "bio_plate" in rdr.fieldnames
        for r in rdr:
            bp = r["bio_plate"] if has_bio else r["plate"].split("_Replicate")[0]
            data[bp][r["Condition"]][r["animal_id"]].append(float(r["pixel_diff"]))

    def animal_means(plate, cond, log=False):
        out = []
        for pix in data[plate][cond].values():
            mu = statistics.mean(pix)
            out.append(math.log2(mu + 1) if log else mu)
        return out

    def ssmd2(numer, dd, dc):  # two-MAD robust SSMD
        dd = 0 if math.isnan(dd) else dd
        dc = 0 if math.isnan(dc) else dc
        den = math.sqrt(dd * dd + dc * dc)
        return numer / den if den else float("nan")

    def glass(numer, dc):      # control-MAD only (robust Glass's delta)
        return numer / dc if dc and not math.isnan(dc) else float("nan")

    # Read-only model estimates (Drug_effect / Rescue) to copy alongside.
    def load_model(plate):
        cand = os.path.join(OUTPUT_DIR, plate, "master_results.csv")
        if not os.path.exists(cand):
            cand = os.path.join(OUTPUT_DIR, "master_results.csv")
        d = {}
        if not os.path.exists(cand):
            return d
        with open(cand) as fh:
            for r in csv.DictReader(fh):
                m = re.match(r"(KO\+\S+) -", r["Passed_Contrasts"])
                if not m:
                    continue
                d.setdefault(m.group(1), {})[r["Test_Family"]] = r
        return d

    cols = [
        "Plate", "drug", "n_animals",
        "mean_drug_raw", "mean_KO_raw", "mean_WT_raw",
        "logmean_drug", "logmean_KO", "logmean_WT",
        "MAD_drug_raw", "MAD_KO_raw", "MAD_WT_raw",
        "MAD_drug_log", "MAD_KO_log", "MAD_WT_log",
        "SSMD_vsKO_raw", "SSMD_vsWT_raw", "SSMD_vsKO_log", "SSMD_vsWT_log",
        "Glass_vsKO_raw", "Glass_vsWT_raw", "Glass_vsKO_log", "Glass_vsWT_log",
        "model_DrugEffect_vsKO_log2", "model_DrugEffect_adjp",
        "model_Rescue_vsWT_log2", "model_Rescue_adjp",
    ]

    rows = []
    for bp in PLATES:
        if bp not in data:
            print(f"  WARNING: {bp} absent from aggregated data; skipping.")
            continue
        mod = load_model(bp)
        ko_raw = animal_means(bp, POS); wt_raw = animal_means(bp, NEG)
        ko_log = animal_means(bp, POS, True); wt_log = animal_means(bp, NEG, True)
        mad_ko_raw, mad_wt_raw = mad(ko_raw), mad(wt_raw)
        mad_ko_log, mad_wt_log = mad(ko_log), mad(wt_log)
        med_ko_raw, med_wt_raw = med(ko_raw), med(wt_raw)
        med_ko_log, med_wt_log = med(ko_log), med(wt_log)
        for cond in sorted(data[bp]):
            if cond in (POS, NEG):
                continue
            d_raw = animal_means(bp, cond); d_log = animal_means(bp, cond, True)
            md_raw, md_log = med(d_raw), med(d_log)
            mad_d_raw, mad_d_log = mad(d_raw), mad(d_log)
            r = {
                "Plate": bp, "drug": cond, "n_animals": len(d_raw),
                "mean_drug_raw": md_raw, "mean_KO_raw": med_ko_raw,
                "mean_WT_raw": med_wt_raw,
                "logmean_drug": md_log, "logmean_KO": med_ko_log,
                "logmean_WT": med_wt_log,
                "MAD_drug_raw": mad_d_raw, "MAD_KO_raw": mad_ko_raw,
                "MAD_WT_raw": mad_wt_raw,
                "MAD_drug_log": mad_d_log, "MAD_KO_log": mad_ko_log,
                "MAD_WT_log": mad_wt_log,
                "SSMD_vsKO_raw": ssmd2(md_raw - med_ko_raw, mad_d_raw, mad_ko_raw),
                "SSMD_vsWT_raw": ssmd2(md_raw - med_wt_raw, mad_d_raw, mad_wt_raw),
                "SSMD_vsKO_log": ssmd2(md_log - med_ko_log, mad_d_log, mad_ko_log),
                "SSMD_vsWT_log": ssmd2(md_log - med_wt_log, mad_d_log, mad_wt_log),
                "Glass_vsKO_raw": glass(md_raw - med_ko_raw, mad_ko_raw),
                "Glass_vsWT_raw": glass(md_raw - med_wt_raw, mad_wt_raw),
                "Glass_vsKO_log": glass(md_log - med_ko_log, mad_ko_log),
                "Glass_vsWT_log": glass(md_log - med_wt_log, mad_wt_log),
            }
            de = mod.get(cond, {}).get("Drug_effect")
            rs = mod.get(cond, {}).get("Rescue")
            r["model_DrugEffect_vsKO_log2"] = float(de["estimate"]) if de else None
            r["model_DrugEffect_adjp"] = float(de["adjusted_pvalue"]) if de else None
            r["model_Rescue_vsWT_log2"] = float(rs["estimate"]) if rs else None
            r["model_Rescue_adjp"] = float(rs["adjusted_pvalue"]) if rs else None
            rows.append(r)

    out = os.path.join(PROD_DIR, "SSMD_vs_model_ALL.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c) for c in cols})
    print(f"Saved: {out}  ({len(rows)} drug rows, {len(cols)} columns)")


if __name__ == "__main__":
    main()
