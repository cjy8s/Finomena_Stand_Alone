#!/usr/bin/env python
"""Consolidate the BAM's per-plate output into the wide production tables.

Reconstructed from the original inline consolidation step (run ad-hoc on
2026-06-11, never previously committed). Reads the per-plate output that
initial_screen_script.r writes to <OUTPUT_DIR>/<Plate N>/ and produces:

  ALL_PLATES_full_results.csv   one row per drug per plate (wide): movement
                                (DrugEffect_vsKO_*), landing (Rescue_vsWT_*),
                                and the posterior P(equivalent to WT/KO).
  full_results_<Plate>.csv      the same, split per plate (compact column names).

This is a pure RESHAPE of the BAM's own committed output — it computes no new
statistics, so it does not touch any frozen quantity. The initial-screen
master_results.csv tags families "Genotype" / "Drug_effect" / "Rescue"; only
the per-drug Drug_effect (movement) and Rescue (landing) rows are pivoted here.

Paths are env-overridable (SageMaker mounts), falling back to the local repo
layout so behaviour is unchanged when run by hand.
"""
import csv
import os
import re

# <OUTPUT_DIR>/<Plate N>/{master_results,posterior_equivalence_summary}.csv
OUTPUT_DIR = os.environ.get(
    "FINOMENA_OUTPUT_DIR",
    "Data/Initial screen/Production Output")
# Where the consolidated tables are written.
PROD_DIR = os.environ.get("FINOMENA_ADDON_PROD_DIR", OUTPUT_DIR)
# Biological plate folder names produced by the per-plate loop.
PLATES = os.environ.get("FINOMENA_PLATES", "Plate 1,Plate 2,Plate 3").split(",")


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def plate_dir(plate):
    """Per-plate output folder. The loop names it by bio_plate (e.g. 'Plate 1');
    a single-pass run writes straight to OUTPUT_DIR."""
    cand = os.path.join(OUTPUT_DIR, plate)
    return cand if os.path.isdir(cand) else OUTPUT_DIR


def load_master(plate):
    """drug -> {Test_Family -> row} from a plate's master_results.csv."""
    path = os.path.join(plate_dir(plate), "master_results.csv")
    mr = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            m = re.match(r"(KO\+\S+) -", r["Passed_Contrasts"])
            if not m:
                continue
            mr.setdefault(m.group(1), {})[r["Test_Family"]] = r
    return mr


def load_pe(plate):
    """(lhs, rhs) -> posterior-equivalence row from a plate's PE summary."""
    path = os.path.join(plate_dir(plate), "posterior_equivalence_summary.csv")
    pe = {}
    if not os.path.exists(path):
        return pe
    with open(path) as fh:
        for r in csv.DictReader(fh):
            pe[(r["lhs"], r["rhs"])] = r
    return pe


def main():
    os.makedirs(PROD_DIR, exist_ok=True)
    POS = os.environ.get("FINOMENA_POS_CONTROL", "KO+DMSO")  # KO+DMSO
    NEG = os.environ.get("FINOMENA_NEG_CONTROL", "WT+DMSO")   # WT+DMSO

    all_cols = [
        "Plate", "drug",
        "DrugEffect_vsKO_est_log2", "DrugEffect_vsKO_SE", "DrugEffect_vsKO_z",
        "DrugEffect_vsKO_raw_p", "DrugEffect_vsKO_adj_p_FDR",
        "Rescue_vsWT_est_log2", "Rescue_vsWT_SE", "Rescue_vsWT_z",
        "Rescue_vsWT_raw_p", "Rescue_vsWT_adj_p_FDR",
        "equiv_M_median_vsKO", "equiv_Pr_vsKO",
        "equiv_M_median_vsWT", "equiv_Pr_vsWT",
    ]
    # Compact per-plate schema (matches the original full_results_PlateN.csv).
    compact_cols = ["drug", "DE_est", "DE_SE", "DE_z", "DE_raw_p", "DE_adj_p",
                    "RS_est", "RS_SE", "RS_z", "RS_raw_p", "RS_adj_p",
                    "M_WT", "PrEq_WT", "PrEq_KO"]

    all_rows = []
    for plate in PLATES:
        mr = load_master(plate)
        pe = load_pe(plate)
        plate_rows = []
        for drug, fam in mr.items():
            de = fam.get("Drug_effect")   # movement: drug vs KO+DMSO
            rs = fam.get("Rescue")        # landing:  drug vs WT+DMSO
            row = {"Plate": plate, "drug": drug}
            for f, tag in [(de, "DrugEffect_vsKO"), (rs, "Rescue_vsWT")]:
                if f:
                    est, se = num(f["estimate"]), num(f["SE"])
                    row[f"{tag}_est_log2"] = est
                    row[f"{tag}_SE"] = se
                    row[f"{tag}_z"] = round(est / se, 4) if (est is not None and se) else None
                    row[f"{tag}_raw_p"] = num(f["raw_pvalue"])
                    row[f"{tag}_adj_p_FDR"] = num(f["adjusted_pvalue"])
            pk = pe.get((drug, POS))   # equivalence vs KO
            pw = pe.get((drug, NEG))   # equivalence vs WT
            row["equiv_M_median_vsWT"] = num(pw["M_median_log2"]) if pw else None
            row["equiv_Pr_vsWT"] = num(pw["Pr_equiv"]) if pw else None
            row["equiv_M_median_vsKO"] = num(pk["M_median_log2"]) if pk else None
            row["equiv_Pr_vsKO"] = num(pk["Pr_equiv"]) if pk else None
            all_rows.append(row)
            plate_rows.append(row)

        # Per-plate compact file.
        safe = plate.replace(" ", "")
        with open(os.path.join(PROD_DIR, f"full_results_{safe}.csv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=compact_cols)
            w.writeheader()
            for r in sorted(plate_rows, key=lambda x: x["drug"]):
                w.writerow({
                    "drug": r["drug"],
                    "DE_est": r.get("DrugEffect_vsKO_est_log2"),
                    "DE_SE": r.get("DrugEffect_vsKO_SE"),
                    "DE_z": r.get("DrugEffect_vsKO_z"),
                    "DE_raw_p": r.get("DrugEffect_vsKO_raw_p"),
                    "DE_adj_p": r.get("DrugEffect_vsKO_adj_p_FDR"),
                    "RS_est": r.get("Rescue_vsWT_est_log2"),
                    "RS_SE": r.get("Rescue_vsWT_SE"),
                    "RS_z": r.get("Rescue_vsWT_z"),
                    "RS_raw_p": r.get("Rescue_vsWT_raw_p"),
                    "RS_adj_p": r.get("Rescue_vsWT_adj_p_FDR"),
                    "M_WT": r.get("equiv_M_median_vsWT"),
                    "PrEq_WT": r.get("equiv_Pr_vsWT"),
                    "PrEq_KO": r.get("equiv_Pr_vsKO"),
                })

    out = os.path.join(PROD_DIR, "ALL_PLATES_full_results.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=all_cols)
        w.writeheader()
        for row in sorted(all_rows, key=lambda x: (x["Plate"], x["drug"])):
            w.writerow({c: row.get(c) for c in all_cols})
    print(f"Saved: {out}  ({len(all_rows)} drug rows across {len(PLATES)} plates)")


if __name__ == "__main__":
    main()
