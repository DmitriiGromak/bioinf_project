import pandas as pd
import subprocess

def line_count(f):
    return int(subprocess.check_output(["wc", "-l", f]).split()[0])

zhunt_total   = line_count("zhunt_merged.bed")
zdnabert_total = line_count("zdnabert_merged.bed")
g4_total = line_count("g4_merged.bed")

print(f"Totals — zhunt: {zhunt_total:,}  zdnabert: {zdnabert_total:,}  g4: {g4_total:,}\n")

FEATURES = ["exons", "introns", "promoters", "downstream", "intergenic"]
t1 = pd.read_csv("table1_counts.txt", sep="\t")
t1["G4_frac"]       = (t1["G4_count"]       / g4_total      ).round(4)
t1["Zhunt_frac"]    = (t1["Zhunt_count"]    / zhunt_total   ).round(4)
t1["ZDNABERT_frac"] = (t1["ZDNABERT_count"] / zdnabert_total).round(4)
t1 = t1[["Feature",
          "G4_count", "G4_frac",
          "Zhunt_count", "Zhunt_frac",
          "ZDNABERT_count", "ZDNABERT_frac"]]
t1.columns = ["Feature",
              "G4 count", "G4 fraction",
              "Zhunt count", "Zhunt fraction",
              "ZDNABERT count", "ZDNABERT fraction"]
t1.to_csv("table1_final.csv", index=False)
print("=== TABLE 1: Prediction regions per feature ===")
print(t1.to_string(index=False))
t2 = pd.read_csv("table2_counts.txt", sep="\t")
t2["G4_frac"]       = (t2["G4_with_hit"]       / t2["Total_intervals"]).round(4)
t2["Zhunt_frac"]    = (t2["Zhunt_with_hit"]    / t2["Total_intervals"]).round(4)
t2["ZDNABERT_frac"] = (t2["ZDNABERT_with_hit"] / t2["Total_intervals"]).round(4)
t2 = t2[["Feature", "Total_intervals",
          "G4_with_hit", "G4_frac",
          "Zhunt_with_hit", "Zhunt_frac",
          "ZDNABERT_with_hit", "ZDNABERT_frac"]]
t2.columns = ["Feature", "Total intervals",
              "G4 with hit", "G4 fraction",
              "Zhunt with hit", "Zhunt fraction",
              "ZDNABERT with hit", "ZDNABERT fraction"]
t2.to_csv("table2_final.csv", index=False)
print("\n=== TABLE 2: Feature intervals containing ≥1 prediction ===")
print(t2.to_string(index=False))
def bp_merged(bed):
    total = 0
    with open(bed) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 3:
                total += int(p[2]) - int(p[1])
    return total
genome_bp = sum(int(line.split()[1]) for line in open("genome.sizes"))
bg = {}
for feat in FEATURES:
    bg[feat] = bp_merged(f"{feat}.bed")

print("\n=== Background (fraction of genome bp, merged) ===")
print(f"{'Feature':<15} {'bp':>15} {'fraction':>10}")
for feat in FEATURES:
    print(f"{feat:<15} {bg[feat]:>15,} {bg[feat]/genome_bp:>10.4f}")
print(f"{'genome total':<15} {genome_bp:>15,} {'1.0000':>10}")

bg_df = pd.DataFrame([
    {"Feature": f, "bp": bg[f], "fraction": round(bg[f]/genome_bp, 4)}
    for f in FEATURES
])
bg_df.to_csv("background_final.csv", index=False)
t3 = t1.copy()
t3 = t3.merge(bg_df[["Feature", "fraction"]].rename(columns={"fraction":"bg_fraction"}),
               on="Feature")
for tool, col in [("G4", "G4 fraction"), ("Zhunt", "Zhunt fraction"), ("ZDNABERT", "ZDNABERT fraction")]:
    t3[f"{tool} fold"] = (t3[col] / t3["bg_fraction"]).round(3)
t3 = t3[["Feature", "bg_fraction",
          "G4 fraction", "G4 fold",
          "Zhunt fraction", "Zhunt fold",
          "ZDNABERT fraction", "ZDNABERT fold"]]
t3.columns = ["Feature", "Background fraction",
              "G4 fraction", "G4 fold-enrichment",
              "Zhunt fraction", "Zhunt fold-enrichment",
              "ZDNABERT fraction", "ZDNABERT fold-enrichment"]
t3.to_csv("table3_enrichment.csv", index=False)
print("\n=== TABLE 3: Fold-enrichment vs genomic background ===")
print(t3.to_string(index=False))
import os, re

def parse_attrs(s):
    return dict(x.split("=", 1) for x in s.strip().split(";") if "=" in x)

GFF_PATH = "genomic.gff"
prot2locus = {}
with open(GFF_PATH) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9: continue
        if p[2] != "CDS": continue
        attrs = parse_attrs(p[8])
        pid = attrs.get("protein_id", "")
        if pid:
            prot2locus[pid] = (p[0], int(p[3])-1, int(p[4]), p[6])
def df_to_md(df):
    cols = list(df.columns)
    header = "| " + " | ".join(cols) + " |"
    sep    = "| " + " | ".join(["---"]*len(cols)) + " |"
    rows   = ["| " + " | ".join(str(v) for v in r) + " |"
              for r in df.itertuples(index=False)]
    return "\n".join([header, sep] + rows)

with open("readme_tables.md", "w") as f:
    f.write("## Epigenetic gene families\n\n")
    f.write("\n\n## Distribution of secondary structures across genomic features\n\n")
    f.write(df_to_md(t1))
    f.write("\n\n## Feature intervals with ≥1 prediction\n\n")
    f.write(df_to_md(t2))
    f.write("\n\n## Fold-enrichment vs genomic background\n\n")
    f.write(df_to_md(t3))
