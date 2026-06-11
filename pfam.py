import re, sys, gzip, shutil, urllib.request
from pathlib import Path
from collections import defaultdict
import pyhmmer
import pandas as pd

FAA_PATH    = "protein.faa"
GFF_PATH    = "genomic.gff"
OUT_HITS    = "epigenetic_hits.tsv"
OUT_SUMMARY = "epigenetic_summary.tsv"
OUT_TABLE   = "epigenetic_table.csv"
OUT_MD      = "epigenetic_table.md"
EVALUE      = 1e-5
CACHE_DIR   = Path("pfam_cache")
SHEET_CSV_URL = (
    "https://docs.google.com/spreadsheets/d/"
    "1d3114LoBXHQACG9iqVDZtQ5pmf-rUYtbHhpEBF6q_Qc"
    "/export?format=csv&gid=1397071643"
)
COL_GENE, COL_FUNCTION, COL_PFAM = 0, 2, 7


def build_prot2locus(gff_path):
    gene_loc = {}
    rna2gene = {}
    prot2rna = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) < 9:
                continue
            feat  = p[2]
            attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
            if feat == "gene":
                gid = attrs.get("ID", "")
                if gid:
                    gene_loc[gid] = (p[0], int(p[3]) - 1, int(p[4]), p[6])
            elif feat in ("mRNA", "transcript"):
                rid = attrs.get("ID", "")
                par = attrs.get("Parent", "")
                if rid and par:
                    rna2gene[rid] = par
            elif feat == "CDS":
                pid = attrs.get("protein_id", "")
                par = attrs.get("Parent", "")
                if pid and par and pid not in prot2rna:
                    prot2rna[pid] = par
    prot2locus = {}
    for pid, rna_id in prot2rna.items():
        gene_id = rna2gene.get(rna_id, "")
        locus   = gene_loc.get(gene_id)
        if locus:
            prot2locus[pid] = locus
    return prot2locus


def fetch_pfam_map():
    import csv, io, requests
    resp = requests.get(SHEET_CSV_URL, timeout=30)
    resp.raise_for_status()
    pf_re    = re.compile(r"(\S+)\s+(PF\d{5})")
    pfam_map = defaultdict(lambda: {"domain_name": "", "functions": set(), "genes": set()})
    reader   = csv.reader(io.StringIO(resp.text))
    next(reader, None)
    for row in reader:
        if len(row) <= max(COL_GENE, COL_FUNCTION, COL_PFAM):
            continue
        gene  = row[COL_GENE].strip()
        func  = row[COL_FUNCTION].strip()
        pfams = row[COL_PFAM].strip()
        if not gene or not pfams:
            continue
        for m in pf_re.finditer(pfams):
            dn, acc = m.group(1), m.group(2)
            pfam_map[acc]["domain_name"] = dn
            pfam_map[acc]["functions"].add(func)
            pfam_map[acc]["genes"].add(gene)
    return dict(pfam_map)


def download_hmm(acc):
    CACHE_DIR.mkdir(exist_ok=True)
    hmm_path = CACHE_DIR / f"{acc}.hmm"
    if hmm_path.exists():
        return hmm_path
    url = f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{acc}/?annotation=hmm"
    gz  = CACHE_DIR / f"{acc}.hmm.gz"
    try:
        urllib.request.urlretrieve(url, gz)
        with gzip.open(gz, "rb") as fi, open(hmm_path, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        gz.unlink(missing_ok=True)
        return hmm_path
    except Exception as e:
        gz.unlink(missing_ok=True)
        return None


def decode(b):
    return b.decode() if isinstance(b, bytes) else str(b)

prot2locus = build_prot2locus(GFF_PATH)
pfam_map = fetch_pfam_map()
hmm_paths = {}
for acc in sorted(pfam_map):
    p = download_hmm(acc)
    if p:
        hmm_paths[acc] = p
profiles = []
for acc, path in hmm_paths.items():
    try:
        with pyhmmer.plan7.HMMFile(str(path)) as hf:
            for hmm in hf:
                hmm.accession = acc.encode()
                profiles.append(hmm)
    except Exception as e:
        print(f"  WARNING: {acc}: {e}")
with pyhmmer.easel.SequenceFile(FAA_PATH, digital=True) as sf:
    sequences = list(sf)
rows = []
for top_hits in pyhmmer.hmmscan(sequences, profiles, cpus=0, E=EVALUE):
    query_id = decode(top_hits.query.name)
    for hit in top_hits:
        if hit.evalue > EVALUE:
            continue
        hmm_acc = decode(hit.accession or hit.name)  # PFxxxxx
        info    = pfam_map.get(hmm_acc, {})
        for dom in hit.domains.included:
            rows.append({
                "protein_id":   query_id,
                "pfam_acc":     hmm_acc,
                "domain_name":  info.get("domain_name", ""),
                "functions":    "; ".join(sorted(info.get("functions", []))),
                "genes_human":  "; ".join(sorted(info.get("genes",     []))[:10]),
                "e_value":      hit.evalue,
                "score":        hit.score,
                "dom_i_evalue": dom.i_evalue,
                "hmm_from":     dom.alignment.hmm_from,
                "hmm_to":       dom.alignment.hmm_to,
                "seq_from":     dom.alignment.target_from,
                "seq_to":       dom.alignment.target_to,
            })

df = pd.DataFrame(rows)
df.to_csv(OUT_HITS, sep="\t", index=False)
summary = (df.sort_values("e_value")
                .groupby(["protein_id", "pfam_acc"], as_index=False)
                .first())
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)
prot2name = {}
with open(GFF_PATH) as f:
    gene_names = {}
    rna2gene_local = {}
    prot2rna_local = {}
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9: continue
        attrs = dict(x.split("=",1) for x in p[8].split(";") if "=" in x)
        if p[2] == "gene":
            gid = attrs.get("ID","")
            if gid:
                gene_names[gid] = attrs.get("Name", attrs.get("gene", gid))
        elif p[2] in ("mRNA","transcript"):
            rid, par = attrs.get("ID",""), attrs.get("Parent","")
            if rid and par: rna2gene_local[rid] = par
        elif p[2] == "CDS":
            pid, par = attrs.get("protein_id",""), attrs.get("Parent","")
            if pid and par and pid not in prot2rna_local:
                prot2rna_local[pid] = par
    for pid, rna in prot2rna_local.items():
        gid = rna2gene_local.get(rna,"")
        if gid: prot2name[pid] = gene_names.get(gid, gid)

out_rows = []
for _, r in summary.iterrows():
    locus  = prot2locus.get(r["protein_id"])
    coords = f"{locus[0]}:{locus[1]}-{locus[2]}({locus[3]})" if locus else "N/A"
    out_rows.append({
        "Family":      r["domain_name"] or r["pfam_acc"],
        "Gene":        prot2name.get(r["protein_id"], r["protein_id"]),
        "Coordinates": coords,
    })
out_df = pd.DataFrame(out_rows).drop_duplicates()
out_df.to_csv(OUT_TABLE, index=False)
cols = list(out_df.columns)
md   = "| " + " | ".join(cols) + " |\n"
md  += "| " + " | ".join(["---"] * len(cols)) + " |\n"
for row in out_df.itertuples(index=False):
    md += "| " + " | ".join(str(v) for v in row) + " |\n"
with open(OUT_MD, "w") as f:
    f.write("## Epigenetic gene families\n\n" + md)
