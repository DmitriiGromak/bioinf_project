from Bio import SeqIO

fasta = "GCF_024542735.1_iyBomHunt1.1_genomic.fna"
zhunt_tsv = "BomHunt_zhunt_unfiltered.tsv"
out_bed   = "zhunt_filtered.bed"
chrom_offsets = []
cumulative = 0
for rec in SeqIO.parse(fasta, "fasta"):
    chrom_len = len(rec.seq)
    chrom_offsets.append((cumulative, cumulative + chrom_len, rec.id))
    cumulative += chrom_len

def abs_to_bed(abs_start, abs_end):
    for cum_start, cum_end, chrom_id in chrom_offsets:
        if cum_start <= abs_start < cum_end:
            return chrom_id, abs_start - cum_start, abs_end - cum_start
    return None, None, None

written = 0
with open(zhunt_tsv) as fin, open(out_bed, "w") as fout:
    for line in fin:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        abs_start, abs_end = int(parts[0]), int(parts[1])
        chrom, start, end = abs_to_bed(abs_start, abs_end)
        if chrom:
            fout.write(f"{chrom}\t{start}\t{end}\n")
            written += 1