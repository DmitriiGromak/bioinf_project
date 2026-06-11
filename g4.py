import re
from Bio import SeqIO
from Bio.Seq import Seq
PATTERN = re.compile(r'(?:G{3,5}[ATGC]{1,7}){3,}G{3,5}', re.IGNORECASE)
fasta = "GCF_024542735.1_iyBomHunt1.1_genomic.fna"
out   = "g4_predictions.bed"
written = 0
with open(out, "w") as fh:
    for rec in SeqIO.parse(fasta, "fasta"):
        seq = str(rec.seq).upper()
        L   = len(seq)
        for m in PATTERN.finditer(seq):
            fh.write(f"{rec.id}\t{m.start()}\t{m.end()}\t.\t40\t+\n")
            written += 1
        rc = str(Seq(seq).reverse_complement())
        for m in PATTERN.finditer(rc):
            fh.write(f"{rec.id}\t{L - m.end()}\t{L - m.start()}\t.\t40\t-\n")
            written += 1