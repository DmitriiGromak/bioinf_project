import argparse
import os
import pickle
import sys
from concurrent.futures import ThreadPoolExecutor
from queue import Queue

import numpy as np
import scipy.ndimage
import torch
import torch.nn.functional as F
from Bio import SeqIO
from tqdm import tqdm
from transformers import BertForTokenClassification, BertTokenizer

def seq2kmer(seq: str, k: int = 6) -> list:
    return [seq[x: x + k] for x in range(len(seq) + 1 - k)]

def split_seq(seq: list, length: int = 512, pad: int = 16) -> list:
    res = []
    for st in range(0, len(seq), length - pad):
        end = min(st + 512, len(seq))
        res.append(seq[st:end])
    return res


def stitch(np_seqs: list, pad: int = 16) -> np.ndarray:
    res = np_seqs[0]
    for seq in np_seqs[1:]:
        res = np.concatenate([res[:-pad], seq])
    return res

def tokenize_pieces(pieces: list, tokenizer) -> list:
    return [
        torch.LongTensor(tokenizer.encode(" ".join(p), add_special_tokens=False))
        for p in pieces
    ]


def predict_chromosome(seq: str, model, tokenizer, device,
                        chunk: int = 512, pad: int = 16,
                        batch_size: int = 64,
                        num_workers: int = 4) -> np.ndarray:
    seq_upper = seq.upper()
    kmer_seq  = seq2kmer(seq_upper, k=6)
    pieces    = split_seq(kmer_seq, length=chunk, pad=pad)
    n_pieces  = len(pieces)
    chunk_size = max(1, n_pieces // num_workers)
    piece_batches = [pieces[i:i+chunk_size] for i in range(0, n_pieces, chunk_size)]
    token_lists: list = [None] * len(piece_batches)
    with ThreadPoolExecutor(max_workers=num_workers) as pool:
        futures = {
            pool.submit(tokenize_pieces, pb, tokenizer): idx
            for idx, pb in enumerate(piece_batches)
        }
        for fut in futures:
            token_lists[futures[fut]] = fut.result()
    all_token_tensors = [t for sublist in token_lists for t in sublist]
    prob_chunks = []
    model.eval()
    with torch.no_grad():
        for batch_start in range(0, n_pieces, batch_size):
            batch_tensors = all_token_tensors[batch_start: batch_start + batch_size]
            max_len = max(t.shape[0] for t in batch_tensors)
            padded     = torch.zeros(len(batch_tensors), max_len, dtype=torch.long)
            attn_mask  = torch.zeros(len(batch_tensors), max_len, dtype=torch.long)
            for i, t in enumerate(batch_tensors):
                padded[i, :t.shape[0]] = t
                attn_mask[i, :t.shape[0]] = 1

            padded    = padded.to(device)
            attn_mask = attn_mask.to(device)
            logits = model(input_ids=padded, attention_mask=attn_mask)[-1]
            probs = F.softmax(logits, dim=-1)[:, :, 1]
            for i in range(probs.shape[0]):
                real_len = int(attn_mask[i].sum().item())
                prob_chunks.append(probs[i, :real_len].cpu().numpy())
    return stitch(prob_chunks, pad=pad)

def call_regions(prob_track: np.ndarray, chrom: str,
                 threshold: float, min_len: int) -> list:
    binary = (prob_track >= threshold).astype(np.int8)
    labeled, n_labels = scipy.ndimage.label(binary)
    regions = []
    for lbl in range(1, n_labels + 1):
        coords = np.where(labeled == lbl)[0]
        if coords.shape[0] >= min_len:
            regions.append((chrom, int(coords[0]), int(coords[-1] + 1)))
    return regions

def split_fasta(fasta_path: str, out_dir: str,
                chromosomes_only: bool = True) -> list:
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if chromosomes_only and not rec.id.startswith("NC_"):
            continue
        out_path = os.path.join(out_dir, f"{rec.id}.fa")
        if not os.path.exists(out_path):
            SeqIO.write(rec, out_path, "fasta")
        written.append((rec.id, out_path))
    return written

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--fasta",       required=True)
    p.add_argument("--model_dir",   required=True)
    p.add_argument("--output",      default="zdna_regions.bed")
    p.add_argument("--threshold",   type=float, default=0.5)
    p.add_argument("--min_len",     type=int,   default=10)
    p.add_argument("--batch_size",  type=int,   default=64,
                   help="Chunks per GPU forward pass. Lower if OOM.")
    p.add_argument("--num_workers", type=int,   default=4,
                   help="CPU threads for parallel tokenization.")
    p.add_argument("--chunk",       type=int,   default=512)
    p.add_argument("--pad",         type=int,   default=16)
    p.add_argument("--all_sequences", action="store_true",
                   help="Include NW_* scaffolds (default: NC_* only)")
    p.add_argument("--save_pkl",    action="store_true")
    p.add_argument("--split_dir",   default=None)
    p.add_argument("--device",      default=None)
    return p.parse_args()


def main():
    args   = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[zdnabert_scan] Device: {device}")

    print(f"[zdnabert_scan] Loading model from: {args.model_dir}")
    tokenizer = BertTokenizer.from_pretrained(args.model_dir)
    model     = BertForTokenClassification.from_pretrained(args.model_dir)
    model.to(device)
    model.eval()
    print(f"[zdnabert_scan] Model loaded ({model.config.num_labels} labels).")

    split_dir = args.split_dir or os.path.join(
        os.path.dirname(os.path.abspath(args.fasta)), "split_chroms"
    )
    chrom_files = split_fasta(args.fasta, split_dir,
                               chromosomes_only=not args.all_sequences)
    print(f"[zdnabert_scan] {len(chrom_files)} sequences to process → {args.output}")
    if not chrom_files:
        print("[zdnabert_scan] ERROR: no sequences found.")
        sys.exit(1)

    all_regions = []
    prob_store  = {} if args.save_pkl else None

    with open(args.output, "w") as bed_fh:
        for seqid, fa_path in tqdm(chrom_files, desc="Chromosomes", unit="chrom"):
            rec     = next(SeqIO.parse(fa_path, "fasta"))
            raw_seq = str(rec.seq)
            n_chunks = max(1, (len(raw_seq) - 5) // (args.chunk - args.pad))
            print(f"\n  {seqid}: {len(raw_seq):,} nt  "
                  f"(~{n_chunks:,} chunks, ~{max(1,n_chunks//args.batch_size):,} GPU batches)")

            prob_track = predict_chromosome(
                raw_seq, model, tokenizer, device,
                chunk=args.chunk, pad=args.pad,
                batch_size=args.batch_size,
                num_workers=args.num_workers,
            )

            if args.save_pkl:
                prob_store[seqid] = prob_track

            regions = call_regions(prob_track, seqid, args.threshold, args.min_len)
            for chrom, start, end in regions:
                bed_fh.write(f"{chrom}\t{start}\t{end}\n")
            all_regions.extend(regions)

    if args.save_pkl:
        pkl_path = args.output.replace(".bed", "_prob_tracks.pkl")
        with open(pkl_path, "wb") as fh:
            pickle.dump(prob_store, fh)


if __name__ == "__main__":
    main()
