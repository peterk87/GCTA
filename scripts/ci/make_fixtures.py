#!/usr/bin/env python3
"""Build deterministic COJO fixtures under tests/fixtures/.

Also writes hand-tweaked edge cases (allele swap, freq mismatch, collinear, empty region).
"""
from __future__ import annotations

import math
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "sandbox/2026-07-13-manc-cojo/scripts"))
from gen_synth_data import (  # noqa: E402
    draw_block_genotypes,
    marginal_stats,
    pack_bed,
)


FIX = ROOT / "tests" / "fixtures"


def write_plink(prefix: Path, G: np.ndarray, bim_rows, freqs, y, n: int) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(prefix.with_suffix(".ma"), "w") as ma:
        ma.write("SNP A1 A2 freq b se p N\n")
        for j in range(G.shape[1]):
            gj = G[:, j].astype(np.float64)
            b, se, p = marginal_stats(gj, y)
            # clamp extreme p for GCTA stability
            p = max(p, 1e-300)
            freq_a1 = float(gj.mean() / 2.0)
            name = bim_rows[j][1]
            ma.write(f"{name} A G {freq_a1:.6g} {b:.6g} {se:.6g} {p:.6g} {n}\n")
    with open(prefix.with_suffix(".bim"), "w") as bim:
        for c, name, cm, bp, a1, a2 in bim_rows:
            bim.write(f"{c}\t{name}\t{cm}\t{bp}\t{a1}\t{a2}\n")
    with open(prefix.with_suffix(".fam"), "w") as fam:
        for i in range(n):
            fam.write(f"FID{i}\tIID{i}\t0\t0\t0\t-9\n")
    with open(prefix.with_suffix(".bed"), "wb") as bed:
        bed.write(pack_bed(G))


def make_pheno(rng, G, causal_idx, betas, h2=0.5):
    n = G.shape[0]
    Gc = G[:, causal_idx].astype(np.float64)
    Gc_std = (Gc - Gc.mean(0)) / (Gc.std(0) + 1e-12)
    genetic = Gc_std @ betas
    genetic *= np.sqrt(h2 / (genetic.var() + 1e-12))
    noise = rng.standard_normal(n)
    noise *= np.sqrt((1 - h2) / (noise.var() + 1e-12))
    return genetic + noise


def case_tiny_le(seed=1):
    """3 SNPs in LE, 1 causal — joint ≈ marginal."""
    rng = np.random.default_rng(seed)
    n, m = 200, 3
    freqs = np.array([0.2, 0.3, 0.4])
    G = np.zeros((n, m), dtype=np.int8)
    for j, f in enumerate(freqs):
        G[:, j] = rng.binomial(2, f, size=n).astype(np.int8)
    # force near-independence by reshuffling columns relative to each other
    for j in range(m):
        G[:, j] = rng.permutation(G[:, j])
    bim = [(1, f"rs_le_{j}", 0, 1_000_000 + j * 5_000_000, "A", "G") for j in range(m)]
    causal = np.array([1])
    betas = np.array([1.2])
    y = make_pheno(rng, G, causal, betas, h2=0.4)
    out = FIX / "tiny_le" / "data"
    write_plink(out, G, bim, freqs, y, n)
    (FIX / "tiny_le" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 --thread-num 1 --maf 0.01\n"
    )


def case_tiny_ld(seed=2):
    """2 SNPs in tight LD, 1 causal — collinearity / conditional path."""
    rng = np.random.default_rng(seed)
    n = 200
    g, freqs = draw_block_genotypes(rng, n, 2, rho=0.98)
    # place within 1Mb so default cojo-wind covers both
    bim = [
        (1, "rs_ld_0", 0, 10_000_000, "A", "G"),
        (1, "rs_ld_1", 0, 10_050_000, "A", "G"),
    ]
    causal = np.array([0])
    betas = np.array([1.5])
    y = make_pheno(rng, g, causal, betas, h2=0.5)
    out = FIX / "tiny_ld" / "data"
    write_plink(out, g, bim, freqs, y, n)
    (FIX / "tiny_ld" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 --thread-num 1 --maf 0.01\n"
    )


def case_small_multi(seed=3):
    """Multi-signal stepwise across a few LD blocks spanning >2Mb (region tests)."""
    rng = np.random.default_rng(seed)
    n = 400
    n_blocks, snps_per = 4, 40
    all_g, bim, freqs_all = [], [], []
    causal_idx, causal_beta = [], []
    snp_i = 0
    bp = 5_000_000
    for blk in range(n_blocks):
        g, freqs = draw_block_genotypes(rng, n, snps_per, rho=0.85)
        pos = np.sort(rng.integers(bp, bp + 400_000, size=snps_per))
        bp += 600_000
        caus = rng.choice(snps_per, size=2, replace=False)
        for cl in caus:
            causal_idx.append(snp_i + int(cl))
            causal_beta.append(rng.normal(0, 1.0))
        for jj in range(snps_per):
            bim.append((1, f"rs_sm_{snp_i+jj}", 0, int(pos[jj]), "A", "G"))
        all_g.append(g)
        freqs_all.append(freqs)
        snp_i += snps_per
    G = np.concatenate(all_g, axis=1)
    y = make_pheno(rng, G, np.array(causal_idx), np.array(causal_beta), h2=0.55)
    out = FIX / "small_multi" / "data"
    write_plink(out, G, bim, np.concatenate(freqs_all), y, n)
    (FIX / "small_multi" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 --thread-num 1 --maf 0.01\n"
    )
    # region mid-point; wind 1000Kb (GCTA minimum for --extract-region-bp)
    mid = (bim[0][3] + bim[len(bim) // 2][3]) // 2
    rdir = FIX / "region_extract"
    rdir.mkdir(parents=True, exist_ok=True)
    for ext in (".bed", ".bim", ".fam", ".ma"):
        shutil.copy2(out.with_suffix(ext), rdir / f"data{ext}")
    (rdir / "cmd").write_text(
        f"--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 "
        f"--thread-num 1 --maf 0.01 --extract-region-bp 1 {mid} 1000\n"
    )


def case_edge_collinear(seed=4):
    """Two identical SNPs — only one should survive collinearity."""
    rng = np.random.default_rng(seed)
    n = 200
    g0 = rng.binomial(2, 0.25, size=n).astype(np.int8)
    G = np.column_stack([g0, g0.copy()])
    bim = [
        (1, "rs_col_0", 0, 20_000_000, "A", "G"),
        (1, "rs_col_1", 0, 20_001_000, "A", "G"),
    ]
    y = make_pheno(rng, G, np.array([0]), np.array([2.0]), h2=0.6)
    out = FIX / "edge_collinear" / "data"
    write_plink(out, G, bim, G.mean(0) / 2, y, n)
    (FIX / "edge_collinear" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 "
        "--cojo-collinear 0.9 --thread-num 1 --maf 0.01\n"
    )


def case_collinear_reject(seed=9):
    """Two independent signals plus a near-copy of A that must not enter the joint model.

    A and B are strong and nearly independent; C ≈ A (~0.92 corr). Stepwise should
    select A and B only. Useful as a multi-insert regression alongside small_multi.
    """
    rng = np.random.default_rng(seed)
    n = 400
    a = rng.binomial(2, 0.30, size=n).astype(np.int8)
    b = rng.binomial(2, 0.35, size=n).astype(np.int8)  # independent of A
    # C: near-copy of A (~0.92 corr) via redrawing ~8% of calls independently.
    c = a.copy()
    flip = rng.random(n) < 0.08
    c[flip] = rng.binomial(2, 0.30, size=int(flip.sum())).astype(np.int8)
    G = np.column_stack([a, b, c]).astype(np.int8)
    bim = [
        (1, "rs_cr_A", 0, 25_000_000, "A", "G"),
        (1, "rs_cr_B", 0, 25_400_000, "A", "G"),
        (1, "rs_cr_C", 0, 25_050_000, "A", "G"),  # near A, inside default window
    ]
    # A,B,C all strong and near-equal so C stays the top conditional candidate after
    # A is picked (C≈0.93·A) and is REJECTED on insert rather than dropped at the
    # conditional stage. A kept marginally strongest so selection order is stable.
    y = make_pheno(rng, G, np.array([0, 1, 2]), np.array([1.25, 1.15, 1.10]), h2=0.6)
    out = FIX / "collinear_reject" / "data"
    write_plink(out, G, bim, G.mean(0) / 2, y, n)
    (FIX / "collinear_reject" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 "
        "--cojo-collinear 0.9 --thread-num 1 --maf 0.01\n"
    )


def case_edge_allele(seed=5):
    """A1/A2 swapped in .ma vs .bim → .badsnps."""
    rng = np.random.default_rng(seed)
    n = 150
    G = rng.binomial(2, 0.3, size=(n, 1)).astype(np.int8)
    bim = [(1, "rs_bad_allele", 0, 30_000_000, "A", "G")]
    y = make_pheno(rng, G, np.array([0]), np.array([1.0]), h2=0.3)
    out = FIX / "edge_allele" / "data"
    write_plink(out, G, bim, G.mean(0) / 2, y, n)
    # rewrite .ma with swapped alleles that don't match bim
    ma_path = out.with_suffix(".ma")
    lines = ma_path.read_text().splitlines()
    hdr, row = lines[0], lines[1].split()
    # Force alleles that match neither order expected... actually GCTA accepts A1==allele2 flip.
    # Use alleles C/T that match neither A/G.
    row[1], row[2] = "C", "T"
    ma_path.write_text(hdr + "\n" + " ".join(row) + "\n")
    (FIX / "edge_allele" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 --thread-num 1\n"
    )
    (FIX / "edge_allele" / "expect_exit").write_text("1\n")  # likely aborts: no matchable SNPs


def case_edge_freqdiff(seed=6):
    """Large freq gap between sumstat and genotype → freq.badsnps."""
    rng = np.random.default_rng(seed)
    n = 150
    G = rng.binomial(2, 0.05, size=(n, 1)).astype(np.int8)
    bim = [(1, "rs_freq", 0, 40_000_000, "A", "G")]
    y = make_pheno(rng, G, np.array([0]), np.array([1.0]), h2=0.2)
    out = FIX / "edge_freqdiff" / "data"
    write_plink(out, G, bim, G.mean(0) / 2, y, n)
    ma_path = out.with_suffix(".ma")
    lines = ma_path.read_text().splitlines()
    hdr, row = lines[0], lines[1].split()
    row[3] = "0.45"  # far from ~0.05 geno freq
    ma_path.write_text(hdr + "\n" + " ".join(row) + "\n")
    (FIX / "edge_freqdiff" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --cojo-p 5e-8 --thread-num 1\n"
    )
    (FIX / "edge_freqdiff" / "expect_exit").write_text("1\n")


def case_edge_region_empty(seed=7):
    """Region with no SNPs — expect exit failure."""
    rng = np.random.default_rng(seed)
    n = 100
    G = rng.binomial(2, 0.2, size=(n, 2)).astype(np.int8)
    bim = [
        (1, "rs_e0", 0, 50_000_000, "A", "G"),
        (1, "rs_e1", 0, 50_100_000, "A", "G"),
    ]
    y = make_pheno(rng, G, np.array([0]), np.array([1.0]), h2=0.2)
    out = FIX / "edge_region_empty" / "data"
    write_plink(out, G, bim, G.mean(0) / 2, y, n)
    # Extract far away (≥1000kb wind around distant bp)
    (FIX / "edge_region_empty" / "cmd").write_text(
        "--bfile data --cojo-file data.ma --cojo-slct --thread-num 1 "
        "--extract-region-bp 1 1 1000\n"
    )
    (FIX / "edge_region_empty" / "expect_exit").write_text("1\n")


def case_forward_and_topn(seed=8):
    """Reuse small_multi inputs with alternate cmds."""
    src = FIX / "small_multi"
    for name, cmd in [
        (
            "forward_slct",
            "--bfile data --cojo-file data.ma --cojo-forward --cojo-p 0.01 "
            "--thread-num 1 --maf 0.01\n",
        ),
        (
            "top_n",
            "--bfile data --cojo-file data.ma --cojo-slct --cojo-top-SNPs 3 "
            "--thread-num 1 --maf 0.01\n",
        ),
        (
            "joint_only",
            "--bfile data --cojo-file data.ma --cojo-joint --thread-num 1 --maf 0.01\n",
        ),
    ]:
        d = FIX / name
        d.mkdir(parents=True, exist_ok=True)
        for ext in (".bed", ".bim", ".fam", ".ma"):
            shutil.copy2(src / f"data{ext}", d / f"data{ext}")
        (d / "cmd").write_text(cmd)


def main():
    FIX.mkdir(parents=True, exist_ok=True)
    case_tiny_le()
    case_tiny_ld()
    case_small_multi()
    case_edge_collinear()
    case_edge_allele()
    case_edge_freqdiff()
    case_edge_region_empty()
    case_forward_and_topn()
    case_collinear_reject()
    print(f"Fixtures written under {FIX}")
    for p in sorted(FIX.iterdir()):
        if p.is_dir() and (p / "cmd").exists():
            print(f"  - {p.name}")


if __name__ == "__main__":
    main()
