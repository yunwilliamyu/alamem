#!/usr/bin/env python3
"""Generate a per-alignment local ANI accuracy scatter plot."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.colors import Normalize


REPO = Path(__file__).resolve().parents[1]
FIGURES = REPO / "scripts" / "figures"
DETAIL_DIR = FIGURES / "local_ani_accuracy_details"
OUTPUT_TSV = FIGURES / "local_ani_accuracy.tsv"
OUTPUT_BALANCED_TSV = FIGURES / "local_ani_accuracy_balanced.tsv"
OUTPUT_PNG = FIGURES / "accuracy_scatterplot.png"
#OUTPUT_SVG = FIGURES / "alamem_local_ani_accuracy.svg"
TARGET_FASTA = REPO / "scripts" / "1mb.fna"
SIMULATOR = REPO / "target" / "release" / "examples" / "ani_simulator"
RUSTUP_HOME = Path("/tmp/alamem-rustup")
RUSTC = RUSTUP_HOME / "toolchains" / "stable-x86_64-unknown-linux-gnu" / "bin" / "rustc"
CARGO = RUSTUP_HOME / "toolchains" / "stable-x86_64-unknown-linux-gnu" / "bin" / "cargo"
CARGO_HOME = Path("/tmp/alamem-cargo")

DEFAULT_ANIS = list(range(90, 101))
DEFAULT_LENGTHS = [50, 75, 100, 150, 250, 400, 650, 1000, 1500, 2000]
BIN_SPECS = [
    ("40-100 bp", 40.0, 100.0),
    ("100-500 bp", 100.0, 500.0),
    ("500-1000 bp", 500.0, 1000.0),
    ("1000+ bp", 1000.0, None),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", type=int, default=500)
    parser.add_argument("--long-bin-reads", type=int)
    parser.add_argument("--anis", type=int, nargs="+", default=DEFAULT_ANIS)
    parser.add_argument("--lengths", type=int, nargs="+", default=DEFAULT_LENGTHS)
    parser.add_argument("--target-per-bin", type=int)
    parser.add_argument("--seed", type=int, default=20260930)
    parser.add_argument("--skip-sim", action="store_true")
    parser.add_argument("--kmer", type=int, default=11)
    return parser.parse_args()


def configure_fonts() -> font_manager.FontProperties:
    plt.rcParams.update(
        {
            "font.family": 'sans-serif',
            "font.sans-serif": ['Helvetica', 'Arial', 'DejaVu Sans', 'Liberation Sans'],
            "font.size": 12,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "grid.linestyle": ":",
            "svg.fonttype": "path",
        }
    )


def build_simulator() -> None:
    env = os.environ.copy()
    env.setdefault("RUSTFLAGS", "-C target-cpu=native")
    cargo = "cargo"
    if CARGO.is_file() and RUSTC.is_file():
        cargo = str(CARGO)
        env["RUSTUP_HOME"] = str(RUSTUP_HOME)
        env["RUSTC"] = str(RUSTC)
        if CARGO_HOME.is_dir():
            env["CARGO_HOME"] = str(CARGO_HOME)
    subprocess.run(
        [cargo, "build", "--example", "ani_simulator", "--release"],
        cwd=REPO,
        env=env,
        check=True,
    )


def run_simulations(
    anis: list[int], lengths: list[int], reads: int, long_bin_reads: int | None
, kmer: int) -> list[Path]:
    DETAIL_DIR.mkdir(parents=True, exist_ok=True)
    detail_paths: list[Path] = []
    for ani in anis:
        for length in lengths:
            read_count = max(reads, long_bin_reads) if long_bin_reads and length >= 2000 else reads
            detail_path = DETAIL_DIR / f"local_ani_ani{ani}_len{length}.tsv"
            detail_paths.append(detail_path)
            print(f"Simulating ANI={ani} length={length} reads={read_count}", flush=True)
            subprocess.run(
                [
                    str(SIMULATOR),
                    str(TARGET_FASTA),
                    str(read_count),
                    str(length),
                    str(ani),
                    "--detail-output",
                    str(detail_path),
                    "-k",
                    str(kmer),
                    "-m",
                    "1000",
                    "--min-ani",
                    "90"
                ],
                cwd=REPO,
                check=True,
                stdout=subprocess.DEVNULL,
            )
    return detail_paths


def combine_details(detail_paths: list[Path]) -> None:
    OUTPUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    with OUTPUT_TSV.open("w", newline="") as out_handle:
        writer = None
        for path in detail_paths:
            with path.open(newline="") as in_handle:
                reader = csv.DictReader(in_handle, delimiter="\t")
                if writer is None:
                    writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames, delimiter="\t")
                if not wrote_header:
                    writer.writeheader()
                    wrote_header = True
                for row in reader:
                    writer.writerow(row)


def read_points() -> dict[str, np.ndarray]:
    columns: dict[str, list[float]] = {
        "true_local_ani": [],
        "estimated_ani": [],
        "local_alignment_len": [],
    }
    with OUTPUT_TSV.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            for key in columns:
                columns[key].append(float(row[key]))
    if not columns["true_local_ani"]:
        raise RuntimeError(f"No rows found in {OUTPUT_TSV}")
    return {key: np.asarray(values, dtype=float) for key, values in columns.items()}


def bin_masks(lengths: np.ndarray) -> list[tuple[str, np.ndarray]]:
    masks = []
    for label, lower, upper in BIN_SPECS:
        if upper is None:
            mask = lengths >= lower
        else:
            mask = (lengths >= lower) & (lengths < upper)
        masks.append((label, mask))
    return masks


def balanced_indices(lengths: np.ndarray, target_per_bin: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    selected = []
    for label, mask in bin_masks(lengths):
        indices = np.flatnonzero(mask)
        if len(indices) < target_per_bin:
            raise RuntimeError(
                f"{label} has {len(indices):,} rows, fewer than requested {target_per_bin:,}"
            )
        selected.append(rng.choice(indices, size=target_per_bin, replace=False))
    return np.concatenate(selected)


def subset_points(points: dict[str, np.ndarray], indices: np.ndarray) -> dict[str, np.ndarray]:
    return {key: values[indices] for key, values in points.items()}


def write_balanced_tsv(indices: np.ndarray) -> None:
    selected = set(int(index) for index in indices)
    with OUTPUT_TSV.open(newline="") as in_handle, OUTPUT_BALANCED_TSV.open("w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        for row_index, row in enumerate(reader):
            if row_index in selected:
                writer.writerow(row)


def plot(points: dict[str, np.ndarray], bold_font: font_manager.FontProperties, kmer) -> None:
    true_ani = points["true_local_ani"]
    estimated_ani = points["estimated_ani"]
    lengths = points["local_alignment_len"]
    norm = Normalize(vmin=float(np.min(lengths)), vmax=float(np.max(lengths)))

    lower = 88.0
    upper = 100.5
    ticks = np.arange(88.0, 101.0, 2.0)

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.4), sharex=True, sharey=True)
    axes_flat = axes.ravel()
    scatter = None
    for ax, (title, mask) in zip(axes_flat, bin_masks(lengths)):
        ax.plot([lower, upper], [lower, upper], color="#6E6E6E", linewidth=1.0, alpha=0.65)
        ax.set_xlim(lower, upper)
        ax.set_ylim(lower, upper)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_title(title, fontproperties=bold_font, pad=8)

        if np.any(mask):
            order = np.argsort(lengths[mask])
            x = true_ani[mask][order]
            y = estimated_ani[mask][order]
            c = lengths[mask][order]
            scatter = ax.scatter(
                x,
                y,
                c=c,
                cmap="viridis",
                norm=norm,
                s=12,
                alpha=0.58,
                linewidths=0,
                rasterized=True,
            )
            pearson = float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 else float("nan")
            label = f"Pearson r = {pearson:.3f}\nn = {len(x):,}"
        else:
            label = "Pearson r = NA\nn = 0"
        ax.text(
            0.04,
            0.96,
            label,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=10,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.72,
            },
        )

    fig.suptitle("alamem local ANI accuracy; k=" + str(kmer), fontproperties=bold_font, fontsize=16, y=0.985)
    fig.supxlabel("True ANI from base-level alignment (%)", y=0.04)
    fig.supylabel("alamem estimated ANI (%)", x=0.045)

    if scatter is None:
        raise RuntimeError("No points available for plotting")
    colorbar = fig.colorbar(
        scatter,
        ax=axes_flat.tolist(),
        pad=0.018,
        fraction=0.038,
    )
    colorbar.set_label("Reported local alignment length (bp)")

    fig.subplots_adjust(left=0.11, right=0.88, bottom=0.1, top=0.92, wspace=0.08, hspace=0.16)
    OUTPUT_PNG_MOD = str(OUTPUT_PNG).rstrip('.png') + "-k" + str(kmer) + ".png"
    fig.savefig(OUTPUT_PNG_MOD, dpi=300)
    #fig.savefig(OUTPUT_SVG)
    print(f"Wrote {OUTPUT_TSV}")
    print(f"Wrote {OUTPUT_PNG_MOD}")
    #print(f"Wrote {OUTPUT_SVG}")


def main() -> None:
    args = parse_args()
    bold_font = configure_fonts()
    if not args.skip_sim:
        build_simulator()
        long_bin_reads = args.long_bin_reads
        if args.target_per_bin and long_bin_reads is None:
            long_bin_reads = 3000
        detail_paths = run_simulations(args.anis, args.lengths, args.reads, long_bin_reads, args.kmer)
        combine_details(detail_paths)
    points = read_points()
    if args.target_per_bin:
        indices = balanced_indices(points["local_alignment_len"], args.target_per_bin, args.seed)
        write_balanced_tsv(indices)
        points = subset_points(points, indices)
        print(f"Wrote {OUTPUT_BALANCED_TSV}")
    plot(points, bold_font, kmer=args.kmer)


if __name__ == "__main__":
    main()

