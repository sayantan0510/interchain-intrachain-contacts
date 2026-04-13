#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from typing import Dict, List, Tuple

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Compute intra-chain contacts from a GROMACS (.xtc or .trr) trajectory, "
            "excluding near-neighbor residue pairs, write a global .dat file, "
            "optionally write chain-wise .csv, and make plots."
             "this code assumes backbone bead names as BB and sidechains as SC1, SC2, SC3, ..."
            "Please change the bead name(s) directly inside the code (lines 107-112) if otherwise"
        )
    )

    p.add_argument("-s", "--structure", required=True, help="Topology file (.tpr, .gro, ...)")
    p.add_argument("-f", "--trajectory", required=True, help="Trajectory file (.xtc, .trr, ...)")
    p.add_argument(
        "--out-prefix",
        required=True,
        help="Prefix for all output files, e.g. intrachain_bb"
    )

    p.add_argument(
        "--contact-mode",
        choices=["bb", "bbsc"],
        default="bb",
        help="bb = backbone only (BB); bbsc = backbone + sidechains (BB + SC*)"
    )

    p.add_argument(
        "--chain-by",
        choices=["molnum", "fragment", "auto", "chainID", "segid"],
        default="molnum",
        help="How to identify chains. For GROMACS proteins, molnum is usually best."
    )

    p.add_argument(
        "--r0-nm",
        type=float,
        required=True,
        help="PLUMED R_0 in nm"
    )
    p.add_argument(
        "--nn",
        type=int,
        default=6,
        help="PLUMED NN exponent (default: 6)"
    )
    p.add_argument(
        "--mm",
        type=int,
        default=12,
        help="PLUMED MM exponent (default: 12)"
    )

    p.add_argument(
        "--exclude-neighbors",
        type=int,
        default=3,
        help=(
            "Exclude atom pairs whose residues are within this sequence separation. "
            "Default: 3, meaning pairs with |i-j| <= 3 are ignored."
        )
    )

    p.add_argument("--tmin-ps", type=float, default=None, help="Minimum time in ps to include")
    p.add_argument("--tmax-ps", type=float, default=None, help="Maximum time in ps to include")
    p.add_argument("--stride", type=int, default=1, help="Analyze every Nth selected frame")

    p.add_argument(
        "--write-chainwise-csv",
        action="store_true",
        help="Also write chain-wise intra-chain contacts as a CSV"
    )
    p.add_argument(
        "--plot-chainwise-hists",
        action="store_true",
        help="Also plot one histogram panel per chain"
    )

    p.add_argument("--bins", type=int, default=50, help="Number of histogram bins (default: 50)")
    p.add_argument(
        "--backend",
        choices=["serial", "OpenMP", "distopia"],
        default="serial",
        help="MDAnalysis distance backend (default: serial)"
    )
    p.add_argument("--verbose", action="store_true", help="Print progress information")

    return p.parse_args()


def build_atom_selection(contact_mode: str) -> str:
    if contact_mode == "bb":
        return "name BB"
    if contact_mode == "bbsc":
        return "name BB or name SC*"
    raise ValueError(f"Unknown contact mode: {contact_mode}")


def plumed_rational_switch(distances: np.ndarray, r0: float, nn: int = 6, mm: int = 12) -> np.ndarray:
    """
    PLUMED-style rational switching function with d0 = 0:

        s(r) = (1 - (r/r0)^nn) / (1 - (r/r0)^mm)

    The removable singularity at r/r0 = 1 is handled with the limit nn/mm.
    """
    x = distances / r0
    numerator = 1.0 - np.power(x, nn)
    denominator = 1.0 - np.power(x, mm)

    out = np.empty_like(distances, dtype=np.float64)
    singular = np.isclose(denominator, 0.0, atol=1e-14, rtol=1e-12)
    nonsingular = ~singular

    out[nonsingular] = numerator[nonsingular] / denominator[nonsingular]
    out[singular] = float(nn) / float(mm)

    return out


def get_chain_labels(ag, mode: str) -> np.ndarray:
    if len(ag) == 0:
        raise ValueError("Selection returned zero atoms.")

    def nonempty_unique(arr) -> int:
        vals = [str(x).strip() for x in arr]
        vals = [x for x in vals if x not in ("", "None")]
        return len(set(vals))

    if mode == "auto":
        if hasattr(ag.atoms, "chainIDs"):
            arr = np.asarray(ag.atoms.chainIDs, dtype=object)
            if nonempty_unique(arr) > 1:
                return np.asarray([str(x) for x in arr], dtype=object)

        if hasattr(ag.atoms, "segids"):
            arr = np.asarray(ag.atoms.segids, dtype=object)
            if nonempty_unique(arr) > 1:
                return np.asarray([str(x) for x in arr], dtype=object)

        if hasattr(ag.atoms, "molnums"):
            arr = np.asarray(ag.atoms.molnums, dtype=object)
            if len(set(arr.tolist())) > 1:
                return np.asarray([str(x) for x in arr], dtype=object)

        mode = "fragment"

    if mode == "chainID":
        if not hasattr(ag.atoms, "chainIDs"):
            raise ValueError("Topology does not contain chainIDs.")
        return np.asarray([str(x) for x in ag.atoms.chainIDs], dtype=object)

    if mode == "segid":
        if not hasattr(ag.atoms, "segids"):
            raise ValueError("Topology does not contain segids.")
        return np.asarray([str(x) for x in ag.atoms.segids], dtype=object)

    if mode == "molnum":
        if not hasattr(ag.atoms, "molnums"):
            raise ValueError("Topology does not contain molnums.")
        return np.asarray([str(x) for x in ag.atoms.molnums], dtype=object)

    if mode == "fragment":
        labels = []
        frag_map: Dict[int, int] = {}
        next_id = 1
        for atom in ag.atoms:
            frag = atom.fragment
            key = int(frag.atoms.indices[0])
            if key not in frag_map:
                frag_map[key] = next_id
                next_id += 1
            labels.append(f"frag{frag_map[key]}")
        return np.asarray(labels, dtype=object)

    raise ValueError(f"Unknown chain definition mode: {mode}")


def build_chain_index_map(labels: np.ndarray) -> Tuple[List[str], List[np.ndarray]]:
    ordered = []
    seen = set()
    for lab in labels:
        if lab not in seen:
            seen.add(lab)
            ordered.append(lab)
    index_arrays = [np.where(labels == lab)[0] for lab in ordered]
    return ordered, index_arrays


def ensure_parent_dir(path_prefix: str) -> None:
    parent = os.path.dirname(os.path.abspath(path_prefix))
    if parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)


def residue_ranks_for_chain(chain_atoms) -> np.ndarray:
    """
    Map each selected atom in a chain to a residue rank 0..Nres-1 within that chain.
    This is safer than using raw resid values directly.
    """
    atom_resindices = np.asarray(chain_atoms.resindices, dtype=int)
    unique_res = []
    seen = set()
    for r in atom_resindices:
        if r not in seen:
            seen.add(r)
            unique_res.append(r)

    rank_map = {residx: rank for rank, residx in enumerate(unique_res)}
    return np.asarray([rank_map[r] for r in atom_resindices], dtype=int)


def intrachain_contact_for_chain(
    chain_atoms,
    residue_ranks: np.ndarray,
    exclusion: int,
    box: np.ndarray,
    r0: float,
    nn: int,
    mm: int,
    backend: str,
) -> float:
    """
    Compute intra-chain contact sum for one chain, excluding pairs with
    residue separation <= exclusion.
    """
    n_atoms = len(chain_atoms)
    if n_atoms < 2:
        return 0.0

    pos = chain_atoms.positions
    dmat = distance_array(pos, pos, box=box, backend=backend)

    iu, ju = np.triu_indices(n_atoms, k=1)
    residue_sep = np.abs(residue_ranks[iu] - residue_ranks[ju])

    mask = residue_sep > exclusion
    if not np.any(mask):
        return 0.0

    distances = dmat[iu[mask], ju[mask]]
    return float(plumed_rational_switch(distances, r0=r0, nn=nn, mm=mm).sum())


def plot_global_timeseries(times_ps: np.ndarray, values: np.ndarray, out_png: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(times_ps, values)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Average intra-chain contacts")
    ax.set_title("Average intra-chain contacts vs time")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_global_hist(values: np.ndarray, bins: int, out_png: str) -> None:
    mean_val = float(np.mean(values))
    std_val = float(np.std(values))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(values, bins=bins)
    ax.axvline(mean_val, linestyle="--", color="orange", linewidth=2)

    textstr = f"Mean = {mean_val:.3f}\nStd = {std_val:.3f}"
    ax.text(
        0.98, 0.98,
        textstr,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7)
    )

    ax.set_xlabel("Average intra-chain contacts")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of average intra-chain contacts")
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def plot_chainwise_hists(chainwise_values: np.ndarray, chain_names: List[str], bins: int, out_png: str) -> None:
    n_chains = chainwise_values.shape[1]
    ncols = 4 if n_chains >= 4 else n_chains
    nrows = math.ceil(n_chains / ncols)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4 * ncols, 3 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for i in range(n_chains):
        ax = axes[i]
        vals = chainwise_values[:, i]
        mean_val = float(np.mean(vals))
        std_val = float(np.std(vals))

        ax.hist(vals, bins=bins)
        ax.axvline(mean_val, linestyle="--", color="orange", linewidth=2)

        textstr = f"Mean = {mean_val:.3f}\nStd = {std_val:.3f}"
        ax.text(
            0.98, 0.98,
            textstr,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7)
        )

        ax.set_title(str(chain_names[i]))
        ax.set_xlabel("Intra-chain contacts")
        ax.set_ylabel("Count")

    for j in range(n_chains, len(axes)):
        axes[j].axis("off")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def main() -> int:
    args = parse_args()

    if args.r0_nm <= 0:
        raise ValueError("R_0 must be > 0.")
    if args.mm <= args.nn:
        print("Warning: MM is usually chosen greater than NN.", file=sys.stderr)
    if args.stride <= 0:
        raise ValueError("--stride must be >= 1")
    if args.exclude_neighbors < 0:
        raise ValueError("--exclude-neighbors must be >= 0")

    ensure_parent_dir(args.out_prefix)

    selection = build_atom_selection(args.contact_mode)
    r0 = args.r0_nm * 10.0  # nm -> Angstrom

    u = mda.Universe(args.structure, args.trajectory)
    sel = u.select_atoms(selection)

    if len(sel) == 0:
        raise ValueError(f'No atoms matched selection: "{selection}"')

    chain_labels = get_chain_labels(sel, args.chain_by)
    chain_names, chain_indices = build_chain_index_map(chain_labels)
    n_chains = len(chain_names)

    if n_chains < 1:
        raise ValueError("No chains were detected.")

    # Pre-build per-chain atom groups and residue-rank arrays
    chain_atomgroups = []
    chain_residue_ranks = []

    for idx in chain_indices:
        chain_ag = sel[idx]
        chain_atomgroups.append(chain_ag)
        chain_residue_ranks.append(residue_ranks_for_chain(chain_ag))

    if args.verbose:
        print(f"Selection: {selection}", file=sys.stderr)
        print(f"Selected atoms: {len(sel)}", file=sys.stderr)
        print(f"Detected chains: {n_chains}", file=sys.stderr)
        print(f"Chain labels: {', '.join(chain_names)}", file=sys.stderr)
        print(f"Excluding residue pairs with |i-j| <= {args.exclude_neighbors}", file=sys.stderr)

    global_dat = args.out_prefix + "_global.dat"
    global_ts_png = args.out_prefix + "_global_timeseries.png"
    global_hist_png = args.out_prefix + "_global_hist.png"
    chainwise_csv = args.out_prefix + "_chainwise.csv"
    chainwise_hist_png = args.out_prefix + "_chainwise_hists.png"

    times_ps = []
    global_values = []
    chainwise_rows = []

    selected_counter = 0

    with open(global_dat, "w") as fout:
        fout.write("# time_ps avg_intrachain_contacts\n")

        for ts in u.trajectory:
            time_ps = float(ts.time)

            if args.tmin_ps is not None and time_ps < args.tmin_ps:
                continue
            if args.tmax_ps is not None and time_ps > args.tmax_ps:
                break

            if selected_counter % args.stride != 0:
                selected_counter += 1
                continue
            selected_counter += 1

            per_chain_values = np.zeros(n_chains, dtype=np.float64)

            for i in range(n_chains):
                per_chain_values[i] = intrachain_contact_for_chain(
                    chain_atoms=chain_atomgroups[i],
                    residue_ranks=chain_residue_ranks[i],
                    exclusion=args.exclude_neighbors,
                    box=ts.dimensions,
                    r0=r0,
                    nn=args.nn,
                    mm=args.mm,
                    backend=args.backend,
                )

            global_avg = float(np.mean(per_chain_values))

            fout.write(f"{time_ps:.6f} {global_avg:.10f}\n")

            times_ps.append(time_ps)
            global_values.append(global_avg)

            if args.write_chainwise_csv or args.plot_chainwise_hists:
                chainwise_rows.append(per_chain_values.copy())

            if args.verbose and len(times_ps) % 10 == 0:
                print(
                    f"Processed frame {ts.frame:8d}  time = {time_ps:.3f} ps  global_avg = {global_avg:.6f}",
                    file=sys.stderr,
                )

    times_ps = np.asarray(times_ps, dtype=np.float64)
    global_values = np.asarray(global_values, dtype=np.float64)

    if len(times_ps) == 0:
        raise ValueError("No frames were analyzed. Check --tmin-ps, --tmax-ps, and --stride.")

    plot_global_timeseries(times_ps, global_values, global_ts_png)
    plot_global_hist(global_values, args.bins, global_hist_png)

    if args.write_chainwise_csv or args.plot_chainwise_hists:
        chainwise_values = np.asarray(chainwise_rows, dtype=np.float64)

        if args.write_chainwise_csv:
            with open(chainwise_csv, "w", newline="") as fc:
                writer = csv.writer(fc)
                writer.writerow(["time_ps"] + chain_names)
                for i in range(len(times_ps)):
                    writer.writerow(
                        [f"{times_ps[i]:.6f}"] + [f"{x:.10f}" for x in chainwise_values[i]]
                    )

        if args.plot_chainwise_hists:
            plot_chainwise_hists(chainwise_values, chain_names, args.bins, chainwise_hist_png)

    if args.verbose:
        print(f"Wrote: {global_dat}", file=sys.stderr)
        print(f"Wrote: {global_ts_png}", file=sys.stderr)
        print(f"Wrote: {global_hist_png}", file=sys.stderr)
        if args.write_chainwise_csv:
            print(f"Wrote: {chainwise_csv}", file=sys.stderr)
        if args.plot_chainwise_hists:
            print(f"Wrote: {chainwise_hist_png}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
