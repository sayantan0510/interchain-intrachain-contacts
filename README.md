# Interchain and Intrachain Contact Analysis  

---

## Overview

This repository provides two Python tools to compute **interchain** and **intrachain** contacts in multi-chain biomolecular systems from GROMACS trajectories.
Contacts are defined using a **PLUMED-style switching function**, which provides a smooth and physically meaningful measure of contacts instead of a hard (Heaviside step function) cutoff.

The scripts support:

- Interchain contact analysis (between different chains)
- Intrachain contact analysis (within the same chain)
- Exclusion of short-range bonded neighbors (for intrachain only)
- Time-trajectory analysis
- Statistical distributions (histograms)
- Backbone-only or (backbone + sidechain) contacts

---

## Scientific Definition

### Switching Function

Contacts are computed using:

s(r) = (1 - (r/r0)^n) / (1 - (r/r0)^m)

Default parameters:
- n = 6
- m = 12

---

### Interchain Contacts

C(t) = (1 / Nchains) * sum_{i<j} C_ij(t)

---

### Intrachain Contacts

C_i(t) = sum over (a < b in chain i, |i - j| > k) of s(r_ab)

where:
- k is the residue exclusion (default = 3)

---

## Installation

Clone the repository:

- git clone https://github.com/sayantan0510/interchain-intrachain-contacts.git

- cd interchain-intrachain-contacts

- pip install -r requirements.txt

---

## Quick Start

Run the provided example:

- cd example

- ./run.sh

---

## Usage and Help options

# INTER-chain

python scripts/inter_chain_contact.py --help

usage: inter_chain_contact.py [-h] -s STRUCTURE -f TRAJECTORY --out-prefix OUT_PREFIX [--contact-mode {bb,bbsc}] [--chain-by {molnum,fragment,auto,chainID,segid}] --r0-nm R0_NM [--nn NN]
                              [--mm MM] [--tmin-ps TMIN_PS] [--tmax-ps TMAX_PS] [--stride STRIDE] [--write-chainwise-csv] [--plot-chainwise-hists] [--bins BINS]
                              [--backend {serial,OpenMP,distopia}] [--verbose]

Compute INTER-chain contacts from a GROMACS (.xtc or .trr) trajectory, write a global .dat file, optionally write chain-wise .csv, and make plots.this code assumes backbone bead names as
BB and sidechains as SC1, SC2, SC3, ... Please change the bead name(s) directly inside the code (lines 123-128) if otherwise

options: <br>
  -h, --help            show this help message and exit <br>
  -s STRUCTURE, --structure STRUCTURE <br>
                        Topology file (.tpr, .gro, ...) <br>
  -f TRAJECTORY, --trajectory TRAJECTORY 
                        Trajectory file (.xtc, .trr, ...)  <br>
  --out-prefix OUT_PREFIX <br>
                        Prefix for all output files, e.g. interchain_bb <br>
  --contact-mode {bb,bbsc} <br>
                        bb = backbone only (BB); bbsc = backbone + sidechains (BB + SC*) <br>
  --chain-by {molnum,fragment,auto,chainID,segid} <br>
                        How to identify chains. For GROMACS proteins, molnum is usually best.<br>
  --r0-nm R0_NM         R_0 in nm <br>
  --nn NN               n exponent (default: 6) <br>
  --mm MM               m exponent (default: 12) <br>
  --tmin-ps TMIN_PS     Minimum time in ps to include <br>
  --tmax-ps TMAX_PS     Maximum time in ps to include <br>
  --stride STRIDE       Analyze every Nth selected frame <br>
  --write-chainwise-csv
                        Also write chain-wise interchain contacts divided by (Nchains - 1) <br>
  --plot-chainwise-hists
                        Also plot one histogram panel per chain <br>
  --bins BINS           Number of bins for histograms (default: 50) <br>
  --backend {serial,OpenMP,distopia} <br>
                        MDAnalysis distance backend (default: serial) <br>
  --verbose             Print progress information <br>


# INTRA-chain

python scripts/intra_chain_contact.py --help<br>

usage: intra_chain_contact.py [-h] -s STRUCTURE -f TRAJECTORY --out-prefix OUT_PREFIX [--contact-mode {bb,bbsc}] [--chain-by {molnum,fragment,auto,chainID,segid}] --r0-nm R0_NM [--nn NN]
                              [--mm MM] [--exclude-neighbors EXCLUDE_NEIGHBORS] [--tmin-ps TMIN_PS] [--tmax-ps TMAX_PS] [--stride STRIDE] [--write-chainwise-csv] [--plot-chainwise-hists]
                              [--bins BINS] [--backend {serial,OpenMP,distopia}] [--verbose]

Compute intra-chain contacts from a GROMACS (.xtc or .trr) trajectory, excluding near-neighbor residue pairs, write a global .dat file, optionally write chain-wise .csv, and make
plots.this code assumes backbone bead names as BB and sidechains as SC1, SC2, SC3, ...Please change the bead name(s) directly inside the code (lines 107-112) if otherwise

options:
  -h, --help            show this help message and exit<br>
  -s STRUCTURE, --structure STRUCTURE 
                        Topology file (.tpr, .gro, ...)<br>
  -f TRAJECTORY, --trajectory TRAJECTORY 
                        Trajectory file (.xtc, .trr, ...)<br>
  --out-prefix OUT_PREFIX 
                        Prefix for all output files, e.g. intrachain_bb <br>
  --contact-mode {bb,bbsc} 
                        bb = backbone only (BB); bbsc = backbone + sidechains (BB + SC*)<br>
  --chain-by {molnum,fragment,auto,chainID,segid} 
                        How to identify chains. For GROMACS proteins, molnum is usually best.<br>
  --r0-nm R0_NM         PLUMED R_0 in nm <br>
  --nn NN               PLUMED NN exponent (default: 6)<br>
  --mm MM               PLUMED MM exponent (default: 12)<br>
  --exclude-neighbors EXCLUDE_NEIGHBORS 
                        Exclude atom pairs whose residues are within this sequence separation. Default: 3, meaning pairs with |i-j| <= 3 are ignored.<br>
  --tmin-ps TMIN_PS     Minimum time in ps to include <br>
  --tmax-ps TMAX_PS     Maximum time in ps to include <br>
  --stride STRIDE       Analyze every Nth selected frame<br>
  --write-chainwise-csv 
                        Also write chain-wise intra-chain contacts as a CSV<br>
  --plot-chainwise-hists 
                        Also plot one histogram panel per chain<br>
  --bins BINS           Number of histogram bins (default: 50)<br>
  --backend {serial,OpenMP,distopia}<br>
                        MDAnalysis distance backend (default: serial)<br>
  --verbose             Print progress information<br>

---

## License

-This project is licensed under the MIT License (see LICENSE file).

--- 

## Author

Sayantan Mondal<br>
Department of Chemistry and Chemical Biology<br>
Harvard University<br>
smondal@fas.harvard.edu<br>
sayantan0510@gmail.com<br>




