# Interchain and Intrachain Contact Analysis (HOMOTYPIC System) 

---

## Overview

This repository provides two Python tools to compute **interchain** and **intrachain** contacts in multi-chain **homotypic** biomolecular systems from GROMACS trajectories (.xtc and .trr).
Contacts are defined using a **PLUMED-style switching function**, which provides a smooth and physically meaningful measure of contacts instead of a hard (Heaviside step function) cutoff.

### Limitation
This particular code is limited to homotypic systems. That means, if you have 100 chains of kind-A, the code can give you inter- and intra-chain contacts. However, if the system is heterotypic, for example, 100 chains of kind-A and 50 chains of kind-B, it cannot give you the interchain contacts between these two kinds. The heterotypic code is available at https://github.com/sayantan0510/interchain-contact-heterotypic

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

$s(r) = \frac{(1 - (r/r0)^n}{(1 - (r/r0)^m}$

Default parameters:
- n = 6
- m = 12

---

### Interchain Contacts

$C(t) = \frac{1}{N_{Chains}} * \sum_{i<j} C_{ij}(t)$

---

### Intrachain Contacts

$C_i(t) = \sum s(r_{ab})$, (a<b in chain i, |i - j| > k)

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

### INTER-chain

python scripts/inter_chain_contact.py --help

Compute INTER-chain contacts from a GROMACS (.xtc or .trr) trajectory, write a global .dat file, optionally write chain-wise .csv, and make plots.this code assumes backbone bead names as
BB and sidechains as SC1, SC2, SC3, ... Please change the bead name(s) directly inside the code (lines 123-128) if otherwise

### INTRA-chain

python scripts/intra_chain_contact.py --help<br>

Compute intra-chain contacts from a GROMACS (.xtc or .trr) trajectory, excluding near-neighbor residue pairs, write a global .dat file, optionally write chain-wise .csv, and make
plots.this code assumes backbone bead names as BB and sidechains as SC1, SC2, SC3, ...Please change the bead name(s) directly inside the code (lines 107-112) if otherwise

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




