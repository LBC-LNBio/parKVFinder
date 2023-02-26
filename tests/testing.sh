#!/bin/bash
cd "$(dirname "$0")"

# Standard mode
../parKVFinder ../input/1FMO.pdb

# Change probe out
../parKVFinder ../input/1FMO.pdb -o 8.0

# Ligand mode
../parKVFinder ../input/1FMO.pdb -L ../input/ligs_1FMO.pdb

# Custom box mode
../parKVFinder ../input/1FMO.pdb -B --custom_box ../input/1FMO.box.KVFinder.in

# Residue box mode
../parKVFinder ../input/1FMO.pdb -B --residues_box ../input/1FMO.residues.KVFinder.in 

# Run from parameters.toml
../parKVFinder
