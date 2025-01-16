import pandas as pd
import csv
import os
import sys

from HeadTailAssign.assigner import Assigner
from HeadTailAssign.gamess import GAMESS
from HeadTailAssign.directory import Directory
from HeadTailAssign.extractor import Extractor
from HeadTailAssign.logger import Logger

directory = Directory()
directory.starting(disable=False)

#Generation of input data
print("########## GENERATING INPUT DATA #############")
df = directory.get_data(args=sys.argv, print_df=True)

name_dir = 'output-directory'
if os.path.isdir(name_dir) == False:
    directory.create_output_dir(name_dir)
    print('')
    print(f'Results will be uploaded at {name_dir}')
else:
    print('')
    print(f'Results will be uploaded at {name_dir}')

# Open logger
f = open(f'{name_dir}/log.out', 'w', encoding="utf-8")
original = sys.stdout
sys.stdout = Logger(sys.stdout, f)

## Define monomer from reaction smiles
print("########## DEFINE MONOMERS FOR ANALYSIS #############\n")
assign = Assigner()
if 'reaction' in df:
    cleanSmiles = assign.clean_atom_mapping(df)
    df['reaction_clean'] = cleanSmiles

    print("User did not provide monomers. Algorithm is searching...")
    monomers = assign.find_monomer(df)
    df['monomers'] = monomers
    print(df)

elif 'monomers_list' in df:
    print("Algorithm is searching for reactants...")
    monomers = assign.get_reactants(df)
    df['monomers'] = monomers
    print(df)

# Create polymer id
createName = assign.create_name(df)
df = createName

# Define frontier orbitals
print("\n########## GENERATING GAMESS INPUT #############\n")
gamess = GAMESS()
gamess.generate_gamess_input(df, name_dir, 'scf', 'sto3g', 'uff', 5000)

# Run GAMESS

print("\n########## RUNNING GAMESS #############\n")
gamess.run_gamess_for_files(name_dir)

# Extraction of head and tail positions

extraction = Extractor()

## Extraction of Mulliken Atomic Overlap Populations from .out file
print("########## EXTRACTION OF MULLIKEN ATOMIC OVERLAP POPULATIONS #############")
population = extraction.op_extractor(name_dir)

## Extraction of the Rx value organized by higher to lower value
print("########## EXTRACTION OF THE RX VALUE #############")
rx = extraction.rx_extractor(name_dir)

## Define only the head and tail data from the Rx value
print("########## DEFINE HEAD AND TAIL #############")

## Get monomer classes
classes = assign.get_class(df, name_dir)
df['classes'] = classes

## Get monomer polymerization mechanism
mechanism = assign.get_mechanism(df, name_dir)
df['mechanism'] = mechanism

## Get head and tail
head_tail = assign.get_head_tail(df, name_dir, notation='usual', structure='oligomer')
df['head_tail'] = head_tail
print(df)

# Save results
print("\nThe results.csv is saved on this directory.\nThanks for using HTA!")
df.to_csv('results.csv', sep=",", index=False, quoting=csv.QUOTE_ALL)