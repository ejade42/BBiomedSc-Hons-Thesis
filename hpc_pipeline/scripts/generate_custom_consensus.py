## Consensus generator from MSA
## Want to make a decision for ambiguous bases, as standard is to present A/T/C/G-only consensus sequences
## Consensus is not necessarily sensible as it may not reflect any single molecule,
## but we will make one anyway for comparability to the literature

## Will take approach similar to that described in Sone 2019:
## - columns with plurality gap will be discard (Sone discarded only if gap was >50%).
## - most common base will be chosen in other columns, even if <50%.

## Reads from MSA file (MAFFT output)
## Requires that all sequences are the same length (should be the case for an alignment)

## Usage: python script.py -i input_filename -o output_filename [-w wrapping_number]


import argparse
from datetime import datetime


## ------------------------------------------------------------------------------------------------------------------
## Returns the current time, formatted as (YYYY-MM-DD HH:MM:SS)
def now():
    time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return f"({time})"


## Takes FASTA multiple sequence alignment, where each sequence might span multiple lines,
## and creates a dictionary where the sequence IDs are the keys and the concatenated sequences are the items
def parse_msa(msa_lines_list):
    msa_dict = {}

    for i in range(len(msa_lines_list)):
        line = msa_lines_list[i].strip()
        if line[0] == ">":
            current_key = line
            msa_dict[current_key] = ""
        else:
            msa_dict[current_key] += line
    
    if verbose:
        print(f"{now()} Parsed MSA input file")
    return(msa_dict)


## Checks that all sequences in MSA are the same length
def verify_msa_lengths(msa_dict):
    lengths = []
    for id in msa_dict:
        lengths.append(len(msa[id]))
    
    if len(set(lengths)) == 1:
        if verbose:
            print(f"{now()} Verified that all sequences in MSA are same length")
        return(lengths[0])
    else:
        print("Error: All sequences must be the same length")
        exit()


## Find most common element of list - ties broken by first appearance
def mode(a_list):
    return(max(a_list, key = a_list.count))


## Creates consensus by simple plurality
def create_consensus(msa, sequence_length):
    consensus = ""
    for i in range(sequence_length):
        bases = []
        for key in msa:
            bases.append(msa[key][i])

        mode_base = mode(bases)
        if mode_base == "-":
            consensus_base = ""
        elif mode_base.upper() in ["A", "C", "G", "T"]:
            consensus_base = mode_base.upper()
        else:
            print("Error: alignment contains a non A/C/G/T/- base")
            exit()
        
        consensus += consensus_base

        if verbose:
            if i % 500 == 0:
                print(f"{now()} Found consensus base for position {i}/{sequence_length} ({round(i/sequence_length*100,2)}% complete)")
    
    return(consensus)
## ------------------------------------------------------------------------------------------------------------------


## Arguments

parser = argparse.ArgumentParser(description = "Plurality consensus generator (from multiple sequence alignment)")
parser.add_argument("-i", "--input", help = "input filename (MSA from MAFFT or similar)", metavar = "filename", required = True)
parser.add_argument("-o", "--output", help = "output filename (FASTA)", metavar = "filename", required = True)
parser.add_argument("-w", "--wrapping", help = "length of lines when printing to FASTA (default: 75)", type = int, metavar = "int", default = 75)
parser.add_argument("-v", "--verbose", help = "print progress/diagnostic information", action = "store_true", default = False)
args = vars(parser.parse_args())

input_filename  = args["input"]
output_filename = args["output"]
output_wrapping = args["wrapping"]
verbose = args["verbose"]
    

## Read in MSA
if verbose: print(f"{now()} Reading in MSA")
input_file = open(input_filename, "r")
input_msa  = input_file.readlines()
input_file.close


## Parse MSA file and generate consensus
if verbose: print(f"{now()} Parsing MSA")
msa = parse_msa(input_msa)
sequence_length = verify_msa_lengths(msa)
consensus = create_consensus(msa, sequence_length)


## Write output, with wrapping specified (e.g. at 75)
if verbose: print(f"{now()} Writing output, wrapping at {output_wrapping}")
output_file = open(output_filename, "w")
output_file.write(">Consensus\n")

num_full_lines   = len(consensus) // output_wrapping
length_last_line = len(consensus)  % output_wrapping
for i in range(num_full_lines):
    output_file.write(consensus[i*output_wrapping : (i+1)*output_wrapping] + "\n")
if length_last_line > 0:
    output_file.write(consensus[num_full_lines*output_wrapping : ] + "\n")
output_file.close()

if verbose: print(f"{now()} Done")

