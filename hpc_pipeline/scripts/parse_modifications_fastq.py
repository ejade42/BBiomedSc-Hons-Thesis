## Mod-BAM FASTQ parser
## BAM files store modification/methylation information in the MM and ML tags
## Using samtools fastq -T MM,ML, these tags can be written to the header row of the fastq, tab-separated
## This script will parse the fastq and create a new similar file containing the sequence, modifications, and modification probabilities.
## Expects all modification types (e.g. "C+h?") to be present in all reads, even if empty


import argparse
from datetime import datetime



## GENERAL HELPER FUNCTIONS
## ------------------------------------------------------------------------------------------------------------------
## Prints a completion message if verbose is set to True, and returns i+1
def print_completion_if_verbose(verb, i, total_number, print_interval, verbose, items):
    if verbose:
        if i % print_interval == 0:
            print(f"{verb} {i} of {total_number} {items} ({round(i/total_number*100,2)}% complete)")
    return i+1

## Returns the current time, formatted as (YYYY-MM-DD HH:MM:SS)
def now():
    time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return f"({time})"
## ------------------------------------------------------------------------------------------------------------------




## FUNCTIONS TO READ IN FILE AND CREATE DICTIONARY OF READ INFORMATION
## ------------------------------------------------------------------------------------------------------------------
## Takes a list of all lines in a fastq file (header-sequence-+-quality) and unpacks them into separate lists
def unpack_fastq(all_lines):
    all_headers   = []
    all_sequences = []
    all_qualities = []

    for i in range(len(all_lines) // 4):
        all_headers.append(all_lines[4*i])
        all_sequences.append(all_lines[4*i + 1])
        all_qualities.append(all_lines[4*i + 3])

    return all_headers, all_sequences, all_qualities


## Takes a list of all FASTQ headers, with the MM and ML tags, and unpacks into three lists of readnames, MM tags, and ML tags
## Should work with additional tags provided MM and ML are first (in that order), but may not
def unpack_headers(all_headers):
    all_readnames = []
    all_MMs = []
    all_MLs = []

    for header in all_headers:
        header_list = header.split("\t")
        all_readnames.append(header_list[0])
        all_MMs.append(header_list[1])
        all_MLs.append(header_list[2])

    return all_readnames, all_MMs, all_MLs


## Takes lists of readnames, sequences, qualities, MMs, and MLs and packages them into a dictionary of readname:readinformation items
## Each item is itself a dictionary, with "sequence", "quality", "MM", and "ML" items
## e.g. reads_dict["readname_x"]["sequence"] will retrieve the sequence associated with read "readname_x"
def create_reads_dict(all_readnames, all_sequences, all_qualities, all_MMs, all_MLs):
    reads_dict = {}
    for i in range(len(all_readnames)):
        reads_dict[all_readnames[i]] = {"read_id":  all_readnames[i],
                                        "sequence": all_sequences[i].replace("\n", ""),
                                        "quality":  all_qualities[i].replace("\n", ""),
                                        "MM": all_MMs[i].replace("\n", ""),
                                        "ML": all_MLs[i].replace("\n", "")}
    return reads_dict


## Packages the previous three functions together to create the reads dictionary from input
def create_reads_dict_from_all_lines(all_lines):
    all_headers, all_sequences, all_qualities = unpack_fastq(all_lines)
    all_readnames, all_MMs, all_MLs = unpack_headers(all_headers)
    reads_dict = create_reads_dict(all_readnames, all_sequences, all_qualities, all_MMs, all_MLs)
    return reads_dict
## ------------------------------------------------------------------------------------------------------------------




## FUNCTIONS TO PARSE READ INFORMATION AND RETURN MODIFICATION PROBABILITIES/POSITIONS
## ------------------------------------------------------------------------------------------------------------------
## Take a list of strings, and converts each item to an integer
## But, if the list is [""], returns an empty list
def convert_list_to_int(my_list):
    if len(my_list) == 1 and my_list[0] == "":
        my_list = []
    else:
        for i in range(len(my_list)):
            my_list[i] = int(my_list[i])
    return my_list


## Takes the "MM" value associated with a read (a string "MM:Z:C+h?,num,num,...,num;C+m?,num,num,...,num;")
## Returns a dictionary with each modification type (e.g. "C+h?") as the keys and the corresponding list of numbers as the items
def create_MM_dict(MM_contents):
    MM_contents = MM_contents[5: ]   # Remove the "MM:Z:" at the start of each MM entry
    MM_list = MM_contents.split(";")[ :-1]   # Split each different modification into its own list entry, removing the empty last entry (as there is a ; after the last modification)

    MM_dict = {}
    MM_types = []
    for modification_info in MM_list:
        modification_list = modification_info.split(",")
        modification_type = modification_list[0]
        modification_contents = modification_list[1: ]
        MM_types.append(modification_type)
        
        MM_dict[modification_type] = convert_list_to_int(modification_contents)

    return MM_dict, MM_types


## Takes the "ML" value associated with a read (a string "MM:B:C,num,num,...,num")
## These are the scaled probabilities for each index in each tag, concatenated together
    # Scaled to an integer 0-255. Each integer value represents the probability space from n/256 to (n+1)/256
## E.g. with "MM:Z:C+h?,2,4;C+m?,2,4;" you might have "MM:B:C,40,36,204,195" meaning the 3rd C base (after 2 Falses) has a 40/256 probability
##      of being hydroxymethylated, a 204/256 probability of being methylated, and the remainder probability of being unmodified.
##      Likewise, the 8th C base (after 2 Falses, True, 4 Falses) has a 36/256 prob of HM-ation and a 195/256 prob of M-ation.
## This function returns a dictionary of the list of probabilities belonging to each methylation, matching up in format with the MM_dict function
## IMPORTANT: modification_types must be passed in the correct order (hence a separate list is used, rather than retrieving keys from MM_dict)
def create_ML_dict(ML_contents, MM_dict, modification_types):
    ML_contents = convert_list_to_int(ML_contents.split(",")[1: ])  ## Split into list, and remove the tag at the front. Everything else should be numbers

    lengths_per_modification_type = []
    for modification_type in modification_types:
        length_this_type = len(MM_dict[modification_type])
        lengths_per_modification_type.append(length_this_type)

    ML_dict = {}
    for i in range(len(modification_types)):
        lengths_list = lengths_per_modification_type[0:i+1]
        prior_length = sum(lengths_list[0:-1])
        total_length = sum(lengths_list)
        ML_dict[modification_types[i]] = ML_contents[prior_length:total_length]

    contents_length = len(ML_contents)  ## Checks that lengths match up
    dict_lengths = []
    for item in ML_dict.keys():
        dict_lengths.append(len(ML_dict[item]))
    if contents_length != sum(dict_lengths):
        print("Error: ML list has not been correctly divided into separate lists for each modification type")
        exit

    return ML_dict


## Takes a sequence and a single base to match, returns a list of all indices where that base is located
def find_matching_indices(sequence, base):
    indicies = []
    for i in range(len(sequence)):
        if sequence[i] == base:
            indicies.append(i)
    return indicies


## Takes a modification series (list of numbers of unmodified/uncalled bases between called bases, per the SAM MM specification)
## Returns a Boolean list listing True for every modified/called base and False for the base in between (the numbers of Falses are given by the modification series)
def expand_modification_series(modification_series, length = 0):
    indices_to_pick = []
    for num in modification_series:
        for i in range(num):
            indices_to_pick.append(False)
        indices_to_pick.append(True)

    if length != 0 and len(indices_to_pick) < length:    ## Can optionally add Falses onto the end to make sure the list reaches a desired length
        length_to_add = length - len(indices_to_pick)
        for i in range(length_to_add):
            indices_to_pick.append(False)
            
    return indices_to_pick


## Takes a list of indices and a True/False list of at least equal length, and returns the indicies that match up with a True in the selection list
def find_modified_called_indices(indices, selection_list):
    modified_or_called_indices = []
    for i in range(len(indices)):
        if selection_list[i] == True:
            modified_or_called_indices.append(indices[i])
    return(modified_or_called_indices)


## Takes input sequence, modification types, and MM dictionary
## Produces a dictionary of called indices (modificaiton type: index_list)
def create_modification_dicts(sequence, modification_types, MM_dict):
    called_indices_dict = {}
    for modification_type in modification_types:
        base_of_interest    = modification_type[0]
        indices_of_interest = find_matching_indices(sequence, base_of_interest)
        indices_to_pick     = expand_modification_series(MM_dict[modification_type], len(indices_of_interest))
        called_indices_dict[modification_type] = find_modified_called_indices(indices_of_interest, indices_to_pick)

    return called_indices_dict


## Takes the dictionary containing all read information for a single read
## Adds list of modification types, and lists of locations and probabilities for each type, to read info dictionary
def parse_modifications(read_information_original):
    read_information = read_information_original.copy()
    sequence = read_information["sequence"]
    MM_dict, modification_types = create_MM_dict(read_information["MM"])
    called_indices_dict = create_modification_dicts(sequence, modification_types, MM_dict)
    ML_dict = create_ML_dict(read_information["ML"], MM_dict, modification_types)

    read_information["modification_types"] = modification_types
    for modification_type in modification_types:
        read_information[f"{modification_type}_locations"] = called_indices_dict[modification_type]
        read_information[f"{modification_type}_probabilities"] = ML_dict[modification_type]

    return read_information
## ------------------------------------------------------------------------------------------------------------------



## FUNCTIONS TO WRITE VARIOUS OUTPUT FILES
## ------------------------------------------------------------------------------------------------------------------

## Takes reads dict (where each key is a read id and each item is a dictionary containing read information)
## Finds all the unique names of keys (i.e. info types) across all reads
def find_all_keys(reads_dict):
    all_keys = []
    for read in reads_dict:
        read_information = reads_dict[read]
        for key in read_information:
            if key not in all_keys:
                all_keys.append(key)
    return all_keys


## Outputs read information as a csv, one row per read
## Contains read id, sequence, quality, and location and probability lists for each modification type
def output_read_information_csv(reads_dict, output_filename):
    all_keys = find_all_keys(reads_dict)
    output_file = open(output_filename, "w")
    output_columns = [x for x in all_keys if x not in ["MM", "ML"]]   ## exclude MM and ML from output as they are parsed to create more useful outputs
    output_file.write('"' + '","'.join(output_columns) + '"\n')       ## prints csv header row, with "" around each item to sanitise

    i = 1
    for read_id in reads_dict:
        read_information = reads_dict[read_id]
        read_output = []
        for colname in output_columns:
            if colname in read_information:
                colname_data = read_information[colname]
                if type(colname_data) == str:
                    read_output.append(colname_data)
                elif type(colname_data) == list:
                    read_output.append(",".join(map(str, colname_data)))
                else:
                    print("Error: unexpected data type in read_information dictionary")
                    exit
            else:
                read_output.append("")
        output_file.write('"' + '","'.join(read_output) + '"\n')
        i = print_completion_if_verbose("Processed", i, len(reads_dict), 100, verbose, "reads")
    
    output_file.close()


## Takes a probability stored as an integer from 0 to 255, where each integer value represents the probability space from n/256 to (n+1)/256
## Returns an integer from 0 to 9, where each integer value represents the probability space from n/10 to n+1/10
## i.e. 0 <= x < 0.1 maps to 0, 0.1 <= x < 0.2 maps to 1, ... 0.9 <= x < 1 maps to 9 (1 is impossible as the max integer is 255).
## accomplished by taking the first digit of the decimal probability and returning as an integer
def convert_prob_256_to_prob_9(probability_256):
    probability_decimal = round(probability_256/256, 4)
    return int(str(probability_decimal)[2])


## Takes a sequence, a letter to represent modifications, a list of all modified/called indices, a list of corresponding probabilities, and optionally a blank character
## Returns a string where uncalled indicies are the blank character and called indices are the specified letter,
##     and a string where uncalled indicies are the blank character and probability of modification for called indices is represented as 0-9.
def make_line_of_dots(sequence, letter, modification_locations, modification_probabilities, blank_pos = "-", blank_prob = " "):
    position_string = ""
    probability_string = ""
    for i in range(len(sequence)):
        if i in modification_locations:
            position_string += letter
            index = modification_locations.index(i)
            probability_string += str(convert_prob_256_to_prob_9(modification_probabilities[index]))
        else:
            position_string += blank_pos
            probability_string += blank_prob
    return position_string, probability_string


## Outputs modification information in a fastq-like file, with a header line per read, then a sequence line,
## followed by a dashed line for each modification type representing bases where modifications are called
## and a line of numbers underneath representing the probability of that modification from low (0) to high (9)
def output_seqs_with_modifications_illustrated(reads_dict, output_filename, blank_pos, blank_prob):
    output_file = open(output_filename, "w")
    
    i = 1
    for read_id in reads_dict:
        read_information = reads_dict[read_id]
        sequence = read_information["sequence"]
        modification_types = read_information["modification_types"]

        output_file.write(read_id + "\n")
        output_file.write(sequence + "\n")

        for modification_type in modification_types:
            letter = modification_type[2].upper()
            modification_locations = read_information[f"{modification_type}_locations"]
            modification_probabilities = read_information[f"{modification_type}_probabilities"]
            position_string, probability_string = make_line_of_dots(sequence, letter, modification_locations, modification_probabilities, blank_pos, blank_prob)
            output_file.write(position_string + "\n")
            output_file.write(probability_string + "\n")
        i = print_completion_if_verbose("Processed", i, len(reads_dict), 100, verbose, "reads")

    output_file.close()
## ------------------------------------------------------------------------------------------------------------------


def main():
    ## Arguments
    parser = argparse.ArgumentParser(description = "DNA modification parser (from modified fastq)")
    parser.add_argument("-i", "--input", help = "input filename (fastq with modifications in header). required", metavar = "filename", required = True)
    parser.add_argument("-o-r", "--output-read-info", help = "output filename for read information csv. if provided, file will be generated", metavar = "filename", default = None)
    parser.add_argument("-o-m", "--output-modification-view", help = "output filename for modification visualiser txt. if provided, file will be generated", metavar = "filename", default = None)
    parser.add_argument("--blank-pos", help = "change blank character for position output lines of modification visualiser (default = '-')", metavar = "character", default = "-")
    parser.add_argument("--blank-prob", help = "change blank character for probability output lines of modification visualiser (default = ' ')", metavar = "character", default = " ")
    parser.add_argument("-v", "--verbose", help = "print progress/diagnostic information", action = "store_true", default = False)
    args = vars(parser.parse_args())
    
    global verbose
    verbose = args["verbose"]


    ## Open input file
    input_filename = args["input"]
    input_file = open(input_filename, "r")
    reads_dict = create_reads_dict_from_all_lines(input_file.readlines())
    input_file.close


    ## Process input file
    if verbose: print(f"{now()} Parsing input file")
    i = 1
    for read in reads_dict:
        reads_dict[read] = parse_modifications(reads_dict[read])
        i = print_completion_if_verbose("Processed", i, len(reads_dict), 100, verbose, "reads")
        
    
    ## If requested, output read information csv
    output_read_info_filename = args["output_read_info"]
    if output_read_info_filename != None:
        if verbose: print(f"\n{now()} Creating read information csv output file")
        output_read_information_csv(reads_dict, output_read_info_filename)
    elif verbose:
        print(f"\n{now()} Skipping read information csv as no output filename provided")


    ## If requested, output modification visualiser txt
    output_modification_filename = args["output_modification_view"]
    blank_pos = args["blank_pos"]
    blank_prob = args["blank_prob"]
    if output_modification_filename != None:
        if verbose: print(f"\n{now()} Creating modification visualiser txt output file")
        output_seqs_with_modifications_illustrated(reads_dict, output_modification_filename, blank_pos, blank_prob)
    elif verbose:
        print(f"\n{now()} Skipping modification visualiser as no output filename provided")


    print(f"\n{now()} Done")


main()
