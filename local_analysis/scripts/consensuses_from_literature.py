## FUNCTIONS
##----------------------------------------------------------------------------------------------------------
# Function to print a long string wrapped (linebreak) at intervals specified by argument, indendented by a tab
from colorama import Fore, Style, Back
import argparse

# Takes a string, and prints it with each character colour-coded
def print_line_with_colour(string_to_print):
    colours =  {"A": Fore.YELLOW,
                "C": Fore.GREEN,
                "G": Fore.BLUE,
                "T": Fore.RED}
    reset = Fore.RESET
    list_to_print = list(string_to_print)
    for character in list_to_print:
        if character in colours.keys():
            print(colours[character] + character + reset, end = "")
        else:
            print(character, end = "")
    print()

# Takes a long string and a line length, and prints the string wrapping at the line length
# Starts the lines indented by 4, with each line numbered
def print_with_wrap(string, line_length, with_colour = False):
    num_full_lines   = len(string) // line_length
    length_last_line = len(string)  % line_length
    for i in range(num_full_lines):
        number = str(i+1)
        n_spaces = 2 - len(number)
        print(" "*n_spaces + number + " "*2, end = "")
        if with_colour == False:
            print(string[i*line_length:(i+1)*line_length])
        elif with_colour == True:
            print_line_with_colour(string[i*line_length:(i+1)*line_length])

    if length_last_line > 0:
        number = str(num_full_lines+1)
        n_spaces = 2 - len(number)
        print(" "*n_spaces + number + " "*2, end = "")
        if with_colour == False:
            print(string[num_full_lines*line_length:])
        elif with_colour == True:
            print_line_with_colour(string[num_full_lines*line_length:])

# Function to print relevant info about a sequence, given relevant wrapping interval
def print_sequence_info(sequence, line_length, sequence_name, with_colour = True):
    if silent == False:
        print(f"Sequence name:\n    {sequence_name}")
        print(f"Length of sequence:\n    {len(sequence)} bases")
        print("Sequence:")
        print_with_wrap(sequence, line_length, with_colour)
        print()


## Reverse complement function
def reverse_sequence(sequence):
    reversed_sequence = ""
    for i in range(len(sequence)):
        base = sequence[-1-i].upper()

        if base == "A":
            new_base = "T"
        elif base == "C":
            new_base = "G"
        elif base == "G":
            new_base = "C"
        elif base == "T":
            new_base = "A"
        else:
            print("Error: can only reverse sequences consisting of A/C/G/T")
        reversed_sequence += new_base
    return(reversed_sequence)
        

## Take arguments to decide whether to be silent
## Arguments
parser = argparse.ArgumentParser(description = "Printing consensuses from literature to csv, and optionally to terminal")
parser.add_argument("-s", "--silent", help = "do not print output to terminal", action = "store_true", default = False)
args = vars(parser.parse_args())
global silent
silent = args["silent"]


all_sequences = {}
##----------------------------------------------------------------------------------------------------------






## 2. Sone et al, 2019
##----------------------------------------------------------------------------------------------------------
## MAIN PAPER
# F1-1: length should be 531, wraps at 75
sone_2019_f1_1 = "GGC"*131 + "GGA"*3 + "GGC"*3 + ("GGA"*3 + "GGC"*2)*7 + "GGC"*5
print_sequence_info(sone_2019_f1_1, 75, "2. Sone et al (2019) F1-1")
all_sequences["sone_2019_f1_1"] = sone_2019_f1_1

# F2-2: length should be 597, wraps at 72
sone_2019_f2_2 = "GGC"*52 + ("GGA"*2 + "GGC"*4)*23 + "GGA"*2 + "GGC"*3 + "GGA"*2 + "GGC"*2
print_sequence_info(sone_2019_f2_2, 72, "2. Sone et al (2019) F2-2")
all_sequences["sone_2019_f2_2"] = sone_2019_f2_2

# F9-1: length should be 430, wraps at 75
sone_2019_f9_1 = "GGC"*39 + "C" + "GGC"*100 + "GGA"*2 + "GGC"*2
print_sequence_info(sone_2019_f9_1, 75, "2. Sone et al (2019) F9-1")
all_sequences["sone_2019_f9_1"] = sone_2019_f9_1


## SUPPLEMENT
# F1-8: length should be 450, wraps at 75
sone_2019_f1_8 = "GGC"*85 + ("GGA"*3 + "GGC"*2)*12 + "GGC"*5
print_sequence_info(sone_2019_f1_8, 75, "2. Sone et al (2019) F1-8")
all_sequences["sone_2019_f1_8"] = sone_2019_f1_8

# F3-1: Length should be 420, wraps at 72
sone_2019_f3_1 = "GGC"*91 + ("GGA"*1 + "GGC"*2)*16 + "GGC"*1
print_sequence_info(sone_2019_f3_1, 72, "2. Sone et al (2019) F3-1")
all_sequences["sone_2019_f3_1"] = sone_2019_f3_1

# F4-2: Length should be 690, wraps at 75
sone_2019_f4_2 = "GGC"*35 + "C" + "GGC"*17 + "G" + "GGC"*25 + "GC" + "GGC"*26 + "GC" + "GGC"*25 + "GTGC" + "GGC"*30 + "C" + "GGC"*23 + "GC" + "GGC"*12 + "GC" + "GGC"*32
print_sequence_info(sone_2019_f4_2, 75, "2. Sone et al (2019) F4-2")
all_sequences["sone_2019_f4_2"] = sone_2019_f4_2

# F7-1: Length should be 392, wraps at 75
sone_2019_f7_1 = "GGC"*30 + "GC" + "GGC"*20 + "C" + "GGC"*19 + "GC" + "GGC"*60
print_sequence_info(sone_2019_f7_1, 75, "2. Sone et al (2019) F7-1")
all_sequences["sone_2019_f7_1"] = sone_2019_f7_1

# F8-1: Length should be 381, wraps at 75
sone_2019_f8_1 = "GGC"*57 + "GGT" + "GGC"*66 + "GGA" + "GGC"*2
print_sequence_info(sone_2019_f8_1, 75, "2. Sone et al (2019) F8-1")
all_sequences["sone_2019_f8_1"] = sone_2019_f8_1

# FD-2: Length should be 279, wraps at 75
sone_2019_fd_2 = "GGC"*93
print_sequence_info(sone_2019_fd_2, 75, "2. Sone et al (2019) FD-2")
all_sequences["sone_2019_fd_2"] = sone_2019_fd_2

# 3619: Length should be 336, wraps at 75
sone_2019_3619 = "GGC"*112
print_sequence_info(sone_2019_3619, 75, "2. Sone et al (2019) 3619")
all_sequences["sone_2019_3619"] = sone_2019_3619

# 3624: Length should be 363, wraps at 75
sone_2019_3624 = "GGC"*121
print_sequence_info(sone_2019_3624, 75, "2. Sone et al (2019) 3624")
all_sequences["sone_2019_3624"] = sone_2019_3624

# 3692: Length should be 405, wraps at 75
sone_2019_3692 = "GGC"*44 + "GGAGGC"*2 + "GGC"*3 + "GGA" + "GGC"*4 + ("GGA"*2 + "GGC"*3)*15 + "GGA"*2 + "GGC"*2 
print_sequence_info(sone_2019_3692, 75, "2. Sone et al (2019) 3692")
all_sequences["sone_2019_3692"] = sone_2019_3692
##----------------------------------------------------------------------------------------------------------




# 11. Yu, Luan et al, 2021
##----------------------------------------------------------------------------------------------------------
# F1-III-5: wraps at 49? Nope, it's inconsistent: sometimes 49, sometimes 50. Thank you authors who published this as an imagine in Times New Roman.
# yu_luan_2021_f1_III_5 = "GGC"*14 + "GCA" + "GGC"*8 + "GC" + "GGC" + "GGT" + "GGC"*10 + "GC" + "GGC"*4 + "GGT" + "GGC"*8 + "GT" + "GGA"*2 + "GAGGGA" + "GGT" + "GGC" + "GGGC" + "GGC"*4 + "AT" + "GGC"*7 + "GGT"*2 + "GGC" + "GGT"*5 + "GGC"*2 + "A" + "GGC" + "GGT"*5
# Not sure this is feasible - extremely long and irregular structure
yu_luan_2021_f1_iii_5 = """GGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGCAGGCG
GCGGCGGCGGCGGCGGCGGCGCGGCGGTGGCGGCGGCGGCGGCGGCGGC
GGCGGCGGCGCGGCGGCGGCGGCGGTGGCGGCGGCGGCGGCGGCGGCGG
CGTGGAGGAGAGGGAGGTGGCGGGCGGCGGCGGCGGCATGGCGGCGGCG
GCGGCGGCGGCGGTGGTGGCGGTGGTGGTGGTGGTGGCGGCAGGCGGTG
GTGGTGGTGGTGGCGGTGGTGGCGGTGGCGGTGCGGTGGCGGCGGTGGTG
GTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGCGGCGGTGGTGGA
GGCAGGTGGTGGTGCAGGGCAGGAGGAGGTGGTGGTGGAGGAGGTGGTG
GCGGAGGAGGTGCGCGGGAGGGAGGTGGCGGTGGAGGAGGTGGTGCGGC
AGGGTGGCGGAGGAGGTGGCGGCGGAGGAGGCGGTGGCGGGAGGGAGG
TGGCGGCGGAGGAGGTGGCGGCGGGAGGAGGCGGTGGCGGAGGAGGCG
GCGGCGGAGGAGGCGGCGGCGGAGGTGGAGGCGGCGGCGGAGGAGGCG
GCGGCGGAGGAGGCGGCGGCGGAGAGGCGGCGGCGGAGGAGGTGGCGG
CGGAGGAGGCGGCGGTGGAGGAGGCGGCG
""".replace("\n", "")

print_sequence_info(yu_luan_2021_f1_iii_5, 49, "11. Yu, Luan et al (2021) F1-III-5")
all_sequences["yu_luan_2021_f1_iii_5"] = yu_luan_2021_f1_iii_5


# F2-IV-5
yu_luan_2021_f2_iv_5 = """GGCAGCAGCGGCGAGCAGCAGCAGCGGCGCTGGCAGCTGGCTGGCGAAG
CAGCAGCAGCAGCGGCGGCGGCAGCAGCGGCAGCAGCGGCAGCGGCAAC
GAGCAACGGCAACGGCAACGGCAGCCTTGGCAGCAGCAGCAGTAGCAGC
GGCGCGGCAGCGGCGGCAGCAACGGCAACGGCAACGGCAACGGCAGCGG
CGGCAGCGGCGGCAACGGCAACGGCGGCGGCCCGAGCGGCAGAGCGGCG
GAGAACGGAGGAACGGCAGAACGGAGGGCGGCGGAGGAGGCAGCGGAG
GAGAGAGGAACGGCGGAGGAAGGAGGGCACGGAGGAGGAGGCGGCGGA
GGAGGAGGCGGCGGAGGAGGAGGCAAGCGGAGGGGCGGCGGAGGAGGA
GGCAGCGGAGGAGGAGGCAGCGGAGGAGGCAGCTTGGAGGGAGAGGCA
GCGGAGGAGGAGGCCCAGCGGAGGAGGAGACTGACGGAGGAGGGAGCG
GCGGAGGAGGCAGCGGAGGAGGAGGCGGCGGAGGAGGAGGCCAGCGGA
GGAGGGCCCAGCGGAGGGAGGAGAAGGCAGCGGAGAGGGAGGCAGCGG
AGGAGGCGAAGCGGAGGAGGCGGCGGAGGCGAGCGGAGGAGGCAGCGG
AGGGAGGGAGGCGAGCGGCCGGAGGAGGAGGCGGCGGAGGAGGAGGCG
CGGAGGAGGAGGCGGCGGAGGAGGAGGCAGCGGAGGAGGGAGAACGGC
GAAGGAGAGGAGGCGGCGGAGGAGGCAGCGGAGGAGGCGGCTGGAGGA
GGCAGTAGCGGAGGGAGGCGGCGGAGGAGGAGGCGGCGGAGGAGGCAG
CGGAGGGAGGCGGCGGAGGGCGGCGGAGGAGGCGGCGGAGGAGGGGCA
GCGGAGGAGGAGGCGGCGGAGGAGGAGGCGGCGGAGAGCGGCGGAGGA
GGCGCGGAGGAGGCAGCGGAGAGGCGGCGGAGGAGGAGGCGGCGGAGA
GGAGGCGGCGGCGGAGGGAGGCGGCG
""".replace("\n", "")

print_sequence_info(yu_luan_2021_f2_iv_5, 49, "11. Yu, Luan et al (2021) F2-IV-5")
all_sequences["yu_luan_2021_f2_iv_5"] = yu_luan_2021_f2_iv_5


# S1
yu_luan_2021_s1 = """GGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCG
GCGGCGGCGGCGGCGGCGGCGGCGGGCGCGCGCGCGCGGCGGCGGCGCG
CGCGGCATGGCGGCGGCGGCGGCGGTGCGCGGGCGCGGCGTGCGGCGGC
GGCGGCGGCGGCGGCGGCGGTGGCGGCGGCGGCGGCGGCGGCGGCGGGC
GGAGGTGGCGGCGGCG
""".replace("\n", "")

print_sequence_info(yu_luan_2021_s1, 49, "11. Yu, Luan et al (2021) S1")
all_sequences["yu_luan_2021_s1"] = yu_luan_2021_s1
##----------------------------------------------------------------------------------------------------------



# 12. Yu, Deng et al, 2021
##----------------------------------------------------------------------------------------------------------
# F1-III-8 - should be 559 bases according to the document
# Adding up characters from the document: there are 193 not highlighted and 122 GCC units, 193 + 3*122 = 559 so that checks out.
yu_deng_2021_f1_iii_8 = ("GGC"*10 + "AGC" + "GGC"*7 + "GCGC" + "GGC"*20 + "AGC" + "GGC"*3 + "AC" + "GGC"*2 + "AGCAGC" + "GGC"*5 + "AGC" + "GGC"*4 + "AC" + "GGC"*1 +
                         "AGC" + "GGC"*1 + "AGC" + "GGC"*6 + "AAC" + "GGC"*2 + "AAC" + "GGC"*1 + "AGCAGCAGCAAC" + "GGC"*1 + "AGCAGC" + "GGC"*1 + "AAAG" +
                         "CAGCAGCGACGAGCACAGCAAC" + "GGC"*1 + "AGCAGCAGCAGCGCT" + "GGC"*1 + "AGCAGC" + "GGC"*3 + "GC" + "GG" +
                         "C"*1 + "AGCAACAGCCAACAGCAACAGCAGCAGCAGCAGC" + "GGC"*1 + "AGC" + "GGC"*3 + "C" + "GGC"*6 + "AGCAGCGC" +
                         "GGC"*4 + "AGC" + "GGC"*9 + "G" + "GGC"*6 + "T" + "GGC"*1 + "GCAGCAGC" + "GGC"*1 + "AGC" + "GGC"*5 + "AGCAGCAGCAGC" + "GGC"*2 + "A" +
                         "GC" + "GGC"*4 + "AC" + "GGC"*1 + "AGC" + "GGC"*1 + "AGC" + "GGC"*5 + "AGCAGC" + "GGC"*3)
print_sequence_info(yu_deng_2021_f1_iii_8, 75, "12. Yu, Deng et al (2021) F1-III-8")
all_sequences["yu_deng_2021_f1_iii_8"] = yu_deng_2021_f1_iii_8


# F1-III-10 - should be 395 bases according to the document
# Adding up characters from the document: there are 61 not highlighted and 113 GGC units, 61 + 3*113 = 400. Did they miscount?
yu_deng_2021_f1_iii_10 = ("GGC"*4 + "GC" + "GGC"*1 + "C" + "GGC"*4 + "GC" + "GGC"*9 + "AGCAAC" + "GGC"*1 + "AAC" + "GGC"*2 + "CCTTC" + "GGC"*14 + "AC" + "GGC"*11 + "AAC" +
                          "GGC"*4 + "AGCAAC" + "GGC"*3 + "AGCAGC" + "GGC"*11 + "AAC" + "GGC"*10 + "AGC" + "GGC"*14 + "GC" + "GGC"*2 + "AACGCAAC" + "GGC"*4 +
                          "AGC" + "GGC"*12 + "AGCAC" + "GGC"*1 + "A" + "GGC"*6)
print_sequence_info(yu_deng_2021_f1_iii_10, 75, "12. Yu, Deng et al (2021) F1-III-10")
all_sequences["yu_deng_2021_f1_iii_10"] = yu_deng_2021_f1_iii_10


# S1 - should be 574 bases according to the document
# Adding up characters from the document: there are 184 not highlighted and 130 GGC units, 184 + 3*130 = 574 so that checks out.
yu_deng_2021_s1 = ("G" + "GGC"*14 + "GC" + "GGC"*3 + "AGC" + "GGC"*4 + "AGCAAC" + "GGC"*1 + "AAC" + "GGC"*1 + "AGC" + "GGC"*1 + "AGCAGCCCAAC" + "GGC"*1 +
                   "AGC" + "GGC"*6 + "AGC" + "GGC"*4 + "AGC" + "GGC"*1 + "AGC" + "GGC"*9 + "GAC" + "GGC"*1 + "AAC" + "GGC"*1 + "AAC" + "GGC"*4 + "AAC" + "GGC"*1 + "A" +
                   "AC" + "GGC"*1 + "AACTGACAGC" + "GGC"*2 + "CTTGCCGC" + "GGC"*1 + "AGC" + "GGC"*2 + "GAGC" + "GGC"*1 + "AGCAGCAAC" + "GGC"*6 +
                   "GC" + "GGC"*1 + "AGC" + "GGC"*1 + "ACCGC" + "GGC"*1 + "AGCAGC" + "GGC"*2 + "TGACAGCAAC" + "GGC"*1 + "AACCTTAT" + "GG" + 
                   "C"*1 + "GCAGCAGC" + "GGC"*1 + "AGCAGCGC" + "GGC"*1 + "AGC" + "GGC"*17 + "AAC" + "GGC"*2 + "AGC" + "GGC"*4 + "AGCTGAC" + "GGC"*1 +
                   "AAC" + "GGC"*3 + "AGC" + "GGC"*2 + "GAGC" + "GGC"*21 + "GCT" + "GGC"*2 + "AGCAGCAGC" + "GGC"*1 + "AGC" + "GGC"*1 + "AAC" + "GGC"*2)
print_sequence_info(yu_deng_2021_s1, 75, "12. Yu, Deng et al (2021) S1")
all_sequences["yu_deng_2021_s1"] = yu_deng_2021_s1


# NIID-1 - should be 424 bases according to the document
yu_deng_2021_niid_1 = ("GGC"*7 + "AGC" + "GGC"*14 + "AAC" + "GGC"*2 + "AGC" + "GGC"*3 + "AGC" + "GGC"*1 + "AGC" + "GGC"*17 + "AGC" + "GGC"*2 + "AGCAGC" + "GGC"*15 +
                       "AGC" + "GGC"*2 + "AGC" + "GGC"*17 + "CGGAAGATGCCCGCCCTGCGCCGCTCACTGTGACGCCGCCGCCA" + 
                       "CGCC" + "GGC"*2 + "AGC" + "GGC"*2 + "T" + "GGC"*1 + "TGCAAC" + "GGC"*3 + "AGC" + "GGC"*4 + "AGC" + "GGC"*3 + "AGC" + "GGC"*1 + "AGCAGC" +
                       "GGC"*2 + "AGC" + "GGC"*8)
print_sequence_info(yu_deng_2021_niid_1, 75, "12. Yu, Deng et al (2021) NIID-1")
all_sequences["yu_deng_2021_niid_1"] = yu_deng_2021_niid_1


# NIID-2 - should be 460 bases according to the document
yu_deng_2021_niid_2 = ("GGC"*14 + "GGTGGT" + "GGC"*12 + "GGT" + "GGC"*13 + "GGT" + "GGC"*3 + "GGT" + "GGC"*3 + "GGT" + "GGC"*2 + "GGT" + "GGC"*16 + "GGT" + "GGC"*8 +
                       "G" + "GGC"*3 + "GGTGC" + "GGC"*2 + "GGT" + "GGC"*1 + "GGTGGT" + "GGC"*3 + "GGGT" + "GGC"*3 + "GGTGA" + "GGC"*1 + "AGC" + "GGC"*9 +
                       "GGAGGA" + "GGC"*1 + "AGGGTGGT" + "GGC"*2 + "GGGAGGAGGTGGT" + "GGC"*1 + "GGAGA" + "GGC"*11 + "GGAGGTGA" + "GGC"*11 +
                       "GGAGGA" + "GGC"*2)
print_sequence_info(yu_deng_2021_niid_2, 75, "12. Yu, Deng et al (2021) NIID-2")
all_sequences["yu_deng_2021_niid_2"] = yu_deng_2021_niid_2


# NIID-3 - should be 612 bases according to the document
yu_deng_2021_niid_3 = ("GATCTGCCCA" + "GGC"*17 + "GC" + "GGC"*5 + "AAC" + "GGC"*2 + "TGACAACAAC" + "GGC"*5 + "GAC" + "GGC"*10 + "AGC" + "GGC"*7 + "AG" +
                       "C" + "GGC"*3 + "GAGCG" + "GGC"*1 + "AACGAC" + "GGC"*9 + "AAC" + "GGC"*1 + "T" + "GGC"*1 + "AC" + "GGC"*1 + "T" + "GGC"*1 + "T" + "GGC"*1 + "AAC" +
                       "GGC"*3 + "AAC" + "GGC"*2 + "AAC" + "GGC"*6 + "AAC" + "GGC"*1 + "TGCAACAACGAC" + "GGC"*7 + "GAC" + "GGC"*7 + "GAC" + "GGC"*4 + "GAC" + 
                       "GGC"*2 + "GAC" + "GGC"*24 + "GAGC" + "GGC"*2 + "AAC" + "GGC"*7 + "GACGAC" + "GGC"*11 + "GC" + "GGC"*6 + "GAC" + "GGC"*2 + "AAC" + "GGC"*1 + "A" +
                       "ACGAC" + "GGC"*1 + "AAC" + "GGC"*1 + "AAC" + "GGC"*4 + "AGC" + "GGC"*1 + "AAC" + "GGC"*3 + "GGAAAGA")
print_sequence_info(yu_deng_2021_niid_3, 75, "Yu, Deng et al (2021) NIID-3")
all_sequences["yu_deng_2021_niid_3"] = yu_deng_2021_niid_3


# NIID-4 - should be 390 bases according to the document
yu_deng_2021_niid_4 = ("GGC"*11 + "AGC" + "GGC"*39 + "AGCGC" + "GGC"*4 + "AGC" + "GGC"*2 + "AGC" + "GGC"*2 + "G" + "GGC"*3 + "AGCAGCAAC" + "GGC"*4 + "GGAG" +
                       "GA" + "GGC"*3 + "AGC" + "GGC"*1 + "AAC" + "GGC"*3 + "GGA" + "GGC"*11 + "GGAAGAGAC" + "GGC"*8 + "AAC" + "GGC"*3 + "GGA" + "GGC"*1 + "T" +
                       "GGC"*1 + "AGC" + "GGC"*1 + "AGCGGAGGA" + "GGC"*3 + "G" + "GGC"*2 + "GGAAGGA" + "GGC"*3)
print_sequence_info(yu_deng_2021_niid_4, 75, "Yu, Deng et al (2021) NIID-4")
all_sequences["yu_deng_2021_niid_4"] = yu_deng_2021_niid_4


# NIID-5 - should be 401 bases according to the document
yu_deng_2021_niid_5 = ("GGC"*6 + "AGC" + "GGC"*11 + "AGCC" + "GGC"*2 + "AGCAGC" + "GGC"*3 + "AGCAGCAGCAAC" + "GGC"*1 + "ACAACAGCAGCAG" +
                       "C" + "GGC"*2 + "AGCGCTGC" + "GGC"*2 + "AGCAGC" + "GGC"*1 + "AGC" + "GGC"*1 + "AGC" + "GGC"*1 + "AAC" + "GGC"*12 + "AGCAGC" + "GGC"*6 +
                       "AGCAAC" + "GGC"*1 + "AAC" + "GGC"*6 + "AGCAGC" + "GGC"*8 + "AGCAAC" + "GGC"*16 + "AGC" + "GGC"*3 + "AGC" + "GGC"*1 + "AGC" +
                       "GGC"*7 + "AGC" + "GGC"*1 + "AAC" + "GGC"*1 + "AGC" + "GGC"*1 + "AGCAGC" + "GGC"*1 + "AGC" + "GGC"*1)
print_sequence_info(yu_deng_2021_niid_5, 75, "Yu, Deng et al (2021) NIID-5")
all_sequences["yu_deng_2021_niid_5"] = yu_deng_2021_niid_5
##----------------------------------------------------------------------------------------------------------



## 13. Fukuda et al, 2021
##----------------------------------------------------------------------------------------------------------
# F1-patient: length should be 300, wraps at 75
fukuda_2021_f1_patient = "GGC"*90 + "GACCGAGAAGATGCCCGCCCTGC" + "GGC"*2 + "G"
print_sequence_info(fukuda_2021_f1_patient, 75, "13. Fukuda et al (2021) F1-patient")
all_sequences["fukuda_2021_f1_patient"] = fukuda_2021_f1_patient

# F1-father: length should be 718, wraps at 75
fukuda_2021_f1_father = "GGC"*239 + "G"
print_sequence_info(fukuda_2021_f1_father, 75, "13. Fukuda et al (2021) F1-father")
all_sequences["fukuda_2021_f1_father"] = fukuda_2021_f1_father

# F2-patient: length should be 358, wraps at 75
fukuda_2021_f2_patient = "GGC"*109 + "GGA" + "GGC"*6 + "GGA" + "GGC"*2 + "G"
print_sequence_info(fukuda_2021_f2_patient, 75, "13. Fukuda et al (2021) F2-patient")
all_sequences["fukuda_2021_f2_patient"] = fukuda_2021_f2_patient

# F2-father: length should be 712, wraps at 75
fukuda_2021_f2_father = "GGC"*227 + "GGA" + "GGC"*6 + "GGA" + "GGC"*2 + "G"
print_sequence_info(fukuda_2021_f2_father, 75, "13. Fukuda et al (2021) F2-father")
all_sequences["fukuda_2021_f2_father"] = fukuda_2021_f2_father

# F3-patient: length should be 496, wraps at 75
fukuda_2021_f3_patient = "GGC"*98 + ("GGA" + "GGC"*4 + "GGA")*10 + "GGC"*3 + "GGA"*2 + "GGC"*2 + "G"
print_sequence_info(fukuda_2021_f3_patient, 75, "13. Fukuda et al (2021) F3-patient")
all_sequences["fukuda_2021_f3_patient"] = fukuda_2021_f3_patient

# F3-father: length should be 1588, wraps at 75
fukuda_2021_f3_father = "GGC"*70 + "GG" + "GGC"*47 + "G" + "GGC"*16 + "G" + "GGC"*30 + "GC" + "GGC"*3 + "G" + "GGC"*170 + "G" + "GGC"*34 + "G" + "GGC"*77 + "GGA" + "GGC"*5 + ("GGA" + "GGC"*4 + "GGA")*11 + "GGC"*3 + "GGA"*2 + "GGC"*2 + "G"
print_sequence_info(fukuda_2021_f3_father, 75, "13. Fukuda et al (2021) F3-father")
all_sequences["fukuda_2021_f3_father"] = fukuda_2021_f3_father

# F4-patient: length should be 448, wraps at 75
fukuda_2021_f4_patient = "GGC"*149 + "G"
print_sequence_info(fukuda_2021_f4_patient, 75, "13. Fukuda et al (2021) F4-patient")
all_sequences["fukuda_2021_f4_patient"] = fukuda_2021_f4_patient

# F4-father: length should be 1044, wraps at 75
fukuda_2021_f4_father = "GGC"*16 + "GG" + "GGC"*331 + "G"
print_sequence_info(fukuda_2021_f4_father, 75, "13. Fukuda et al (2021) F4-father")
all_sequences["fukuda_2021_f4_father"] = fukuda_2021_f4_father
##----------------------------------------------------------------------------------------------------------




## 16. Liu et al, 2022
##----------------------------------------------------------------------------------------------------------
# Patient 1: Doesn't give numbers, but block appears to be 13*GGC wide (i.e. wraps at 39) and has 10 rows with 25 extra. 39*10 + 25 = 415
liu_2022_patient_1 = "GGC"*27 + "GC" + "GGC"*35 + "GC" + "GGC"*10 + "GCG" + "GGC"*61 + "GGA" + "GGC"*2
print_sequence_info(liu_2022_patient_1, 39, "16. P Liu et al (2022) Patient 1")
all_sequences["liu_2022_patient_1"] = liu_2022_patient_1

# Patient 3: Wraps at 39, has 7 complete rows with 17 extra. 39*7 + 17 = 290
liu_2022_patient_3 = "GGC"*91 + "GC" + "GGC"*5
print_sequence_info(liu_2022_patient_3, 39, "16. P Liu et al (2022) Patient 3")
all_sequences["liu_2022_patient_3"] = liu_2022_patient_3
##----------------------------------------------------------------------------------------------------------




## 19. Kameyama et al, 2022
##----------------------------------------------------------------------------------------------------------
# F1-II-3 (unphased) - 285 bp, wraps at 75
kameyama_2022_f1_ii_3 = "GGC"*95
print_sequence_info(kameyama_2022_f1_ii_3, 75, "19. Kameyama et al (2022) F1-II-3")
all_sequences["kameyama_2022_f1_ii_3"] = kameyama_2022_f1_ii_3

# F2-II-2 (allele 1) - 207 bp, wraps at 75
kameyama_2022_f2_ii_2_allele_1 = "GGC"*65 + "GGA"*2 + "GGC"*2
print_sequence_info(kameyama_2022_f2_ii_2_allele_1, 75, "19. Kameyama et al (2022) F2-II-2, allele 1")
all_sequences["kameyama_2022_f2_ii_2_allele_1"] = kameyama_2022_f2_ii_2_allele_1

# F2-II-2 (allele 2) - 309 bp, wraps at 75
kameyama_2022_f2_ii_2_allele_2 = "GGC"*99 + "GGA"*2 + "GGC"*2
print_sequence_info(kameyama_2022_f2_ii_2_allele_2, 75, "19. Kameyama et al (2022) F2-II-2, allele 2")
all_sequences["kameyama_2022_f2_ii_2_allele_2"] = kameyama_2022_f2_ii_2_allele_2
##----------------------------------------------------------------------------------------------------------




## 20. Fitrah et al, 2023
##----------------------------------------------------------------------------------------------------------
## MAIN PAPER
# Patient 7 Allele 1 - length 328, wrapping inconsistent
fitrah_2023_patient_7_allele_1 = ("GGC"*5 + "GC" + "GGC"*6 + "G" + "GGC"*12 + "GTG" + 
                                  "C" + "GGC"*12 + "G" + "GGC"*12 + 
                                  "GGC"*25 +
                                  "GGC"*17 + "GC" + "GGC"*7 + 
                                  "GGC"*6 + "GGA"*2 + "GGC"*2)
print_sequence_info(fitrah_2023_patient_7_allele_1, 75, "20. Fitrah et al (2023) Patient 7, allele 1")
all_sequences["fitrah_2023_patient_7_allele_1"] = fitrah_2023_patient_7_allele_1

# Patient 7 Allele 2 - length 351, wraps at 75
fitrah_2023_patient_7_allele_2 = "GGC"*69 + "AGCGGC"*2 + "GGC"*23 + "AGC" + "GGC"*16 + "GGA"*2 + "GGC"*2
print_sequence_info(fitrah_2023_patient_7_allele_2, 75, "20. Fitrah et al (2023) Patient 7, allele 2")
all_sequences["fitrah_2023_patient_7_allele_2"] = fitrah_2023_patient_7_allele_2

# Patient 14 Allele 1 - says length is 405, wrapping inconsistent. But, counting the last row starting from 375, ends at 404.
# I think they forgot there's a missing character at the end of the second to last row. If that row was complete, the last position would be 405.
fitrah_2023_patient_14_allele_1 = ("GGC"*25 + 
                                   "GGC"*3 + "GG" + "GGC"*11 + "GG" + "GGC"*3 + "GG" + "GGC"*6 + 
                                   "GG" + "GGC"*2 + "GG" + "GGC"*5 + "GGGGGGCGG" + "GGC"*4 + "G"*18 + "GGC"*2 + "G"*5 +
                                   "GGC" + "G"*6 + "GGC"*3 + "GGGGCAG" + "GGC"*6 + "GG" + "GGC"*3 + "GG" + "GGC" + "G"*4 + "GGC"*4 +
                                   "GGC"*3 + "GC" + "GGC"*4 + "GGGGC"*2 + "GGC"*2 + "GC" + "GGC"*11 +
                                   "GGC"*10)
print_sequence_info(fitrah_2023_patient_14_allele_1, 75, "20. Fitrah et al (2023) Patient 14, allele 1")
all_sequences["fitrah_2023_patient_14_allele_1"] = fitrah_2023_patient_14_allele_1

# Patient 14 Allele 2 - length 1083, wrapping inconsistent. Their base count is wrong e.g. goes from 672 to 745 when 75 bases are added, which should take it to 747
fitrah_2023_patient_14_allele_2 = ("TTA"*2 + "CTA" + "TTA"*6 + "CTA" + "TTA"*2 + "CTA" + "TTA" + "ATA" + "TTA"*2 + "TTTAATATA" + "TTA"*2 + "TCA" + "TTA"*2 +
                                   "TTTATTATTATCATTATTATCATTATTTATCATTATTATTATCATTATCATTTACTATTACCATTATTATCATTA" +
                                   "TACATTAATATTATCATTATTTATTTATTATTACTATTATCATTACTATTATTATTATTTATTATTATTATTTAT" +
                                   "TTACTATTACTTACTATTTATTTATTATTTATTTACGATGATACGATGATGATGATGATGATGATGATGATGAT" + #"." +
                                   "GAT"*7 + "ATTAT" + "ATG"*3 + "AG" + "ATG"*12 + "AT" +
                                   "GAG" + "GAT"*24 +
                                   "GAT"*25 +
                                   "GAT"*9 + "GAC" + "GAT"*4 + "GC" + "GAT"*5 + "GAC" + "GAT"*4 + #"." +
                                   "GAT"*10 + "GACATCACC" + "GAT" + "TCACC" + "GAT" + "GAC"*4 + "GAT" + "GAC" + "GAT"*2 + #"." +
                                   "GAC"*3 + "GAT" + "GCTGACCACAGATTAATGACGACGACCACTAATTAATTATCAGGATTACCACCACTACCTGAC" +
                                   "CAGATTATTTGACGACCATTACCACTATTTTTAATTATGACACACCACTACCACCACTACTATTAATTAATACCA" +
                                   "CTTTTCACTACTTCACTACTACTACCACCACTAATTAATACCACCACCACCACTACATTACATGATTATCATACT" +
                                   "ATTTGTTTACATGATGACATTATTACTATTACTATTTACTTACTATTATTATGCGGCGGCGGCGGCGGCGGCGGC" +
                                   "GGC"*11 + "GC" + "GGC" + "GGGTGC" + "GGC"*2 + "G" + "GGC" + "TGGCA" + "GGC"*5 + #"." +
                                   "GGC"*16)
print_sequence_info(fitrah_2023_patient_14_allele_2, 75, "20. Fitrah et al (2023) Patient 14, allele 2")
all_sequences["fitrah_2023_patient_14_allele_2"] = fitrah_2023_patient_14_allele_2


## SUPPLEMENT
# Patient 1 - length claimed 594, wraps at 75
fitrah_2023_patient_1 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G" +
                         "G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G" +
                         "G C G G C G G C G G C G G C G G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C A G C G G C G G C G G C G G C G G C G G C A G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G A G G C G G A G G C G G C G G A G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_1, 75, "20. Fitrah et al (2023) Patient 1")
all_sequences["fitrah_2023_patient_1"] = fitrah_2023_patient_1

# Patient 2 - length claimed 313, wrapping inconsistent
fitrah_2023_patient_2 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G G C G G C" + #".." +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_2, 75, "20. Fitrah et al (2023) Patient 2")
all_sequences["fitrah_2023_patient_2"] = fitrah_2023_patient_2

# Patient 3 - length claimed 475, wrapping inconsistent
fitrah_2023_patient_3 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" + #".." +
                         "G G C G G C G G C A G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C A G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C"
                         ).replace(" ", "")
print_sequence_info(fitrah_2023_patient_3, 75, "20. Fitrah et al (2023) Patient 3")
all_sequences["fitrah_2023_patient_3"] = fitrah_2023_patient_3

# Patient 4 - length claimed 306, wraps at 75
fitrah_2023_patient_4 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G A" +
                         "G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_4, 75, "20. Fitrah et al (2023) Patient 4")
all_sequences["fitrah_2023_patient_4"] = fitrah_2023_patient_4

# Patient 5 - length claimed 303, wraps at 75
fitrah_2023_patient_5 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G A G G A G G C" +
                         "G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_5, 75, "20. Fitrah et al (2023) Patient 5")
all_sequences["fitrah_2023_patient_5"] = fitrah_2023_patient_5

# Patient 6 - length claimed 320, wrapping inconsistent
fitrah_2023_patient_6 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G G G C" + #"." +
                         "G G C G G C G G C G G C G G C G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_6, 75, "20. Fitrah et al (2023) Patient 6")
all_sequences["fitrah_2023_patient_6"] = fitrah_2023_patient_6

# Patient 8 - length claimed 305, wrapping inconsistent
fitrah_2023_patient_8 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G A C G G C G G C" + 
                         "G G C G G C G G C A G C G G C G G C G G C G G C G G C A G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C A G C" +
                         "G G C A G C G G C G G C G G C G G C G G C G G C G G C G G C G G C A G G C G C G G C G G C G G C G G C G C G G C G A C G G C G G C T G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G A G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G G C G G C G G A G G A" + #"." +
                         "G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_8, 75, "20. Fitrah et al (2023) Patient 8")
all_sequences["fitrah_2023_patient_8"] = fitrah_2023_patient_8

# Patient 9 - length claimed 327, wraps at 75
fitrah_2023_patient_9 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G A A C A C C C A C C C C C C C G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                         "G G C G G C G G C G G C G G C G G A G G A G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_9, 75, "20. Fitrah et al (2023) Patient 9")
all_sequences["fitrah_2023_patient_9"] = fitrah_2023_patient_9

# Patient 10 - length claimed 339, wraps at 75
fitrah_2023_patient_10 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G A G G A G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_10, 75, "20. Fitrah et al (2023) Patient 10")
all_sequences["fitrah_2023_patient_10"] = fitrah_2023_patient_10

# Patient 11 - length claimed 524, wraps at 75
fitrah_2023_patient_11 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G G G A A" +
                          "A C G G C A G C A G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C A G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C A G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_11, 75, "20. Fitrah et al (2023) Patient 11")
all_sequences["fitrah_2023_patient_11"] = fitrah_2023_patient_11

# Patient 12 - length claimed 423, wraps at 75
fitrah_2023_patient_12 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G A G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_12, 75, "20. Fitrah et al (2023) Patient 12")
all_sequences["fitrah_2023_patient_12"] = fitrah_2023_patient_12

# Patient 13 - length claimed 283, wrapping inconsistent
fitrah_2023_patient_13 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G G C G G C G G C G G C G G C" + #".." +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G A G G A G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_13, 75, "20. Fitrah et al (2023) Patient 13")
all_sequences["fitrah_2023_patient_13"] = fitrah_2023_patient_13

# Patient 15 - length claimed 290, wrapping inconsistent
fitrah_2023_patient_15 = ("G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" + #"." +
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C" + 
                          "G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G C G G A G G A G G C G G C").replace(" ", "")
print_sequence_info(fitrah_2023_patient_15, 75, "20. Fitrah et al (2023) Patient 15")
all_sequences["fitrah_2023_patient_15"] = fitrah_2023_patient_15
##----------------------------------------------------------------------------------------------------------




## 27. Watanabe et al, 2024
##----------------------------------------------------------------------------------------------------------
watanabe_2024_ii_7 = ("GGC"*21 + "G" +
                      "GC" + "GGC"*23 + "GC" + "GGC" + "GGGC" + "GGGGC" + "GGC"*3 + "G" +
                      "GC" + "GGC"*4 + "GGGGGGC" + "GGC"*3 + "GCGGGCGGGGC" + "GGC"*13 + "GGGGC" + "GGC"*3 + "G" +
                      "GC" + "GGC"*17 + "GC" + "GGC"*6 + "GC" + "GGC"*2 + "GGGGC" + "GGC" + "GGGC" + "GG" +
                      "C" + "GGC"*5 + "GGGGGC"*2 + "GGC"*20 + "GC" + "GGC" + "GG" + 
                      "C" + "GGC"*31 + "G" +
                      "GC" + "GGC"*9 + "GGGC" + "GGC"*8 + "GGGC" + "GGC"*2 + "GCGC" + "GGC" + "GC" + "GGC"*4 + "GGGC" + "GGG" +
                      "C" + "GGC"*7 + "GGGCGGGCGC" + "GGC"*4 + "GGGC"*2 + "GGC" + "GGGC" + "GGC"*12 + 
                      "GGCGGGC"*2 + "GGC"*10 + "GGGCGGC" + "GGGC"*2 + "GGC"*4 + "GC" + "GGC"*7 + "G" +
                      "GC" + "GGC"*3 + "GGGC" + "GGC"*3 + "GC"*3 + "GGGGC" + "GGC"*4 + "GC" + "GGC"*15 + "G" +
                      "GC" + "GGC"*31 +
                      "GGC"*30 + "GGGCG" +
                      "GC" + "GGC"*3 + "GGAGA" + ("GGC"*4 + "GGA"*2)*4 + "GGC"*2 + "G" +
                      "GCGGC" + "AGGA" + "GGC"*4 + "GGAGA" + "GGA"*3 + ("GGC"*4 + "GGA"*2)*3 + "GGC"*2 +
                      "GGC"*2 + "GGA"*2 + ("GGC"*4 + "GGA"*2)*4 + "GGCGGC").replace("\n", "")
print_sequence_info(watanabe_2024_ii_7, 95, "27. Watanabe et al (2024) II-7")
all_sequences["watanabe_2024_ii_7"] = watanabe_2024_ii_7

watanabe_2024_iii_12 = ("GGC"*21 + "GG" +
                        "C" + "GGC"*31 + "GG" +
                        "C" + "GGC"*31 + "GG" + 
                        "C" + "GGC"*18 + ("GGA"*2 + "GGC"*4)*2 + "GGAGG" +
                        "A" + ("GGC"*4 + "GGA"*2)*5 + "GGCGG" + 
                        "C" + "GGC"*2 + "GGA"*2 + "GGC"*2).replace("\n", "")
print_sequence_info(watanabe_2024_iii_12, 96, "27. Watanabe et al (2024) III-12")
all_sequences["watanabe_2024_iii_12"] = watanabe_2024_iii_12

watanabe_2024_iii_13 = ("GGC"*21 + "GG" +
                        "C" + "GGC"*31 + "GG" +
                        "C" + "GGC"*31 + "GG" +
                        "C" + "GGC"*11 + "GGAGGCA" + "GGC"*3 + ("GGA"*2 + "GGC"*4)*2 + "GGAGGAGGCG" +
                        "GC" + "GGC"*2 + ("GGA"*2 + "GGC"*4)*4 + "GGAGGAGGCGGCGGCG" +
                        "GC" + ("GGA"*2 + "GGC"*4)*5 + "GGAG" +
                        "GA" + ("GGC"*4 + "GGA"*2)*4 + "GGC"*2).replace("\n", "")
print_sequence_info(watanabe_2024_iii_13, 96, "27. Watanabe et al (2024) III-13")
all_sequences["watanabe_2024_iii_13"] = watanabe_2024_iii_13
##----------------------------------------------------------------------------------------------------------




## Write output to csv
output_file = open("intermediate_files/python_consensuses_literature.csv", "w")
output_file.write("consensus_id,sequence\n")
for sequence_id in all_sequences.keys():
    sequence_content = all_sequences[sequence_id]
    output_file.write(sequence_id + "," + sequence_content + "\n")
output_file.close()
