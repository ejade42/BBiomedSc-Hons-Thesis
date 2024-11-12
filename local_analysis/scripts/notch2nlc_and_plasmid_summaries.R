## Create N2C transcript and plasmid files used for Figures 1-2 and 2-1.

library(tidyverse)
library(RGenetics)
source("scripts/common_functions.R")

## HELPER FUNCTIONS
## -------------------------------------------------------------------------------------------------------
## Takes a DNA sequence and translates the whole thing to amino acids,
## starting with the first three bases as the first codon
translate_sequence <- function(sequence) {
    num_codons <- nchar(sequence) %/% 3
    aa_seq     <- NULL
    for (i in 1:num_codons) {
        codon <- substr(sequence, i*3-2, i*3)
        aa_seq <- c(aa_seq, codonToAAone(codon))
    }
    return(paste(aa_seq, collapse = ""))
}

## Takes an uppercase amino acid sequence (with start codons represented by a specifiable character),
## and changes any bases to lowercase if they are outside an M-----* ORF 
make_lowercase_if_outside_orf <- function(aa_sequence, stop_codon) {
    active_orf <- FALSE
    for (i in 1:nchar(aa_sequence)) {
        letter = substr(aa_sequence, i, i) 
        if (letter == "M") {
            active_orf <- TRUE
        } else if (letter == stop_codon) {
            active_orf <- FALSE
        }
        
        if (active_orf == FALSE) {
            substr(aa_sequence, i, i) <- tolower(letter)
        }
    }
    return(aa_sequence)
}


## Takes a sequence, and inserts some number of spaces at the start, then lists each character of the sequence
## followed by a specified number of spaces (optional whether spaces are listed after the final character)
insert_spaces_between_sequence <- function(sequence, spaces_between, initial_spaces = 0, end_with_spaces = FALSE, blank = " ") {
    output_sequence <- rep(blank, initial_spaces)
    for (i in 1:(nchar(sequence)-1)) {
        output_sequence <- c(output_sequence, substr(sequence, i, i), rep(blank, spaces_between))
    }
    i <- nchar(sequence)
    output_sequence <- c(output_sequence, substr(sequence, i, i))
    if (end_with_spaces == TRUE) {
        output_sequence <- c(output_sequence, rep(blank, spaces_between))
    }
    return(paste(output_sequence, collapse = ""))
}
## -------------------------------------------------------------------------------------------------------



## PROCESS TRANSCRIPT VARIANT 1
## -------------------------------------------------------------------------------------------------------
## Read in transcript from FASTA, combine into single string, and validate length
## from: https://www.ncbi.nlm.nih.gov/nuccore/NM_001364012.2?report=fasta
notch2nlc_transcript_1 <- readLines("input/notch2nlc_transcript_variant_1.fasta")
notch2nlc_transcript_1 <- notch2nlc_transcript_1[2:length(notch2nlc_transcript_1)]
notch2nlc_transcript_1 <- paste(notch2nlc_transcript_1, collapse = "")

if (nchar(notch2nlc_transcript_1) != 8737) {
    print("Warning: length is incorrect", quote = F)
} else {
    print("Length is correct", quote = F)
}

## Translate
stop_codon <- "*"
reading_frame_1 <- translate_sequence(notch2nlc_transcript_1) %>% 
    gsub("Stop", stop_codon, .) %>%
    make_lowercase_if_outside_orf(., stop_codon)
reading_frame_2 <- translate_sequence(substr(notch2nlc_transcript_1, 2, nchar(notch2nlc_transcript_1))) %>% 
    gsub("Stop", stop_codon, .) %>%
    make_lowercase_if_outside_orf(., stop_codon)
reading_frame_3 <- translate_sequence(substr(notch2nlc_transcript_1, 3, nchar(notch2nlc_transcript_1))) %>% 
    gsub("Stop", stop_codon, .) %>%
    make_lowercase_if_outside_orf(., stop_codon)

## Process exons
exon_lengths <- c(302, 82, 260, 336, 7757)
exon_string <- NULL
blank = "-"
for (i in 1:length(exon_lengths)) {
    this_exon <- paste0(i, paste(rep(blank, exon_lengths[i]-2), collapse = ""), "|")
    exon_string <- paste0(exon_string, this_exon)
}
exon_string

## Output
output_exons <- c(">NOTCH2NLC transcript variant 1 (NCBI GenBank: NM_001364012.2) sequence with forward reading frames",
                  exon_string,
                  notch2nlc_transcript_1,
                  insert_spaces_between_sequence(reading_frame_3, spaces_between = 2, initial_spaces = 3),
                  insert_spaces_between_sequence(reading_frame_1, spaces_between = 2, initial_spaces = 1),
                  insert_spaces_between_sequence(reading_frame_2, spaces_between = 2, initial_spaces = 2))

writeLines(output_exons, "intermediate_files/notch2nlc_transcript_variant_1_reading_frames_exons.txt")
## -------------------------------------------------------------------------------------------------------





## PROCESS PLASMID
## -------------------------------------------------------------------------------------------------------
plasmid_GGC98_fastq <- readLines("input/n2c_uorf_3tag_plasmid_sequence.fastq")
plasmid_GGC98 <- reverse_sequence(plasmid_GGC98_fastq[2])


## Translate
stop_codon <- "*"
reading_frame_1 <- translate_sequence(plasmid_GGC98) %>% 
    gsub("Stop", stop_codon, .)
reading_frame_2 <- translate_sequence(substr(plasmid_GGC98, 2, nchar(plasmid_GGC98))) %>% 
    gsub("Stop", stop_codon, .)
reading_frame_3 <- translate_sequence(substr(plasmid_GGC98, 3, nchar(plasmid_GGC98))) %>% 
    gsub("Stop", stop_codon, .)


## Annotate promoter
cmv_promoter <- "CGTTACATAACTTACGGTAAATGGCCCGCCTGGCTGACCGCCCAACGACCCCCGCCCATTGACGTCAATAATGACGTATGTTCCCATAGTAACGCCAATAGGGACTTTCCATTGACGTCAATGGGTGGAGTATTTACGGTAAACTGCCCACTTGGCAGTACATCAAGTGTATCATATGCCAAGTACGCCCCCTATTGACGTCAATGACGGTAAATGGCCCGCCTGGCATTATGCCCAGTACATGACCTTATGGGACTTTCCTACTTGGCAGTACATCTACGTATTAGTCATCGCTATTACCATGGTGATGCGGTTTTGGCAGTACATCAATGGGCGTGGATAGCGGTTTGACTCACGGGGATTTCCAAGTCTCCACCCCATTGACGTCAATGGGAGTTTGTTTTGGCACCAAAATCAACGGGACTTTCCAAAATGTCGTAACAACTCCGCCCCATTGACGCAAATGGGCGGTAGGCGTGTACGGTGGGAGGTCTATATAAGCAGAGCT"
promoter_location <- str_locate(plasmid_GGC98, cmv_promoter)

exon_1_start <- "CCAAACTTCGGGCGGCGGCTGAG"
exon_1_start_location <- str_locate(plasmid_GGC98, exon_1_start)[1]

uORF_start <- "ATGTGGATCTGCCCA"
uORF_start_location <- str_locate(plasmid_GGC98, uORF_start)[1]
uORF_end <- "GATGCCCGCCCTGCGCCGCTCTGCTGTGGGCGCTGCTGGCGCTCTGGCTGTGCTGCGCGACCCCCGCGCA"
uORF_end_location <- str_locate(plasmid_GGC98, uORF_end)[2]

annotation <- rep(" ", nchar(plasmid_GGC98))

annotation[promoter_location[1]] <- "P"
annotation[(promoter_location[1]+1):(promoter_location[2]-1)] <- "-"
annotation[promoter_location[2]] <- "|"

annotation[exon_1_start_location] <- "1"

annotation[uORF_start_location] <- "U"
annotation[(uORF_start_location+1):uORF_end_location] <- "-"



## Output file with header, sequence, and 3 reading frames
output_data <- c(">N2C uORF plasmid (GGC/GGA x97) sequence with forward reading frames",
                 paste(annotation, collapse = ""),
                 plasmid_GGC98,
                 insert_spaces_between_sequence(reading_frame_3, spaces_between = 2, initial_spaces = 3),
                 insert_spaces_between_sequence(reading_frame_1, spaces_between = 2, initial_spaces = 1),
                 insert_spaces_between_sequence(reading_frame_2, spaces_between = 2, initial_spaces = 2))

writeLines(output_data, "intermediate_files/n2c_uorf_3tag_plasmid_reading_frames.txt")
## -------------------------------------------------------------------------------------------------------
