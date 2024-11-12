## Runs in the individual participants' output directories
## Reads the clipped fastq file, and makes a csv of participant ID, sequence ID, sequence, and length

# Read participant ID from argument supplied when starting script
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2 || length(arguments) > 3) {stop("Please provide the participant ID and allele threshold when calling the script, and optionally a second hyper-expansion threshold")}

participant_id <- as.character(arguments[1])
threshold      <- as.numeric(arguments[2])
if (length(arguments) == 3) {
    hyper_threshold <- as.numeric(arguments[3])
} else {
    hyper_threshold <- Inf
}


# Function that takes a DNA sequence (read) and generates the reverse complement
reverse_sequence <- function(sequence) {
    sequence_vector     <- strsplit(sequence, split = "")[[1]]
    reversed_vector     <- rev(sequence_vector)
    new_sequence_vector <- rep(NA, length(reversed_vector))
    
    for (i in 1:length(reversed_vector)) {
        if (reversed_vector[i] == "A") {
            new_sequence_vector[i] <- "T"
        } else if (reversed_vector[i] == "C") {
            new_sequence_vector[i] <- "G"
        } else if (reversed_vector[i] == "G") {
            new_sequence_vector[i] <- "C"
        } else if (reversed_vector[i] == "T") {
            new_sequence_vector[i] <- "A"
        } else {
            stop("Cannot reverse read for non-A/C/G/T")
        }
    }
    
    new_sequence <- paste(new_sequence_vector, collapse = "")
    return(new_sequence)
}


# Read in forward and reverse reads
fastq_forward <- readLines(paste0(participant_id, "_4_NOTCH2NLC_clipped_reads_forward.fastq"))
fastq_reverse <- readLines(paste0(participant_id, "_4_NOTCH2NLC_clipped_reads_reverse.fastq"))
fastq_combined <- c(fastq_forward, fastq_reverse)


# Reformat, calculate length, and write forward sequence
blank <- rep(NA, length(fastq_combined) / 4)
participant_data <- data.frame(participant_id = participant_id, read_id = blank, sequence = blank, forward = blank, length = blank, forward_seq = blank)

for (i in 1:(length(fastq_combined)/4)) {
    participant_data[i, "read_id"] <- fastq_combined[4*i - 3]
    sequence <- fastq_combined[4*i - 2]
    participant_data[i, "sequence"] <- sequence
    participant_data[i, "length"]   <- nchar(sequence)
    
    if (nchar(sequence) >= hyper_threshold) {
        participant_data[i, "allele"] <- "hyper"
    } else if (nchar(sequence) >= threshold) {
        participant_data[i, "allele"] <- "long"
    } else {
        participant_data[i, "allele"] <- "short"
    }
    
    if (i <= length(fastq_forward)/4) {
        participant_data[i, "forward"]     <- "forward"
        participant_data[i, "forward_seq"] <- sequence
    } else {
        participant_data[i, "forward"]     <- "reverse"
        participant_data[i, "forward_seq"] <- reverse_sequence(sequence)
    }
}

write.csv(participant_data, paste0(participant_id, "_read_information.csv"), row.names = F)