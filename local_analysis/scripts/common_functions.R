## Function to reverse complement a DNA sequence
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


## Basic string and vector operations
string_to_vector  <- function(string) {as.numeric(unlist(strsplit(string, split = ",")))}
vector_to_string  <- function(vector) {paste(vector, collapse = ",")}


## Function to update all participant IDs using the new publication IDs
assign_publication_ids <- function(all_participant_data) {
    for (id in unique(all_participant_data$participant_id)) {
        if (id == "22NZ3355R") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F1-II-3"
        } else if (id == "20JK4316D") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F1-III-1"
        } else if (id == "21LE9817N") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F1-III-3"
        } else if (id %in% c("327", "327_2", "327_merged")) {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F2-I-1"
        } else if (id %in% c("313", "313_2", "313_merged")) {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F2-II-1"
        } else if (id %in% c("317", "317_2", "317_merged")) {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F2-II-2"
        } else if (id %in% c("318", "318_2", "318_merged")) {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F2-II-3"
        } else if (id %in% c("311", "311_2", "311_merged")) {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F3-II-1"
        } else if (id %in% c("319", "319_2", "319_merged")) {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "F4-II-1"
        } else if (id == "20JK4322C") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "MND1"
        } else if (id == "20JO4660C") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "MND2"
        } else if (id == "22PP0924B") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "MND3"
        } else if (id == "22RR5603Q") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "MND4"
        } else if (id == "22RS5734T") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "MND5"
        } else if (id == "GA5817C") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "MND6"
        } else if (id == "GQ1307H") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "MND7"
        } else if (id == "ishiura_2019_f9193_ii_5_corrected") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Ishiura 2019 F9193-II-5"
        } else if (id == "ishiura_2019_f9193_ii_5_raw") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Ishiura 2019 F9193-II-5 (raw)"
        } else if (id == "podar_2023_21073LRa006") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Podar 2023"
        } else if (id == "tian_2019_f1_iv_15") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Tian 2019 F1-IV-15"
        } else if (id == "tian_2019_f1_iv_7") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Tian 2019 F1-IV-7"
        } else if (id == "tian_2019_f2_ii_3") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Tian 2019 F2-II-3"
        } else if (id == "tian_2019_f4_ii_2") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Tian 2019 F4-II-2"
        } else if (id == "tian_2019_f5_ii_1") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Tian 2019 F5-II-1"
        } else if (id == "tian_2019_f5_ii_4") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Tian 2019 F5-II-4"
        } else if (id == "tian_2019_f9_ii_6") {
            all_participant_data[all_participant_data$participant_id == id, "publication_id"] <- "Tian 2019 F9-II-6"
        } else {
            print(paste("Error: participant", id, "has no publication ID"))
        }
    }

    return(all_participant_data$publication_id)
}




## Function to assign a batch and data origin to each participants
assign_participant_batches <- function(all_participant_data) {
    batch_1  <- c("20JK4316D", "21LE9817N", "22NZ3355R")
    batch_2  <- c("GA5817C", "GQ1307H", "20JK4322C", "20JO4660C", "22PP0924B", "22RR5603Q", "22RS5734T")
    batch_3_4  <- c("311", "313", "317", "318", "319", "327", "311_2", "313_2", "317_2", "318_2", "319_2", "327_2", "311_merged", "313_merged", "317_merged", "318_merged", "319_merged", "327_merged")
    ishiura  <- c("ishiura_2019_f9193_ii_5_raw", "ishiura_2019_f9193_ii_5_corrected")
    podar    <- c("podar_2023_21073LRa006")
    tian     <- c("tian_2019_f1_iv_7", "tian_2019_f1_iv_15", "tian_2019_f2_ii_3", "tian_2019_f4_ii_2", "tian_2019_f5_ii_1", "tian_2019_f5_ii_4", "tian_2019_f9_ii_6")
    
    
    for (id in unique(all_participant_data$participant_id)) {
        
        if (id %in% batch_1) {
            all_participant_data[all_participant_data$participant_id == id, "data_source"] <- "Novel (ours)"
            all_participant_data[all_participant_data$participant_id == id, "batch"] <- "Batch 1"
        } else if (id %in% batch_2) {
            all_participant_data[all_participant_data$participant_id == id, "data_source"] <- "Novel (ours)"
            all_participant_data[all_participant_data$participant_id == id, "batch"] <- "Batch 2"
        } else if (id %in% batch_3_4) {
            all_participant_data[all_participant_data$participant_id == id, "data_source"] <- "Novel (ours)"
            all_participant_data[all_participant_data$participant_id == id, "batch"] <- "Batch 3/4"
        } else if (id %in% ishiura) {
            all_participant_data[all_participant_data$participant_id == id, "data_source"] <- "Literature"
            all_participant_data[all_participant_data$participant_id == id, "batch"] <- "Ishiura"
        } else if (id %in% podar) {
            all_participant_data[all_participant_data$participant_id == id, "data_source"] <- "Literature"
            all_participant_data[all_participant_data$participant_id == id, "batch"] <- "Podar"
        } else if (id %in% tian) {
            all_participant_data[all_participant_data$participant_id == id, "data_source"] <- "Literature"
            all_participant_data[all_participant_data$participant_id == id, "batch"] <- "Tian"
        } else {
            print("Error: participant not in any group")
        }
    }
    return(list(data_source = all_participant_data$data_source, batch = all_participant_data$batch))
}



## Calculate sequence statistics
calculate_sequence_statistics <- function(sequence_data) {
    
    # Total length of sequence
    sequence_data$length_bp <- nchar(sequence_data$sequence)
    
    # Single-nucleotide proportions
    nucleotides_to_get_proportion <- c("A", "T", "C", "G")
    for (nucleotide in nucleotides_to_get_proportion) {
        var_name <- paste0("proportion_", nucleotide)
        sequence_data[[var_name]] <- str_count(sequence_data$sequence, nucleotide) / sequence_data$length_bp
    }
    
    # Number of each repeat unit (without worrying about frame)
    # GGG excluded as then some bases might be double-counted e.g. GGGT would count as GGG and GGT
    seqs_to_count <- c("GGC", "GGA", "GGT", "AGC")
    for (seq_of_interest in seqs_to_count) {
        var_name <- paste0("count_", seq_of_interest)
        sequence_data[[var_name]] <- str_count(sequence_data$sequence, seq_of_interest)
    }
    ## Check for GGAGC overlap
    sequence_data[["count_GGAGC"]] <- str_count(sequence_data$sequence, "GGAGC")
    
    # Total number of repeat units defined above
    sequence_data[["total_repeat_units"]] <- rowSums(sequence_data %>% dplyr::select(paste0("count_", seqs_to_count)))
    # Proportion of the total number of repeat units represented by each distinct repeat unit
    for (seq_of_interest in seqs_to_count) {
        proportion_name <- paste0("proportion_", seq_of_interest)
        count_name      <- paste0("count_", seq_of_interest)
        sequence_data[[proportion_name]] <- sequence_data[[count_name]] / sequence_data$total_repeat_units
    }
    # Proportion of the whole sequence represented by each repeat unit
    for (seq_of_interest in seqs_to_count) {
        proportion_name <- paste0("whole_seq_proportion_", seq_of_interest)
        count_name      <- paste0("count_", seq_of_interest)
        sequence_data[[proportion_name]] <- sequence_data[[count_name]] * 3 / sequence_data$length_bp
    }
    sequence_data[["whole_seq_proportion_total"]] <- rowSums(sequence_data %>% dplyr::select(paste0("whole_seq_proportion_", seqs_to_count)))
    
    return(sequence_data)
}



## Calculate family-wide equivalent of test-wise p value given number of comparisons
## Uses mpfr to avoid losing very small P values
family_wide_p <- function(p_value, comparisons, family = "Bonferroni") {
    if (family == "Bonferroni") {
        family_wide_p <- comparisons * mpfr(as.character(p_value), 256)
    } else if (family == "Sidak") {
        family_wide_p <- 1 - (1 - mpfr(as.character(p_value), 256)) ^ comparisons
    } else {print("Error: Family not recognised")}
    
    return(as.numeric(family_wide_p))
}


## Calculate test-wise significance threshold of chosen alpha given number of comparisons
test_wise_alpha <- function(alpha_value, comparisons, family = "Bonferroni") {
    if (family == "Bonferroni") {
        test_wise_alpha <- alpha_value / comparisons
    } else if (family == "Sidak") {
        test_wise_alpha <- 1 - (1 - alpha_value) ^ (1 / comparisons)
    } else {print("Error: Family not recognised")}
    
    return(test_wise_alpha)
}