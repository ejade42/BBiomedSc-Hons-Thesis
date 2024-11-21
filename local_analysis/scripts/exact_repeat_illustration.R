## Creates visualisation of exact trimmed repeats (Figure 3-8a)

library(tidyverse)
library(raster)
source("scripts/common_functions.R")

## A = 1, C = 2, G = 3, T = 4, blank = 0
convert_base_to_number <- function(base) {
    if (base == "A") {
        number <- 1
    } else if (base == "C") {
        number <- 2
    } else if (base == "G") {
        number <- 3
    } else if (base == "T") {
        number <- 4
    } else {
        print("Error: base must be A/T/C/G")
    }
    return(number)
}


convert_sequence_to_numbers <- function(sequence, length) {
    numerical_vector <- NULL
    for (i in 1:length) {
        if (i <= nchar(sequence)) {
            numerical_vector[i] <- convert_base_to_number(substr(sequence, i, i))
        } else {
            numerical_vector[i] <- 0
        }
    }
    return(numerical_vector)
}

plot_image <- function(sequences, colour_values) {
    max_length <- max(nchar(sequences))
    image_matrix <- matrix(NA, nrow = length(sequences), ncol = max_length)
    for (i in 1:length(sequences)) {
        numeric_sequence_representation <- convert_sequence_to_numbers(sequences[i], max_length)
        image_matrix[i, ] <- numeric_sequence_representation
    }
    image_data <- as.data.frame(raster(image_matrix), xy = TRUE)
    
    plot <- ggplot(image_data, aes(x = x, y = y, fill = as.character(layer))) +
        geom_raster() +
        scale_fill_manual(values = colour_values) +
        coord_flip(expand = FALSE) +
        guides(x = "none", y = "none", fill = "none") +
        theme(axis.title = element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"))
    
    return(plot)
}



#colour_values <- c("0" = "white", "1" = "yellow", "2" = "green", "3" = "blue", "4" = "red")
#colour_values <- c("0" = "white", "1" = "#F07800", "2" = "#008C00", "3" = "#006EDC", "4" = "#DC1E1E")
colour_values <- c("0" = "white", "1" = "#FFAA00", "2" = "#00BC00", "3" = "#0000DC", "4" = "#FF1E1E")


sequence_data_unsorted <- read.csv("intermediate_files/read_information_all_participants.csv") %>%
    filter(participant_id != "ishiura_2019_f9193_ii_5_raw")
    


sequence_data_unsorted$publication_id <- assign_publication_ids(sequence_data_unsorted)
sequence_data_unsorted$publication_id <- factor(sequence_data_unsorted$publication_id,
                                                c("F1-II-3", "F1-III-1", "F1-III-3",
                                                   "F2-I-1", "F2-II-1", "F2-II-2", "F2-II-3",
                                                   "F3-II-1", "F4-II-1", paste0("MND", 1:7),
                                                   "Ishiura 2019 F9193-II-5", "Podar 2023",
                                                   "Tian 2019 F1-IV-7", "Tian 2019 F1-IV-15", "Tian 2019 F2-II-3",
                                                   "Tian 2019 F4-II-2", "Tian 2019 F5-II-1", "Tian 2019 F5-II-4", "Tian 2019 F9-II-6"))




sequence_data <- sequence_data_unsorted %>%
    as.data.frame() %>%
    arrange(allele, publication_id, desc(exact_length))

participant_gap <- 2
allele_gap      <- 7
sequences <- NULL
for (this_allele in unique(sequence_data$allele)) {
    for (this_participant_id in unique(sequence_data[sequence_data$allele == this_allele, ]$participant_id)) {
        sequences <- c(sequences, sequence_data %>% filter(allele == this_allele, participant_id == this_participant_id) %>% pull(exact_repeat), rep("", participant_gap))
    }
    sequences <- c(sequences, rep("", allele_gap))
}
sequences <- sequences[1:(length(sequences)-allele_gap-participant_gap)]

plot_image(sequences, colour_values)

ggsave("output_figures/Figure 3-7a - Repeat illustration.png", dpi = 10, width = max(nchar(sequences)), height = length(sequences), limitsize = FALSE)
