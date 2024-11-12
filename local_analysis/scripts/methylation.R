## Creates qualitative (F3-10) and quantitative (F3-11) methylation plots

library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(raster)
library(ggnewscale)
library(ggtext)
library(ggpubr)
library(ggsignif)
library(gglm)
library(Rmpfr)
source("scripts/common_functions.r")

call_threshold <- 0.70
hyper_allele_threshold <- 800
long_allele_threshold  <- 200
methylation_data <- read.csv("input/methylation_information_all_participants.csv")
methylation_data$direction <- ifelse(methylation_data$direction == "forward", "Forward", "Reverse")

methylation_data$publication_id <- assign_publication_ids(methylation_data)



## CREATE METHYLATION VARIABLES
calculate_decimal_probability <- function(score_255, mode = "centre") {
    ## Choose numerator based on specified mode
    ## Per samtools specification, the probability space corresponding to integer N is
    ## from N/256 to (N+1)/256. This specifies where in that range to set the decimal probability
    if (mode == "lower") {numerator <- score_255}
    else if (mode == "centre") {numerator <- score_255 + 0.5}
    else if (mode == "upper")  {numerator <- score_255 + 1}
    
    return(numerator/256)
}

for (i in 1:nrow(methylation_data)) {
    methylation_probabilities_255        <- string_to_vector(methylation_data[i, "C_m_probabilities"])
    methylation_probabilities            <- calculate_decimal_probability(methylation_probabilities_255)
    hydroxymethylation_probabilities_255 <- string_to_vector(methylation_data[i, "C_hm_probabilities"])
    hydroxymethylation_probabilities     <- calculate_decimal_probability(hydroxymethylation_probabilities_255)
    any_modification_probabilities       <- calculate_decimal_probability(hydroxymethylation_probabilities_255 + methylation_probabilities_255)
    
    CpG_sites              <- length(methylation_probabilities)
    methylation_calls        <- sum(methylation_probabilities >= call_threshold)
    hydroxymethylation_calls <- sum(hydroxymethylation_probabilities >= call_threshold)
    any_modification_calls   <- sum(any_modification_probabilities >= call_threshold)
    
    length <- methylation_data[i, "length"]
    if (length >= hyper_allele_threshold) {methylation_data[i, "allele"] <- "Hyper"}
    else if (length >= long_allele_threshold) {methylation_data[i, "allele"] <- "Long"}
    else {methylation_data[i, "allele"] <- "Short"}
    
    # Call information
    methylation_data[i, "CpG_sites"]              <- CpG_sites    
    methylation_data[i, "call_threshold"]           <- call_threshold
    methylation_data[i, "methylation_calls"]        <- methylation_calls
    methylation_data[i, "hydroxymethylation_calls"] <- hydroxymethylation_calls
    methylation_data[i, "modification_calls"]       <- any_modification_calls
    methylation_data[i, "methylation_proportion"]        <- methylation_calls / CpG_sites
    methylation_data[i, "hydroxymethylation_proportion"] <- hydroxymethylation_calls / CpG_sites
    methylation_data[i, "modification_proportion"]       <- any_modification_calls / CpG_sites
    
    
    # Probability information
    methylation_data[i, "locations"]                        <- methylation_data[i, "C_m_locations"]
    methylation_data[i, "methylation_probabilities"]        <- vector_to_string(methylation_probabilities)
    methylation_data[i, "hydroxymethylation_probabilities"] <- vector_to_string(hydroxymethylation_probabilities)
    methylation_data[i, "modification_probabilities"]       <- vector_to_string(any_modification_probabilities)
    
    methylation_data[i, paste0("methylation_probabilities_", c("min", "q1", "median", "mean", "q3", "max"))] <- as.vector(summary(methylation_probabilities))
    methylation_data[i, paste0("hydroxymethylation_probabilities_", c("min", "q1", "median", "mean", "q3", "max"))] <- as.vector(summary(hydroxymethylation_probabilities))
    methylation_data[i, paste0("modification_probabilities_", c("min", "q1", "median", "mean", "q3", "max"))] <- as.vector(summary(any_modification_probabilities))
    
    methylation_data[i, "methylation_probabilities_sum"] <- sum(methylation_probabilities)
    methylation_data[i, "hydroxymethylation_probabilities_sum"] <- sum(hydroxymethylation_probabilities)
    methylation_data[i, "modification_probabilities_sum"] <- sum(any_modification_probabilities)
}


## Merge in phenotypes
phenotype_information <- read.csv("input/phenotype_information.csv")[, c("participant_id", "family_id", "condition", "phenotype")]
phenotype_information <- phenotype_information[!duplicated(phenotype_information), ]
phenotype_information$publication_id <- phenotype_information$participant_id
phenotype_information$participant_id <- NULL

merged_data <- merge(methylation_data, phenotype_information, by = "publication_id", all.x = T)
merged_data[merged_data$publication_id %in% c(paste0("MND", 1:7)), c("condition", "phenotype")] <- "MND"
merged_data[merged_data$phenotype == "Asymptomatic", "condition"] <- "Asymptomatic"
for (i in 1:7) {
    merged_data[merged_data$publication_id %in% c(paste0("MND", i)), "family_id"] <- paste0("MND", i)
}

write.csv(merged_data, "intermediate_files/methylation_data.csv", row.names = F)
merged_data$condition <- factor(merged_data$condition, c("NIID", "OPDM", "MND", "Asymptomatic"))
merged_data$phenotype <- factor(merged_data$phenotype, c("Asymptomatic", "OPDM", "NIID: Weakness", "NIID: Other/unknown", "MND"))





## Filter based on number of islands - graphs to justify
hist <- ggplot(merged_data, aes(x = CpG_sites)) +
    geom_histogram(binwidth = 1, boundary = 0, closed = "left") +
    geom_vline(xintercept = 6, col = "red", alpha = 0.5, linetype = "dashed") +
    labs(x = "CpG sites per read", y = "Read count") +
    theme_bw() + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())

scat <- ggplot(merged_data, aes(x = CpG_sites, y = length)) +
    geom_point() +
    geom_vline(xintercept = 6, col = "red", alpha = 0.5, linetype = "dashed") +
    labs(x = "CpG sites per read", y = "Trimmed read length (bp)") +
    theme_bw()

qc_plot <- egg::ggarrange(hist, scat, nrow = 2)
qc_plot
ggsave("output_figures/Appendix Figure 2-2 - Methylation QC.png", qc_plot, dpi = 600, width = 7, height = 5)
    
merged_data_filtered <- merged_data %>% filter(CpG_sites > 5)





## METHYLATION ILLUSTRATION
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
convert_modification_to_number_vector <- function(modification_locations, modification_probabilities, max_length) {
    locations <- string_to_vector(modification_locations)
    probabilities <- string_to_vector(modification_probabilities)
    
    output_vector <- rep(0, max_length)
    for (i in 1:max_length) {
        if (i %in% locations) {
            output_vector[i] <- probabilities[which(locations == i)]
        }
    }
    
    return(output_vector)
}

create_plot_masks <- function(modification_locations, modification_probabilities, max_length, sequence_length) {
    locations <- string_to_vector(modification_locations)
    probabilities <- string_to_vector(modification_probabilities)
    
    ## 0 = background
    ## 1 = non-CpG
    ## 2 = CpG
    output_vector <- rep(0, max_length)
    for (i in 1:max_length) {
        if (i %in% locations) {
            output_vector[i] <- 2
        } else if (i <= sequence_length) {
            output_vector[i] <- 1
        }
    }
    
    return(output_vector)
}

plot_image <- function(modification_locations_col, modification_probabilities_col, sequence_lengths, low_col = "blue", high_col = "red", background = "white", other_bases = "grey") {
    max_length <- max(sequence_lengths)
    image_matrix <- matrix(NA, nrow = length(modification_locations_col), ncol = max_length)
    mask_matrix  <- matrix(NA, nrow = length(modification_locations_col), ncol = max_length)
    for (i in 1:length(modification_locations_col)) {
        numeric_sequence_representation <- convert_modification_to_number_vector(modification_locations_col[i], modification_probabilities_col[i], max_length)
        image_matrix[i, ] <- numeric_sequence_representation
        
        mask <- create_plot_masks(modification_locations_col[i], modification_probabilities_col[i], max_length, sequence_lengths[i])
        mask_matrix[i, ] <- mask
    }
    image_data <- as.data.frame(raster(image_matrix), xy = TRUE)
    mask_data  <- as.data.frame(raster(mask_matrix), xy = TRUE)
    
    plot <- ggplot(mapping = aes(x = x, y = y)) +
        geom_raster(data = image_data, aes(fill = layer)) +
        scale_fill_gradient(low = low_col, high= high_col) +
        guides(fill = "none") +
        new_scale_fill() +
        geom_raster(data = mask_data, aes(fill = as.character(layer))) +
        scale_fill_manual(values = c("0" = background, "1" = other_bases, "2" = alpha("white", 0))) +
        coord_flip(expand = FALSE) +
        guides(x = "none", y = "none", fill = "none") +
        theme(axis.title = element_blank(), plot.margin=grid::unit(c(0,0,0,0), "mm"))
    
    return(plot)
}

construct_vector <- function(sequence_data, col_name, participant_gap = 2, allele_gap = 7) {
    sequences <- NULL
    for (this_allele in unique(sequence_data$allele)) {
        for (this_participant_id in unique(sequence_data[sequence_data$allele == this_allele, ]$participant_id)) {
            sequences <- c(sequences, sequence_data %>% filter(allele == this_allele, participant_id == this_participant_id) %>% pull(col_name), rep("", participant_gap))
        }
        sequences <- c(sequences, rep("", allele_gap))
    }
    sequences <- sequences[1:(length(sequences)-allele_gap-participant_gap)]
}


## Reverse sequences if direction == "Reverse"
reverse_sequence_if_needed <- function(sequence_vector, reverse_vector) {
    if (length(sequence_vector) != length(reverse_vector)) {
        print("Error: vectors need to be same length")
        return(NA)
    }
    
    new_sequence_vector <- rep(NA, length(sequence_vector))
    for (i in 1:length(sequence_vector)) {
        if (reverse_vector[i] == "Forward") {
            new_sequence_vector[i] <- sequence_vector[i]
        } else if (reverse_vector[i] == "Reverse") {
            new_sequence_vector[i] <- reverse_sequence(sequence_vector[i])
        }
    }
    return(new_sequence_vector)
}

## Reverse locations if direction == "Reverse"
reverse_locations_if_needed <- function(locations_vector, reverse_vector, length_vector) {
    if (length(locations_vector) != length(reverse_vector) || length(locations_vector) != length(length_vector)) {
        print("Error: vectors need to be same length")
        return(NA)
    }
    
    new_locations_vector <- rep(NA, length(locations_vector))
    for (i in 1:length(locations_vector)) {
        if (reverse_vector[i] == "Forward") {
            new_locations_vector[i] <- locations_vector[i]
        } else if (reverse_vector[i] == "Reverse") {
            length <- length_vector[i]
            reverse_positions <- as.numeric(strsplit(locations_vector[i], ",")[[1]])
            new_positions <- rev((length + 1) - reverse_positions)
            new_locations_vector[i] <- paste(new_positions, collapse = ",")
        }
    }
    return(new_locations_vector)
}

## Reverse probabilities if direction == "Reverse"
reverse_probabilities_if_needed <- function(probabilities_vector, reverse_vector, length_vector) {
    if (length(probabilities_vector) != length(reverse_vector) || length(probabilities_vector) != length(length_vector)) {
        print("Error: vectors need to be same length")
        return(NA)
    }
    
    new_probabilities_vector <- rep(NA, length(probabilities_vector))
    for (i in 1:length(probabilities_vector)) {
        if (reverse_vector[i] == "Forward") {
            new_probabilities_vector[i] <- probabilities_vector[i]
        } else if (reverse_vector[i] == "Reverse") {
            length <- length_vector[i]
            probabilities_to_reverse <- strsplit(probabilities_vector[i], ",")[[1]]
            new_probabilities_vector[i] <- paste(rev(probabilities_to_reverse), collapse = ",")
        }
    }
    return(new_probabilities_vector)
}

methylation_sorted <- merged_data_filtered %>%
    arrange(allele, publication_id, desc(length))

methylation_sorted$sequence_forward <- reverse_sequence_if_needed(methylation_sorted$sequence, methylation_sorted$direction)
methylation_sorted$C_m_locations_forward <- reverse_locations_if_needed(methylation_sorted$C_m_locations, methylation_sorted$direction, methylation_sorted$length)
methylation_sorted$C_m_probabilities_forward <- reverse_probabilities_if_needed(methylation_sorted$C_m_probabilities, methylation_sorted$direction, methylation_sorted$length)


locations     <- construct_vector(methylation_sorted, "C_m_locations_forward")
probabilities <- construct_vector(methylation_sorted, "C_m_probabilities_forward")
lengths       <- as.numeric(construct_vector(methylation_sorted, "length")) %>%
    replace_na(0)


plot_image(locations, probabilities, lengths, "blue", "red", "white", "lightgrey")
ggsave("output_figures/Figure 3-10a - Methylation illustration.png", dpi = 10, width = max(lengths), height = length(lengths), limitsize = FALSE)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------



## QUALITATIVE GRAPHS
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
point_alpha <- 0.85

## Methylation probabilities
mean_combined <- ggplot(merged_data_filtered, aes(x = length, y = methylation_probabilities_mean, col = phenotype, shape = direction)) +
    geom_point(alpha = point_alpha) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Trimmed read length over <i>NOTCH2NLC</i> repeat region (bp)", y = "Mean CpG methylation probability",
         col = "Phenotype", shape = "Read direction") +
    theme_bw() +
    scale_x_continuous(breaks = seq(0,2200,400)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    geom_vline(xintercept = 800, col = "red", linetype = "dashed") +
    geom_vline(xintercept = 200, col = "red", linetype = "dashed") +
    #ggtitle("(b) Mean CpG methylation probability") +
    theme_bw() + theme(axis.title.x = element_markdown())
mean_combined

mean_faceted <- ggplot(merged_data_filtered, aes(x = length, y = methylation_probabilities_mean, col = phenotype, shape = direction)) +
    geom_point(alpha = point_alpha) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Trimmed read length over <i>NOTCH2NLC</i> repeat region (bp)", y = "Mean CpG methylation probability",
         col = "Phenotype", shape = "Read direction") +
    scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
    geom_vline(xintercept = 800, col = "red", linetype = "dashed", linewidth = 0.25) +
    geom_vline(xintercept = 200, col = "red", linetype = "dashed", linewidth = 0.25) +
    facet_wrap(~publication_id) +
    #ggtitle("(c) Mean CpG methylation probability by participant") +
    theme_bw() + theme(axis.title.x = element_markdown())

sum_combined <- ggplot(merged_data_filtered, aes(x = length, y = methylation_probabilities_sum, col = phenotype, shape = direction)) +
    geom_point(alpha = point_alpha) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Trimmed read length over <i>NOTCH2NLC</i> repeat region (bp)", y = "Sum CpG methylation probability",
         col = "Phenotype", shape = "Read direction") +
    scale_x_continuous(breaks = seq(0,2200,400)) +
    ylim(0, 700) +
    geom_vline(xintercept = 800, col = "red", linetype = "dashed") +
    geom_vline(xintercept = 200, col = "red", linetype = "dashed") +
    #ggtitle("(d) Sum CpG methylation probability") +
    theme_bw() + theme(axis.title.x = element_markdown())
    
sum_faceted <- ggplot(merged_data_filtered, aes(x = length, y = methylation_probabilities_sum, col = phenotype, shape = direction)) +
    geom_point(alpha = point_alpha) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Trimmed read length over <i>NOTCH2NLC</i> repeat region (bp)", y = "Sum CpG methylation probability",
         col = "Phenotype", shape = "Read direction") +
    ylim(0, 700) +
    geom_vline(xintercept = 800, col = "red", linetype = "dashed", linewidth = 0.25) +
    geom_vline(xintercept = 200, col = "red", linetype = "dashed", linewidth = 0.25) +
    facet_wrap(~publication_id) +
    #ggtitle("(e) Sum CpG methylation probability per participant") +
    theme_bw() + theme(axis.title.x = element_markdown())
    
ggarrange(mean_combined, NULL, mean_faceted, NULL, NULL, NULL, sum_combined, NULL, sum_faceted, 
          nrow = 3, ncol = 3, heights = c(1, 0.1, 1), widths = c(1, 0.06, 1), common.legend = TRUE, legend = "right")
ggsave("output_figures/Figure 3-10bcde - Qualitative methylation graphs.png", dpi = 1400, width = 16, height = 10)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------













## PARAMETRIC (LMM) MODELS
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
merged_data_filtered$allele_releveled <- factor(merged_data_filtered$allele, levels = c("Short", "Long", "Hyper"))


## Mean model of allele and direction. No fixed effect of participant ID
methylation_mean_model_direction <- lmer(methylation_probabilities_mean ~ allele * direction + (1 | participant_id) + (1 | family_id), data = merged_data_filtered)
methylation_mean_model_direction_releveled <- lmer(methylation_probabilities_mean ~ allele_releveled * direction + (1 | participant_id) + (1 | family_id), data = merged_data_filtered)

## Sum model of allele and direction. No fixed effect of participant ID
methylation_sum_model_direction <- lmer(methylation_probabilities_sum ~ allele * direction + (1 | participant_id) + (1 | family_id), data = merged_data_filtered)
methylation_sum_model_direction_releveled <- lmer(methylation_probabilities_sum ~ allele_releveled * direction + (1 | participant_id) + (1 | family_id), data = merged_data_filtered)



## Mean model with phenotypes (now need to exclude suspicious barcode)
merged_data_filtered$suspicious_barcoding <- ifelse((merged_data_filtered$publication_id == "F2-I-1" & merged_data_filtered$allele == "Long") |
                                                        (merged_data_filtered$publication_id == "F2-II-3" & merged_data_filtered$allele == "Hyper"), TRUE, FALSE)
merged_data_unsuspicious <- merged_data_filtered %>% filter(suspicious_barcoding == FALSE)

methylation_mean_model_phenotype_long <- lmer(methylation_probabilities_mean ~ phenotype + (1 | direction) + (1 | participant_id) + (1 | family_id), data = merged_data_unsuspicious %>% filter(allele == "Long"))
methylation_mean_model_phenotype_short <- lmer(methylation_probabilities_mean ~ phenotype + (1 | direction) + (1 | participant_id) + (1 | family_id), data = merged_data_unsuspicious %>% filter(allele == "Short"))

## Sum model with phenotypes
methylation_sum_model_phenotype_long <- lmer(methylation_probabilities_sum ~ phenotype + (1 | direction) + (1 | participant_id) + (1 | family_id), data = merged_data_unsuspicious %>% filter(allele == "Long"))
methylation_sum_model_phenotype_short <- lmer(methylation_probabilities_sum ~ phenotype + (1 | direction) + (1 | participant_id) + (1 | family_id), data = merged_data_unsuspicious %>% filter(allele == "Short"))


## Save model outputs
methylation_model_anovas <- data.frame(response = rep(c("Mean", "Sum"), each = 5),
                                       explanatory = rep(c(rep("Allele/direction", 3), "Phenotype (long allele)", "Phenotype (short allele)"), times = 2),
                                       reads = rep(NA, 10), variable = rep(NA, 10), chisq = rep(NA, 10), df = rep(NA, 10), p_value = rep(NA, 10))

methylation_model_anovas[c(1:3, 6:8), "reads"] <- nrow(merged_data_filtered)


methylation_model_anovas[1:3, "variable"] <- rownames(Anova(methylation_mean_model_direction))
methylation_model_anovas[1:3, c("chisq", "df", "p_value")] <- Anova(methylation_mean_model_direction)

methylation_model_anovas[4, "variable"] <- rownames(Anova(methylation_mean_model_phenotype_long))
methylation_model_anovas[4, "reads"] <- nrow(merged_data_unsuspicious %>% filter(allele == "Long"))
methylation_model_anovas[4, c("chisq", "df", "p_value")] <- Anova(methylation_mean_model_phenotype_long)

methylation_model_anovas[5, "variable"] <- rownames(Anova(methylation_mean_model_phenotype_short))
methylation_model_anovas[5, "reads"] <- nrow(merged_data_unsuspicious %>% filter(allele == "Short"))
methylation_model_anovas[5, c("chisq", "df", "p_value")] <- Anova(methylation_mean_model_phenotype_short)

methylation_model_anovas[6:8, "variable"] <- rownames(Anova(methylation_sum_model_direction))
methylation_model_anovas[6:8, c("chisq", "df", "p_value")] <- Anova(methylation_sum_model_direction)

methylation_model_anovas[9, "variable"] <- rownames(Anova(methylation_mean_model_phenotype_long))
methylation_model_anovas[9, "reads"] <- nrow(merged_data_unsuspicious %>% filter(allele == "Long"))
methylation_model_anovas[9, c("chisq", "df", "p_value")] <- Anova(methylation_sum_model_phenotype_long)

methylation_model_anovas[10, "variable"] <- rownames(Anova(methylation_sum_model_phenotype_short))
methylation_model_anovas[10, "reads"] <- nrow(merged_data_unsuspicious %>% filter(allele == "Short"))
methylation_model_anovas[10, c("chisq", "df", "p_value")] <- Anova(methylation_sum_model_phenotype_short)

methylation_model_anovas[, c("chisq", "df", "p_value")] <- signif(methylation_model_anovas[, c("chisq", "df", "p_value")], 3)
write.table(methylation_model_anovas, "intermediate_files/methylation_model_anovas.csv", row.names = F)


## Assumptions
ggarrange(NULL, NULL, NULL,
          gglm(methylation_mean_model_direction, theme = ggplot2::theme_bw()), 
          gglm(methylation_mean_model_phenotype_long, theme = ggplot2::theme_bw()),
          gglm(methylation_mean_model_phenotype_short, theme = ggplot2::theme_bw()),  NULL, NULL, NULL,
          gglm(methylation_sum_model_direction, theme = ggplot2::theme_bw()),
          gglm(methylation_sum_model_phenotype_long, theme = ggplot2::theme_bw()),
          gglm(methylation_sum_model_phenotype_short, theme = ggplot2::theme_bw()),
          labels = c(rep(NA, 3), "(a) Mean ~ direction * allele", "(b) Mean ~ phenotype, long only", "(c) Mean ~ phenotype, short only",  
                     rep(NA, 3), "(d) Sum ~ direction * allele", "(e) Sum ~ phenotype, long only", "(f) Sum ~ phenotype, short only"),
          nrow = 4, ncol = 3, heights = c(0.1, 1, 0.1, 1, 0.1, 1), label.y = 1.05, label.x = 0)
ggsave("output_figures/Appendix Figure 2-3 - Methylation model assumptions.png", dpi = 600, width = 15, height = 10)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------





## MANN-WHITNEY
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
do_wilcox_comparisons <- function(data, response_var, group_var) {
    ## Determine set of combinations to be tested
    all_groups <- sort(unique(data[, group_var]))
    all_combinations_matrix <- combn(all_groups, 2)
    all_combinations <- list()
    for (i in 1:ncol(all_combinations_matrix)) {
        all_combinations[[i]] <- c(all_combinations_matrix[1, i], all_combinations_matrix[2, i])
    }
    
    output_data <- as.data.frame(matrix(NA, 0, 4))
    colnames(output_data) <- c("Group 1", "Group 2", "n_1", "n_2")
    
    ## Perform each comparison
    for (i in 1:length(all_combinations)) {
        comparison <- all_combinations[[i]]
        wilcox_results <- wilcox.test(data[data[, group_var] == comparison[1], response_var],
                                      data[data[, group_var] == comparison[2], response_var],
                                      exact = TRUE, paired = FALSE, alternative = "two.sided")
        
        output_data[i, "Group 1"] <- comparison[1]
        output_data[i, "Group 2"] <- comparison[2]
        output_data[i, "n_1"]     <- nrow(data[data[, group_var] == comparison[1], ])
        output_data[i, "n_2"]     <- nrow(data[data[, group_var] == comparison[2], ])
        output_data[i, "p_value"] <- wilcox_results$p.value
    }
    
    return(output_data)
}


merged_data_filtered$phenotype_character <- as.character(merged_data_filtered$phenotype)

comparisons_per_response <- 19
wilcox_p_values <- data.frame(response = rep(c("Mean", "Sum"), each = comparisons_per_response),
                              data = rep(NA, comparisons_per_response * 2),
                              group_1 = rep(NA, comparisons_per_response * 2), 
                              group_2 = rep(NA, comparisons_per_response * 2),
                              n_1 = rep(NA, comparisons_per_response * 2),
                              n_2 = rep(NA, comparisons_per_response * 2),
                              p_value = rep(NA, comparisons_per_response * 2))

wilcox_p_values[c(1:3, 1:3 + comparisons_per_response), "data"] <- c("Short reads, all QC'd", "Long reads, all QC'd", "Hyper reads, all QC'd")
wilcox_p_values[1, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Short"), "methylation_probabilities_mean", "direction")
wilcox_p_values[2, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Long"), "methylation_probabilities_mean", "direction")
wilcox_p_values[3, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Hyper"), "methylation_probabilities_mean", "direction")
wilcox_p_values[1 + comparisons_per_response, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Short"), "methylation_probabilities_sum", "direction")
wilcox_p_values[2 + comparisons_per_response, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Long"), "methylation_probabilities_sum", "direction")
wilcox_p_values[3 + comparisons_per_response, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Hyper"), "methylation_probabilities_sum", "direction")

wilcox_p_values[c(4:6, 4:6 + comparisons_per_response), "data"] <- "All QC'd reads"
wilcox_p_values[4:6, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(merged_data_filtered, "methylation_probabilities_mean", "allele")
wilcox_p_values[4:6 + comparisons_per_response, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(merged_data_filtered, "methylation_probabilities_sum", "allele")

wilcox_p_values[c(7:16, 7:16 + comparisons_per_response), "data"] <- "Short reads, suspicious barcoding removed"
wilcox_p_values[7:16, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Short", suspicious_barcoding == FALSE), "methylation_probabilities_mean", "phenotype_character")
wilcox_p_values[7:16 + comparisons_per_response, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Short", suspicious_barcoding == FALSE), "methylation_probabilities_sum", "phenotype_character")

wilcox_p_values[c(17:19, 17:19 + comparisons_per_response), "data"] <- "Long reads, suspicious barcoding removed"
wilcox_p_values[17:19, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Long", suspicious_barcoding == FALSE), "methylation_probabilities_mean", "phenotype_character")
wilcox_p_values[17:19 + comparisons_per_response, c("group_1", "group_2", "n_1", "n_2", "p_value")] <- do_wilcox_comparisons(filter(merged_data_filtered, allele == "Long", suspicious_barcoding == FALSE), "methylation_probabilities_sum", "phenotype_character")

wilcox_p_values$family_wide_p <- family_wide_p(wilcox_p_values$p_value, comparisons_per_response, "Sidak")
wilcox_p_values[, c("p_value", "family_wide_p")] <- signif(wilcox_p_values[, c("p_value", "family_wide_p")], 3)
write.csv(wilcox_p_values, "output_tables/Appendix Table 3-2 and 3-3 - Methylation Mann-Whitney P values.csv", row.names = F)


## Wilcox plots
adjusted_thresholds <- test_wise_alpha(c("***" = 0.001, "**" = 0.01, "*" = 0.05), comparisons_per_response, "Sidak")
point_alpha <- 0.85
set.seed(1234)

mean_direction_short <- ggplot(merged_data_filtered %>% filter(allele == "Short"), aes(x = direction, y = methylation_probabilities_mean)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    geom_signif(comparisons = list(c("Forward", "Reverse")),
                map_signif_level = adjusted_thresholds, vjust = 0.5,
                y_position = 0.65) +
    scale_y_continuous(limits = c(0, 1.175), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Read direction", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(i) Direction (wildtype reads)") +
    theme_bw()
mean_direction_short

mean_direction_long <- ggplot(merged_data_filtered %>% filter(allele == "Long"), aes(x = direction, y = methylation_probabilities_mean)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    scale_y_continuous(limits = c(0, 1.175), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Read direction", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(ii) Direction (expanded reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
mean_direction_long

mean_direction_hyper <- ggplot(merged_data_filtered %>% filter(allele == "Hyper"), aes(x = direction, y = methylation_probabilities_mean)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    scale_y_continuous(limits = c(0, 1.175), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Read direction", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(iii) Direction (hyperexpanded reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
mean_direction_hyper

mean_allele <- ggplot(merged_data_filtered, aes(x = factor(allele, levels = c("Short", "Long", "Hyper")), y = methylation_probabilities_mean)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    geom_signif(comparisons = list(c("Short", "Long"),
                                   c("Long", "Hyper"),
                                   c("Hyper", "Short")),
                map_signif_level = adjusted_thresholds, vjust = 0.5,
                y_position = seq(1, 1.2, 0.0625)) +
    scale_x_discrete(labels = c("Wildtype", "Expanded", "Hyperexpanded")) +
    scale_y_continuous(limits = c(0, 1.175), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Allele", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(iv) Allele (all reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
mean_allele

mean_phenotype_short <- ggplot(merged_data_filtered %>% filter(allele == "Short", suspicious_barcoding == FALSE), aes(x = phenotype, y = methylation_probabilities_mean)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    scale_x_discrete(labels = c("Asymptomatic", "ODPM", "NIID:\nWeakness", "NIID:\nOther/unknown", "MND")) +
    scale_y_continuous(limits = c(0, 1.175), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Phenotype", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(v) Phenotype (wildtype, non-suspicious reads)") +
    theme_bw()
mean_phenotype_short

mean_phenotype_long <- ggplot(merged_data_filtered %>% filter(allele == "Long", suspicious_barcoding == FALSE), aes(x = phenotype, y = methylation_probabilities_mean)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    geom_signif(comparisons = list(c("NIID: Weakness", "OPDM"),
                                   c("NIID: Weakness", "NIID: Other/unknown")),
                map_signif_level = adjusted_thresholds, vjust = 0.5,
                y_position = seq(1, 1.2, 0.0625)) +
    scale_x_discrete(labels = c("ODPM", "NIID:\nWeakness", "NIID:\nOther/unknown")) +
    scale_y_continuous(limits = c(0, 1.175), breaks = seq(0, 1, 0.25)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Phenotype", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(vi) Phenotype (expanded, non-suspicious reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
mean_phenotype_long

mean_direction_plots <- ggarrange(mean_direction_short, mean_direction_long, mean_direction_hyper, nrow = 1, common.legend = TRUE, legend = "none")
ggarrange(mean_direction_plots, mean_allele, mean_phenotype_short, mean_phenotype_long, widths = c(2, 1), nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("output_figures/Figure 3-11a - Mean methylation boxplots.png", dpi = 600, width = 17.5, height = 8)




set.seed(1234)

sum_direction_short <- ggplot(merged_data_filtered %>% filter(allele == "Short"), aes(x = direction, y = methylation_probabilities_sum)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 10)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Read direction", y = "Sum methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(i) Direction (wildtype reads)") +
    theme_bw()
sum_direction_short

sum_direction_long <- ggplot(merged_data_filtered %>% filter(allele == "Long"), aes(x = direction, y = methylation_probabilities_sum)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 1000, 20)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Read direction", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(ii) Direction (expanded reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
sum_direction_long

sum_direction_hyper <- ggplot(merged_data_filtered %>% filter(allele == "Hyper"), aes(x = direction, y = methylation_probabilities_sum)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 700, 200)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Read direction", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(iii) Direction (hyperexpanded reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
sum_direction_hyper

sum_allele <- ggplot(merged_data_filtered, aes(x = factor(allele, levels = c("Short", "Long", "Hyper")), y = methylation_probabilities_sum)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    geom_signif(comparisons = list(c("Short", "Long"),
                                   c("Long", "Hyper"),
                                   c("Hyper", "Short")),
                map_signif_level = adjusted_thresholds, vjust = 0.5,
                y_position = seq(600, 750, 50)) +
    scale_x_discrete(labels = c("Wildtype", "Expanded", "Hyperexpanded")) +
    scale_y_continuous(limits = c(0, 750), breaks = seq(0, 700, 200)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Allele", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(iv) Allele (all reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
sum_allele

signif_y_pos_start <- 18.5
sum_phenotype_short <- ggplot(merged_data_filtered %>% filter(allele == "Short", suspicious_barcoding == FALSE), aes(x = phenotype, y = methylation_probabilities_sum)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    geom_signif(comparisons = list(c("OPDM", "NIID: Other/unknown"),
                                   c("OPDM", "NIID: Weakness"),
                                   c("NIID: Other/unknown", "MND"),
                                   c("Asymptomatic", "OPDM"),
                                   c("NIID: Weakness", "MND"),
                                   c("OPDM", "MND"),
                                   c("Asymptomatic", "MND")),
                map_signif_level = adjusted_thresholds, vjust = 0.5,
                y_position = signif_y_pos_start + c(0, 2, 2, 4, 4, 6, 8)) +
    scale_x_discrete(labels = c("Asymptomatic", "ODPM", "NIID:\nWeakness", "NIID:\nOther/unknown", "MND")) +
    scale_y_continuous(limits = c(0, 28), breaks = seq(0, 20, 10)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Phenotype", y = "Sum methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(v) Phenotype (wildtype, non-suspicious reads)") +
    theme_bw()
sum_phenotype_short

sum_phenotype_long <- ggplot(merged_data_filtered %>% filter(allele == "Long", suspicious_barcoding == FALSE), aes(x = phenotype, y = methylation_probabilities_sum)) +
    geom_jitter(aes(col = phenotype, shape = direction), width = 0.3, alpha = point_alpha) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    geom_signif(comparisons = list(c("NIID: Weakness", "NIID: Other/unknown")),
                map_signif_level = adjusted_thresholds, vjust = 0.5,
                y_position = 95) +
    scale_x_discrete(labels = c("ODPM", "NIID:\nWeakness", "NIID:\nOther/unknown")) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "NIID: Weakness" = "#008C00", "NIID: Other/unknown" = "#781EF0" , "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    labs(x = "Phenotype", y = "Mean methylation probability", col = "Phenotype", shape = "Read direction") +
    ggtitle("(vi) Phenotype (expanded, non-suspicious reads)") +
    theme_bw() + theme(axis.title.y = element_blank())
sum_phenotype_long

sum_direction_plots <- ggarrange(sum_direction_short, sum_direction_long, sum_direction_hyper, nrow = 1, common.legend = TRUE, legend = "none")
ggarrange(sum_direction_plots, sum_allele, sum_phenotype_short, sum_phenotype_long, widths = c(2, 1), nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("output_figures/Figure 3-11b - Sum methylation boxplots.png", dpi = 600, width = 17.5, height = 8)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------






## PHENOTYPE CLASSES TABLE
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
methylation_phenotype_classes_summary_table <- merged_data_filtered %>%
    filter(suspicious_barcoding == FALSE) %>%
    group_by(allele, phenotype) %>%
    summarise(reads = n(), participants = length(unique(publication_id))) %>%
    arrange(desc(allele))

write.csv(methylation_phenotype_classes_summary_table, "output_tables/Appendix Table 2-4 - Methylation phenotype categories.csv", row.names = F)
## -----------------------------------------------------------------------------------------------------------------------------------------------------------
