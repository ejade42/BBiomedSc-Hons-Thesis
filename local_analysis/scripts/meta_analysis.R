## Performs meta-analysis of genotype-phenotype relationships

library(tidyverse)
library(ggsignif)
library(ggtext)
library(RGenetics)
library(ggnewscale)
library(lme4)
library(ggpubr)
library(car)
library(gglm)
library(Rmpfr)
source("scripts/common_functions.R")
set.seed(1234)

## Reading in data
##--------------------------------------------------------------------------------------------------
# Read in manually entered information about sequences
manual_info <- read_csv("input/phenotype_information.csv")

# Read in sequence information from Python output
python_sequences  <- rbind(read_csv("intermediate_files/python_consensuses_literature.csv"),
                           read_csv("intermediate_files/python_consensuses_pipeline.csv"))
colnames(python_sequences)[1] <- "sequence_id"

# Merge together
sequence_data <- as_tibble(merge(manual_info, python_sequences, all.x = FALSE, all.y = FALSE))
rm(manual_info, python_sequences)

# Remove participants marked for exclusion (i.e. duplicates)
sequence_data <- sequence_data %>% 
    filter(include == TRUE) %>% 
    replace_na(list(ancestry_group = "Unknown"))

## Create summary variables
sequence_data <- calculate_sequence_statistics(sequence_data)

## Classify phenotypes as major or minor based on n >= 5
for (phenotype in sequence_data$phenotype) {
    if (nrow(sequence_data[sequence_data$phenotype == phenotype, ]) >= 5) {
        sequence_data[sequence_data$phenotype == phenotype, "major_class"] <- "Major"
    } else {
        sequence_data[sequence_data$phenotype == phenotype, "major_class"] <- "Minor"
    }
}

## Create factors with levels in specified order
sequence_data$phenotype_factor <- factor(sequence_data$phenotype, c("Asymptomatic",
                                                                    "OPDM",
                                                                    "HDMNM",
                                                                    "Parkinson's",
                                                                    "NIID: Parkinsonism",
                                                                    "NIID: Renal",
                                                                    "NIID: Weakness",
                                                                    "NIID: Dementia", 
                                                                    "NIID: Paroxysmal",
                                                                    "NIID: Movement disorder",
                                                                    "NIID: Other/unknown"))

sequence_data$condition_factor <- factor(sequence_data$condition, c("OPDM", "HDMNM", "PD", "NIID"))
##--------------------------------------------------------------------------------------------------



## MANN-WHITNEY
##--------------------------------------------------------------------------------------------------
## Calculate weight as reciprocal of number of rows for each level of the chosen column
weight_by_column <- function(data_to_weight, grouping_variable) {
    for (group in unique(pull(data_to_weight, grouping_variable))) {
        #print(group)
        count <- nrow(data_to_weight[pull(data_to_weight, grouping_variable) == group, ])
        #print(count)
        weight <- 1 / count
        #print(weight)
        data_to_weight[pull(data_to_weight, grouping_variable) == group, "weight"] <- weight
    }
    return(data_to_weight)
}


sequence_data_filtered <- sequence_data %>% 
    filter(suspicious_barcoding == FALSE) %>%
    weight_by_column(., "participant_id")

phenotype_summary_table <- sequence_data_filtered %>%
    group_by(phenotype_factor) %>%
    summarise(alleles = n(), participants = sum(weight))
write.csv(phenotype_summary_table, "output_tables/Appendix Table 2-5 - Meta-analysis phenotype categories.csv", row.names = F)





## Function to calculate P values for all pairs, for a given column
all_wilcoxon_tests <- function(data, all_combinations, column) {
    p_values   <- rep(NA, length(all_combinations))
    statistics <- rep(NA, length(all_combinations))
    
    for (i in 1:length(all_combinations)) {
        comparison <- all_combinations[[i]]
        wilcox_results <- wilcox.test(filter(data, phenotype_factor == comparison[1]) %>% pull(column),
                                      filter(data, phenotype_factor == comparison[2]) %>% pull(column),
                                      exact = TRUE, paired = FALSE, alternative = "two.sided")
        
        p_values[i]   <- wilcox_results$p.value
        statistics[i] <- wilcox_results$statistic
    }
    
    return(list(p_values = p_values, statistics = statistics))
}


## Calculating all P-values
major_phenotypes <- c("Asymptomatic", "OPDM", "NIID: Weakness", "NIID: Dementia", "NIID: Paroxysmal", "NIID: Other/unknown")
all_combinations_matrix <- combn(major_phenotypes, 2)
all_combinations <- list()
for (i in 1:ncol(all_combinations_matrix)) {
    all_combinations[[i]] <- c(all_combinations_matrix[1, i], all_combinations_matrix[2, i])
}


all_p_values <- data.frame(phenotype_1 = all_combinations_matrix[1, ], phenotype_2 = all_combinations_matrix[2, ])
columns_to_test <- c(c("length_bp", "whole_seq_proportion_GGA", "whole_seq_proportion_GGC", "whole_seq_proportion_GGT", "whole_seq_proportion_AGC", "whole_seq_proportion_total", "proportion_A", "proportion_C", "proportion_G", "proportion_T"))
for (column in columns_to_test) {
    all_p_values[, column] <- all_wilcoxon_tests(sequence_data_filtered, all_combinations, column)[["p_values"]]
}


## Calculating adjusted significance
original_thresholds <- c("***" = 0.001, "**" = 0.01, "*" = 0.05)

adjusted_thresholds_15 <- test_wise_alpha(original_thresholds, 15, "Sidak")
all_significances_15 <- all_p_values
for (i in 1:nrow(all_p_values)) {
    for (j in 3:ncol(all_p_values)) {
        if (!is.na(all_p_values[i, j])) {
            if (all_p_values[i, j] <= adjusted_thresholds_15[1]) {all_significances_15[i, j] <- names(adjusted_thresholds_15)[1]}
            else if (all_p_values[i, j] <= adjusted_thresholds_15[2]) {all_significances_15[i, j] <- names(adjusted_thresholds_15)[2]}
            else if (all_p_values[i, j] <= adjusted_thresholds_15[3]) {all_significances_15[i, j] <- names(adjusted_thresholds_15)[3]}
            else {all_significances_15[i, j] <- ""}
        }
    }
}

adjusted_p_values <- all_p_values[, 1:2]
for (i in 3:ncol(all_p_values)) {
    adjusted_p_values[, i] <- family_wide_p(all_p_values[, i], 15, "Sidak")
}
colnames(adjusted_p_values) <- colnames(all_p_values)

adjusted_p_values_rounded <- cbind(adjusted_p_values[, 1:2], signif(adjusted_p_values[, 3:ncol(adjusted_p_values)], 3))
all_p_values_rounded <- cbind(all_p_values[, 1:2], signif(all_p_values[, 3:ncol(all_p_values)], 3))


write.csv(all_p_values_rounded, "output_tables/Appendix Table 3-4 - Unadjusted meta-analysis comparisons P values.csv", row.names = F)
write.csv(adjusted_p_values_rounded, "output_tables/Appendix Table 3-5 - Adjusted meta-analysis comparisons P values.csv", row.names = F)
##--------------------------------------------------------------------------------------------------



## BOXPLOTS
##--------------------------------------------------------------------------------------------------
point_alpha <- 0.8
thick_stroke <- 1.15

## Length
length <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = length_bp)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    geom_signif(comparisons = list(c("Asymptomatic", "OPDM"), 
                                   c("Asymptomatic", "NIID: Weakness"),
                                   c("Asymptomatic", "NIID: Dementia"),
                                   c("Asymptomatic", "NIID: Paroxysmal"),
                                   c("Asymptomatic", "NIID: Other/unknown")),
                map_signif_level = adjusted_thresholds_15, vjust = 0.5,
                y_position = seq(1700, 2200, 125)) +
    geom_signif(comparisons = list(c("NIID: Weakness", "NIID: Paroxysmal")),
                map_signif_level = adjusted_thresholds_15, vjust = 0.5,
                y_position = 1575) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "Repeat length (bp)", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    ylim(0, 2300) +
    theme_bw() + theme(axis.title.x = element_blank())
length


## GGA proportion
gga <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = whole_seq_proportion_GGA)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    geom_signif(comparisons = list(c("NIID: Weakness", "NIID: Dementia"), 
                                   c("NIID: Weakness", "NIID: Paroxysmal"),
                                   c("NIID: Weakness", "NIID: Other/unknown"),
                                   c("NIID: Weakness", "OPDM"),
                                   c("NIID: Weakness", "Asymptomatic")),
                map_signif_level = adjusted_thresholds_15, vjust = 0.5,
                y_position = seq(0.3, 0.4, 0.025)) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "GGA proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    ylim(-0.01, 0.42) +
    theme_bw() + theme(axis.title.x = element_blank())
gga


## GGC proportion
ggc <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = whole_seq_proportion_GGC)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    geom_signif(comparisons = list(c("NIID: Weakness", "NIID: Dementia")),
                map_signif_level = adjusted_thresholds_15, vjust = 0.5,
                y_position = 1.01) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "GGC proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(limits = c(-0.01, 1.06), breaks = seq(0, 1, 0.2)) +
    theme_bw() + theme(axis.title.x = element_blank())
ggc


### GGT proportion
ggt <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = whole_seq_proportion_GGT)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "GGT proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(breaks = seq(0, 0.6, 0.2)) +
    theme_bw() + theme(axis.title.x = element_blank())
ggt


### AGC proportion
agc <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = whole_seq_proportion_AGC)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "AGC proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(limits = c(-0.01, 0.3), breaks = seq(0, 0.3, 0.1)) +
    theme_bw() + theme(axis.title.x = element_blank())
agc


### Total unit proportion
units <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = whole_seq_proportion_total)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "GGA+GGC+GGT+AGC proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(limits = c(-0.01, 1.01), breaks = seq(0, 1, 0.2)) +
    theme_bw() + theme(axis.title.x = element_blank())
units


## A proportion
a <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = proportion_A)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    geom_signif(comparisons = list(c("NIID: Weakness", "Asymptomatic"), 
                                   c("NIID: Weakness", "NIID: Dementia"),
                                   c("NIID: Weakness", "NIID: Paroxysmal")),
                map_signif_level = adjusted_thresholds_15, vjust = 0.5,
                y_position = seq(0.25, 0.4, 0.025)) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "A proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(limits = c(-0.01, 0.32), breaks = seq(0, 0.4, 0.1)) +
    theme_bw() + theme(axis.title.x = element_blank())
a

## C proportion
c <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = proportion_C)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    geom_signif(comparisons = list(c("NIID: Weakness", "NIID: Dementia"), 
                                   c("NIID: Weakness", "NIID: Paroxysmal"),
                                   c("NIID: Weakness", "NIID: Other/unknown"),
                                   c("NIID: Weakness", "OPDM"),
                                   c("NIID: Weakness", "Asymptomatic")),
                map_signif_level = adjusted_thresholds_15, vjust = 0.5,
                y_position = seq(0.375, 0.5, 0.025)) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "C proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(limits = c(-0.01, 0.5), breaks = seq(0, 0.5, 0.1)) +
    theme_bw() + theme(axis.title.x = element_blank())
c

## G proportion
g <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = proportion_G)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "G proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
    theme_bw() + theme(axis.title.x = element_blank())
g

## T proportion
t <- ggplot(sequence_data_filtered, aes(x = phenotype_factor, y = proportion_T)) +
    geom_boxplot(aes(col = major_class), alpha = 0, outlier.shape = NA) +
    scale_colour_manual(values = c("Major" = "black", "Minor" = "#888")) +
    labs(col = "Phenotype\ncategory (n ≥ 5)") +
    new_scale_colour() +
    geom_jitter(aes(col = condition_factor, shape = ancestry_group, stroke = ancestry_group), width = 0.3, alpha = point_alpha) +
    scale_shape_manual(values = c("East Asian" = 19, "Unknown" = 1, "Polynesian" = 4, "European" = 6)) +
    scale_discrete_manual(aesthetics = "stroke", values = c("East Asian" = 1, "Unknown" = thick_stroke, "Polynesian" = thick_stroke, "European" = thick_stroke)) +
    scale_colour_manual(values = c("OPDM" = "#DC1E1E", "HDMNM" = "#006EDC", "PD" = "#DC00DC", "NIID" = "#008C00")) +
    scale_x_discrete(labels = c("Asymptomatic", "OPDM", "HDMNM", "Parkinson's", "NIID:\nParkinsonism", "NIID:\nRenal", "NIID:\nWeakness", "NIID:\nDementia", "NIID:\nParoxysmal", "NIID:\nMovement\ndisorder", "NIID:\nOther/\nunknown")) +
    labs(x = "Phenotype", y = "T proportion", col = "Disease", shape = "Ancestry", stroke = "Ancestry") +
    scale_y_continuous(limits = c(-0.01, 0.35), breaks = seq(0, 0.4, 0.1)) +
    theme_bw() + theme(axis.title.x = element_blank())
t

set.seed(1234)
#ggarrange(length, units, gga, ggc, agc, ggt, a, c, g, t, nrow = 5, ncol = 2, common.legend = T, legend = "right")
#ggsave("Graphs/241029_repeat_structure_all.png", dpi = 600, width = 17.5, height = 20)

ggarrange(length + ggtitle("(a) Repeat length"), agc + ggtitle("(b) AGC proportion"), 
          gga + ggtitle("(c) GGA proportion"), ggc + ggtitle("(d) GGC proportion"), 
          a + ggtitle("(e) A proportion"), c + ggtitle("(f) C proportion"), 
          nrow = 3, ncol = 2, common.legend = T, legend = "right")
ggsave("output_figures/Figure 3-12 - Meta-analysis six variables.png", dpi = 600, width = 17.5, height = 13.5)

ggarrange(units + ggtitle("(a) Repeat unit proportion"), ggt + ggtitle("(b) GGT proportion"), 
          g + ggtitle("(c) G proportion"), t + ggtitle("(d) T proportion"), 
          nrow = 2, ncol = 2, common.legend = T, legend = "right")
ggsave("output_figures/Figure 3-12 - Meta-analysis four variables.png", dpi = 600, width = 17.5, height = 9)

ggarrange(length + ggtitle("(a) Repeat length"), units + ggtitle("(b) Repeat unit proportion"), 
          gga + ggtitle("(c) GGA proportion"), ggc + ggtitle("(d) GGC proportion"), 
          agc + ggtitle("(e) AGC proportion"), ggt + ggtitle ("(f) GGT proportion"),
          a + ggtitle("(g) A proportion"), c + ggtitle("(h) C proportion"), 
          g + ggtitle("(i) G proportion"), t + ggtitle("(j) T proportion"),
          nrow = 5, ncol = 2, common.legend = T, legend = "right")
ggsave("output_figures/Figure 3-12 - Meta-analysis ten variables.png", dpi = 600, width = 17.5, height = 20)
##--------------------------------------------------------------------------------------------------




## WEIGHTED MIXED-EFFECTS LINEAR MODELS
##--------------------------------------------------------------------------------------------------
test_wise_alpha(original_thresholds, 30, "Sidak")

## Fit weighted mixed effects models with ancestry
ancestry_model_length <- lmer(length_bp ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_gga <- lmer(whole_seq_proportion_GGA ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_ggc <- lmer(whole_seq_proportion_GGC ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
## GGT fails to converge
ancestry_model_ggt <- lmer(whole_seq_proportion_GGT ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_agc <- lmer(whole_seq_proportion_AGC ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_units <- lmer(whole_seq_proportion_total ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_a <- lmer(proportion_A ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_c <- lmer(proportion_C ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_g <- lmer(proportion_G ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)
ancestry_model_t <- lmer(proportion_T ~ phenotype_factor * ancestry_group + (1 | family_id) + (1 | participant_id), weights = weight, data = sequence_data_filtered)



ancestry_model_output <- data.frame(Model = rep(c("Repeat length", "GGA proportion", "GGC proportion", "AGC proportion", "Repeat unit proportion", "A proportion", "C proportion", "G proportion", "T proportion"), each = 3),
                                    Variable = rep(c("Phenotype", "Ancestry", "Interaction"), times = 9),
                                    Chisq = rep(NA, 27),
                                    Df = rep(NA, 27),
                                    P = rep(NA, 27))
ancestry_model_output[1:3, 3:5] <- Anova(ancestry_model_length)
ancestry_model_output[4:6, 3:5] <- Anova(ancestry_model_gga)
ancestry_model_output[7:9, 3:5] <- Anova(ancestry_model_ggc)
ancestry_model_output[10:12, 3:5] <- Anova(ancestry_model_agc)
ancestry_model_output[13:15, 3:5] <- Anova(ancestry_model_units)
ancestry_model_output[16:18, 3:5] <- Anova(ancestry_model_a)
ancestry_model_output[19:21, 3:5] <- Anova(ancestry_model_c)
ancestry_model_output[22:24, 3:5] <- Anova(ancestry_model_g)
ancestry_model_output[25:27, 3:5] <- Anova(ancestry_model_t)

ancestry_model_output$Family_Wide_P <- family_wide_p(ancestry_model_output$P, 30, "Sidak")
ancestry_model_output_rounded <- cbind(ancestry_model_output[1:2], signif(ancestry_model_output[3:6], 3))
write.csv(ancestry_model_output_rounded, "output_tables/Table 3-1 - Meta-analysis linear mixed model outputs.csv", row.names = F)


## Check assumptions
ggarrange(NULL, NULL, NULL,
          gglm(ancestry_model_length, theme = ggplot2::theme_bw()), 
          gglm(ancestry_model_gga, theme = ggplot2::theme_bw()),
          gglm(ancestry_model_ggc, theme = ggplot2::theme_bw()),  NULL, NULL, NULL,
          gglm(ancestry_model_agc, theme = ggplot2::theme_bw()),
          gglm(ancestry_model_units , theme = ggplot2::theme_bw()),
          gglm(ancestry_model_a, theme = ggplot2::theme_bw()),  NULL, NULL, NULL,
          gglm(ancestry_model_c, theme = ggplot2::theme_bw()), 
          gglm(ancestry_model_g, theme = ggplot2::theme_bw()),
          gglm(ancestry_model_t, theme = ggplot2::theme_bw()),
          labels = c(rep(NA, 3), "(a) Repeat length", "(b) GGA proportion", "(c) GGC proportion",  rep(NA, 3),
                     "(d) AGC proportion", "(e) Repeat unit proportion", "(f) A proportion",  rep(NA, 3),
                     "(g) C proportion", "(h) G proportion", "(i) T proportion"),
          nrow = 6, ncol = 3, heights = c(0.1, 1, 0.1, 1, 0.1, 1), label.y = 1.05)
ggsave("output_figures/Appendix Figure 2-4 - Meta-analysis model assumptions.png", dpi = 600, width = 20, height = 15)
##--------------------------------------------------------------------------------------------------

