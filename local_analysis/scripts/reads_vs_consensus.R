## Creates Appendix Tables 2-3 (quantitative consensus method comparison) and 3-1 (read depth and consensus success summary),
## and Figure 2-7 (dimensionality reduction to compare consensus methods)

library(tidyverse)
library(ggnewscale)
library(ggpubr)
library(tsne)
source("scripts/common_functions.R")

read_data <- read.csv("intermediate_files/read_information_all_participants.csv")

read_data$sequence <- read_data$exact_repeat
read_data$sequence_id <- read_data$read_id
read_data$type <- "Read"
read_data$method <- "None"

consensus_data <- read.csv("intermediate_files/python_consensuses_pipeline.csv")
consensus_info <- read.csv("input/manual_pipeline_consensus_info.csv")
consensus_data_final <- merge(consensus_data, consensus_info)
consensus_data_final$sequence_id <- consensus_data_final$consensus_id

sequence_data <- merge(read_data, consensus_data_final, all = T)[ , c("participant_id", "sequence_id", "sequence", "type", "method", "allele")]
sequence_data <- sequence_data %>% filter(method != "Garvan", participant_id != "ishiura_2019_f9193_ii_5_raw")
rm(consensus_data, consensus_info, consensus_data_final, read_data)

sequence_data$publication_id <- assign_publication_ids(sequence_data)
sequence_data$publication_id <- factor(sequence_data$publication_id,
                                       c("F1-II-3", "F1-III-1", "F1-III-3",
                                         "F2-I-1", "F2-II-1", "F2-II-2", "F2-II-3",
                                         "F3-II-1", "F4-II-1", paste0("MND", 1:7),
                                         "Ishiura 2019 F9193-II-5", "Podar 2023",
                                         "Tian 2019 F1-IV-7", "Tian 2019 F1-IV-15", "Tian 2019 F2-II-3",
                                         "Tian 2019 F4-II-2", "Tian 2019 F5-II-1", "Tian 2019 F5-II-4", "Tian 2019 F9-II-6"))

sequence_data$allele_original <- sequence_data$allele
sequence_data$allele <- convert_allele_names(sequence_data$allele_original)
sequence_data$allele <- factor(sequence_data$allele, c("Wildtype", "Expanded", "Hyperexpanded"))

sequence_data <- calculate_sequence_statistics(sequence_data)
write.csv(sequence_data, "intermediate_files/consensus_read_data_combined.csv", row.names = F)

sequence_data_for_pca <- sequence_data[ , c("length_bp", "proportion_A", "proportion_C", "proportion_G", "proportion_T", "whole_seq_proportion_GGA", "whole_seq_proportion_GGC", "whole_seq_proportion_GGT", "whole_seq_proportion_AGC", "whole_seq_proportion_total")] %>%
    scale()


## Checking whether GGAGC double-counting is an issue:
print(paste("Number of GGAGC:", sum(sequence_data$count_GGAGC)), quote = F)
print(paste("Sequences with GGAGC:", nrow(sequence_data %>% filter(count_GGAGC != 0))), quote = F)

## Compare to GGC:
print(paste("Number of GGC:", sum(sequence_data$count_GGC)), quote = F)
print(paste("Sequences with GGC:", nrow(sequence_data %>% filter(count_GGC != 0))), quote = F)



## Summary table of read depth and consensus success
consensus_summary_table <- sequence_data %>%
    group_by(publication_id, allele, method) %>%
    summarise(number = n()) %>%
    pivot_wider(id_cols = c(publication_id, allele), names_from = method, values_from = number) %>%
    relocate(c(publication_id, allele, None, Local, Global, Flye)) %>%
    arrange(publication_id, desc(allele))


read_depth_means <- c(colMeans(consensus_summary_table[consensus_summary_table$allele == "Wildtype", "None"], na.rm = T),
                      colMeans(consensus_summary_table[consensus_summary_table$allele == "Expanded", "None"], na.rm = T),
                      colMeans(consensus_summary_table[consensus_summary_table$allele == "Hyperexpanded", "None"], na.rm = T),
                      sum(consensus_summary_table$None) / length(unique(consensus_summary_table$publication_id)))
names(read_depth_means) <- c("Wildtype", "Expanded", "Hyperexpanded", "Total")
print("Mean read depth by allele:", quote = F)
print(read_depth_means, quote = F)

new_rows <- data.frame("publication_id" = rep(c("Total"), each = 4), 
                       "allele" = rep(c("Wildtype", "Expanded", "Hyperexpanded", "Total"), times = 1))
new_rows[1, c("None", "Local", "Global", "Flye")] <- colSums(consensus_summary_table[consensus_summary_table$allele == "Wildtype", c("None", "Local", "Global", "Flye")], na.rm = T)
new_rows[2, c("None", "Local", "Global", "Flye")] <- colSums(consensus_summary_table[consensus_summary_table$allele == "Expanded", c("None", "Local", "Global", "Flye")], na.rm = T)
new_rows[3, c("None", "Local", "Global", "Flye")] <- colSums(consensus_summary_table[consensus_summary_table$allele == "Hyperexpanded", c("None", "Local", "Global", "Flye")], na.rm = T)
new_rows[4, c("None", "Local", "Global", "Flye")] <- colSums(consensus_summary_table[, c("None", "Local", "Global", "Flye")], na.rm = T)

consensus_summary_table <- rbind(consensus_summary_table, new_rows)
colnames(consensus_summary_table) <- c("Participant Id", "Allele", "Reads", "Flye consensus", "Global consensus", "Local consensus")
write.csv(consensus_summary_table, "output_tables/Appendix Table 3-1 - Read depth and consensus success summary.csv", row.names = F)
rm(consensus_summary_table, new_rows)


## Define colours
allele_colours <- c("Wildtype" = "#800030", "Expanded" = "darkblue", "Hyperexpanded" = "darkgreen")
method_colours <- c("Flye" = "red", "Global" = "orange", "Local" = "#FB61D7")


## GGA and repeat length plot
medians <- sequence_data %>%
    filter(type == "Read") %>%
    group_by(publication_id, allele) %>%
    summarise(length_bp = median(length_bp), 
              whole_seq_proportion_AGC = median(whole_seq_proportion_AGC),
              whole_seq_proportion_GGC = median(whole_seq_proportion_GGC),
              whole_seq_proportion_GGA = median(whole_seq_proportion_GGA),
              whole_seq_proportion_GGT = median(whole_seq_proportion_GGT))


bp_gga <- ggplot(sequence_data, aes(x = length_bp, y = whole_seq_proportion_GGA)) +
    geom_vline(data = medians, aes(xintercept = length_bp, col = allele),
               linetype = "dashed", alpha = 0.3) +
    geom_hline(data = medians, aes(yintercept = whole_seq_proportion_GGA, col = allele),
               linetype = "dashed", alpha = 0.3) +
    geom_point(data = sequence_data %>% filter(type == "Read"), aes(col = allele), 
               size = 1.5, alpha = 0.7, stroke = 0.75) +
    scale_colour_manual(values = allele_colours) +
    labs(col = "Allele type") +
    new_scale_colour() +
    geom_point(data = sequence_data %>% filter(type == "Consensus") %>% arrange(desc(method)),
               aes(col = method, shape = method), size = 2, stroke = 1) +
    scale_colour_manual(values = method_colours) +
    scale_shape_manual(values = c("Flye" = 4, "Global" = 10, "Local" = 2)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(0,2200,500), limits = c(0,2200)) + ylim(0, 0.32) +
    facet_wrap(~publication_id, scales = "fixed") +
    guides(alpha = "none") +
    ggtitle("(a) Repeat length vs GGA") +
    labs(x = "Total repeat length (bp)", y = "GGA proportion of whole repeat", 
         col = "Consensus\nestimation\nmethod", shape = "Consensus\nestimation\nmethod", alpha = "Consensus\nestimation\nmethod")
bp_gga



## TSNE
set.seed(1234)
tsne_data <- tsne(sequence_data_for_pca, k = 2, max_iter = 1000, perplexity = 50) %>%
    as.data.frame %>%
    cbind(sequence_data)

tsne_medians <- tsne_data %>%
    filter(type == "Read") %>%
    group_by(publication_id, allele) %>%
    summarise(V1 = median(V1), 
              V2 = median(V2))

tsne <- ggplot(tsne_data, aes(x = V1, y = V2)) +
    geom_vline(data = tsne_medians, aes(xintercept = V1, col = allele),
               linetype = "dashed", alpha = 0.3) +
    geom_hline(data = tsne_medians, aes(yintercept = V2, col = allele),
               linetype = "dashed", alpha = 0.3) +
    geom_point(data = tsne_data %>% filter(type == "Read"), aes(col = allele), 
               size = 1.5, alpha = 0.7, stroke = 0.75) +
    scale_colour_manual(values = allele_colours) +
    labs(col = "Allele type") +
    new_scale_colour() +
    geom_point(data = tsne_data %>% filter(type == "Consensus") %>% arrange(desc(method)),
               aes(col = method, shape = method), size = 2, stroke = 1) +
    scale_colour_manual(values = method_colours) +
    scale_shape_manual(values = c("Flye" = 4, "Global" = 10, "Local" = 2)) +
    theme_bw() +
    #xlim(0, 700) + ylim(0, 0.3) +
    facet_wrap(~publication_id, scales = "fixed") +
    guides(alpha = "none") +
    ggtitle("(b) t-SNE") +
    labs(x = "t-SNE V1", y = "t-SNE V2", 
         col = "Consensus\nestimation\nmethod", shape = "Consensus\nestimation\nmethod", alpha = "Consensus\nestimation\nmethod")
tsne




## PCA
pca_scaled <- prcomp(sequence_data_for_pca, center = FALSE, scale. = FALSE)
pca_data <- cbind(sequence_data, pca_scaled$x)


pca_medians <- pca_data %>%
    filter(type == "Read") %>%
    group_by(publication_id, allele) %>%
    summarise(PC1 = median(PC1), 
              PC2 = median(PC2))

pca_faceted <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_vline(data = pca_medians, aes(xintercept = PC1, col = allele),
               linetype = "dashed", alpha = 0.3) +
    geom_hline(data = pca_medians, aes(yintercept = PC2, col = allele),
               linetype = "dashed", alpha = 0.3) +
    geom_point(data = pca_data %>% filter(type == "Read"), aes(col = allele), 
               size = 1.5, alpha = 0.7, stroke = 0.75) +
    scale_colour_manual(values = allele_colours) +
    labs(col = "Allele type") +
    new_scale_colour() +
    geom_point(data = pca_data %>% filter(type == "Consensus") %>% arrange(desc(method)),
               aes(col = method, shape = method), size = 2, stroke = 1, alpha = 0.9) +
    scale_colour_manual(values = method_colours) +
    scale_shape_manual(values = c("Flye" = 4, "Global" = 10, "Local" = 2)) +
    theme_bw() +
    ggtitle("(c) PCA") +
    labs(col = "Consensus\nestimation\nmethod", shape = "Consensus\nestimation\nmethod", alpha = "Consensus\nestimation\nmethod") +
    facet_wrap(~publication_id, scales = "fixed") +
    guides(alpha = "none")
pca_faceted



## MAY OR MAY NOT BE INCLUDED IN FIGURE, KEPT HERE FOR OPTIONS
rotation_data <- data.frame(pca_scaled$rotation, variable = rownames(pca_scaled$rotation))
arrow_scale <- 15
pca_combined <- ggplot(rotation_data, aes(x = PC1, y = PC2)) +
    geom_segment(aes(x = 0, xend = PC1*arrow_scale, y = 0, yend = PC2*arrow_scale), arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
    geom_text(aes(x = PC1*arrow_scale, y = PC2*arrow_scale, label = variable), size = 2.5, col = "#444") +
    geom_point(data = pca_data %>% filter(type == "Read"), aes(col = allele), 
               size = 1.5, alpha = 0.7, stroke = 0.75) +
    scale_colour_manual(values = allele_colours) +
    labs(col = "Allele type") +
    new_scale_colour() +
    geom_point(data = pca_data %>% filter(type == "Consensus") %>% arrange(desc(method)),
               aes(col = method, shape = method), size = 1.5, stroke = 1, alpha = 0.9) +
    scale_colour_manual(values = method_colours) +
    scale_shape_manual(values = c("Flye" = 4, "Global" = 10, "Local" = 2)) +
    scale_x_continuous(limits = c(-12,8.5), breaks = seq(-12,8,4)) +
    ggtitle("(d) PCA loadings") +
    labs(col = "Consensus\nestimation\nmethod", shape = "Consensus\nestimation\nmethod", alpha = "Consensus\nestimation\nmethod") +
    theme_bw()
pca_combined




ggarrange(bp_gga, tsne, pca_faceted, pca_combined, nrow = 2, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("output_figures/Figure 2-7 - Dimensionality reduction to compare consensus construction methods.png", dpi = 600, width = 18, height = 14)

#ggarrange(bp_gga, tsne, pca_faceted, nrow = 3, ncol = 1, common.legend = TRUE, legend = "right")
#ggsave("output_figures/Figure 2-7 - Dimensionality reduction to compare consensus construction methods.png", dpi = 600, width = 8, height = 18)





## quantitative analysis


find_distances <- function(dataframe, colname_1, colname_2) {
    distance_data <- data.frame(participant_id = NA, allele = NA, read_id = NA, consensus_method = NA, consensus_x = NA, consensus_y = NA, read_x = NA, read_y = NA, distance = NA)
    
    row <- 1
    for (participant in unique(dataframe$publication_id)) {
        participant_data <- dataframe %>% filter(publication_id == participant)
        for (each_allele in unique(participant_data$allele)) {
            allele_data <- participant_data %>% filter(allele == each_allele)
            locations <- allele_data[allele_data$type == "Read", c(colname_1, colname_2, "sequence_id")]
            for (consensus_method in unique(allele_data[allele_data$type == "Consensus", "method"])) {
                consensus_location <- allele_data[allele_data$method == consensus_method , c(colname_1, colname_2)]
                for (i in 1:nrow(locations)) {
                    consensus_x <- consensus_location[1, 1]
                    consensus_y <- consensus_location[1, 2]
                    read_x <- locations[i, 1]
                    read_y <- locations[i, 2]
                    
                    distance <- sqrt((consensus_x - read_x)^2 + (consensus_y - read_y)^2)
    
                    distance_data[row, "participant_id"] <- participant
                    distance_data[row, "allele"] <- each_allele
                    distance_data[row, "read_id"] <- locations[i, 3]
                    distance_data[row, "consensus_method"] <- consensus_method
                    distance_data[row, "consensus_x"] <- consensus_x
                    distance_data[row, "consensus_y"] <- consensus_y
                    distance_data[row, "read_x"] <- read_x
                    distance_data[row, "read_y"] <- read_y
                    distance_data[row, "distance"] <- distance
                    
                    row <- row + 1
                }
            }
        }
    }
    
    return(distance_data)
}


find_distances_medians <- function(dataframe, medians, colname_1, colname_2) {
    distance_data <- data.frame(participant_id = NA, allele = NA, consensus_method = NA, consensus_x = NA, consensus_y = NA, median_x = NA, median_y = NA, distance = NA)
    row <- 1
    for (participant in unique(dataframe$publication_id)) {
        participant_data <- dataframe %>% filter(publication_id == participant)
        for (each_allele in unique(participant_data$allele)) {
            allele_data <- participant_data %>% filter(allele == each_allele)
            location <- medians[medians$publication_id == participant & medians$allele == each_allele, c(colname_1, colname_2)]

            for (consensus_method in unique(allele_data[allele_data$type == "Consensus", "method"])) {
                consensus_location <- allele_data[allele_data$method == consensus_method , c(colname_1, colname_2)]

                consensus_x <- consensus_location[1, 1]
                consensus_y <- consensus_location[1, 2]
                median_x <- location[1, 1]
                median_y <- location[1, 2]
                
                distance <- sqrt((consensus_x - median_x)^2 + (consensus_y - median_y)^2)

                distance_data[row, "participant_id"] <- participant
                distance_data[row, "allele"] <- each_allele
                distance_data[row, "consensus_method"] <- consensus_method
                distance_data[row, "consensus_x"] <- consensus_x
                distance_data[row, "consensus_y"] <- consensus_y
                distance_data[row, "median_x"] <- median_x
                distance_data[row, "median_y"] <- median_y
                distance_data[row, "distance"] <- distance
                
                row <- row + 1
            }
        }
    }
    
    return(distance_data)
}


## Reads
## Not meaningful to do length/GGA unscaled - then it is essentially only looking at length
scaled_data <- cbind(sequence_data[ , c(1:7)], sequence_data_for_pca)
length_gga_distances <- find_distances(scaled_data, "length_bp", "whole_seq_proportion_GGA")
length_gga_results <- length_gga_distances %>% 
    group_by(consensus_method) %>%
    summarise(mean = mean(distance), median = median(distance), count = n())

tsne_distances <- find_distances(tsne_data, "V1", "V2")
tsne_results <- tsne_distances %>% 
    group_by(consensus_method) %>%
    summarise(mean = mean(distance), median = median(distance), count = n())

pca_distances <- find_distances(pca_data, "PC1", "PC2")
pca_results <- pca_distances %>% 
    group_by(consensus_method) %>%
    summarise(mean = mean(distance), median = median(distance), count = n())


## Medians
scaled_medians <- scaled_data %>%
    filter(type == "Read") %>%
    group_by(publication_id, allele) %>%
    summarise(length_bp = median(length_bp), 
              whole_seq_proportion_AGC = median(whole_seq_proportion_AGC),
              whole_seq_proportion_GGC = median(whole_seq_proportion_GGC),
              whole_seq_proportion_GGA = median(whole_seq_proportion_GGA),
              whole_seq_proportion_GGT = median(whole_seq_proportion_GGT))

length_gga_distances_medians <- find_distances_medians(scaled_data, scaled_medians, "length_bp", "whole_seq_proportion_GGA")
length_gga_results_medians <- length_gga_distances_medians %>% 
    group_by(consensus_method) %>%
    summarise(mean = mean(distance), median = median(distance), count = n())

tsne_distances_medians <- find_distances_medians(tsne_data, tsne_medians, "V1", "V2")
tsne_results_medians <- tsne_distances_medians %>% 
    group_by(consensus_method) %>%
    summarise(mean = mean(distance), median = median(distance), count = n())

pca_distances_medians <- find_distances_medians(pca_data, pca_medians, "PC1", "PC2")
pca_results_medians <- pca_distances_medians %>% 
    group_by(consensus_method) %>%
    summarise(mean = mean(distance), median = median(distance), count = n())



reformat_results_dataframe <- function(results_dataframe, observations_level, dimensionality_reduction) {
    colnames(results_dataframe) <- c("consensus_method", "mean_distance", "median_distance", "number_of_observations")
    results_dataframe$observations_level <- observations_level
    results_dataframe$dimensionality_reduction <- dimensionality_reduction
    return(results_dataframe[ , c(5, 6, 1, 4, 2, 3)])
}

output_results_table <- rbind(reformat_results_dataframe(length_gga_results, "Reads", "Length/GGA"),
                              reformat_results_dataframe(tsne_results, "Reads", "t-SNE"),
                              reformat_results_dataframe(pca_results, "Reads", "PCA"),
                              reformat_results_dataframe(length_gga_results_medians, "Allele medians", "Length/GGA"),
                              reformat_results_dataframe(tsne_results_medians, "Allele medians", "t-SNE"),
                              reformat_results_dataframe(pca_results_medians, "Allele medians", "PCA"))
output_results_table[, c("mean_distance", "median_distance")] <- round(output_results_table[, c("mean_distance", "median_distance")], 3)
write.csv(output_results_table, "output_tables/Appendix Table 2-3 - Quantitative consensus method comparison.csv", row.names = F)
