## Creates exact read length histogram (Figure 3-8b/c) and the alpha scaling diagram (Appendix Fig 3-1).

library(tidyverse)
library(ggtext)
library(ggpattern)
library(ggpubr)
source("scripts/common_functions.R")

all_participant_data <- read.csv("intermediate_files/read_information_all_participants.csv")
all_participant_data$participant_id <- as.character(all_participant_data$participant_id)

all_participant_data$publication_id <- assign_publication_ids(all_participant_data)
all_participant_data$batch          <- assign_participant_batches(all_participant_data)$batch
all_participant_data$data_source    <- assign_participant_batches(all_participant_data)$data_source

all_participant_data <- all_participant_data %>% 
    filter(participant_id != "ishiura_2019_f9193_ii_5_raw")


## Merge in phenotypes
phenotype_information <- read.csv("input/phenotype_information.csv")[, c("participant_id", "condition", "phenotype")]
phenotype_information <- phenotype_information[!duplicated(phenotype_information), ]
phenotype_information$publication_id <- phenotype_information$participant_id
phenotype_information$participant_id <- NULL

merged_data <- merge(all_participant_data, phenotype_information, by = "publication_id", all.x = T)
merged_data[merged_data$publication_id %in% c(paste0("MND", 1:7)), c("condition", "phenotype")] <- "MND"
merged_data[merged_data$phenotype == "Asymptomatic", "condition"] <- "Asymptomatic"


## Relevel factors
merged_data$publication_id <- factor(merged_data$publication_id,
                                     c("F1-II-3", "F1-III-1", "F1-III-3",
                                       "F2-I-1", "F2-II-1", "F2-II-2", "F2-II-3",
                                       "F3-II-1", "F4-II-1", paste0("MND", 1:7),
                                       "Ishiura 2019 F9193-II-5", "Podar 2023",
                                       "Tian 2019 F1-IV-7", "Tian 2019 F1-IV-15", "Tian 2019 F2-II-3",
                                       "Tian 2019 F4-II-2", "Tian 2019 F5-II-1", "Tian 2019 F5-II-4", "Tian 2019 F9-II-6"))
merged_data$condition <- factor(merged_data$condition, c("NIID", "OPDM", "MND", "Asymptomatic"))



## Calculating flanking lengths for fuzzy allele boundaries
all_participant_data$flanking_lengths <- all_participant_data$flanked_length - all_participant_data$exact_length
print("Flanking sequence length in each trimmed repeat:", quote = F)
summary(all_participant_data$flanking_lengths)
table(all_participant_data$flanking_lengths)
print(paste("Mean:", mean(all_participant_data$flanking_lengths)), quote = F)
print(paste("SD:", sd(all_participant_data$flanking_lengths)), quote = F)



## Main figures
full <- ggplot(merged_data, aes(x = exact_length, fill = condition, pattern = data_source)) +
    scale_x_continuous(breaks = seq(0,2200,200), limits = c(0,2200)) +
    scale_pattern_manual(values = c("Literature" = "stripe", "Novel (ours)" = "none")) +
    scale_fill_manual(values = c("OPDM" = "#DC1E1E", "NIID" = "#008C00", "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    geom_histogram_pattern(col = alpha("#333", 1), linewidth = 0.3, binwidth = 25, boundary = 0, closed = "left",
                           pattern_angle = 45, pattern_colour = NA, pattern_fill = "black", pattern_spacing = 0.01, pattern_density = 0.2) +
    geom_vline(xintercept = 173, col = "red", linetype = "dashed") +
    geom_vline(xintercept = 773, col = "red", linetype = "dashed") +
    ylim(0, 150) +
    theme_bw() +
    labs(x = "Exact <i>NOTCH2NLC</i> repeat length (bp)", y = "Read count", pattern = "Data source", fill = "Phenotype") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),  order = 1,
           fill    = guide_legend(override.aes = list(pattern = "none"), order = 2)) +
    theme(title = element_markdown()) 



faceted <- ggplot(merged_data, aes(x = exact_length, fill = condition)) +
    scale_x_continuous(breaks = seq(0,2200,400), limits = c(0,2200)) +
    scale_y_continuous(breaks = seq(0,50,25), limits = c(0,50)) +
    scale_fill_manual(values = c("OPDM" = "#DC1E1E", "NIID" = "#008C00", "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    scale_pattern_manual(values = c("Literature" = "stripe", "Novel" = "none")) +
    geom_histogram_pattern(aes(pattern = data_source), col = alpha("#333", 1), linewidth = 0.3, binwidth = 25, boundary = 0, closed = "left",
                           pattern_angle = 45, pattern_colour = NA, pattern_fill = "black", pattern_spacing = 0.01, pattern_density = 0.2) +
    geom_histogram(col = alpha("#333", 1), linewidth = 0.3, binwidth = 25, boundary = 0, closed = "left") +
    geom_vline(xintercept = 173, col = "red", linetype = "dashed") +
    geom_vline(xintercept = 773, col = "red", linetype = "dashed") +
    theme_bw() +
    labs(x = "Exact <i>NOTCH2NLC</i> repeat length (bp)", y = "Read count", pattern = "Data source", fill = "Phenotype") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),  order = 1,
           fill    = guide_legend(override.aes = list(pattern = "none"), order = 2),
           alpha   = "none") +
    facet_wrap(~publication_id) +
    theme(title = element_markdown()) 

height_diff <- 2.6
left <- ggarrange(NULL, full, nrow = 2, ncol = 1, heights = c(height_diff, 6.25), common.legend = TRUE, legend = "none")
combined <- ggarrange(left, faceted, nrow = 1, ncol = 2, widths = c(7.5, 10.5), common.legend = T, legend = "right")

ggsave("output_figures/Figure 3-7bc - Exact repeat sizes histogram.png", plot = combined, dpi = 1200, width = 18, height = 6.25+height_diff)



### Test bin closing
bins_data <- data.frame(x = c(0, 0.1, 4, 5, 6, 10, 11))
ggplot(bins_data, aes(x = x)) +
    geom_histogram(binwidth = 5, boundary = 0, closed = "left") +
    ylim(0, 5)

