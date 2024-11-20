## Creates Figure 2-6: Histogram of trimmed read lengths over NOTCH2NLC repeat region

library(tidyverse)
library(ggtext)
library(ggpubr)
source("scripts/common_functions.R")

all_participant_data <- read.csv("intermediate_files/read_information_all_participants.csv")

all_participant_data$batch <- assign_participant_batches(all_participant_data)$batch
all_participant_data$data_source <- assign_participant_batches(all_participant_data)$data_source
all_participant_data$publication_id <- assign_publication_ids(all_participant_data)
all_participant_data$publication_id <- factor(all_participant_data$publication_id,
                                              c("F1-II-3", "F1-III-1", "F1-III-3",
                                                "F2-I-1", "F2-II-1", "F2-II-2", "F2-II-3",
                                                "F3-II-1", "F4-II-1", paste0("MND", 1:7),
                                                "Ishiura 2019 F9193-II-5", "Podar 2023",
                                                "Tian 2019 F1-IV-7", "Tian 2019 F1-IV-15", "Tian 2019 F2-II-3",
                                                "Tian 2019 F4-II-2", "Tian 2019 F5-II-1", "Tian 2019 F5-II-4", "Tian 2019 F9-II-6"))


main_histogram <- ggplot(all_participant_data, aes(x = flanked_length, fill = batch, pattern = data_source)) +
    scale_x_continuous(breaks = seq(0,2200,200), limits = c(0,2200)) +
    scale_fill_manual(values = c("Batch 1" = "#DC1E1E", "Batch 2" = "#F07800", "Batch 3/4" = "#008C00", "Ishiura" = "#006EDC", "Podar" = "#781EF0", "Tian" = "#DC00DC")) +
    geom_histogram(col = alpha("#333", 1), linewidth = 0.3, binwidth = 25, boundary = 0, closed = "left") +
    geom_vline(xintercept = 200, col = "red", linetype = "dashed") +
    geom_vline(xintercept = 800, col = "red", linetype = "dashed") +
    theme_bw() +
    labs(x = "Trimmed read length over <i>NOTCH2NLC</i> repeat region (bp)", y = "Read count", pattern = "Data source", fill = "Participant group") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),  order = 1,
           fill    = guide_legend(override.aes = list(pattern = "none"), order = 2)) +
    theme(title = element_markdown())

subset_histogram <- ggplot(all_participant_data %>% filter(participant_id %in% c("313_merged", "327_merged")), aes(x = flanked_length, fill = batch, pattern = data_source)) +
    scale_x_continuous(breaks = seq(0,2200,200), limits = c(0,2200)) +
    scale_y_continuous(breaks = seq(0,15,5)) +
    scale_fill_manual(values = c("Batch 1" = "#DC1E1E", "Batch 2" = "#F07800", "Batch 3/4" = "#008C00", "Ishiura" = "#006EDC", "Podar" = "#781EF0", "Tian" = "#DC00DC")) +
    geom_histogram(col = alpha("#333", 1), linewidth = 0.3, binwidth = 25, boundary = 0, closed = "left") +
    geom_vline(xintercept = 200, col = "red", linetype = "dashed") +
    geom_vline(xintercept = 800, col = "red", linetype = "dashed") +
    theme_bw() +
    labs(x = "Trimmed read length over <i>NOTCH2NLC</i> repeat region (bp)", y = "Read count", pattern = "Data source", fill = "Participant group") +
    guides(fill = guide_legend(override.aes = list(pattern = "none"), order = 2)) +
    facet_wrap(~publication_id, nrow = 2) +
    theme(title = element_markdown())

ggarrange(main_histogram, NULL, subset_histogram, ncol = 1, common.legend = TRUE, legend = "right", heights = c(5, 0.1, 3.5), labels = c("(a)", "", "(b)"))
ggsave("output_figures/Figure 2-6 - Trimmed read length histogram.png", dpi = 600, width = 8, height = 9)
