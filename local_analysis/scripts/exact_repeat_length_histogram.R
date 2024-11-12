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

create_rectangle_parameters <- function(original_position, diff_vector, alpha_transform = 0.5) {
    diff_table <- as.data.frame(table(diff_vector))
    all_layers <- vector(mode = "list", length = nrow(diff_table))
    alpha_interval <- 1 / max(diff_table[, 2])
    parameters <- data.frame(matrix(NA, nrow(diff_table), 4))
    
    for (i in 1:nrow(diff_table)) {
        position_shift <- as.numeric(as.character((diff_table[i, 1])))
        xmin_position <- original_position - position_shift
        xmax_position <- original_position - position_shift + 1
        alpha_value <- (alpha_interval * diff_table[i, 2]) ^ (alpha_transform)

        parameters[i, 1] <- position_shift
        parameters[i, 2] <- xmin_position
        parameters[i, 3] <- xmax_position
        parameters[i, 4] <- alpha_value
    }
    
    colnames(parameters) <- c("shift", "xmin", "xmax", "alpha")
    return(parameters)
}


all_participant_data$flanking_lengths <- all_participant_data$flanked_length - all_participant_data$exact_length
print("Flanking sequence length in each trimmed repeat:", quote = F)
summary(all_participant_data$flanking_lengths)
table(all_participant_data$flanking_lengths)
print(paste("Mean:", mean(all_participant_data$flanking_lengths)), quote = F)
print(paste("SD:", sd(all_participant_data$flanking_lengths)), quote = F)


parameters_200 <- create_rectangle_parameters(200, all_participant_data$flanking_lengths, 1/3)
parameters_800 <- create_rectangle_parameters(800, all_participant_data$flanking_lengths, 1/3)



## Make brief figure explaining alpha scaling for appendix
scaling_data <- data.frame(x_values = c(seq(0, 424, 0.1), seq(0, 30, 0.1)))
scaling_data$y_values_cube <- (scaling_data$x_values / 424) ^ (1/3)
scaling_data$y_values_line <- (scaling_data$x_values / 424)
scaling_data$range <- c(rep("Whole range", length(seq(0, 424, 0.1))),
                        rep("Zoomed", length(seq(0, 30, 0.1))))

ggplot(scaling_data, aes(x = x_values)) +
    scale_colour_manual(values = c("Linear" = "#006EDC", "Cube root" = "#DC1E1E")) +
    geom_line(aes(y = y_values_cube, alpha = y_values_cube, col = "Cube root"), linewidth = 2) +
    geom_line(aes(y = y_values_line, alpha = y_values_line, col = "Linear"), linewidth = 2) +
    scale_alpha_continuous(range = c(0,1)) +
    guides(alpha = "none") +
    labs(x = "Occurences", y = "Alpha") +
    facet_wrap(~range, scales = "free_x") +
    theme_bw() + theme(legend.title = element_blank())

ggsave("output_figures/Appendix Figure 3-1 - Alpha scaling.png", dpi = 600, width = 8, height = 4)






## Main figures
full <- ggplot(merged_data, aes(x = exact_length, fill = condition, pattern = data_source)) +
    scale_x_continuous(breaks = seq(0,2200,200), limits = c(0,2200)) +
    scale_alpha_continuous(range = c(0,1)) +
    scale_pattern_manual(values = c("Literature" = "stripe", "Novel (ours)" = "none")) +
    scale_fill_manual(values = c("OPDM" = "#DC1E1E", "NIID" = "#008C00", "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    geom_histogram_pattern(col = alpha("#333", 1), linewidth = 0.3, binwidth = 25, boundary = 0, closed = "left",
                           pattern_angle = 45, pattern_colour = NA, pattern_fill = "black", pattern_spacing = 0.01, pattern_density = 0.2) +
    geom_rect(data = parameters_200, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = alpha), inherit.aes = F, fill = "red") +
    geom_rect(data = parameters_800, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = alpha), inherit.aes = F, fill = "red") +
    theme_bw() +
    labs(x = "Exact <i>NOTCH2NLC</i> repeat length (bp)", y = "Read count", pattern = "Data source", fill = "Phenotype") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),  order = 1,
           fill    = guide_legend(override.aes = list(pattern = "none"), order = 2),
           alpha   = "none") +
    theme(title = element_markdown()) 


faceted <- ggplot(merged_data, aes(x = exact_length, fill = condition)) +
    scale_x_continuous(breaks = seq(0,2200,500), limits = c(0,2200)) +
    scale_y_continuous(breaks = seq(0,50,25), limits = c(0,50)) +
    scale_alpha_continuous(range = c(0,1)) +
    scale_fill_manual(values = c("OPDM" = "#DC1E1E", "NIID" = "#008C00", "MND" = "#F07800", "Asymptomatic" = "#006EDC")) +
    geom_histogram(col = alpha("#333", 1), linewidth = 0.3, binwidth = 25, boundary = 0, closed = "left") +
    geom_rect(data = parameters_200, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = alpha), inherit.aes = F, fill = "red") +
    geom_rect(data = parameters_800, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = alpha), inherit.aes = F, fill = "red") +
    theme_bw() +
    labs(x = "Exact <i>NOTCH2NLC</i> repeat length (bp)", y = "Read count", pattern = "Data source", fill = "Phenotype") +
    guides(pattern = guide_legend(override.aes = list(fill = "white")),  order = 1,
           fill    = guide_legend(override.aes = list(pattern = "none"), order = 2),
           alpha   = "none") +
    facet_wrap(~publication_id) +
    theme(title = element_markdown()) 

ggarrange(full, faceted, nrow = 1, ncol = 2, heights = 0.37, widths = c(0.5, 0.5), common.legend = T, legend = "right")
ggsave("output_figures/Figure 3-8bc - Exact repeat sizes histogram.png", dpi = 1200, width = 18, height = 6.25)





### Test bin closing
bins_data <- data.frame(x = c(0, 0.1, 4, 5, 6, 10, 11))
ggplot(bins_data, aes(x = x)) +
    geom_histogram(binwidth = 5, boundary = 0, closed = "left") +
    ylim(0, 5)

