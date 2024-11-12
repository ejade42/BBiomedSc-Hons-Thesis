## Runs in output directory, compiles all the individual read information csvs.
## Takes all participant IDs as arguments

arguments <- as.character(commandArgs(trailingOnly = TRUE))
all_participant_ids <- arguments[1:length(arguments)-1]
output_filename     <- arguments[length(arguments)]

all_participant_data <- NULL
for (participant_id in all_participant_ids) {
    new_data <- read.csv(paste0(participant_id, "/", participant_id, "_methylation_information.csv"))
    all_participant_data <- rbind(all_participant_data, new_data)
}

write.csv(all_participant_data, output_filename, row.names = F)
