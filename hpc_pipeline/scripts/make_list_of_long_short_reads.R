## Runs in the individual participants' output directories
## Reads the csv of read information, and creates a file of long and short allele read IDs


arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1) {stop("Please provide the participant ID when calling the script")}

participant_id <- as.character(arguments[1])

participant_data <- read.csv(paste0(participant_id, "_read_information.csv"))

for (i in 1:nrow(participant_data)) {
    participant_data[i, "read_id"] <- substring(participant_data[i, "read_id"], 2)
}


short_reads <- as.data.frame(participant_data[participant_data$allele == "short", "read_id"])
long_reads  <- as.data.frame(participant_data[participant_data$allele == "long", "read_id"])
hyper_reads <- as.data.frame(participant_data[participant_data$allele == "hyper", "read_id"])


write.table(short_reads, paste0(participant_id, "_short_reads.txt"), row.names = F, col.names = F, quote = F)
write.table(long_reads, paste0(participant_id, "_long_reads.txt"), row.names = F, col.names = F, quote = F)
write.table(hyper_reads, paste0(participant_id, "_hyper_reads.txt"), row.names = F, col.names = F, quote = F)
