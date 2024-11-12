## Takes read information csvs (output of pipeline, plus "repeat_seq" column which is the same as "forward_seq"
## but manually trimmed exactly to the repeat region), outputs a combined csv for downstream analysis.

source("scripts/common_functions.R")
all_data <- rbind(read.csv("input/batch_1_read_information_all_participants_with_trimming.csv"),
                  read.csv("input/batch_2_read_information_all_participants_with_trimming.csv"),
                  read.csv("input/batch_3_4_read_information_all_participants_with_trimming.csv"))

all_data$check_backwards <- lapply(all_data$repeat_seq, reverse_sequence)

output_data <- all_data[ , c("participant_id", "read_id", "allele", "forward_seq", "length", "repeat_seq")]
colnames(output_data)[4:6] <- c("flanked_repeat", "flanked_length", "exact_repeat")
output_data$exact_length <- nchar(output_data$exact_repeat)

## Make allele all caps
output_data[output_data$allele == "short", "allele"] <- "Short"
output_data[output_data$allele == "long", "allele"]  <- "Long"
output_data[output_data$allele == "hyper", "allele"] <- "Hyper"

write.csv(output_data, "intermediate_files/read_information_all_participants.csv", row.names = F)
