## Runs in dorado/$participant_id for each participant
## Analyses methylation information


## ARGUMENTS
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3) {stop("Please provide the participant id, input filename, and output filename when calling the script")}
participant_id  <- as.character(arguments[1])
filename_base   <- as.character(arguments[2])
output_filename <- as.character(arguments[3])


## FUNCTIONS
string_to_vector  <- function(string) {as.numeric(unlist(strsplit(string, split = ",")))}
vector_to_string  <- function(vector) {paste(vector, collapse = ",")}


## MAIN PROCESSING
methylation_data <- read.csv(paste0(filename_base, "_information.csv"))
colnames(methylation_data) <- c("read_id", "sequence", "quality", "modification_types", "C_hm_locations", "C_hm_probabilities", "C_m_locations", "C_m_probabilities")
methylation_data$participant_id <- participant_id

methylation_data$length <- nchar(methylation_data$sequence)

forward_reads <- readLines(paste0(filename_base, "_forward_ids.txt"))
reverse_reads <- readLines(paste0(filename_base, "_reverse_ids.txt"))
for (i in 1:nrow(methylation_data)) {
    if (methylation_data[i, "read_id"] %in% forward_reads) {methylation_data[i, "direction"] <- "forward"}
    else if (methylation_data[i, "read_id"] %in% reverse_reads) {methylation_data[i, "direction"] <- "reverse"}
    else {methylation_data[i, "direction"] <- NA}
}

write.csv(methylation_data[ , c("participant_id", "read_id", "direction", "length", "sequence", "quality", "modification_types", "C_hm_locations", "C_hm_probabilities", "C_m_locations", "C_m_probabilities")], output_filename, row.names = F)
