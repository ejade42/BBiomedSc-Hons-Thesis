input_filenames=("intermediate_files/notch2nlc_transcript_variant_1_reading_frames_exons.txt"
				 "input/n2c_uorf_3tag_plasmid_sequence.fastq"
				 "intermediate_files/n2c_uorf_3tag_plasmid_reading_frames.txt"
				 "input/phenotype_information.csv"
				 "intermediate_files/all_local_consensuses.fasta")

output_filenames=("1 - NOTCH2NLC transcript variant 1 cDNA with exon annotations and translations.txt"
				  "2 - pcDNA3.1 N2C uORF-3Tag plasmid sequence.fastq"
				  "3 - N2C uORF-3Tag plasmid with annotations and translations.txt"
				  "6 - Phenotype information.csv"
				  "8 - All local consensuses.fasta")

output_location="/Volumes/ScotterLab/Evelyn Jade/BBiomedSci (Hons) Assignments/240821 thesis with git/appendix/Supplementary files/"


for i in {1..${#input_filenames[*]}}
do
	old_file=$input_filenames[$i]
	new_file=$output_filenames[$i]
	echo $old_file" > "$new_file; echo
	
	cp $old_file $output_location$new_file
done
