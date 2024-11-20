#!/bin/bash -e
echo \($(date +"%Y-%m-%d %T")\) "start"

echo \($(date +"%Y-%m-%d %T")\) "exporting consensus files from python scripts"
python3 scripts/consensuses_from_pipeline.py
python3 scripts/consensuses_from_literature.py -s

echo \($(date +"%Y-%m-%d %T")\) "creating read information summary file"
Rscript scripts/process_read_information_from_pipeline.R

echo \($(date +"%Y-%m-%d %T")\) "creating N2C sequence/translation reference files for Figures 1-2 and 2-1"
Rscript scripts/notch2nlc_and_plasmid_summaries.R

echo; echo \($(date +"%Y-%m-%d %T")\) "creating Figure 2-6: Flanked read sizes histogram"
Rscript scripts/trimmed_read_length_histogram.R

echo; echo \($(date +"%Y-%m-%d %T")\) "creating all read vs consensus outputs (Appendix Tables 2-3 and 3-1; Figure 2-7)"
Rscript scripts/reads_vs_consensus.R

echo; echo \($(date +"%Y-%m-%d %T")\) "creating exact repeat illustration and histograms (Figure 3-8)"
Rscript scripts/exact_repeat_illustration.R
Rscript scripts/exact_repeat_length_histogram.R

echo; echo \($(date +"%Y-%m-%d %T")\) "creating selected consensus sequence visualisations (Figure 3-9)"
Rscript scripts/consensus_illustrations.R

echo; echo \($(date +"%Y-%m-%d %T")\) "creating methylation figures (Figures 3-10 and 3-11)"
Rscript scripts/methylation.R

echo; echo \($(date +"%Y-%m-%d %T")\) "creating meta-analysis figures (Figure 3-12)"
Rscript scripts/meta_analysis.R


echo \($(date +"%Y-%m-%d %T")\) "done"
