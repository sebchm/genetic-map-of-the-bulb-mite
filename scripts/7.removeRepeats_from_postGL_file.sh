#!/bin/bash
#### get positions of markers from FamLim=3 and MissingLim=0.25:

cd ~/LM_new_genome/results/5.RemovingRepeats/famInfLim3
cut -f1,3 snps_data.Filtering_dataTol0.05_familyInfLimit3_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions200bp.bed > marker_positions_dataTol0.05_familyInfLimit3_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions200bp.txt

#### select only markers from contigs which were missing in FamLim4:
        # these markers have been filtered using Rmd script: ~/LM_new_genome/results/4.Filtering2/filtering_markers_dataTol0.05_minMAF0.15_flank200bp.Rmd
grep -f contigs_to_complement_from_FamLim3_dataTol0.05_minMAF0.05_flank200bp.txt marker_positions_dataTol0.05_familyInfLimit3_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions200bp.txt > markers_from_FamLim3_to_complement_FamLim4_data.Filtering_dataTol0.05_familyInfLimit3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp.txt

#### merge markers from FamLim4 and complemented FamLim3:
        MARKERS_FAMLIM3_TO_COMPLEMENT=~/LM_new_genome/results/5.RemovingRepeats/famInfLim3/markers_from_FamLim3_to_complement_FamLim4_data.Filtering_dataTol0.05_familyInfLimit3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp.txt

cut -f 1,3 ~/LM_new_genome/results/5.RemovingRepeats/famInfLim4/snps_data.Filtering_dataTol0.05_familyInfLimit4_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions100bp.bed > ~/LM_new_genome/results/5.RemovingRepeats/famInfLim4/marker_positions_data.Filtering_dataTol0.05_familyInfLimit4_missingLimit0.25_RepeatsRemovedFlankingRegions200bp.txt
        MARKERS_FAMLIM4=~/LM_new_genome/results/5.RemovingRepeats/famInfLim4/marker_positions_data.Filtering_dataTol0.05_familyInfLimit4_missingLimit0.25_RepeatsRemovedFlankingRegions200bp.txt

        cat $MARKERS_FAMLIM3_TO_COMPLEMENT $MARKERS_FAMLIM4 | sort | uniq > ~/LM_new_genome/results/5.RemovingRepeats/famInfLim4_complementedFamLim3/marker_positions_data.Filtering_dataTol0.05_familyInfLimit4ComplFamInfLimit3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15.txt


cd ~/LM_new_genome/results/5.RemovingRepeats/famInfLim4_complementedFamLim3/
MARKERS_FAMLIM4_COMPLEMENTED_FAMLIM3=marker_positions_data.Filtering_dataTol0.05_familyInfLimit4ComplFamInfLimit3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15.txt

# in the post file, keep only markers which were retained after removing repeats and complementing markers from FamLim3:
    awk 'NR==FNR {
    # read file with filtered markers and store the markers in an array
    markers[$1"\t"$2] = 1
    next
}
{
    # for each line in the post file, copy the second field
    pos = $2 
    # remove any trailing asterisk from the temporary variable (asterisks denot sex markers)
    sub(/\*$/, "", pos) 
    # check if the cleaned-up key (CHR and cleaned POS) exists in the markers array
    if (($1"\t"pos) in markers) { 
        # if yes, print the original line (with asterisk included)
        print 
    }}' $MARKERS_FAMLIM4_COMPLEMENTED_FAMLIM3 <(zcat ~/LM_new_genome/results/4.Filtering2/famInfLim3/data.Filtering_dataTol0.05_familyInfLimit3_missingLimit0.25.gz) \
> data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15.txt

# prepare header
zcat ~/LM_new_genome/results/4.Filtering2/famInfLim3/data.Filtering_dataTol0.05_familyInfLimit3_missingLimit0.25.gz | head -n 7 > header.txt

# add header
cat header.txt data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15.txt | gzip > data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15.gz
