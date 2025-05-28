#!/bin/bash
        #------------------------------------------------
        # 6.1 prepare BED files based on marker positions
        #------------------------------------------------

        # expected output: contig start end

for FAM_INF_LIM in 3 4
        do
PATH_FILTERED_GZ=~/LM_new_genome/results/4.Filtering2/famInfLim${FAM_INF_LIM}
zcat ${PATH_FILTERED_GZ}/data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25.gz | tail -n+8 | cut -f1,2 > ${PATH_FILTERED_GZ}/snps_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25.txt
awk -F'\t' '{gsub("\\*", "", $2); print $1, $2 - 1, $2 }' OFS='\t' ${PATH_FILTERED_GZ}/snps_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25.txt > ${PATH_FILTERED_GZ}/snps_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25.bed

        #----------------------------------------
        # 6.2 Remove repeats + flanking sequences
        #----------------------------------------
REPEATS_BED=~/LM_new_genome/data/repeats/before_anchoring/rhizoglyphusRobini.filteredRepeats.bed
GENOME_FILE=~/assembly/23.EarlGrey_TE_annotation/host_data/IW24_minLen5kb_minQ10_contaminantsRemoved_polished_Racon3_medaka1_Pilon3_pugedDups.fasta.fai
NAME_PREFIX=data.Filtering_dataTol0.05
MISS_LIM=0.25
        for FLANK in 100 200
                do
# flank repeats by 100 and 200 bp to be more conservative during repeat removal:

bedtools slop -i ${REPEATS_BED} -g ${GENOME_FILE} -b ${FLANK} > ~/LM_new_genome/data/repeats/rhizoglyphusRobini.filteredRepeats_AnnotationExtended${FLANK}bp.bed


PATH_FILTERED_GZ=~/LM_new_genome/results/5.RemovingRepeats/famInfLim${FAM_INF_LIM}

cd ${PATH_FILTERED_GZ}

REPEATS=~/LM_new_genome/data/repeats/before_anchoring/rhizoglyphusRobini.filteredRepeats_AnnotationExtended${FLANK}bp.bed
MARKERS_BED=~/LM_new_genome/results/4.Filtering2/famInfLim${FAM_INF_LIM}/snps_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25.bed

bedtools intersect \
    -a ${MARKERS_BED} \
    -b ${REPEATS} \
    -v > ${PATH_FILTERED_GZ}/snps_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions${FLANK}bp.bed

# count number of markers in contigs
cut -f1 ${PATH_FILTERED_GZ}/snps_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions${FLANK}bp.bed | sort | uniq -c > ${PATH_FILTERED_GZ}/marker_count_per_contigs_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions${FLANK}bp.bed
# count number of markers
wc -l ${PATH_FILTERED_GZ}/marker_count_per_contigs_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions${FLANK}bp.bed > number_of_contigs_data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_MAF0.15_missingLimit0.25_RepeatsRemovedFlankingRegions${FLANK}bp.bed

    done
done
