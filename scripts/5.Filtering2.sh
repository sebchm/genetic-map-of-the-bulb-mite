#!/bin/bash
for MISS_LIM in 0.5 0.25 0.75
    do
    for FAM_INF_LIM in 3 4
        do
                        echo LM3 ParentCall2 with missing limit=${MISS_LIM} and famInformativeLimit=${FAM_INF_LIM} started at `date`
zcat ${RESULTS_DIR_PC2}/data.ParentCall_q20_Q20.gz | java -cp ~/software/Lep_Map3_0.5/bin Filtering2 data=- dataTolerance=0.05 MAFLimit=0.05 missingLimit=${MISS_LIM} removeNonInformative=1 familyInformativeLimit=${FAM_INF_LIM} | gzip > ${RESULTS_DIR_F2}/data.Filtering_dataTol0.05_familyInfLimit${FAM_INF_LIM}_missingLimit${MISS_LIM}.gz
    done
done

# use R script proportion_of_genome_size_covered_by_contigsWithMarkers.R
