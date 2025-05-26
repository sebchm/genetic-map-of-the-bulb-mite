#!/bin/bash
### Prepare intervals
                # copy the intervals to the LA directory
                MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21

                mkdir -p ~/LM_new_genome/results/8.LA/${MAP_NAME}/intervals
                cp ~/LM_new_genome/results/7.OM2/${MAP_NAME}/infMask123/intervals/*int ~/LM_new_genome/results/8.LA/${MAP_NAME}/intervals
                    # sex chromosome- LG7
                cp ~/LM_new_genome/results/7.OM2/${MAP_NAME}/infMask2/intervals/*int ~/LM_new_genome/results/8.LA/${MAP_NAME}/intervals

                # change dir to interval dir:
                cd ~/LM_new_genome/results/8.LA/${MAP_NAME}/intervals

                # "name" of the SNP file (2 columns: contig contig_physical_position)
                SNP_FILE=~/LM_new_genome/results/5.RemovingRepeats/famInfLim4_complementedFamLim3/snps_data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15.txt

                INTERVAL_NAME_PREFIX=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_LOD21_infMask123

                # "name" the markers: add contig and contig position
                for CHR in 1 2 3 4 5 6
                        do
                # assign the contigs and physical position to each marker
                        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
                        ${INTERVAL_NAME_PREFIX}_LG${CHR}_int > ${INTERVAL_NAME_PREFIX}_LG${CHR}_int_named.txt
                                done

                # special care for LG8:
                CHR=8
                        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
                        ${INTERVAL_NAME_PREFIX}_LG${CHR}_60iterations_rmEdgeMarkers_int > ${INTERVAL_NAME_PREFIX}_LG${CHR}_int_named.txt

                #add CHR column
                awk '$2 = $2 FS "1"' ${INTERVAL_NAME_PREFIX}_LG1_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG1_int_named_POC_input.txt
                awk '$2 = $2 FS "2"' ${INTERVAL_NAME_PREFIX}_LG2_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG2_int_named_POC_input.txt
                awk '$2 = $2 FS "3"' ${INTERVAL_NAME_PREFIX}_LG3_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG3_int_named_POC_input.txt
                awk '$2 = $2 FS "4"' ${INTERVAL_NAME_PREFIX}_LG4_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG4_int_named_POC_input.txt
                awk '$2 = $2 FS "5"' ${INTERVAL_NAME_PREFIX}_LG5_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG5_int_named_POC_input.txt
                awk '$2 = $2 FS "6"' ${INTERVAL_NAME_PREFIX}_LG6_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG6_int_named_POC_input.txt
                awk '$2 = $2 FS "8"' ${INTERVAL_NAME_PREFIX}_LG8_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG8_int_named_POC_input.txt

                # special care for the sex chromosome:
                    # assign contig and contig position
                INTERVAL_NAME_PREFIX=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_LOD21_infMask2
                CHR=7
                         awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
                        ${INTERVAL_NAME_PREFIX}_LG${CHR}_int > ${INTERVAL_NAME_PREFIX}_LG${CHR}_int_named.txt
                    # add chromosome number
                awk '$2 = $2 FS "7"' ${INTERVAL_NAME_PREFIX}_LG${CHR}_int_named.txt | sed "s/ /\t/g" > ${INTERVAL_NAME_PREFIX}_LG${CHR}_int_named_POC_input.txt


                # merge all intervals files into one file
                cat *POC* > ${MAP_NAME}_allLGs_int_named_POC_input.txt

### Run Lep-Anchor
MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21
REF_GENOME=~/LM_new_genome/data/genome/IW24_minLen5kb_minQ10_contaminantsRemoved_polished_Racon3_medaka1_Pilon3_pugedDups.fasta
INTERVALS=~/LM_new_genome/results/8.LA/${MAP_NAME}/intervals/${MAP_NAME}_allLGs_int_named_POC_input.txt

~/software/Lep-Anchor/lepanchor_wrapper2_20Runs.sh -t 40 -T 20 -f ${REF_GENOME} -n 8 -m ${INTERVALS}
