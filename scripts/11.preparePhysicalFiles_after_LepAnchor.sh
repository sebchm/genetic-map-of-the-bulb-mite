#!/bin/bash
LA_PATH=~/software/Lep-Anchor
MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15
RESULTS_PATH=~/software/Lep-Anchor/new_genome/${MAP_NAME}_lodSC21
SNP_FILE=~/LM_new_genome/results/5.RemovingRepeats/famInfLim4_complementedFamLim3/snps_${MAP_NAME}.txt
INTERVALS_PATH=~/LM_new_genome/results/8.LA/${MAP_NAME}_lodSC21/intervals

cd ${RESULTS_PATH}

for LG in 1 2 3 4 5 6 8
        do

awk -f ${LA_PATH}/liftover.awk ${RESULTS_PATH}/chr${LG}.agp ${INTERVALS_PATH}/${MAP_NAME}_LOD21_infMask123_LG${LG}_int_named_POC_input.txt |sort -V | grep LG > ${RESULTS_PATH}/order${LG}.liftover
awk -vinverse=1 -f ${LA_PATH}/liftover.awk ${RESULTS_PATH}/chr${LG}.agp ${RESULTS_PATH}/order${LG}.liftover|awk '(NR==FNR){m[$1"\t"($2+0)]=NR-1}(NR!=FNR){print m[$1"\t"($2+0)]}' ${SNP_FILE} - > ${RESULTS_PATH}/order${LG}.phys
        done

# LG7 (sex chromosome)
LG=7
awk -f ${LA_PATH}/liftover.awk ${RESULTS_PATH}/chr${LG}.agp ${INTERVALS_PATH}/${MAP_NAME}_LOD21_infMask2_LG${LG}_int_named_POC_input.txt |sort -V | grep LG > ${RESULTS_PATH}/order${LG}.liftover
awk -vinverse=1 -f ${LA_PATH}/liftover.awk ${RESULTS_PATH}/chr${LG}.agp ${RESULTS_PATH}/order${LG}.liftover|awk '(NR==FNR){m[$1"\t"($2+0)]=NR-1}(NR!=FNR){print m[$1"\t"($2+0)]}' ${SNP_FILE} - > ${RESULTS_PATH}/order${LG}.phys

# merge phys and liftover files:
for LG in {1..8}
                do
paste ${RESULTS_PATH}/order${LG}.phys ${RESULTS_PATH}/order${LG}.liftover > ${RESULTS_PATH}/order${LG}_phys_liftover.txt
                done
# merge all LGs into one file:
cat ${RESULTS_PATH}/order*_phys_liftover.txt > ${RESULTS_PATH}/order_allLGs_phys_liftover.txt

# change spaces to tabs and remove unused columns. The output columns are: 1st column: in 1st row is a 1st marker in 1st LG and the value is a line number from the snps_lifover file, then LG, position on LG and LG again
sed "s/ /\t/g" ${RESULTS_PATH}/order_allLGs_phys_liftover.txt | cut -f1,3,4 > ${RESULTS_PATH}/order_allLGs_phys_liftover.txt2
# change name of the file:
mv ${RESULTS_PATH}/order_allLGs_phys_liftover.txt2 ${RESULTS_PATH}/order_allLGs_phys_liftover.txt

     awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
        order_allLGs_phys_liftover.txt > order_allLGs_phys_liftover_named.txt
