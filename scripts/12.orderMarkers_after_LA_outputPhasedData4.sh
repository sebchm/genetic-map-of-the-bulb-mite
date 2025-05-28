#! /bin/bash

#### MASTER SCRIPT FOR ORDERING AFTER ANCHORING
# OUTPUT PHASED DATA4

LM=$(echo java -cp ~/software/Lep_Map3_0.51/bin)
MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15
DATA=~/LM_new_genome/results/5.RemovingRepeats/famInfLim4_complementedFamLim3/${MAP_NAME}.gz
#--------------
# OrderMarkers2 after anchoring: sex-Avg map
#--------------

### FINAL ORDERING OF THE MAPS
MAP=~/LM_new_genome/results/6.SC2/dataTol0.05/raw_maps/${MAP_NAME}_lodSC21.txt
SNP_FILE=~/LM_new_genome/results/5.RemovingRepeats/famInfLim4_complementedFamLim3/snps_${MAP_NAME}.txt
PHYS_FILE_DIR=~/software/Lep-Anchor/new_genome/${MAP_NAME}_lodSC21

        ##################
        ### INFMASK123 ###
        ##################

INFMASK=123
RESULTS_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_lodSC21/infMask${INFMASK}
mkdir -p ${RESULTS_DIR}/intervals

for CHR in 1 2 3 4 5 6 8
          do
echo "gunzip -c ${DATA} | ${LM}  OrderMarkers2 data=- map=${MAP} useMorgan=1 chromosome=${CHR} numThreads=40 \
        outputPhasedData=4 \
        numMergeIterations=40 \
        grandparentPhase=1 \
        sexAveraged=1 \
        informativeMask=${INFMASK} \
        improveOrder=0 \
        evaluateOrder=${PHYS_FILE_DIR}/order${CHR}.phys \
        calculateIntervals=${RESULTS_DIR}/intervals/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_int_physical \
        > ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical.txt \
        2> ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_log.txt"
        done|parallel --jobs 30


# name the markers
for CHR in 1 2 3 4 5 6 8
                do
        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical.txt > \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_named.txt
                done


        #############################
        ### INFMASK13 & INFMASK23 ###
        #############################
for INFMASK in 13 23
        do
RESULTS_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_lodSC21/infMask${INFMASK}
mkdir -p ${RESULTS_DIR}/intervals
for CHR in 1 2 3 4 5 6 8
          do
echo "gunzip -c ${DATA} | ${LM}  OrderMarkers2 data=- map=${MAP} useMorgan=1 chromosome=${CHR} numThreads=40 \
        outputPhasedData=4 \
        numMergeIterations=40 \
        grandparentPhase=1 \
        improveOrder=0 \
        informativeMask=${INFMASK} \
        evaluateOrder=${PHYS_FILE_DIR}/order${CHR}.phys \
        calculateIntervals=${RESULTS_DIR}/intervals/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_int_physical \
        > ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical.txt \
        2> ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_log.txt"
        done|parallel --jobs 30

# name the markers
for CHR in 1 2 3 4 5 6 8
                do
        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical.txt > \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_named.txt
                done

        done

################################################
#### LG7: sex chromosome. InfMask=2 only #####
################################################
LM=$(echo java -cp ~/software/Lep_Map3_0.51/bin)

INFMASK=2
RESULTS_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_lodSC21/infMask${INFMASK}
mkdir -p ${RESULTS_DIR}/intervals

CHR=7
gunzip -c ${DATA} | ${LM}  OrderMarkers2 data=- map=${MAP} useMorgan=1 chromosome=${CHR} numThreads=40 \
        outputPhasedData=4 \
        numMergeIterations=40 \
        grandparentPhase=1 \
        informativeMask=${INFMASK} \
        improveOrder=0 \
        evaluateOrder=${PHYS_FILE_DIR}/order${CHR}.phys \
        calculateIntervals=${RESULTS_DIR}/intervals/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_int_physical \
        > ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical.txt \
        2> ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_log.txt

# name the markers
        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical.txt > \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_named.txt

##########################
#### PREPARE SUMMARY #####
##########################
# prepare map summary- sexAvg map
cd ~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_lodSC21
rm -f summary_all_maps.txt

for INFMASK in 123 23
        do
      for CHR in 1 2 3 4 5 6 8
                do
        tail -n+4 infMask${INFMASK}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_named.txt |  \
        cut -f1,2,4 | sed -e "s/$/\t${CHR}\tInfMask${INFMASK}/" >> ~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_lodSC21/summary_all_maps.txt
                done

done
INFMASK=13
      for CHR in 1 2 3 4 5 6 8
                do
        tail -n+4 infMask${INFMASK}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_named.txt |  \
        cut -f1,2,3 | sed -e "s/$/\t${CHR}\tInfMask${INFMASK}/" >> ~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_lodSC21/summary_all_maps.txt
                done
        
# sex chromosome with InfMask2
INFMASK=2
CHR=7
        tail -n+4 infMask${INFMASK}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_named.txt |  \
        cut -f1,2,4 | sed -e "s/$/\t${CHR}\tInfMask${INFMASK}/" >> ~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_lodSC21/summary_all_maps.txt
