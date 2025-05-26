#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40      
#SBATCH --mem=65gb
#SBATCH --time=10:00:00           
#SBATCH --job-name="OM2"
#SBATCH --mail-type=ALL,TIME_LIMIT
  
module load openjdk/17.0.11_9-gcc-14.2.0

MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15
DATA=/mnt/storage_3/home/sebchm/pl0288-01/project_data/SC/data/${MAP_NAME}.gz
LM=$(echo java -cp /mnt/storage_3/home/sebchm/pl0288-01/project_data/SC/software/LepMap3/bin)
lod_SC=21
OUTPUT_DIR=../../results/OM2/${MAP_NAME}_lodSC${lod_SC}
MAP=../../results/SC2/${MAP_NAME}/raw_maps/${MAP_NAME}_lodSC${lod_SC}.txt
SNP_FILE=/mnt/storage_3/home/sebchm/pl0288-01/project_data/SC/data/snps_${MAP_NAME}.txt

#--------------
# OrderMarkers2
#--------------
CHR=8
INFMASK=123

RESULTS_DIR=${OUTPUT_DIR}/infMask${INFMASK}
mkdir -p ${RESULTS_DIR}
mkdir -p ${RESULTS_DIR}/intervals

gunzip -c ${DATA} | ${LM}  OrderMarkers2 data=- map=${MAP} useMorgan=1 chromosome=$CHR numThreads=40 \
        outputPhasedData=1 \
        numMergeIterations=60 \
        grandparentPhase=1 \
        sexAveraged=1 \
        informativeMask=${INFMASK} \
        calculateIntervals=${RESULTS_DIR}/intervals/${MAP_NAME}_LOD${lod_SC}_infMask${INFMASK}_LG${CHR}_60iterations_int \
        > ${RESULTS_DIR}/${MAP_NAME}_LOD${lod_SC}_infMask${INFMASK}_LG${CHR}_60iterations.txt \
        2> ${RESULTS_DIR}/${MAP_NAME}_LOD${lod_SC}_infMask${INFMASK}_LG${CHR}_60iterations_log.txt

        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
        ${RESULTS_DIR}/${MAP_NAME}_LOD${lod_SC}_infMask${INFMASK}_LG${CHR}_60iterations.txt > \
        ${RESULTS_DIR}/${MAP_NAME}_LOD${lod_SC}_infMask${INFMASK}_LG${CHR}_60iterations_named.txt
