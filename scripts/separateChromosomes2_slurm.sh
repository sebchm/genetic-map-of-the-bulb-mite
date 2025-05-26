#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40      
#SBATCH --mem=20gb                
#SBATCH --time=9:00:00     
#SBATCH --job-name=LM_SC_maps  
#SBATCH --mail-type=ALL,TIME_LIMIT

module load java8/jdk1.8.0_40

MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15
DATA=${MAP_NAME}.gz
LM=$(echo java -cp software/LepMap3/bin)

mkdir ../../results/SC2/${MAP_NAME}/raw_maps
mkdir ../../results/SC2/${MAP_NAME}/sorted_maps

lod_SC=$1
zcat ${DATA}| ${LM} SeparateChromosomes2 data=- lodLimit=${lod_SC} numThreads=40 > ../../results/SC2/${MAP_NAME}/raw_maps/${MAP_NAME}_lodSC${lod_SC}.txt

    # sort the maps
sort ../../results/SC2/${MAP_NAME}/raw_maps/${MAP_NAME}_lodSC${lod_SC}.txt | uniq -c | sort -n > ../../results/SC2/${MAP_NAME}/sorted_maps/${MAP_NAME}_lodSC${lod_SC}_sorted.txt
