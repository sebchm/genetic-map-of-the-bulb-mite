#!/bin/bash
set -e
set -u
set -o pipefail

# this pipeline was used to create the genetic map of the bulb mite, Rhizoglyphus robini. 
# author: Sebastian Chmielewski, sebchm@amu.edu.pl

#=======================================================================================================================================
#== 1. trim the Illumina reads, map them to the genome, mark duplicates, remove ambigously mapping reads and run QC with qualimap ======
#=======================================================================================================================================
while read r1 r2 r3 ; do source 1.rr_map_removeDuplicates.sh ${r1} ${r2} ${r3} ; done < ~/LM_new_genome/data/sample_name_R1_R2_usedForMapping.txt

#==========================================================
#== 2. calculate stats for each sample after mapping ======
#==========================================================
bash 2.calculate_stats_for_each_sample_after_mapping.sh

#===========================================
#== 3. calculate genotype likelihoods ======
#===========================================
bash 3.calculate_genotype_likelihoods.sh

#========================
#== 4. run ParentCall2 ==
#========================
bash 4.ParentCall2.sh

#========================
#== 5. run Filtering 2 ==
#========================
bash 5.Filtering2.sh

#============================
#== 6. remove repeats =======
#============================
bash 6.removeRepeats.sh

#===============================================
#== 7. prepare post file without repeats ==
#===============================================
bash 7.removeRepeats_from_postGL_file.sh

#====================================
#== 8. SeparateChromosomes (slurm) ==
#====================================

for LOD in {18..30}
        do
                sbatch separateChromosomes2_slurm.sh
        done

#==============================
#== 9. OrderMarkers2 (slurm) ==
#==============================

# order autosomes
for LG in 1 2 3 4 5 6 8;
  do
  sbatch orderMarkers_autosomes_infMask123_LOD21_10hours.sh ${LG}
  done
  
# order sex chromosome (infMask2)
srun orderMarkers_LG7_infMask2_LOD21_10hours.sh

# order LG8 (it had a misjoin- corrected by increasing number of iterations)
srun orderMarkers_LG8_infMask123_LOD21_10hours.sh

#=========================
#== 10. Lep-Anchor =======
#=========================
bash 10.LepAnchor.sh

#=================================================
#== 11. prepare physical files after Lep-Anchor ==
#=================================================
bash 11.preparePhysicalFiles_after_LepAnchor.sh

#==============================================================
#== 12. Order markers after anchoring with OutputPhasedData4 ==
#==============================================================
12.orderMarkers_after_LA_outputPhasedData4.sh

#==============================================================
#== 13. Order markers after anchoring with OutputPhasedData1 ==
#==============================================================
        # only for diagnostics to draw LMPlots
13.orderMarkers_after_LA_outputPhasedData1.sh

#=========================
#== 14. QUALITY CONTROL ==
#=========================
        # draw LMPlots
14_1.draw_LMPlots.sh
        # genotyping errors tend to be moved to map ends and they inflate its length
14_2.check_EdgeMarkers_MapInflation.R 
        # count number of crossovers for eachc chromosome- high numbers indicate errors
14_3.count_recombinations_per_chromosome_individual.sh

#========================
#== 15. Final ordering ==
#========================
        # problematic LGs were ordered with ProximityScale=100
        # this script outputs final maps
15.FINAL_SCRIPT_order_markers_after_LA_outputPhasedData4.sh

#======================
#== 16. Map analysis ==
#======================
map_qc_and_summary_pipeline.R


