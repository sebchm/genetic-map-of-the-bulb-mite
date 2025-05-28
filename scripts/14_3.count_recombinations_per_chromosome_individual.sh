
MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp
OUTPUT=~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_minMAF0.15_lodSC21/crossover_count_per_chromosome.txt
rm -f ${OUTPUT}

# for INFMASK in 123 13 23

INFMASK=13
        RESULTS_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_minMAF0.15_lodSC21/infMask${INFMASK}
            for CHR in 1 2 3 4 8
                do
                grep "recombines" ${RESULTS_DIR}/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData4/g" >> ${OUTPUT}

                grep "recombines" ${RESULTS_DIR}/outputPhasedData1/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_outputPhasedData1_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData1/g">> ${OUTPUT}
            done

            for CHR in 5 6; do
                grep "recombines" ${RESULTS_DIR}/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_proximityScale100_physical_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData4/g">> ${OUTPUT}
                grep "recombines" ${RESULTS_DIR}/outputPhasedData1/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_outputPhasedData1_proximityScale100_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData1/g" >> ${OUTPUT}
                done

INFMASK=23
        RESULTS_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_minMAF0.15_lodSC21/infMask${INFMASK}
            for CHR in 1 2 3 4 5 8
                do
                grep "recombines" ${RESULTS_DIR}/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData4/g" >> ${OUTPUT}
                grep "recombines" ${RESULTS_DIR}/outputPhasedData1/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_outputPhasedData1_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData1/g" >> ${OUTPUT}
            done

        CHR=6
                grep "recombines" ${RESULTS_DIR}/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_proximityScale100_physical_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData4/g" >> ${OUTPUT}
                grep "recombines" ${RESULTS_DIR}/outputPhasedData1/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_outputPhasedData1_proximityScale100_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData1/g" >> ${OUTPUT}
                

INFMASK=2
        RESULTS_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_NAME}_minMAF0.15_lodSC21/infMask${INFMASK}
        CHR=7
                grep "recombines" ${RESULTS_DIR}/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData4/g" >> ${OUTPUT}
                grep "recombines" ${RESULTS_DIR}/outputPhasedData1/${MAP_NAME}_minMAF0.15_infMask${INFMASK}_LG${CHR}_physical_outputPhasedData1_log.txt | cut -f3,5 | sed  "s/$/\t$CHR\tinfMask$INFMASK\tOutputPhasedData1/g" >> ${OUTPUT}

 
