MAP_PREFIX=data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15

for INFMASK in 123 13 23
        do
MAP_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_PREFIX}_lodSC21/infMask${INFMASK}/outputPhasedData1
mkdir -p ${MAP_DIR}/LMPlots

        for LG in 1 2 3 4 5 6 8
                do
java -cp ~/software/Lep_Map3_0.5/bin LMPlot ${MAP_DIR}/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.txt > ${MAP_DIR}/LMPlots/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.txt.dot
dot -Tpng ${MAP_DIR}/LMPlots/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.txt.dot > ${MAP_DIR}/LMPlots/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.txt.png
            done
        done

INFMASK=2
LG=7
MAP_DIR=~/LM_new_genome/results/10.OM2_after_LA/${MAP_PREFIX}_lodSC21/infMask${INFMASK}/outputPhasedData1
mkdir -p ${MAP_DIR}/LMPlots
java -cp ~/software/Lep_Map3_0.5/bin LMPlot ${MAP_DIR}/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.txt > ${MAP_DIR}/LMPlots/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.dot
dot -Tpng ${MAP_DIR}/LMPlots/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.dot > ${MAP_DIR}/LMPlots/${MAP_PREFIX}_infMask${INFMASK}_LG${LG}_physical_outputPhasedData1.png

