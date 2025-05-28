#! /bin/bash
# pipeline used for constructing genetic map of bulb mite and anchoring the genome assembly 
# author: Sebastian Chmielewski, Evolutionary Biology Group, Adam Mickiewicz University, Poland
  # schmielewski_a_amu.edu.pl
  
#-------------------------------
# calculate genotype likelihoods 
#-------------------------------

UNMASKED_REGIONS=(...)/Rhizoglyphus_robini_Genome_mod_UNMASKED_PART_GENOME.bed
RESULTS_DIR=(...)
LM=$(echo java -cp (...)/Lep-Map3/bin/)

for i in $(cat contig_names.txt)
	do
echo "samtools mpileup -r \"$i\" -l $UNMASKED_REGIONS -q 20 -Q 10 -s \$(cat sorted_bams) | ${LM} Pileup2Likelihoods |gzip >\"$i\".post.gz"
	done > SNP_calling.txt
parallel --jobs 30 < SNP_calling.txt
zcat *.post.gz | awk '(NR==1 || ($1!="CHR"))'|gzip > all_post.gz

#------------
# ParentCall2 
#------------
	echo LM3 ParentCall2 started at `date`
zcat ${DATA_GL} | cut -f95 --complement | ${LM} ParentCall2 data=pedigree.txt posteriorFile=- removeNonInformative=1 XLimit=2 | gzip > ${RESULTS_DIR}/data.ParentCall.gz 2> ${RESULTS_DIR}/data.ParentCall.err

#--------------
# Filtering2
#--------------

	 echo LM3 Filtering2 started at `date`
zcat ${RESULTS_DIR}/data.ParentCall.gz | ${LM} Filtering2 data=- dataTolerance=0.05 MAFLimit=0.15 removeNonInformative=1 familyInformativeLimit=4 missingLimit=0.1 |gzip> ${RESULTS_DIR}/data.Filtering_dataTol0.05_familyInfLimit4_MAF0.15_missingLim0.1.gz
zcat ${RESULTS_DIR}/data.ParentCall.gz | ${LM} Filtering2 data=- dataTolerance=0.05 MAFLimit=0.15 removeNonInformative=1 familyInformativeLimit=3 missingLimit=0.1 |gzip> ${RESULTS_DIR}/data.Filtering_dataTol0.05_familyInfLimit3_MAF0.15_missingLim0.1.gz

	# using a custom R script, only markers informative in 4 families were selected. In a case when contig had less than 5 markers, markers informative in at least 3 families were selected. 

#---------------------
# SeparateChromosomes2
#---------------------

for lod_SC in {15..20}
        do
zcat ${RESULTS_DIR}/data_dataTol0.05_famInfLimit4_complfamInfLim3_MAF0.15_missingLim0.1.gz | ${LM} SeparateChromosomes2 data=- lodLimit=${lod_SC} sizeLimit=10 > ${RESULTS_DIR}/data.Filtering_dataTol0.05_familyInfLimit3_MAF0.15_missingLim0.1_lodSC${lod_SC}_sizeLim10.txt
        # sort raw maps
sort ${RESULTS_DIR}/data.Filtering_dataTol0.05_familyInfLimit3_MAF0.15_missingLim0.1_lodSC${lod_SC}_sizeLim10.txt | uniq -c | sort -n > ${RESULTS_DIR}/data.Filtering_dataTol0.05_familyInfLimit3_MAF0.15_missingLim0.1_lodSC${lod_SC}_sizeLim10_sorted.txt
        done

#--------------
# OrderMarkers2
#--------------
MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit3_MAF0.15_missingLim0.1_lodSC18_sizeLim10
SNP_FILE=snps_${MAP_NAME}.txt

for CHR in {1..8}
        do
                if [ "$CHR" -eq 4 ]; then # skip the sex chromosome
                continue 
                fi
echo "gunzip -c ${DATA} | ${LM} OrderMarkers2 data=- map=${MAP} useKosambi=1 chromosome=${CHR} numThreads=60 \
        outputPhasedData=1 \
        grandparentPhase=1 \
        sexAveraged=1 \
        numMergeIterations=60 \
        informativeMask=123 \ # set to 2 for the sex chromosome
        calculateIntervals=${RESULTS_DIR}/intervals/${MAP_NAME}_LG${CHR}_int \
        > ${RESULTS_DIR}/${MAP_NAME}_LG${CHR}.txt \
        2> ${RESULTS_DIR}/${MAP_NAME}_LG${CHR}_log.txt"
        done|parallel --jobs 50

        # retrieve marker position and contig
        for CHR in {1..8}
                        if [ "$CHR" -eq 4 ]; then # skip the sex chromosome
                        continue 
                        fi
        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
        ${RESULTS_DIR}/${MAP_NAME}_LG${CHR}.txt > \
        ${RESULTS_DIR}/${MAP_NAME}_LG${CHR}_contigPos.txt

			# preparemap summary for R input
	       tail -n+4 $${RESULTS_DIR}/${MAP_NAME}_LG${CHR}_contigPos.txt |  \
	       cut -f1,2,3,4 | sed -e "s/$/\t${CHR}\tInfMask123/" >> ${RESULTS_DIR}/map_infMask123_summary.txt
	                done

#----------------
# GenomeAnchoring
#----------------

MAP_NAME=data.Filtering_dataTol0.05_familyInfLimit3_MAF0.15_missingLim0.1_lodSC18_sizeLim10
REF_GENOME=(...)/Rhizoglyphus_robini_Genome_mod.fasta.gz
CHAIN_FILE=all.chain.gz
INTERVALS=${MAP_NAME}_infMask123_allLGs_int_named_POC_input.txt # marker intervals with columns: contig position_on_contig chromosome interval_start interval_end
PAF_FILE=alnONT.paf

~/software/Lep-Anchor/lepanchor_wrapper2_20Runs.sh -t 40 -T 20 -f ${REF_GENOME} -n 8 -p ${PAF_FILE} -c ${CHAIN_FILE} -m ${INTERVALS}

#-------------------------------------------------------------
# Prepare liftover file with marker coordinates on chromosomes
#-------------------------------------------------------------

#! /bin/bash
LA_PATH=~/software/Lep-Anchor/
MAP_NAME=dataTol0.05_familyInfLimit4_complemented_famInfLim3_MAF0.15_missingLim0.1_MG235rm_infMask123
RESULTS_PATH=~/software/Lep-Anchor/${MAP_NAME}
SNP_FILE=snps_${MAP_NAME}.txt

cd ${RESULTS_PATH}
for LG in $(seq 1 8)
        do
awk -f ${LA_PATH}/liftover.awk ${RESULTS_PATH}/chr${LG}.agp (...)/intervals/${MAP_NAME}_LG${LG}_int_named_POC_input.txt |sort -V| grep LG > ${RESULTS_PATH}/order${LG}.liftover
awk -vinverse=1 -f ${LA_PATH}/liftover.awk ${RESULTS_PATH}/chr${LG}.agp ${RESULTS_PATH}/order${LG}.liftover|awk '(NR==FNR){m[$1"\t"($2+0)]=NR-1}(NR!=FNR){print m[$1"\t"($2+0)]}' ${SNP_FILE} - > ${RESULTS_PATH}/order${LG}.phys 
        done

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

#-----------------------------------------------
# OrderMarkers2 after anchoring: sex-specific map
#-----------------------------------------------
SNP_FILE=snps_${MAP_NAME}.txt
PHYS_FILE_DIR=(...)
MAP=${MAP_NAME}/${MAP_NAME}_lodSC18_sizeLim10.txt


for INFMASK in 13 23
        do
RESULTS_DIR=infMask${INFMASK}
mkdir -p ${RESULTS_DIR}
mkdir -p ${RESULTS_DIR}/intervals

for CHR in {1..8}
        do
                if [ "$CHR" -eq 4 ]; then # skip the sex chromosome
                continue 
                fi

echo "gunzip -c ${DATA} | ${LM}  OrderMarkers2 data=- map=${MAP} useKosambi=1 chromosome=${CHR} numThreads=40 \
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

# assign marker coordinates
for CHR in {1..8}
        do
                if [ "$CHR" -eq 4 ]; then # skip the sex chromosome
                continue 
                fi
        awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${SNP_FILE} \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical.txt > \
        ${RESULTS_DIR}/${MAP_NAME}_infMask${INFMASK}_LG${CHR}_physical_named.txt
                done
        done
