#!/bin/bash
# pipeline for trimming reads with Trimmomatic, mapping the reads with bwa mem, masking and removing duplicates with samtools
# uploaded by Sebastian Chmielewski (schmielewski_a_amu.edu.pl), AMU, Poland

# run this script with command:
  # while read r1 r2 r3 ; do source RhRob_TrimReads_MapReads_RmDuplicates_RmAmbReads.sh ${r1} ${r2} ${r3} ; done < sample_name_R1_R2_all_samples.txt,
  # where sample_name_R1_R2_all_samples.txt is a 3-column text file with sample name, path to FASTQ with forward reads and reverse reads, respectively


##prepare refernce genome before running (e.g. index) - to avoid doing it each time for each sample

FQPATH=(...)/fastq_files                                #<--- PATH TO FASTQ FILES
REFSEQ=(...)/Rhizoglyphus_robini_Genome_mod.fasta      #<--- INDEXED (SAMTOOLS), GATK DICTIONARY, BWA INDEX

ADAPTERS=TruSeq3-PE.fa
ADAPTERS_PATH=(...)/Trimmomatic-0.39/adapters

OUTDIR=(...)

###INFORMATION PROVIDED FROM COMMAND LINE:

IDNAME=$1     #INDIVIDUAL ID    $1
FASTQ_1=$2    #FASTQ R1 FILE    $2
FASTQ_2=$3    #FASTQ R2 FILE    $3


##COMPUTATIONAL OPTIONS (MEM, CPUs etc.)

NTHREADS=30

##LINKS TO SOFTWARE NOT IN $PATH
TRIMMOMATIC=(...)/Trimmomatic-0.39/trimmomatic-0.39.jar
#PICARD=~/software/picard/build/libs/picard.jar


##!!!!!!!!!!!!!
#########################
#######END OF EDITING####
#########################


echo Starting at `date`

#=======DIR FOR SAMPLE===============
##create tmp working directory and copy fastq files
if [ ! -d "${OUTDIR}/${IDNAME}" ]
then
    mkdir ${OUTDIR}/${IDNAME}
    mkdir ${OUTDIR}/${IDNAME}/tmp_${IDNAME}
fi

cd ${OUTDIR}/${IDNAME}


#====================================
#====TRIM READS WITH TRIMMOMATIC=====
#====================================
echo Starting trimming reads at `date`

cp ${ADAPTERS_PATH}/${ADAPTERS} .
java -jar ${TRIMMOMATIC} PE -threads ${NTHREADS} ${FASTQ_1} ${FASTQ_2} ${FASTQ_1}_trimmed.fastq_P.gz ${FASTQ_1}_trimmed.fastq_U.gz \
   ${FASTQ_2}_trimmed.fastq_P.gz ${FASTQ_2}_trimmed.fastq_U.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

# earlier version:
#cp ${ADAPTERS_PATH}/${ADAPTERS} .
#java -jar ${TRIMMOMATIC} PE -threads ${NTHREADS} ${FQPATH}/${FASTQ_1} ${FQPATH}/${FASTQ_2} ${FASTQ_1}_trimmed.fastq_P.gz ${FASTQ_1}_trimmed.fastq_U.gz \
#   ${FASTQ_2}_trimmed.fastq_P.gz ${FASTQ_2}_trimmed.fastq_U.gz ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 


echo Finished trimming reads at `date`

#====================================
#===========MAPPING==================
#====================================

echo Starting mapping at `date`

R1=${FASTQ_1}_trimmed.fastq_P.gz
R2=${FASTQ_2}_trimmed.fastq_P.gz

bwa mem -t ${NTHREADS} -R "@RG\tID:${IDNAME}\tSM:${IDNAME}\tLB:library1" ${REFSEQ} ${R1} ${R2} | samtools sort -@${NTHREADS} -n -O BAM -o ${IDNAME}.out_sorted.bam -
rm -f ./*.gz  #remove it or not
wc -l ${IDNAME}.out_sorted.bam > ${IDNAME}_n_lines_sorted_bam

#====================================
#=======ADD TAGS USING FIXMATE ======
#====================================
echo Starting adding tags for duplicate marking at `date`

samtools fixmate --threads ${NTHREADS} -m ${IDNAME}.out_sorted.bam ${IDNAME}.out_sorted.fixmate.bam
samtools sort ${IDNAME}.out_sorted.fixmate.bam --threads ${NTHREADS} -O BAM -o ${IDNAME}.out_sorted.fixmate_sorted.bam
rm ${IDNAME}.out_sorted.fixmate.bam

#====================================
#=======MARK DUPLICATES =============
#====================================
echo Starting marking duplicates at `date`

samtools markdup --threads ${NTHREADS} ${IDNAME}.out_sorted.fixmate_sorted.bam ${IDNAME}.dupmarked.bam
wc -l ${IDNAME}.dupmarked.bam > ${IDNAME}.dupmarked_n_lines

#==========================================
#=======REMOVE AMBIGOUS READS =============
#==========================================
echo Starting removing ambigous reads at `date`

samtools view --threads ${NTHREADS} -q 20 -bS ${IDNAME}.dupmarked.bam | samtools sort --threads ${NTHREADS} -o ${IDNAME}.dupmarked.AmbigRm.bam
rm ${IDNAME}.dupmarked.bam
wc -l MG* > ${IDNAME}.dupmarked.AmbigRm_n_lines

#======================
#=======INDEX BAM======
#======================

samtools index ${IDNAME}.dupmarked.AmbigRm.bam

echo Finished mapping at `date`

cd ~

