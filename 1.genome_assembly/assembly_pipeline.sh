#!/bin/bash
set -euo pipefail

# author: Sebastian Chmielewski, sebchm_@_amu.edu.pl
##############################
# Configuration
##############################

# Base directories
DATA_DIR=~/assembly/data/long_reads/IW24
MAP_DIR=~/assembly/data/long_reads/IW24/mapped_reads_contig_genome/withoutSecondaryReads
FLYE_DIR=~/assembly/flye/IW24_minLen5kb_minQ10_contaminantsRemoved
RACON_DIR=~/assembly/20.racon/IW24_minLen5kb_minQ10_contaminantsRemoved
PILON_BASE=~/assembly/22.pilon/IW24_minLen5kb_minQ10_contaminantsRemoved_polished
PURGE_DIR=~/assembly/18.purge_dups/IW24_minLen5kb_minQ10_contaminantsRemoved_polished
TE_DIR=~/assembly/23.EarlGrey_TE_annotation/host_data

# Common parameters
THREADS=50
MEMORY=100G
GENOME_SIZE=293m

# Input files
RAW_READS="${DATA_DIR}/IW24_K.fastq.gz"
PREV_ASSEMBLY=~/LM_new_bams/bin/1.creating_new_bams/genome/Rhizoglyphus_robini_Genome.fasta
TRIMMED_R1=~/assembly/data/short_reads/short_reads_trimmed/IW24/IW24_1_trimmed.fastq_P.gz
TRIMMED_R2=~/assembly/data/short_reads/short_reads_trimmed/IW24/IW24_2_trimmed.fastq_P.gz

# Output files
FILTERED_READS="${DATA_DIR}/IW24_minLen5kb_minQ10.fastq.gz"
MAPPED_BAM="${MAP_DIR}/IW24_minLen5kb_minQ10_mapped.bam"


#=============================================================
#== 1) filter out reads shorter than 5kb and meanQ < 10 ======
#=============================================================
filtlong --min_length 5000 --min_mean_q 10 ${RAW_READS} | gzip > ${FILTERED_READS}


#=======================================================
#== 2) map reads to previous version of the genome =====
#=======================================================
       cd ${MAP_DIR} 
       minimap2 -ax map-ont ${PREV_ASSEMBLY} ${FILTERED_READS} --secondary=no -t ${THREADS} | samtools view -b | samtools sort -o ${MAPPED_BAM}
       samtools index ${MAPPED_BAM}

       #== run qualimap
       ~/software/qualimap_v2.2.1/qualimap bamqc -bam ${MAPPED_BAM} -nt ${THREADS} --java-mem-size=${MEMORY}

       #== extract mapped reads
       samtools fastq ${MAPPED_BAM} | gzip > ${DATA_DIR}/ReadsThatMappedContigLevelAssembly_noSecondary_IW24_minLen5kb_minQ10.fastq.gz

#===========================================
#== 3) assemble the genome using Flye ======
#===========================================
  echo Running Flye at `date`

cd ${FLYE_DIR}
flye --nano-raw ${DATA_DIR}/ReadsThatMappedContigLevelAssembly_noSecondary_IW24_minLen5kb_minQ10.fastq.gz \
       --genome-size ${GENOME_SIZE} \
       --out-dir ${FLYE_DIR} \
       --threads ${THREADS}

#=================================
#== 4) polishing with Racon ======
#=================================
  echo Polishing with Racon `date`
READS=${DATA_DIR}/ReadsThatMappedContigLevelAssembly_noSecondary_IW24_minLen5kb_minQ10.fastq.gz
FASTA_PATH=${FLYE_DIR}/assembly.fasta
mkdir -p ${RACON_DIR}/{round1,round2,round3}

       # round 1
       cd ${RACON_DIR}/round1
       minimap2 -t ${THREADS} -x map-ont ${FASTA_PATH} ${READS} > raw_to_owl1.paf
       ~/software/racon/build/bin/racon -t ${THREADS} ${READS} raw_to_owl1.paf ${FASTA_PATH} > polished_racon1.fasta

       # round 2
       cd ${RACON_DIR}/round2
       minimap2 -t ${THREADS} -x map-ont ${RACON_DIR}/round1/polished_racon1.fasta ${READS} > racon1_mapped.paf
       ~/software/racon/build/bin/racon -t ${THREADS} ${READS} racon1_mapped.paf ${RACON_DIR}/round1/polished_racon1.fasta > polished_racon2.fasta

       # round 3
       cd ${RACON_DIR}/round3
       minimap2 -t ${THREADS} -x map-ont ${RACON_DIR}/round2/polished_racon2.fasta ${READS} > racon2_mapped.paf
       ~/software/racon/build/bin/racon -t ${THREADS} ${READS} racon2_mapped.paf ${RACON_DIR}/round2/polished_racon2.fasta > polished_racon3.fasta


#===================
#== 5) Medaka ======
#===================
  echo Running Medaka at `date`
DRAFT_GENOME=${RACON_DIR}/round3/polished_racon3.fasta
medaka_consensus -i ${READS} -d ${DRAFT_GENOME} -o "${RACON_DIR}/medaka_output" -t ${THREADS}


#==================
#== 6) Pilon ======
#==================
  echo Running Pilon at `date`
# Round 1
       cd ${PILON_BASE}/prepare_short_reads_for_pilon
       REFSEQ=${PILON_BASE}/../21.medaka/consensus.fasta
       NTHREADS=${THREADS}

       # map to the genome after 3 rounds of racon and 1 medaka
       bwa mem -t ${NTHREADS} ${REFSEQ} ${TRIMMED_R1} ${TRIMMED_R2} | samtools sort -@ ${NTHREADS} -n -O BAM -o pilot_round1.unsorted.bam -
       # add tags using fixmate
       samtools fixmate --threads ${NTHREADS} -m pilot_round1.unsorted.bam pilot_round1.fixmate.bam
       # sort bam
       samtools sort --threads ${NTHREADS} -O BAM -o pilot_round1.fixmate_sorted.bam pilot_round1.fixmate.bam
       # mark duplicates
       samtools markdup --threads ${NTHREADS} pilot_round1.fixmate_sorted.bam pilot_round1.dedup.bam
       # index bam
       samtools index pilot_round1.dedup.bam

       java -Xmx${MEMORY} -jar ~/software/pilon/pilon-1.24.jar --genome ${REFSEQ} --fix all --changes --frags pilot_round1.dedup.bam --threads ${NTHREADS} --output ${PILON_BASE}/pilon_round1

# Round 2
       REFSEQ=${PILON_BASE}/pilon_round1/pilon_round1.fasta
       cd ${PILON_BASE}/pilon_round2

       # map to the genome after 3 rounds of racon and 1 medaka and 1 pilon
       bwa mem -t ${NTHREADS} ${REFSEQ} ${TRIMMED_R1} ${TRIMMED_R2} | samtools sort -@ ${NTHREADS} -n -O BAM -o pilot_round2.unsorted.bam -
       # add tags using fixmate
       samtools fixmate --threads ${NTHREADS} -m pilot_round2.unsorted.bam pilot_round2.fixmate.bam
       # sort bam
       samtools sort --threads ${NTHREADS} -O BAM -o pilot_round2.fixmate_sorted.bam pilot_round2.fixmate.bam
       # mark duplicates
       samtools markdup --threads ${NTHREADS} pilot_round2.fixmate_sorted.bam pilot_round2.dedup.bam
       # index bam
       samtools index pilot_round2.dedup.bam

       java -Xmx${MEMORY} -jar ~/software/pilon/pilon-1.24.jar --genome ${REFSEQ} --fix all --changes --frags pilot_round2.dedup.bam --threads ${NTHREADS} --output ${PILON_BASE}/pilon_round2 --tracks

# Round 3
       REFSEQ=${PILON_BASE}/pilon_round2/pilon_round2.fasta
       cd ${PILON_BASE}/pilon_round3

       # map to the genome after 3 rounds of racon and 1 medaka and 2 pilon
       bwa mem -t ${NTHREADS} ${REFSEQ} ${TRIMMED_R1} ${TRIMMED_R2} | samtools sort -@ ${NTHREADS} -n -O BAM -o pilot_round3.unsorted.bam -
       # add tags using fixmate
       samtools fixmate --threads ${NTHREADS} -m pilot_round3.unsorted.bam pilot_round3.fixmate.bam
       # sort bam
       samtools sort --threads ${NTHREADS} -O BAM -o pilot_round3.fixmate_sorted.bam pilot_round3.fixmate.bam
       samtools markdup --threads ${NTHREADS} pilot_round3.fixmate_sorted.bam pilot_round3.dedup.bam
       # index bam
       samtools index pilot_round3.dedup.bam

       java -Xmx${MEMORY} -jar ~/software/pilon/pilon-1.24.jar --genome ${REFSEQ} --fix all --changes --frags pilot_round3.dedup.bam --threads ${NTHREADS} --output ${PILON_BASE}/pilon_round3 --tracks

#======================
#== 7) PURGEDUPS ======
#======================
  echo Running purgeDups at `date`
conda activate rr_map_env
ASSEMBLY_NAME=IW24_minLen5kb_minQ10_contaminantsRemoved_polished
mkdir -p ${PURGE_DIR}/${ASSEMBLY_NAME} && cd ${PURGE_DIR}/${ASSEMBLY_NAME}
mv ${PILON_BASE}/pilon_round3/pilon_round3.fasta ${PILON_BASE}/pilon_round3/${ASSEMBLY_NAME}.fasta # rename pilon output
ASSEMBLY_FASTA=${PILON_BASE}/pilon_round3/${ASSEMBLY_NAME}.fasta

       # 7.1 map reads to the flye assembly
       minimap2 -x map-ont -t ${READS} ${ASSEMBLY_FASTA} ${DATA_DIR}/ReadsThatMappedContigLevelAssembly_noSecondary_IW24_minLen5kb_minQ10.fastq.gz gzip -c > ont_to_asm.paf.gz
       
       # 7.2 calculate stats
       pbcstat ont_to_asm.paf.gz
       calcuts PB.stat > cutoffs 2> calcults.log
       
       # 7.3 plit the assembly
       split_fa ${ASSEMBLY_FASTA} > assembly_split.fasta
       minimap2 -t ${READS} -x asm5 -DP assembly_split.fasta assembly_split.fasta | gzip -c > assembly_split.self.paf.gz

       # 7.4 Purge haplotigs and overlaps
       purge_dups -2 -T cutoffs -c PB.base.cov assembly_split.self.paf.gz > dups.bed 2> purge_dups.log

       get_seqs -e dups.bed ${ASSEMBLY_FASTA}
       python3 ~/software/purge_dups/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png

#=======================
#== 8) EarlGreyTE ======
#=======================
  echo Running EarlGrey annotation at `date`
conda activate singularity_env
cd "${TE_DIR}"
singularity shell -C -H "$(pwd):/work" --writable-tmpfs -u earlgrey.sif
earlGrey -g ${PILON_BASE}/pilon_round3/${ASSEMBLY_NAME}.fasta -s rhizoglyphusRobini -t 45 -o ${TE_DIR}/IW24_minLen5kb_minQ10_contaminantsRemoved_polished
