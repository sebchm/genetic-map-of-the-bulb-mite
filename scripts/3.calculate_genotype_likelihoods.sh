#!/bin/bash
# calculate genotype likelihoods 

# prepare list of sorted bams
        # remove the list of bams if it exists
        rm -rf ~/LM_new_genome/bin/2.genotypeLikelihood/bam_list.txt

        # make an actual list of bams
for i in $(cat ~/LM_new_genome/data/sample_names_used_for_mapping.txt)
        do
        echo "/media/raid/home/schmielewski/LM_new_genome/results/1.mapping/${i}/${i}.dupmarked.bam" 
        done > ~/LM_new_genome/bin/2.genotypeLikelihood/bam_list.txt

# prepare list of contigs from the reference genome
grep ">" ~/LM_new_genome/data/genome/IW24_minLen5kb_minQ10_contaminantsRemoved_polished_Racon3_medaka1_Pilon3_pugedDups.fasta | awk '{print substr($1,2)}' > ~/LM_new_genome/bin/2.genotypeLikelihood/contigs.txt

BIN_PATH=~/LM_new_genome/bin/2.genotypeLikelihood
BAM_LIST=${BIN_PATH}/bam_list.txt
CONTIG_LIST=${BIN_PATH}/contigs.txt
PEDIGREE=~/LM_new_genome/data/pedigree_LM_new_genome.txt

cd ~/LM_new_genome/results/2.genotypeLikelihood/separate_contigs

for i in $(cat ${CONTIG_LIST})
do
        echo "samtools mpileup -r \"$i\" -q 20 -Q 20 -s \$(cat ${BAM_LIST}) | java -cp ~/software/Lep_Map3_0.5/bin Pileup2Likelihoods minAlleleFreq=0.05 mappingFile=~/LM_new_genome/data/sample_names_used_for_mapping.txt | gzip >\"$i\"_q20_Q20.mpileup.gz"

done > ${BIN_PATH}/SNP_calling_commands.txt
parallel --jobs 50 < ${BIN_PATH}/SNP_calling_commands.txt
zcat *.post.gz | awk '(NR==1 || ($1!="CHR"))'|gzip > all_contigs_q20_Q20_post.gz


DATA_GL=~/LM_new_genome/results/2.genotypeLikelihood/all_contigs_q20_Q20_post.gz
RESULTS_DIR_PC2=/media/raid/home/schmielewski/LM_new_genome/results/3.ParentCall2
RESULTS_DIR_F2=/media/raid/home/schmielewski/LM_new_genome/results/4.Filtering2
PEDIGREE=~/LM_new_genome/data/pedigree_LM_new_genome.txt
