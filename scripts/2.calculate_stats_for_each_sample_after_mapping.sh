#!/bin/bash
# after mapping, calculate per sample: 1) duplication rate; 2) mean mapping quality; 3) number of mapped reads

rm -f ~/LM_new_genome/bin/1.mapping/number_mapped_reads_per_sample.txt
for i in $(cat ~/LM_new_genome/data/sample_names_used_for_mapping.txt)
         do 
                grep "number of mapped reads" ~/LM_new_genome/results/1.mapping/${i}/${i}.coordinate_sorted_stats/genome_results.txt | sed "s/^/${i}\t/g" >> ~/LM_new_genome/bin/1.mapping/number_mapped_reads_per_sample.txt
         done 

### calculate duplication rate per sample
rm -f ~/LM_new_genome/bin/1.mapping/duplicationRate_per_sample.txt
for i in $(cat ~/LM_new_genome/data/sample_names_used_for_mapping.txt)
         do 
                grep "duplication rate" ~/LM_new_genome/results/1.mapping/${i}/${i}.coordinate_sorted_stats/genome_results.txt | sed "s/^/${i}\t/g" >> ~/LM_new_genome/bin/1.mapping/duplicationRate_per_sample.txt
         done

### calculate mean mapping quality per sample
rm -f ~/LM_new_genome/bin/1.mapping/meanMappingQuality_per_sample.txt
for i in $(cat ~/LM_new_genome/data/sample_names_used_for_mapping.txt)
         do 
                grep "mean mapping quality" ~/LM_new_genome/results/1.mapping/${i}/${i}.coordinate_sorted_stats/genome_results.txt | sed "s/^/${i}\t/g" >> ~/LM_new_genome/bin/1.mapping/meanMappingQuality_per_sample.txt
         done

### calculate number of reads per sample
rm -f ~/LM_new_genome/bin/1.mapping/numberOfReads_per_sample.txt
for i in $(cat ~/LM_new_genome/data/sample_names_used_for_mapping.txt)
         do 
                grep "number of reads = " ~/LM_new_genome/results/1.mapping/${i}/${i}.coordinate_sorted_stats/genome_results.txt | sed "s/^/${i}\t/g" >> ~/LM_new_genome/bin/1.mapping/numberOfReads_per_sample.txt
         done

### calculate coverage per sample
rm -f ~/LM_new_genome/bin/1.mapping/meanCoverage_per_sample.txt
for i in $(cat ~/LM_new_genome/data/sample_names_used_for_mapping.txt)
         do 
                grep "mean coverageData" ~/LM_new_genome/results/1.mapping/${i}/${i}.coordinate_sorted_stats/genome_results.txt | sed "s/^/${i}\t/g" >> ~/LM_new_genome/bin/1.mapping/meanCoverage_per_sample.txt
         done 
