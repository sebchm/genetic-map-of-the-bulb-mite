
### check what proportion of genome is covered by contigs with markers:
contigs_markers <- read.table("~/LM_new_genome/bin/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15/number_of_markers_per_contig.txt", col.names = c("n_markers", "contig"))
contig_size <- read.table("~/LM_new_genome/data/genome/IW24_minLen5kb_minQ10_contaminantsRemoved_polished_Racon3_medaka1_Pilon3_pugedDups.fasta.fai", col.names = c("contig", "contig_length", "C", "D", "E"))[,1:2]

genome_size <- sum(contig_size$contig_length)

contigs <- left_join(contigs_markers, contig_size, by = "contig")

size_contigs_with_markers <- sum(contigs$contig_length)

proportion_contigsWithMarkers <- size_contigs_with_markers/genome_size
proportion_contigsWithMarkers
