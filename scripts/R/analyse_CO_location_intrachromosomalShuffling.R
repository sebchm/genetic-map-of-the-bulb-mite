library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(viridis)

# this script allows to extract positions of crossovers based on OrderMarkers2 (Lep-Map3) output. Also, it calculates values of intrachromosomal shuffling, following Veller et al. 2019.

#=================================
#== 1. paternal CO analysis ======
#=================================

 # 1) prepare names of F2 individuals from each family:
   # P2 family
   TEMP <- scan(text="MG051     MG052   MG053   MG054   MG055   MG056   MG057   MG058   MG060   MG061   MG062   MG063   MG064   MG065   MG066   MG067   MG068   MG069   MG070   MG071   MG072   MG073   MG074   MG075   MG076   MG077   MG078   MG079   MG080   MG081   MG082", quiet = TRUE, what="")
   P2_individuals_names <- c(paste(dput(TEMP), "father", sep = "_"), paste(dput(TEMP), "mother", sep = "_"))

   # P10 family
    TEMP <- scan(text="MG110    MG111   MG112   MG113   MG114   MG115   MG116   MG117   MG118   MG119   MG120   MG121   MG122   MG123   MG124   MG125   MG126   MG127   MG128   MG129   MG130   MG131   MG132   MG133   MG134   MG135   MG136   MG137   MG138   MG139   MG140   MG141", quiet = TRUE, what="")
   P10_individuals_names <-  c(paste(dput(TEMP), "father", sep = "_"), paste(dput(TEMP), "mother", sep = "_"))
   
   # P6 family
    TEMP <- scan(text="MG222    MG223   MG224   MG225   MG226   MG227   MG228   MG229   MG230   MG231   MG232   MG233   MG234   MG236   MG237   MG238   MG239   MG240   MG241   MG242   MG243   MG244   MG245   MG246   MG247   MG248   MG249   MG250   MG251   MG252   MG253   MG254   MG255   MG256   MG257   MG258   MG259   MG260   MG261   MG262   MG263   MG264   MG265   MG266   MG267   MG268   MG269   MG270   MG271   MG272   MG273   MG274   MG275   MG276   MG277   MG278   MG279   MG280   MG281", what="", quiet = TRUE)
   P6_individuals_names <- c(paste(dput(TEMP), "father", sep = "_"), paste(dput(TEMP), "mother", sep = "_"))

   #P8_2 family:
    TEMP <- scan(text="MG282    MG283   MG284   MG285   MG286   MG287   MG288   MG289   MG290   MG291   MG292   MG293   MG294   MG296   MG297   MG299   MG300   MG301   MG302   MG303   MG304   MG305   MG306   MG307   MG308   MG309   MG310   MG311   MG312   MG313   MG314   MG315   MG316   MG317   MG318   MG319   MG320   MG321   MG322   MG323   MG324   MG325   MG326   MG327   MG328   MG329   MG330   MG331   MG332   MG333   MG334   MG335   MG336   MG337   MG338   MG339   MG340", what="", quiet = TRUE)
   P8_2_individuals_names <- c(paste(dput(TEMP), "father", sep = "_"), paste(dput(TEMP), "mother", sep = "_"))

   # initialise dfs: 
datalist = data.frame()
haplotypes_i = data.frame()
number_of_recombinations_i = data.frame()

for (chr in c(1,2,3,4,8)){ # exclude LG5 and LG6 (with proximityScale=100) and LG7 (sexChr) 
chr_temp <- read.table(paste0("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/infMask13/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_infMask13_LG",chr,"_physical_named.txt"), colClasses=c("factor", "factor", "numeric", "numeric", "factor", "factor", "factor",rep(c("character", "numeric", "numeric"),4)))[,-c(3,5:7)]

colnames(chr_temp) <- c("contig", "position_contig", "genetic_position", "P2_individuals", "P2_paternal_recombinations", "P2_maternal_recombinations", "P10_individuals", "P10_paternal_recombinations", "P10_maternal_recombinations", "P6_individuals", "P6_paternal_recombinations", "P6_maternal_recombinations", "P8_2_individuals", "P8_2_paternal_recombinations", "P8_2_maternal_recombinations")

chr_temp$position_contig <- sub('\\*$', '', chr_temp$position_contig)
chr_temp$position_contig <- as.numeric(chr_temp$position_contig)
  
#### 4) split every two digits into separate columns:
## 4.1) add underscore after every second digit:
chr_temp$P2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P2_individuals, perl = TRUE)
chr_temp$P10_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P10_individuals, perl = TRUE)
chr_temp$P6_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P6_individuals, perl = TRUE)
chr_temp$P8_2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P8_2_individuals, perl = TRUE)

# split columns to separate individuals
chr_temp <-  chr_temp %>%
   separate(P2_individuals, sep = "_", into = P2_individuals_names) %>%
   separate(P10_individuals, sep = "_", into = P10_individuals_names) %>%
   separate(P6_individuals, sep = "_", into = P6_individuals_names) %>%
   separate(P8_2_individuals, sep = "_", into = P8_2_individuals_names) %>%
  mutate(chromosome = chr, 
         physical_pos_no = row_number()) 

# chnge format to longer
chr_temp_longer <- chr_temp %>%
  pivot_longer(cols  = starts_with("MG"),
              names_to = c("individual", "parent"),
              names_sep = "_",
              values_to = "phase") %>%
  select(-ends_with("recombinations")) %>%
  filter(parent == "father")

# there are some uncretain haplotypes: replace them with adjacent values:
chr_temp_longer_noUncertain <- chr_temp_longer %>%
  mutate(phase = ifelse(phase == "-", NA, phase)) %>%
  group_by(chromosome, individual) %>%
  fill(phase, .direction = "updown")

indv_list <- unique(chr_temp_longer_noUncertain$individual)
haplotypes <- data.frame()

for (i in indv_list){
  #for (p in c("father", "mother")) {
  p="father"
  chr_temp_longer_noUncertain_i <- chr_temp_longer_noUncertain[chr_temp_longer_noUncertain$individual == i,]
  
  rle <- rle(chr_temp_longer_noUncertain_i$phase)
  
  t <- data.frame(lengths = rle$lengths, 
                  values = rle$values, 
                 individual = i, 
                 parent = p, 
                 LG = chr)
  
  haplotypes <- bind_rows(haplotypes, t)
  
}

number_of_recombinations <- haplotypes %>%
  group_by(individual, parent) %>%
  tally() %>%
  mutate(chromosome = chr, 
         n = n - 1)

haplotypes_i <- bind_rows(haplotypes_i, haplotypes)
number_of_recombinations_i <- bind_rows(number_of_recombinations_i,  number_of_recombinations)
# datalist[[chr]] <- chr_temp_longer
datalist = bind_rows(chr_temp_longer_noUncertain, datalist)

}


### ADD LG5 and LG6:
for (chr in c(5,6)){  
#chr=7
chr_temp <- read.table(paste0("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/infMask13/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_infMask13_LG",chr,"_proximityScale100_physical_named.txt"), colClasses=c("factor", "factor", "numeric", "numeric", "factor", "factor", "factor",rep(c("character", "numeric", "numeric"),4)))[,-c(3,5:7)]

colnames(chr_temp) <- c("contig", "position_contig", "genetic_position", "P2_individuals", "P2_paternal_recombinations", "P2_maternal_recombinations", "P10_individuals", "P10_paternal_recombinations", "P10_maternal_recombinations", "P6_individuals", "P6_paternal_recombinations", "P6_maternal_recombinations", "P8_2_individuals", "P8_2_paternal_recombinations", "P8_2_maternal_recombinations")

chr_temp$position_contig <- sub('\\*$', '', chr_temp$position_contig)
chr_temp$position_contig <- as.numeric(chr_temp$position_contig)
  
#### 4) split every two digits into separate columns:
## 4.1) add underscore after every second digit:
chr_temp$P2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P2_individuals, perl = TRUE)
chr_temp$P10_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P10_individuals, perl = TRUE)
chr_temp$P6_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P6_individuals, perl = TRUE)
chr_temp$P8_2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P8_2_individuals, perl = TRUE)

# split columns to separate individuals
chr_temp <-  chr_temp %>%
   separate(P2_individuals, sep = "_", into = P2_individuals_names) %>%
   separate(P10_individuals, sep = "_", into = P10_individuals_names) %>%
   separate(P6_individuals, sep = "_", into = P6_individuals_names) %>%
   separate(P8_2_individuals, sep = "_", into = P8_2_individuals_names) %>%
  mutate(chromosome = chr, 
         physical_pos_no = row_number()) 

# chnge format to longer
chr_temp_longer <- chr_temp %>%
  pivot_longer(cols  = starts_with("MG"),
              names_to = c("individual", "parent"),
              names_sep = "_",
              values_to = "phase") %>%
  select(-ends_with("recombinations")) %>%
  filter(parent == "father")

# there are some uncretain haplotypes: replace them with adjacent values:
chr_temp_longer_noUncertain <- chr_temp_longer  %>%
   mutate(phase = ifelse(phase == "-", NA, phase)) %>%
   group_by(chromosome, individual) %>%
   fill(phase, .direction = "updown")

indv_list <- unique(chr_temp_longer_noUncertain$individual)
haplotypes <- data.frame()

for (i in indv_list){
  #for (p in c("father", "mother")) {
  p="father"
  chr_temp_longer_noUncertain_i <- chr_temp_longer_noUncertain[chr_temp_longer_noUncertain$individual == i,]
  
  rle <- rle(chr_temp_longer_noUncertain_i$phase)
  
  t <- data.frame(lengths = rle$lengths, 
                  values = rle$values, 
                 individual = i, 
                 parent = p, 
                 LG = chr)
  
  haplotypes <- bind_rows(haplotypes, t)
  
}

number_of_recombinations <- haplotypes %>%
  group_by(individual, parent) %>%
  tally() %>%
  mutate(chromosome = chr, 
         n = n - 1)

haplotypes_i <- bind_rows(haplotypes_i, haplotypes)
number_of_recombinations_i <- bind_rows(number_of_recombinations_i,  number_of_recombinations)
# datalist[[chr]] <- chr_temp_longer
datalist = bind_rows(chr_temp_longer_noUncertain, datalist)

}

# end of loading files for paternal CO analysis


#=================================
#== 1. maternal CO analysis ======
#=================================

for (chr in c(1,2,3,4,5,8)){ # exclude LG6 (proximityScale=100) and LG7 (sex chromosome) 
# chr=2
chr_temp <- read.table(paste0("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/infMask23/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_infMask23_LG",chr,"_physical_named.txt"), colClasses=c("factor", "factor", "numeric", "numeric", "factor", "factor", "factor",rep(c("character", "numeric", "numeric"),4)))[,-c(3,5:7)]

colnames(chr_temp) <- c("contig", "position_contig", "genetic_position", "P2_individuals", "P2_paternal_recombinations", "P2_maternal_recombinations", "P10_individuals", "P10_paternal_recombinations", "P10_maternal_recombinations", "P6_individuals", "P6_paternal_recombinations", "P6_maternal_recombinations", "P8_2_individuals", "P8_2_paternal_recombinations", "P8_2_maternal_recombinations")

chr_temp$position_contig <- sub('\\*$', '', chr_temp$position_contig)
chr_temp$position_contig <- as.numeric(chr_temp$position_contig)

  
#### 4) split every two digits into separate columns:
## 4.1) add underscore after every second digit:
chr_temp$P2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P2_individuals, perl = TRUE)
chr_temp$P10_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P10_individuals, perl = TRUE)
chr_temp$P6_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P6_individuals, perl = TRUE)
chr_temp$P8_2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P8_2_individuals, perl = TRUE)

# split columns to separate individuals
chr_temp <-  chr_temp %>%
   separate(P2_individuals, sep = "_", into = P2_individuals_names) %>%
   separate(P10_individuals, sep = "_", into = P10_individuals_names) %>%
   separate(P6_individuals, sep = "_", into = P6_individuals_names) %>%
   separate(P8_2_individuals, sep = "_", into = P8_2_individuals_names) %>%
  mutate(chromosome = chr, 
         physical_pos_no = row_number()) 

# chnge format to longer
chr_temp_longer <- chr_temp %>%
  pivot_longer(cols  = starts_with("MG"),
              names_to = c("individual", "parent"),
              names_sep = "_",
              values_to = "phase") %>%
  select(-ends_with("recombinations")) %>%
  filter(parent == "mother")

# there are some uncretain haplotypes: replace them with adjacent values:
chr_temp_longer_noUncertain <- chr_temp_longer %>%
  mutate(phase = ifelse(phase == "-", NA, phase)) %>%
  group_by(chromosome, individual) %>%
  fill(phase, .direction = "updown")

indv_list <- unique(chr_temp_longer_noUncertain$individual)
haplotypes <- data.frame()

for (i in indv_list){
  p="mother"
  chr_temp_longer_noUncertain_i <- chr_temp_longer_noUncertain[chr_temp_longer_noUncertain$individual == i,]
  
  rle <- rle(chr_temp_longer_noUncertain_i$phase)
  
  t <- data.frame(lengths = rle$lengths, 
                  values = rle$values, 
                 individual = i, 
                 parent = p, 
                 LG = chr)
  
  haplotypes <- bind_rows(haplotypes, t)
  
}

number_of_recombinations <- haplotypes %>%
  group_by(individual, parent) %>%
  tally() %>%
  mutate(chromosome = chr, 
         n = n - 1)

haplotypes_i <- bind_rows(haplotypes_i, haplotypes)
number_of_recombinations_i <- bind_rows(number_of_recombinations_i, number_of_recombinations )
# datalist[[chr]] <- chr_temp_longer
datalist = bind_rows(chr_temp_longer_noUncertain, datalist)


}

##### ADD LG7- sex chromosome ####
chr=7
chr_temp <- read.table(paste0("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/infMask2/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_infMask2_LG7_physical_named.txt"), colClasses=c("factor", "factor", "numeric", "numeric", "factor", "factor", "factor",rep(c("character", "numeric", "numeric"),4)))[,-c(3,5:7)]

colnames(chr_temp) <- c("contig", "position_contig", "genetic_position", "P2_individuals", "P2_paternal_recombinations", "P2_maternal_recombinations", "P10_individuals", "P10_paternal_recombinations", "P10_maternal_recombinations", "P6_individuals", "P6_paternal_recombinations", "P6_maternal_recombinations", "P8_2_individuals", "P8_2_paternal_recombinations", "P8_2_maternal_recombinations")

chr_temp$position_contig <- sub('\\*$', '', chr_temp$position_contig)
chr_temp$position_contig <- as.numeric(chr_temp$position_contig)

#### 4) split every two digits into separate columns:
## 4.1) add underscore after every second digit:
chr_temp$P2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P2_individuals, perl = TRUE)
chr_temp$P10_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P10_individuals, perl = TRUE)
chr_temp$P6_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P6_individuals, perl = TRUE)
chr_temp$P8_2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P8_2_individuals, perl = TRUE)

# split columns to separate individuals
chr_temp <-  chr_temp %>%
   separate(P2_individuals, sep = "_", into = P2_individuals_names) %>%
   separate(P10_individuals, sep = "_", into = P10_individuals_names) %>%
   separate(P6_individuals, sep = "_", into = P6_individuals_names) %>%
   separate(P8_2_individuals, sep = "_", into = P8_2_individuals_names) %>%
  mutate(chromosome = chr, 
         physical_pos_no = row_number()) 

# chnge format to longer
chr_temp_longer <- chr_temp %>%
  pivot_longer(cols  = starts_with("MG"),
              names_to = c("individual", "parent"),
              names_sep = "_",
              values_to = "phase") %>%
  select(-ends_with("recombinations")) %>%
  filter(parent == "mother")

# there are some uncretain haplotypes: replace them with adjacent values:
chr_temp_longer_noUncertain <- chr_temp_longer %>%
   mutate(phase = ifelse(phase == "-", NA, phase)) %>%
   group_by(chromosome, individual) %>%
   fill(phase, .direction = "updown")

indv_list <- unique(chr_temp_longer_noUncertain$individual)
haplotypes <- data.frame()

for (i in indv_list){
  #for (p in c("father", "mother")) {
  p="mother"
  chr_temp_longer_noUncertain_i <- chr_temp_longer_noUncertain[chr_temp_longer_noUncertain$individual == i,]
  
  rle <- rle(chr_temp_longer_noUncertain_i$phase)
  
  t <- data.frame(lengths = rle$lengths, 
                  values = rle$values, 
                 individual = i, 
                 parent = p, 
                 LG = chr)
  
  haplotypes <- bind_rows(haplotypes, t)
  
}

number_of_recombinations <- haplotypes %>%
  group_by(individual, parent) %>%
  tally() %>%
  mutate(chromosome = chr, 
         n = n - 1)

haplotypes_i <- bind_rows(haplotypes_i, haplotypes)
number_of_recombinations_i <- bind_rows(number_of_recombinations_i, number_of_recombinations )
datalist = bind_rows(chr_temp_longer_noUncertain, datalist)
### end of section: add chromosome 7


##### ADD chromosome 6 ####
chr=6
chr_temp <- read.table("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/infMask23/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_infMask23_LG6_proximityScale100_physical_named.txt", colClasses=c("factor", "factor", "numeric", "numeric", "factor", "factor", "factor",rep(c("character", "numeric", "numeric"),4)))[,-c(3,5:7)]

colnames(chr_temp) <- c("contig", "position_contig", "genetic_position", "P2_individuals", "P2_paternal_recombinations", "P2_maternal_recombinations", "P10_individuals", "P10_paternal_recombinations", "P10_maternal_recombinations", "P6_individuals", "P6_paternal_recombinations", "P6_maternal_recombinations", "P8_2_individuals", "P8_2_paternal_recombinations", "P8_2_maternal_recombinations")

chr_temp$position_contig <- sub('\\*$', '', chr_temp$position_contig)
chr_temp$position_contig <- as.numeric(chr_temp$position_contig)

#### 4) split every two digits into separate columns:
## 4.1) add underscore after every second digit:
chr_temp$P2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P2_individuals, perl = TRUE)
chr_temp$P10_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P10_individuals, perl = TRUE)
chr_temp$P6_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P6_individuals, perl = TRUE)
chr_temp$P8_2_individuals <- gsub("(.)(?!$)", "\\1_", chr_temp$P8_2_individuals, perl = TRUE)


# split columns to separate individuals
chr_temp <-  chr_temp %>%
   separate(P2_individuals, sep = "_", into = P2_individuals_names) %>%
   separate(P10_individuals, sep = "_", into = P10_individuals_names) %>%
   separate(P6_individuals, sep = "_", into = P6_individuals_names) %>%
   separate(P8_2_individuals, sep = "_", into = P8_2_individuals_names) %>%
  mutate(chromosome = chr, 
         physical_pos_no = row_number()) 

# chnge format to longer
chr_temp_longer <- chr_temp %>%
  pivot_longer(cols  = starts_with("MG"),
              names_to = c("individual", "parent"),
              names_sep = "_",
              values_to = "phase") %>%
  select(-ends_with("recombinations")) %>%
  filter(parent == "mother")

# there are some uncretain haplotypes: replace them with adjacent values:
chr_temp_longer_noUncertain <- chr_temp_longer %>%
  mutate(phase = ifelse(phase == "-", NA, phase)) %>%
  group_by(chromosome, individual) %>%
  fill(phase, .direction = "updown")

indv_list <- unique(chr_temp_longer_noUncertain$individual)
haplotypes <- data.frame()

for (i in indv_list){
  #for (p in c("father", "mother")) {
  p="mother"
  chr_temp_longer_noUncertain_i <- chr_temp_longer_noUncertain[chr_temp_longer_noUncertain$individual == i,]
  
  rle <- rle(chr_temp_longer_noUncertain_i$phase)
  
  t <- data.frame(lengths = rle$lengths, 
                  values = rle$values, 
                 individual = i, 
                 parent = p, 
                 LG = chr)
  
  haplotypes <- bind_rows(haplotypes, t)
  
}

number_of_recombinations <- haplotypes %>%
  group_by(individual, parent) %>%
  tally() %>%
  mutate(chromosome = chr, 
         n = n - 1)

haplotypes_i <- bind_rows(haplotypes_i, haplotypes)
number_of_recombinations_i <- bind_rows(number_of_recombinations_i, number_of_recombinations )
# datalist[[chr]] <- chr_temp_longer
datalist = bind_rows(chr_temp_longer_noUncertain, datalist)
### end of loading files for mother CO analysis

#===================================
#== 3. add family information ======
#===================================

sample_info <- read.table("~/LM_new_genome/data/sample_info_with_morph_based_on_perigree_P2_P10_P6_P8_2_corrected.txt", h = T)

# sample info stores information about: family; individual_name; generation; sex; morph (binary phenotypic trait)
  # family individual generation sex morph
  # P2 MG001 P female NA
  # P2 MG002 P male fighter
  # P2 MG015 F1 female NA

fam_indv <- sample_info %>%
  select(individual, family)

datalist <- datalist %>%
      left_join(fam_indv, by = "individual")

# change LG labels to chromosomes (linkage groups are ordered by marker number and chromosomes by their physical size)
datalist <- datalist %>%
  mutate(chromosome = recode(as.character(chromosome), # adjust LGs to match sizes:
                        `6` = "chr1",
                        `3` = "chr2",
                        `8` = "chr3",
                        `2` = "chr4",
                        `1` = "chr5",
                        `7` = "chr6",
                        `4` = "chr7",
                        `5` = "chr8"),
    chromosome = factor(chromosome, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")))

number_of_recombinations_i <- number_of_recombinations_i %>%
  mutate(chromosome = recode(as.character(chromosome), # adjust LGs to match sizes:
                        `6` = "chr1",
                        `3` = "chr2",
                        `8` = "chr3",
                        `2` = "chr4",
                        `1` = "chr5",
                        `7` = "chr6",
                        `4` = "chr7",
                        `5` = "chr8"),
    chromosome = factor(chromosome, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")))



#================================
#== 4. chromosome painting ======
#================================

for (chr in 1:8)
p1 <- datalist %>%
    filter(chromosome == chr) %>%
    ggplot(aes(physical_pos_no, individual, color = phase)) +
    geom_tile() +
    facet_grid(family~parent, scales = "free_y") +
    scale_color_viridis(discrete = TRUE) +
    ggtitle(paste0("no uncertain haplotypes, chromosome LG " , chr)) +
    theme(axis.ticks = element_blank())

plot(p1)
}



sample_info <- read.table("~/LM_new_genome/data/sample_info_with_morph_based_on_perigree_P2_P10_P6_P8_2_corrected.txt", h = T)

# count number of CO per each parent per chromosome:
number_of_recombinations_all <- number_of_recombinations_i %>%
  left_join(sample_info, by = "individual") %>%
  select(-generation)

# sum number of COs per each parent:
number_of_recombinations_summedSex <- number_of_recombinations_all %>%
  group_by(family, parent, individual) %>%
  summarise(sum_CO = sum(n))

number_of_recombinations_summedSex <- number_of_recombinations_summedSex %>%
  group_by(parent, sum_CO) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(height = ifelse(parent == "father", -count, count))

# plot number of COs per individual (meiosis) for each parent (n = 8)
CO_per_individual_plot <- ggplot(number_of_recombinations_summedSex, aes(x = sum_CO, y = height, fill = parent)) +
  geom_bar(stat = "identity", position = "identity", width = 1, color = "black") +
  scale_y_continuous(labels = abs) +
  labs(y = "count", x = "crossover count") +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#2ab1ca", "#eb6e50"))

# group CO count by family
N_CO_per_meiosis_plot <- number_of_recombinations_all %>%
  group_by(family, parent, individual) %>%
  summarise(sum_CO = sum(n)) %>%
  ggplot(aes(sum_CO)) +
  geom_histogram(binwidth = 1) +
  facet_grid(parent ~ family) +
  ggtitle("Number of CO per meiosis (N=8 individuals)")

N_CO_per_chr_plot <- number_of_recombinations_all %>%
  group_by(family, parent, individual, chromosome) %>%
  summarise(sum_CO = sum(n)) %>%
  ggplot(aes(chromosome, sum_CO)) +
  geom_col(binwidth = 1) +
  facet_grid(parent ~ family) +
  ggtitle("Number of CO per meiosis (N=8 individuals)")

# add number of COs per chromosome, even no CO occured
  # make a df with all individuals and chromosomes
expand_table <- expand_grid(individual = indv_list,
                           parent = c("mother", "father"),
                           chromosome = paste0("chr", seq(1:8)))


number_of_recombinations_i <- number_of_recombinations_i %>%
  full_join(expand_table, by = c("individual", "parent", "chromosome")) %>%
  select(1:4) %>%
  left_join(sample_info, by = "individual") %>%
  select(-generation) %>%
  mutate(parent_id = paste(parent, family, sep = "_")) %>%
  mutate(n = replace_na(n, 0))

# plot number of recombinations per chromosome for each individual

N_COs_per_indv_plot <- number_of_recombinations_i %>%
  ggplot(aes(individual, n, fill = parent)) +
  geom_col(position = "dodge") +
  facet_grid(chromosome ~ family, scales = "free_x")

number_of_recombinations_i %>%
  group_by(parent) %>%
  summarise(sum(n))

### sum number of COs per individual
umber_of_recombinations_i %>%
  group_by(individual, parent) %>%
  summarise(sum(n))

#===================================
#== 5. obtain crossover positions ==
#===================================
# add liftovered position of a marker
snp_lifotver <- read.table("~/software/Lep-Anchor/new_genome/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/order_allLGs_phys_liftover_named.txt", col.names = c("contig", "position_contig", "position_chromosome", "chromosome"))

snp_lifotver$position_contig <- sub('\\*$', '', snp_lifotver$position_contig)
snp_lifotver$position_contig <- as.numeric(snp_lifotver$position_contig)

# change LG labels to chromosomes:
snp_lifotver <- snp_lifotver %>%
  mutate(chromosome = recode(as.character(chromosome), # adjust LGs to match sizes:
                        `6` = "chr1",
                        `3` = "chr2",
                        `8` = "chr3",
                        `2` = "chr4",
                        `1` = "chr5",
                        `7` = "chr6",
                        `4` = "chr7",
                        `5` = "chr8"),
    chromosome = factor(chromosome, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")))

# these data contains location of crossovers. However note that the precise location of CO sometimes cannot be resolved, because outputPhasedData was set to 1. outputPhasedData=4 masks uncertain phases, but in such cases it's not possible to localise COs.
  # note that the CO happened between the marker in the data and previous marker. However, the CO correlations were analysed in big windows (1Mb - 2 Mb), so it wouldn't affect the results.

datalist_COs <- datalist %>%
  left_join(snp_lifotver, by = c("contig", "position_contig", "chromosome")) %>% # add the liftovered position of the markers
  arrange(individual, parent, chromosome, position_chromosome) %>% # sorting
  group_by(individual, parent, chromosome) %>%
  mutate(phase = as.integer(phase), # change from character to int
    phase_lagged = ifelse( # lag the values of phase. If the parent or chromosome changes, don't take the phase from the previous parent/chromosome
      lag(chromosome, default = chromosome[1]) != chromosome |
      lag(parent, default = parent[1]) != parent,
      NA,
      lag(phase))) %>%
  fill(phase_lagged, .direction = "up") %>% # for the first lagged values: fill with next phase value
  mutate(phase_diff = phase_lagged - phase) %>% # variable to detect CO
  filter(phase_diff != 0)

# count crossovers in each chromosome and parent
datalist_COs %>% 
  group_by(individual, parent, chromosome) %>%
  tally()


CO_location <- datalist_COs %>%
  select(individual, parent, chromosome, position_chromosome) %>%
  arrange(chromosome, position_chromosome)

CO_location %>%
  group_by(chromosome, position_chromosome) %>%
  tally(name = "n_CO") %>%
  ggplot(aes(position_chromosome, n_CO)) +
  geom_point() +
  facet_wrap(~chromosome, scales = "free_x")

# save CO locations:
# write.table(CO_location, "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21//map_analysis/CO_location/all_chr_CO_location.txt", row.names = F, quote = F, sep = "\t")

# plot distribution of COs
CO_location %>%
  group_by(chromosome, position_chromosome, parent) %>%
  tally(name = "n_CO") %>%
  ggplot(aes(position_chromosome, n_CO, fill = parent)) +
  geom_col(width = 400000) +
  facet_wrap(~chromosome, scales = "free_x") +
  ylab("Number of crossovers")
  
n_COs_mother <- CO_location %>%
  filter(parent == "mother") %>%
  ungroup() %>%
  tally() %>%
  pull()

n_COs_father <- CO_location %>%
  filter(parent == "father") %>%
  ungroup() %>%
  tally() %>%
  pull()

paste0(n_COs_mother, " COs in females and ", n_COs_father, " COs in males were detected.")

#=============================================
#== 6. calculate intrachromosomal shuffling ==
#=============================================

# Chromosomal shuffling was calculated following the method outlined by Veller et al. (2019). Initially, chromosomes were "painted" to assess the proportion of parental haplotype. Subsequently, intrachromosomal shuffling was computed for each chromosome and averaged across each sex. To determine whether differences in intrachromosomal (IA) shuffling between sexes were attributable to variations in CO localization or to differences in overall recombination rates, only chromosomes with at least one CO were considered. Finally, the IA shuffling value was divided by the number of crossovers to normalize the data.

# load chromosome length:
chromosome_physical_length <- read.table("~/software/Lep-Anchor/new_genome/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21//REF_LA.fa.fai")[1:8,1:2]
colnames(chromosome_physical_length) <- c("LG", "physical_length")
chromosome_physical_length$LG <- as.integer(gsub("LG", "", chromosome_physical_length$LG))

chromosome_physical_length <- chromosome_physical_length %>%
  mutate(LG = recode(as.character(LG), # adjust LGs to match sizes:
                        `6` = "chr1",
                        `3` = "chr2",
                        `8` = "chr3",
                        `2` = "chr4",
                        `1` = "chr5",
                        `7` = "chr6",
                        `4` = "chr7",
                        `5` = "chr8"),
    LG = factor(LG, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")))
# write.table(datalist, "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21//map_analysis/CO_location/datalist_CO_location.txt", quote = F, row.names = F, sep = "\t")

first_last_phase_markers <- datalist %>%
  left_join(snp_lifotver, by = c("contig", "position_contig", "chromosome")) %>% # add chromosome position of a marker 
  select(-c(contig, position_contig, family, physical_pos_no)) %>% # remove unused columns
  ungroup() %>% # reset grouping
  group_by(individual, chromosome, parent) %>% 
  mutate(change = phase != lag(phase, default = first(phase))) %>% # detect phase change
  mutate(group_id = cumsum(change)) %>% # each phase gets individual number 
  group_by(individual, chromosome, parent, group_id, phase) %>%
  summarise( # calculate start and end of each phase within a chromosome
    phase_start = min(position_chromosome),
    phase_end = max(position_chromosome),
    phase_length = phase_end - phase_start + 1,
    .groups = 'drop') %>%
  left_join(chromosome_physical_length, by = c("chromosome" = "LG")) %>% # add chromosome length
  rename(chromosome_length = physical_length) %>%
  group_by(individual, chromosome, parent, group_id, phase) %>%
  mutate(phase_proportion = phase_length/chromosome_length) # calculate proportion of phase in relation to chromosome length

# write.table(first_last_phase_markers, "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21//map_analysis/CO_location/first_last_phase_markers.txt", quote = F, row.names = F, sep = "\t")

# visualisation to check correctness:
  # nb. that chromosome parts without markers (where the phase is unknown) are not included (see Veller et al. 2019, fig. 5, black chromosome parts)

# supplementary figure with chromosome paiting

painted_chromosomes_plot <- first_last_phase_markers %>%
ggplot(aes(xmin = phase_start/1e6, xmax = phase_end/1e6, y = individual, color = phase)) +
  geom_linerange(linewidth = 1.1) +
  facet_grid(parent ~ chromosome, scales = "free_x") +
  scale_color_manual(values = c("#83878b", "#eb6e50", "#2ab1ca")) +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = "none", 
        text = element_text(size = 20)) +
  labs(x = "physical position [Mb]")

#saveRDS(painted_chromosomes_plot, "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/map_analysis/CO_location/painted_chromosomes_plot.rds")
# ggsave("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/map_analysis/CO_location/plots/painted_chromosomes_plot.png", IA_shuffling_plot)


# prepare a sex-specific histogram with CO distance to chromosome end
CO_dist_to_chr_end <- CO_location %>%
  left_join(chromosome_physical_length, by = c("chromosome" = "chr")) %>% # add chromosome size
  mutate(chromosome_middle = physical_length / 2, 
         distance_chr_end = ifelse(position_chromosome < chromosome_middle, position_chromosome, physical_length - position_chromosome), 
         chromosome = ifelse(chromosome == "chr6", "chr6 (sex)", chromosome), 
         parent = ifelse(parent == "mother", "female map", "male map"))
  
plot_CO_dist_chr_end <- ggplot(CO_dist_to_chr_end, aes(distance_chr_end/1e6, fill = parent)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~chromosome, nrow = 2) +
  theme(legend.position = "top", 
        legend.title = element_blank()) +
  scale_fill_manual(values = c("#2ab1ca", "#eb6e50")) +
  labs(y = "number of COs", x = "distance to chromosome end [Mb]")

#saveRDS(plot_CO_dist_chr_end, "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/map_analysis/CO_location/plot_CO_dist_chr_end.rds")
# ggsave("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/map_analysis/CO_location/plots/plot_CO_dist_chr_end.png", plot_CO_dist_chr_end, width = 7, height = 5)


# calculate genome size:
genome_size <- sum(chromosome_physical_length$physical_length)

# calculate p (proportion of phase p). Keep in mind that the proportion of alternative phase is not 1-p, because the phase is unknown for some chromosome parts (at the beginning and end, or between phase changes)
phase_length <- first_last_phase_markers %>%
  group_by(individual, chromosome, parent, phase) %>%
  summarise(sum_phase_proportion = sum(phase_proportion)) %>%
  pivot_wider(id_cols = c(individual, chromosome, parent), 
              names_from = phase, 
              values_from = sum_phase_proportion, 
              names_prefix = "p_phase_", 
              values_fill = 0) %>%
  left_join(chromosome_physical_length, by = c("chromosome" = "LG")) %>% # add chromosome length
  rename(chromosome_length = physical_length) %>%
  mutate(L_sq = (chromosome_length/genome_size)^2, # squared proportion of a chromosome within a genome
         intra_shuffling_chr = 2 * p_phase_0 * p_phase_1 * L_sq) # first term of eq. 2;, Veller et al. 2019

# intrachromosomal shuffling: first term of eq. 2
intrachromosomal_shuffling_indv <- phase_length %>%
 group_by(individual, parent) %>%
  summarise(intra_shuffling = sum(intra_shuffling_chr)) # summed; first term of eq. 2

  # average within sex:
intrachromosomal_shuffling_sex <- intrachromosomal_shuffling_indv %>%
  group_by(parent) %>%
  summarise(mean_intra_shuffling = mean(intra_shuffling))

# interchromosomal shuffling: second term of eq. 2
interchromosomal_shuffling <- chromosome_physical_length %>%
  mutate(L = physical_length / genome_size) %>%
  summarise(0.5 * (1 - sum(L^2))) %>%
  pull()
  
#### calculate intrachromosomal shuffling for each chromosome:
IA_shuffling_chrom <- phase_length %>%
  ungroup() %>%
  group_by(chromosome, parent) %>%
  summarise(mean_intra_shuffling_chr = mean(intra_shuffling_chr, na.rm = TRUE), .groups = 'drop') 

IA_shuffling_plot <- IA_shuffling_chrom %>%
  mutate(chromosome = gsub("chr", "", chromosome)) %>%
ggplot(aes(chromosome, mean_intra_shuffling_chr, fill = parent)) +
    geom_col(position = position_dodge2(preserve = "single", padding = 0))+
  scale_fill_manual(values = c("#eb6e50", "#2ab1ca")) +
  theme_bw() +
  theme(legend.position = "none", 
        text = element_text(size = 10), 
          plot.tag.position = c(0.03, 1.05), 
        panel.grid.minor = element_blank()) +
  labs(tag = "B", x = "chromosome", y = "intrachromosomal shuffling") 

# sex-averaged IA shuffling:

IA_shuffling_chrom %>%
  summarise(mean(mean_intra_shuffling_chr))

IA_shuffling_chrom %>%
  group_by(parent) %>%
  summarise(mean(mean_intra_shuffling_chr))

# remove sex chromosome:
IA_shuffling_chrom %>%
  filter(chromosome != "chr6") %>%
  group_by(parent) %>%
  summarise(mean(mean_intra_shuffling_chr))


# scale IA shuffling by number of COs: first version (with removing chromosomes without COs)
  # get number of COs:
n_COs_IA_shuff <- first_last_phase_markers %>% 
  group_by(individual, chromosome, parent) %>% 
  tally(name = "n_COs") %>%
  mutate(n_COs = n_COs - 1)

  # remove chromosomes without CO:
std_IA_shuffling <- phase_length %>%
  filter(intra_shuffling_chr != 0) %>%
  left_join(n_COs_IA_shuff, by = c("individual", "chromosome", "parent")) %>%
  mutate(std_intra_shuffling_chr = intra_shuffling_chr / n_COs)

std_IA_shuffling %>%
  group_by(parent) %>%
  summarise(mean(std_intra_shuffling_chr, na.rm = T))

# scale IA shuffling by number of COs: second version (without removing chromosomes without COs)
n_COs_sex <- CO_location %>%
  group_by(parent) %>%
  summarise(sum_COs = sum(n()))

intrachromosomal_shuffling_sex_nCOs <- n_COs_sex %>% 
  left_join(intrachromosomal_shuffling_sex, by = "parent") %>%
  mutate(mean_intra_shuffling_scaled_by_COs = mean_intra_shuffling/sum_COs)


# check correlation between IA shuffling and chromosome length. 
  # expected is a pattern where longer chromosomes recombine less in their middle parts 
chromosome_physical_length2 <- chromosome_physical_length %>%
  mutate(LG = recode(as.character(LG), # adjust LGs to match sizes:
                        `3` = "chr1",
                        `1` = "chr2",
                        `5` = "chr3",
                        `4` = "chr4",
                        `2` = "chr5",
                        `8` = "chr6 (sex)",
                        `7` = "chr7",
                        `6` = "chr8"),
    LG = factor(LG, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6 (sex)", "chr7", "chr8"))) 


# IA VS chromosome size:
IA_shuffling_chrom_l <- IA_shuffling_chrom %>%
  #mutate(chromosome = as.integer(chromosome)) %>%
    left_join(chromosome_physical_length, by = c("chromosome" = "LG")) 

IA_shuffling_chrom_l %>%
  ggplot(aes(physical_length, mean_intra_shuffling_chr, color = parent)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)


LM3 <- lm(mean_intra_shuffling_chr ~ physical_length * parent + I(physical_length^2), data = IA_shuffling_chrom_l)
summary(LM3)
drop1(LM3, test = "F") # remove interactio 

LM3.2 <- lm(mean_intra_shuffling_chr ~ physical_length + parent + I(physical_length^2), data = IA_shuffling_chrom_l)
summary(LM3.2)
plot(LM3.2)

drop1(LM3.2, test = "F")
anova(LM3.2)


# wilcoxon test
IA_shuffling_chrom_wider <- IA_shuffling_chrom_l %>%
  pivot_wider(id_cols = chromosome, 
              names_from = parent, 
              values_from = mean_intra_shuffling_chr)

wilcox.test(IA_shuffling_chrom_wider$father, IA_shuffling_chrom_wider$mother, paired = TRUE)


ggplot(IA_shuffling_chrom_l, aes(physical_length, mean_intra_shuffling_chr, color = parent)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  ylim(0,NA) +
    theme(legend.position = "top", text = element_text(size = 15)) +
    theme_bw() +
        scale_color_manual(values = c( "#2ab1ca","#eb6e50")) +
    labs(x = "chromosome length [Mb]", y = "intrachromosomal shuffling") 
