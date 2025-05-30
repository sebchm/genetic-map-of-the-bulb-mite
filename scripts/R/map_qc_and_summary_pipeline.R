library(tidyverse) 
library(cowplot)
library(gridExtra)
library(ggpubr)

# aim: quality control, calculate recombination rates, periphery bias, permutation tests and linear models

##==============
## 1. load maps
##==============

  # infMask = 13 is the male map and infMask = 23 is the female map. LG7 is the sex chromosome. LG6 in both sexes and LG5 in males were calculated with proximityScale=100. 


map <- data.frame()

for (infMask in c("13", "23")) {
  for (chr in 1:8) {
    if (chr == 7) next
    if (chr == 6 & infMask == "23") next # this was done with ProximityScale=100
   if (chr %in% c(5,6) & infMask == "13") next # these were done with ProximityScale=100

    file_path <- paste0("(...)/results/(...)/infMask", infMask, "/(...)_infMask",infMask,"_LG",chr,"_physical.txt")

        map_i <- read.table(file_path)  

# select columns based on infMask. In male maps (infMask=13), genetic position is in the 3rd col, and for females in 4th. 
    if (infMask == "13") {
      map_i <- map_i %>% select(V1, V2, V3)
    } else if (infMask == "23") {
      map_i <- map_i %>% select(V1, V2,V4)
    }
colnames(map_i) <- c("contig", "position_contig", "position_genetic")

# modify the data
    map_i <- map_i %>%
      mutate(snp_no = row_number(), 
             LG = paste0("LG",chr), 
             infMask = infMask)

    map <- rbind(map, map_i)
  }
}
  # load LG7 (sex chromosome)
chr = 7
infMask = 2
  file_path <- paste0("(...)/results/(...)/infMask2/(...)_infMask2_LG7_physical.txt")

    map_i <- read.table(file_path)  

# select columns based on infMask
map_i <- map_i %>% select(V1, V2,V4)
colnames(map_i) <- c("contig", "position_contig", "position_genetic")

# modify the data
    map_i <- map_i %>%
      mutate(snp_no = row_number(), 
             LG = paste0("LG",chr), 
             infMask = infMask)

  map <- rbind(map, map_i)
  
  # load LG6 with infMask23 
chr = 6
infMask = 23
  file_path <- paste0("(...)/results/(...)/infMask23/(...)_infMask23_LG6_proximityScale100_physical.txt")

    map_i <- read.table(file_path)  

# celect columns based on infMask
map_i <- map_i %>% select(V1, V2,V4)
colnames(map_i) <- c("contig", "position_contig", "position_genetic")

# modify the data
    map_i <- map_i %>%
      mutate(snp_no = row_number(), 
             LG = paste0("LG",chr), 
             infMask = infMask)

  map <- rbind(map, map_i)
  
  # load LG6 and LG8  with infMask13 
for (chr in c(5,6)) {
infMask = 13
  file_path <- paste0("(...)/results/(...)/infMask13/(...)_infMask13_LG",chr,"_proximityScale100_physical.txt")

    map_i <- read.table(file_path)  

# celect columns based on infMask
map_i <- map_i %>% select(V1, V2, V3)
colnames(map_i) <- c("contig", "position_contig", "position_genetic")

# modify the data
    map_i <- map_i %>%
      mutate(snp_no = row_number(), 
             LG = paste0("LG",chr), 
             infMask = infMask)

  map <- rbind(map, map_i)
}
# change physical position to integer
map$position_contig <- sub('\\*$', '', map$position_contig)
map$position_contig <- as.numeric(map$position_contig)


##==============================================
## 2. assign markers to their physical positions
##==============================================

# read the file with columns: order (in relation to snps.txt), position on LG, LG
marker_position_on_LG <- read.table("(...)/order_allLGs_phys_liftover_named.txt", col.names = c("contig", "position_contig", "position_LG", "chromosome"))
marker_position_on_LG$position_contig <- sub('\\*$', '', marker_position_on_LG$position_contig)
marker_position_on_LG$position_contig <- as.numeric(marker_position_on_LG$position_contig)
marker_position_on_LG$chromosome <- paste0("LG", marker_position_on_LG$chromosome)

  # add SNP position on chromosome:
map <-  map %>%
  left_join(marker_position_on_LG, by = c("LG" = "chromosome", "position_contig", "contig"))

  # rename the chromosomes to match their physical size
map <- map  %>%
    mutate(LG = recode(LG,
                         "LG6" = "chr1",
                         "LG3" = "chr2",
                         "LG8" = "chr3",
                         "LG2" = "chr4",
                         "LG1" = "chr5",
                         "LG7" = "chr6",
                         "LG4" = "chr7",
                         "LG5" = "chr8")) %>%
  mutate(LG = factor(LG, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8"))) %>% 
  arrange(LG)
  # load chromosome sizes:
physical_length_chromosomes <- read.table("(...)/REF_LA.fa.fai")[1:8, 1:2]

colnames(physical_length_chromosomes) = c("LG", "physical_length")

# rename the chromosomes to match their physical size
physical_length_chromosomes <- physical_length_chromosomes  %>%
    mutate(LG = recode(LG,
                         "LG6" = "chr1",
                         "LG3" = "chr2",
                         "LG8" = "chr3",
                         "LG2" = "chr4",
                         "LG1" = "chr5",
                         "LG7" = "chr6",
                         "LG4" = "chr7",
                         "LG5" = "chr8")) %>%
  mutate(LG = factor(LG, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8"))) %>% 
  arrange(LG)


##===================================
## 3. calculate inter-marker distance
##===================================

mean_intermarker_distance_females_cM <-  map %>%
  filter(infMask != "13") %>%
      group_by(LG) %>%
  mutate(position_genetic_lagged = lag(position_genetic, default = 0), 
         intermarker_distance = position_genetic - position_genetic_lagged) %>%
    ungroup() %>%
  summarise(mean_intermarker_distance = mean(intermarker_distance)) %>%
  pull()

mean_intermarker_distance_males_cM <-  map %>%
  filter(infMask == "13") %>%
      group_by(LG) %>%
  mutate(position_genetic_lagged = lag(position_genetic, default = 0), 
         intermarker_distance = position_genetic - position_genetic_lagged) %>%
  ungroup() %>%
  summarise(mean_intermarker_distance = mean(intermarker_distance)) %>%
  pull()

paste0("Mean inter-marker distance in female map is ", mean_intermarker_distance_females_cM, " cM, and in male map ", mean_intermarker_distance_males_cM, " cM")

mean_intermarker_distance_females_bp <-  map %>%
  filter(infMask == "23" | infMask == "2") %>%
    group_by(LG) %>%
  mutate(position_physical_lagged = lag(position_LG, default = 0), 
         intermarker_distance = position_LG - position_physical_lagged) %>%
    ungroup() %>%
  summarise(mean_intermarker_distance = mean(intermarker_distance)) %>%
  pull()

mean_intermarker_distance_males_bp <-  map %>%
  filter(infMask == "13") %>%
  group_by(LG) %>%
  mutate(position_physical_lagged = lag(position_LG, default = 0), 
         intermarker_distance = position_LG - position_physical_lagged) %>%
    ungroup() %>%
  summarise(mean_intermarker_distance = mean(intermarker_distance)) %>%
  pull()

paste0("Mean inter-marker distance in female map is ", mean_intermarker_distance_females_bp, " bp, and in male map ", mean_intermarker_distance_males_bp, " bp")


##=================
## 4. count markers
##=================

# count markers on autosomes
map %>%
  filter(LG != "chr6") %>%
  group_by(infMask) %>%
  tally(name = "number_of_markers")

# count markers on the sex chromosome
n_markers_sex_chromosome <- map %>%
  filter(LG == "chr6") %>%
  tally()

# check genetic length of the sex chromosome
genetic_length_sex_chromosome <- map %>%
  filter(LG == "chr6") %>%
  summarise(max(position_genetic))

# calculate recombination rate
  # get size of the sex chromosome
physical_length_sex_chromosome <- subset(physical_length_chromosomes, chr == "chr6")[,2]/1e6
  # calculate recombination rate
rec_rate_sex_chromosome <- round(genetic_length_sex_chromosome / physical_length_sex_chromosome,2)

paste0("The sex chromosome, which is hemizygous in males, contained ",n_markers_sex_chromosome, " markers across ", genetic_length_sex_chromosome, " cM.")

##=======================
## 5. prepare map summary
##=======================

map_summary_per_LG <- map %>%
  mutate(sex = ifelse(infMask == 13, "male", "female")) %>%
  group_by(sex, LG) %>%
  summarise(genetic_length = round(max(position_genetic),1), 
            number_of_markers = n()) %>%
  left_join(physical_length_chromosomes, by = "LG") %>%
  mutate(chromosome_size = physical_length / 1e6, 
         rec_rate = round(genetic_length/chromosome_size, 1)) %>%
  pivot_wider(id_cols = c(LG), 
              names_from = sex,
              values_from = c(genetic_length, number_of_markers, rec_rate)) %>%
  mutate(F_M_length_ratio = round(genetic_length_female/genetic_length_male,1))

general_map_summary <- map %>%
  mutate(sex = ifelse(infMask == "13", "male", "female")) %>%
  group_by(sex, LG) %>%
  select(-infMask) %>%
  summarise(genetic_length = max(position_genetic), 
            n_markers = n()) %>%
  left_join(physical_length_chromosomes, by = "LG") %>%
  mutate(physical_length = physical_length/1e6,
         rec_rate = genetic_length/physical_length) %>%
  pivot_wider(id_cols = "LG", names_from = "sex", values_from = c("genetic_length", "n_markers", "rec_rate", "physical_length")) %>%
  select(-physical_length_male) %>%
  mutate(F_M_length_ratio = genetic_length_female/genetic_length_male) %>%
  rename("physical_length" = "physical_length_female") %>%
  mutate(across(where(is.numeric), ~ round(., 2)))


map_length_male_female <- map %>%
  group_by(LG, infMask) %>%
  filter(LG != "chr6") %>%
  summarise(map_length = max(position_genetic)) %>%
  group_by(infMask) %>%
  summarise(sum_length_map = sum(map_length))

female_map_length <- map_length_male_female[2,2] %>% pull()
male_map_length <- map_length_male_female[1,2] %>% pull()

# calculate % longer is the female map than the male map
female_map_longer <- round((female_map_length - male_map_length) * 100 / male_map_length,0)

paste0("The female autosomal map spanned ", round(female_map_length,2)," cM, being ", female_map_longer , "% longer than the male autosomal map, which had ", male_map_length, " cM")

write.table(general_map_summary, "general_map_summary.txt", quote = F, row.names = FALSE)

### genome-wide sex-specific recombination rate:
map_physical_length <- sum(physical_length_chromosomes$physical_length[physical_length_chromosomes$LG != 6])/1e6

male_GW_RR <- sum(general_map_summary$genetic_length_male, na.rm = T)/map_physical_length
female_GW_RR <- sum(general_map_summary$genetic_length_female[general_map_summary$LG != "chr6"])/map_physical_length

female_map_RR_higher <- round((female_GW_RR - male_GW_RR) * 100 / male_GW_RR,0)


paste0("Genome-wide recombination rate is ", round(female_GW_RR,2)," and ", round(male_GW_RR,2), " cM/Mb in female and male map, respectively, making female recombination rate ", round(female_map_RR_higher), "% higher.")

##============================================
## 6. count crossovers in male and female maps
##============================================

map %>%
  #filter(infMask == 23) %>%
  filter(infMask == 23 | infMask == 2) %>%
  distinct(position_genetic) %>%
  tally()

map %>%
  filter(infMask == 13) %>%
  distinct(position_genetic) %>%
  tally()

##===========================
## 7. analyse genome assembly
##===========================

agp <- read.table("(...)/REF_LA.agp", col.names = c("object", "object_beg", "object_end", "part_number", "component_type", "component_id", "component_beg", "component_end", "orientation"))

# remove gaps from the agp:
agp <- agp %>%
  filter(component_type != "U") %>%
  mutate(component_beg = as.integer(as.character(component_beg)), 
         component_end = as.integer(as.character(component_end)),
         component_length = component_end - component_beg + 1)

agp_chromosomes <- agp %>%
  filter(str_detect(object, "LG"))

paste(nrow(agp_chromosomes), " contigs have been incorporated to the new assembly")

# length of contigs in chromosomes
length_in_chromosomes <- agp %>%
  filter(str_detect(object, "LG")) %>%
  summarise(sum_contig_length = sum(component_length))

paste(round(length_in_chromosomes/1e6,3), " Mb of the contigs have been incorporated to chromosomes")

# total length of the summary
all_contig_sizes <- read.table("(...)/REF_LA.fa.fai")[,1:2]
new_assembly_length <- sum(all_contig_sizes$V2)

paste0("New assembly has ", round(new_assembly_length/1e6, 2), " Mb")

paste0(round(length_in_chromosomes*100/new_assembly_length,1), "% of the assembly length has been incorporated to the map.")
### contigs not used 
not_used_contigs <- read.table("(...)/not_used.agp")[,c(1,3)]
colnames(not_used_contigs) <- c("contig", "contig_length")
not_used_length <- sum(not_used_contigs$contig_length)
paste0("Contigs not used in assembly have ", round(not_used_length/1e6, 2), " Mb")

# rename the chromosomes to match their physical size
agp <- agp %>%
  mutate(object = if_else(str_detect(object, "^sq"), 
                          object, 
                          recode(object,
                         "LG6" = "chr1",
                         "LG3" = "chr2",
                         "LG8" = "chr3",
                         "LG2" = "chr4",
                         "LG1" = "chr5",
                         "LG7" = "chr6",
                         "LG4" = "chr7",
                         "LG5" = "chr8")),
         object = factor(object)) %>%
  mutate(object = forcats::fct_relevel(object, c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8"))) %>%
  arrange(object)

# check number of contigs included in the map
n_contigs_in_map <- agp %>%
  distinct(component_id) %>%
  tally() %>%
  pull()

##==============================================
## 8. check for the presence of chimeric contigs
##==============================================
  # none --> great!
map %>% group_by(LG) %>%
  distinct(contig) %>%
  group_by(contig) %>%
  tally(name = "n_chromosomes") %>%
  arrange(desc(n_chromosomes)) %>%
  filter(n_chromosomes != 1)


##============================
## 9.plot figure 1A: Marey map
##============================

# MDR - morph-determining region
# based on the genome alignment add morph-determining region (MDR)
MDR_start_pos <- 39620728
MDR_end_pos <- 42000261

MDR <- data.frame(LG = "chr2",  xmin = MDR_start_pos / 1e6, xmax = MDR_end_pos / 1e6,
  ymin = 0, ymax = 50)

fig_1A <- map %>%
  mutate(infMask = ifelse(infMask == "13", "male map", "female map"), 
         LG = as.character(LG),
         LG = ifelse(LG == "chr6", "chr6 (sex)", LG)) %>%
ggplot(aes(position_LG/1e6, position_genetic, color = infMask)) +
    geom_rect(data = MDR, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightgrey", alpha = 0.8) +
  geom_point(size = .65, alpha = .5) +
  scale_color_manual(values = c("#2ab1ca", "#eb6e50")) +
  facet_wrap(~LG, ncol = 1, strip.position="left") +
  theme_bw() +
  theme(plot.tag = element_text(size = 15), 
        plot.tag.position = c(0.01, 1),
        plot.background = element_rect(fill = "transparent"),
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = "transparent"), 
        strip.placement = "outside", 
        legend.position = "none") +
  labs(x = "physical position [Mb]", y = "genetic position [cM]", tag = "A")

plot(fig_1A)


##============================================================================
## 10.plot figure 1B: physical length of chromosome VS mean recombination rate
##============================================================================

fig_1B <- recombination_rate_longer %>%
    mutate(LG = gsub("chr", "", LG)) %>%
ggplot(aes(physical_length, rec_rate, color = sex, label = factor(LG))) +
  geom_smooth(method = "lm", alpha = .1, linewidth = .5) +
  geom_point(alpha = .8, size = 2.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#2ab1ca", "#eb6e50")) +
    labs(tag = "B", x = "chromosome size [Mb]", y = "rec. rate\n[cM/Mb]") +
     theme(plot.tag.position = c(0.01, .95), 
           plot.tag = element_text(size = 15), 
        #panel.grid.minor = element_blank(), 
        panel.grid = element_blank()) +
  annotate("text", x = 49, y = 1.2, label = "P == 0.1", parse = TRUE, size = 5) +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) 

plot(fig_1B)

# LM: is recombination rate associated with physical length?
  # full model
LM1.1 <- lm(rec_rate ~ physical_length * sex, data = recombination_rate_longer)
summary(LM1.1)
drop1(LM1.1, test = "F") # interaction is not significant- drop it

LM1 <- lm(rec_rate ~ physical_length + sex, data = recombination_rate_longer)
summary(LM1)
# length is not significant. 
drop1(LM1, test = "F")
anova(LM1)
summary(LM1)

  # exclude chromosome 1 -- no change
LM1.2 <- lm(rec_rate ~ physical_length * sex, data = recombination_rate_longer[recombination_rate_longer$LG != "chr1",])
drop1(LM1.2, test = "F") # drop interaction
  # exclude chromosome 1 
LM1.3 <- lm(rec_rate ~ physical_length + sex, data = recombination_rate_longer[recombination_rate_longer$LG != "chr1",])
drop1(LM1.3, test = "F") # length is not significant


##=================================================
## 11.plot figure 1C: chromosome size VS map length
##=================================================

fig_1C <- map %>%
  mutate(sex = ifelse(infMask == "13", "male", "female")) %>% 
  group_by(sex, LG) %>%
  summarise(genetic_length = max(position_genetic)) %>%
  left_join(physical_length_chromosomes, by = "LG") %>%
  mutate(physical_length = physical_length/1e6, 
         sex = ifelse(str_detect(sex, "male"), paste0(sex, " map"), sex)) %>%
  mutate(LG = gsub("chr", "", LG)) %>%
  ggplot(aes(y = genetic_length, x = physical_length, color = sex)) +
  geom_smooth(method = "lm", alpha = .1, linewidth = .5) +
  geom_point(alpha = .8, size = 2.5) +
  theme_classic() +
  labs(x = "chromosome size [Mb]", y = "genetic length \n[cM]" , color=NULL, text = NULL) +
  scale_color_manual(values = c("#2ab1ca", "#eb6e50" )) +
  scale_fill_manual(values = c("#2ab1ca", "#eb6e50"))+
  labs(tag = "C") +
     theme(plot.tag.position = c(0.01, .95), 
           plot.tag = element_text(size = 15), 
           legend.box = "horizontal", 
          legend.direction = "horizontal",  
           legend.position = "top",
           legend.title = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate("text", x = 50, y = 35, label = "P == 0.3", parse = TRUE, size = 5) +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) 

plot(fig_1C)

### is map length of associated with the physical size of the chromosome? 
  # prepare data
phys_gen_length_longer <- general_map_summary %>%
  select(c(LG, physical_length, genetic_length_female, genetic_length_male)) %>%
  pivot_longer(cols = ends_with("male"), 
               names_to = "sex", 
               values_to = "genetic_length", 
               names_prefix = "genetic_length_")

# full model
LM2.1 <- lm(genetic_length ~ physical_length * sex, data = phys_gen_length_longer)
summary(LM2.1)
drop1(LM2.1, test = "F") # drop interaction

# model with no interaction - no differences
LM2.2 <- lm(genetic_length ~ physical_length + sex, data = phys_gen_length_longer)
summary(LM2.2)
drop1(LM2.2, test = "F")

##=============================================
## 12.plot figure 1D: male VS female map length
##=============================================

fig_1D <- map %>%
  mutate(sex = ifelse(infMask == "13", "male", "female")) %>%
  group_by(sex, LG) %>%
  summarise(genetic_length = max(position_genetic)) %>%
  pivot_wider(id_cols = LG, 
              names_from = sex, 
              values_from = genetic_length) %>%
  filter(LG != 6) %>%
  mutate(LG = gsub("chr", "", LG)) %>%
  ggplot(aes(male, female, label = factor(LG))) +
  geom_smooth(method = "lm", se = F, color = "black", linewidth = 0.5, fullrange=TRUE) +
  geom_point(alpha = .8) +
  theme_classic() +
  ggrepel::geom_text_repel(size = 5, show.legend = FALSE) +
  #lims(x = c(0,30), y = c(0,50)) +
  geom_abline(linetype = "dashed") +
  labs(tag = "D", x = "male map length [cM]", y = "female map length\n[cM]") +
     theme(plot.tag = element_text(size = 15),
           plot.tag.position = c(0.01, .95), 
        panel.grid.minor = element_blank()) +
  annotate("text", x = 26, y = 13, label = "P == 0.14", parse = TRUE, size = 5)

### LM: male VS female map length
LM3 <- lm(genetic_length_female ~ genetic_length_male, data = subset(general_map_summary, LG != "chr6"))
summary(LM3)
drop1(LM3, test = "F")


##==================================================
## 13.combine panels for figure 1 into a single plot
##==================================================

# get a common legend for fig. 1
common_legend <- get_legend(fig_1C)
fig_1C <- fig_1C +
    theme(legend.position = "none") 

# arrange plots 1B-1D 
right_column_plots <- grid.arrange(
  fig_1B,
  fig_1C,
  fig_1D,
  ncol = 1, 
  heights = c(1, 1, 1), 
  layout_matrix = rbind(1, 2, 3))

# arrange the marey_plot with the common legend
left_plot_with_legend <- arrangeGrob(
  common_legend,
  fig_1A,
  ncol = 1,
  heights = c(.3, 14))

# arrange all fig. 1 plots together
fig1 <- grid.arrange(
  left_plot_with_legend,
  right_column_plots,
  ncol = 2,
  widths = c(3, 2))

##===============================================
## 14.save files which will be uploaded to GitHub
##===============================================

# write.table(map[,c(5,7,3,6)], "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/map_analysis/filesToUploadGitHub/Rhrob_map.txt", quote = F, row.names = F, sep = "\t")
# write.table(physical_length_chromosomes, "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21//map_analysis/filesToUploadGitHub/chromosome_lengths.txt", quote = F, row.names = F, sep = "\t")
# write.table(agp, "~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21//map_analysis/filesToUploadGitHub/Rhrob.agp", quote = F, row.names = F, sep = "\t")


##===============================================
## 15. find regions with highest difference in RR
##===============================================
  # morph-determining region shows striking difference in recombination rate between sexes- male recombination is effecively suppressed. I want to formally test it. 

# load 2.38 Mb windows:
mb2_windows <- read.table("~/LM_new_genome/data/genome/chromosome_scale/2.38Mb_windows/REF_LA_chrNames_sorted_2.38Mb_100kbOverlap_windows.txt", col.names = c("LG", "window_start", "window_end"))

mb2_windows <- mb2_windows %>%
  filter(chr != "chr6") %>% # remove sex chromosome
  filter(str_detect(LG, "chr")) %>% # remove contigs not included in the map
  mutate(window_index = row_number()) 

  # add variable: infMask
mb2_windows <- expand_grid(mb2_windows, infMask = c("13", "23"))
mb2_windows$infMask <- as.factor(mb2_windows$infMask)

map_RR_diff <- map %>%
  mutate(infMask = as.factor(infMask)) %>%
  filter(LG != "chr6") %>% # remove sex chromosome 
  left_join(mb2_windows, by = c("LG", "infMask"), relationship = "many-to-many") %>% # add windows
  filter(position_LG >= window_start & position_LG < window_end)

RR_diff <- tibble()

# calculate difference between female and male recombination rate for each window
for (i in unique(map_RR_diff$window_index)){
print(i)

RR_diff_temp <- map_RR_diff %>%
  filter(window_index == i) %>% 
  group_by(infMask) %>%
  summarise(RR = (max(position_genetic) - min(position_genetic)) / # genetic length
              ((max(position_LG) - min(position_LG)) / 1e6), # ...divided by physical length
            min_pos_physical = min(position_LG)/1e6, 
            max_pos_physical = max(position_LG)/1e6, 
            window_size = max_pos_physical - min_pos_physical, 
            chr = unique(LG)) %>% 
  pivot_wider(names_from = infMask, values_from = RR, names_prefix = "RR_") %>%
  mutate(RR_rate = RR_23 - RR_13) %>%
  mutate(window_index = i)
RR_diff <- bind_rows(RR_diff, RR_diff_temp)
}

RR_diff <- RR_diff %>%
  left_join(mb2_windows[,-5], by = "window_index")

colSums(is.na(RR_diff))
RR_diff <- na.omit(RR_diff)

hist(RR_diff$RR_rate)

MDR <- data.frame(LG = "chr2",  xmin = 39620728 / 1e6, xmax = 42000261 / 1e6,
  ymin = 0, ymax = 5)

RR_diff_min2Mb <- RR_diff %>%
  filter(window_size > 2)

# assign windows to MDR/non-MDR
RR_diff_min2Mb <- RR_diff_min2Mb %>%
  mutate(is_MDR = ifelse(chr == "chr2" & window_start >= MDR_start_pos & window_start <= MDR_end_pos | 
                            chr == "chr2" & window_end >= MDR_start_pos & window_end <= MDR_end_pos, TRUE, FALSE))


RR_diff_min2Mb %>%
ggplot(aes(xmin = min_pos_physical, xmax = max_pos_physical, y = RR_rate, color = is_MDR)) +
  facet_wrap(~LG) +
  geom_rect(data = MDR, inherit.aes = FALSE, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "steelblue", alpha = 0.3) +
    geom_errorbar() +
  ylim(-10, 10)


cutoff <- quantile(abs(RR_diff_min2Mb$RR_rate), 0.95, na.rm = TRUE)

### permutation test:

# check if MDR windows have been assigned properly:
table(RR_diff_min2Mb$is_MDR)
t <- RR_diff_min2Mb %>% filter(is_MDR == TRUE)

# calculate observed difference in means
obs_diff_RR <- mean(RR_diff_min2Mb$RR_rate[RR_diff_min2Mb$is_MDR == TRUE]) - mean(RR_diff_min2Mb$RR_rate[RR_diff_min2Mb$is_MDR == FALSE])

# vector for differences in means
perm_diffs_RR <- vector("numeric", 1000)

for(i in 1:1000) {
  # shuffle MDR
  perm_groups <- sample(RR_diff_min2Mb$is_MDR)
  
  # calculate difference in means
  perm_diffs_RR[i] <- mean(RR_diff_min2Mb$RR_rate[perm_groups == TRUE], na.rm = TRUE) - mean(RR_diff_min2Mb$RR_rate[perm_groups == FALSE], na.rm = TRUE)
}
# calculate p-value (proportion of perm_diffs >= obs_diff). Add 1 to nominator to avoid P-value of 1 (Phipson and Smyth 2010)
p_value <- sum(perm_diffs_RR >= obs_diff_RR) / length(perm_diffs_RR)

print(p_value)


# calculate y and x corrdinates for the p-value:
x_mid <- mean(range(perm_diffs_RR))
hist_data <- hist(perm_diffs_RR, breaks = 100, plot = FALSE)

plot_perm_test <- ggplot() +
  geom_histogram(aes(perm_diffs_RR), bins = 100) +
  geom_vline(xintercept = obs_diff_RR, linetype = "dashed") +
  #ggtitle("GC content") +
  labs(y = "", x="null distribution") +
  theme_bw() +
  annotate("text", x = x_mid, y = 1000, size = 5, label = paste("P =", round(p_value, 3))) +
  # labs(tag = "H", fontface = "bold") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.05, 0.95), 
          plot.tag = element_text(size = 20, face="bold" ))

##======================================================================
## 16. calculate sex-specific recombination rates in 1Mb and 2Mb windows
##======================================================================
### 1 Mb
mb_windows <- read.table("(...)/REF_LA_chrNames_sorted_1Mb_windows.txt", col.names = c("LG", "window_start", "window_end"))

# remove contigs not incorporated to chromosomes
mb_windows <- mb_windows %>%
  filter(str_detect(LG, "chr"))

mb_windows <- expand_grid(mb_windows, infMask = c("13", "23"))

# chromosome 6 is a sex chromosome- only infMask2
mb_windows <- mb_windows %>%
  filter(!(LG == "chr6" & infMask == "13")) %>%
  mutate(infMask = ifelse(LG == "chr6", 2, infMask))

### 2 Mb
mb2_windows <- read.table("(...)/REF_LA_chrNames_sorted_2Mb_windows.txt", col.names = c("LG", "window_start", "window_end"))

# remove contigs
mb2_windows <- mb2_windows %>%
  filter(str_detect(LG, "chr"))

mb2_windows <- expand_grid(mb2_windows, infMask = c("13", "23"))

# chromosome 6 is a sex chromosome- only infMask2
mb2_windows <- mb2_windows %>%
  filter(!(LG == "chr6" & infMask == "13")) %>%
  mutate(infMask = ifelse(LG == "chr6", 2, infMask))

  # recombination rate in 1Mb windows

rec_rate_based_on_map <- map %>% 
  full_join(mb_windows, by = c("LG", "infMask"), relationship = "many-to-many") %>%
  filter(position_LG > window_start & position_LG <= window_end) %>%
  group_by(LG, infMask, window_start, window_end) %>%
  mutate(genetic_dist = max(position_genetic) - min(position_genetic), 
         rec_rate = genetic_dist/((window_end - window_start)/1e6)) %>%
  group_by(LG, infMask, window_start, window_end) %>%
  distinct(rec_rate) %>%
  full_join(mb_windows, by = c("LG", "infMask", "window_start", "window_end")) %>%
  mutate(infMask = ifelse(infMask == "13", "male map", "female map"), 
         rec_rate = replace_na(rec_rate, 0)) 

rec_rate_based_on_map_plot <- rec_rate_based_on_map %>%
  ggplot(aes(window_start/1e6, rec_rate, color = infMask)) +
  geom_smooth(span = 0.3, se = F, alpha = .8) +
  facet_wrap(~LG, nrow = 2) +
  ylim(0,NA) +
    theme_minimal() +
  theme(legend.position = "top", 
        plot.tag = element_text(size = 17, face="bold")) +
    labs(y = "recombination rate\n[cM/Mb]", x = "window position [Mb]", color=NULL) +
  scale_color_manual(values = c("#eb6e50", "#2ab1ca")) 

plot(rec_rate_based_on_map_plot)


  # recombination rate in 2Mb windows

rec_rate_based_on_map_2mb <- map %>% 
  full_join(mb2_windows, by = c("LG", "infMask"), relationship = "many-to-many") %>%
  filter(position_LG > window_start & position_LG <= window_end) %>%
  group_by(LG, infMask, window_start, window_end) %>%
  mutate(genetic_dist = max(position_genetic) - min(position_genetic), 
         rec_rate = genetic_dist/((window_end - window_start)/1e6)) %>%
  group_by(LG, infMask, window_start, window_end) %>%
  distinct(rec_rate) %>%
  full_join(mb2_windows, by = c("LG", "infMask", "window_start", "window_end")) %>%
  mutate(infMask = ifelse(infMask == "13", "male map", "female map"), 
         rec_rate = replace_na(rec_rate, 0), 
         infMask = ifelse(infMask == "male map", "male_map", 
                          ifelse(infMask == "female map", "female_map", "error")))

rec_rate_based_on_map_plot_2mb <- rec_rate_based_on_map_2mb %>%
  ggplot(aes(window_start/1e6, rec_rate, color = infMask)) +
  geom_smooth(span = 0.3, se = F, alpha = .8) +
  facet_wrap(~LG, nrow = 2) +
  ylim(0,NA) +
    theme_minimal() +
  theme(legend.position = "top", 
        plot.tag = element_text(size = 17, face="bold")) +
    labs(y = "recombination rate\n[cM/Mb]", x = "window position [Mb]", color=NULL) +
  scale_color_manual(values = c("#eb6e50", "#2ab1ca")) 

RR <- rec_rate_based_on_map

rec_rate_based_on_map <- rec_rate_based_on_map %>%
  left_join(physical_length_chromosomes, by = "LG") %>%
    mutate(relative_position = (window_start + window_end)/2/physical_length) # calculate relative window position

##===============================================
## 17. prepare Fig 2: sex-specific landscape plot
##===============================================
  # NB: no sex chromosome!
  ### scaling by the autosomal-wide recombination rate (rec_rate/female_GW_RR + rec_rate/male_GW_RR)
rec_landscape_plot <- rec_rate_based_on_map %>%
  mutate(rec_rate_scaled = ifelse(infMask == "female map", rec_rate/female_GW_RR, rec_rate/male_GW_RR)) %>%
  filter(LG != "chr6") %>%
  ggplot(aes(x = relative_position, y = rec_rate_scaled, color = infMask)) +
  geom_smooth(alpha = .2, span = 0.75) +
  scale_y_continuous(limits = c(0, NA), oob = scales::squish) +
  labs(tag = "A", x = "relative position", y = "scaled recombination rate [cM/Mb]") +
    theme_bw() +
  scale_color_manual(values = c("#2ab1ca", "#eb6e50")) +
  theme(legend.box = "horizontal", 
        legend.position = "top",
        legend.title = element_blank(),
        text = element_text(size = 10),
        plot.tag.position = c(0.02, 1.05), 
        panel.grid.minor = element_blank()) 


# prepare a legend
legend_heterochiasmy <- get_legend(rec_landscape_plot) 

# re-make the plot
rec_landscape_plot <- rec_landscape_plot +
  theme(legend.position = "none") 

# load plot with intrachromosomal shuffling
IA_shuffling_plot <- readRDS(".../IA_shuffling_plot.rds")
IA_shuffling_grob <- ggplotGrob(IA_shuffling_plot)

Fig2_nolegend <- arrangeGrob(rec_landscape_plot, IA_shuffling_plot, ncol = 2)
Fig2 <- arrangeGrob(legend_heterochiasmy, Fig2_nolegend, ncol = 1, nrow = 2, heights = c(2, 10))

plot(Fig2)

rec_rate_based_on_map %>%
  filter(LG != 6) %>%
ggplot(aes(relative_position, rec_rate, color = LG)) +
  geom_smooth(se = F) + 
  facet_wrap(~infMask) +
  ylim(min = 0, NA)

rec_landscape_plot_sepChr <- rec_rate_based_on_map %>%
  mutate(LG = ifelse(LG == "chr6", "chr 6 (sex chr)", LG)) %>%
  mutate(rec_rate_scaled = ifelse(infMask == "female map", rec_rate/female_GW_RR, rec_rate/male_GW_RR)) %>% # scale RR by RR
  ggplot(aes(relative_position, rec_rate_scaled, color = infMask)) +
  geom_smooth(se = F) + 
  facet_wrap(~LG, ncol = 4, nrow = 2) +
  ylim(min = 0, NA) +
  scale_color_manual(values = c("#2ab1ca", "#eb6e50"), name = "") +
        theme_bw() +
  theme(legend.position = "top", 
        text = element_text(size = 15)) +
    labs(x = "Relative position", y = "Scaled ecombination rate (cM/Mb)") 

rec_landscape_plot_sepChr

##============================
## 18. periphery bias analysis
##============================

# calculate distance to chromosome end and periphery bias
  # periphery bias: 10% at each extremity of each chromosome divided by the mean recombination rate

PB=0.1
periphery_bias <- rec_rate_based_on_map %>%
mutate(window_middle = (window_end + window_start) / 2, 
         chromosome_middle = physical_length / 2, 
         window_closer_to = ifelse(window_middle < chromosome_middle, 1, physical_length),
         distance_to_chromosome_end = abs(window_middle - window_closer_to + 1)) %>%
  group_by(LG) %>%
  mutate(is_periphery = ifelse((window_start+window_end)/2 < PB * physical_length, TRUE, ifelse((window_start+window_end)/2 > (1-PB) * physical_length, TRUE, FALSE))) %>%
  group_by(LG, infMask, is_periphery) %>%
  summarise(mean_RR_periphery = mean(rec_rate)) %>%
  filter(is_periphery == TRUE) 

mean_rec_rate_chromosome <- rec_rate_based_on_map %>%
  group_by(LG, infMask) %>%
  summarise(mean_chr_rec_rate = mean(rec_rate))

# add mean chromosome RR and chromosome size
  periphery_bias <- periphery_bias %>%
    left_join(mean_rec_rate_chromosome, by = c("LG", "infMask")) %>%
    mutate(periphery_bias = round(mean_RR_periphery/mean_chr_rec_rate,2)) %>%
    select(-is_periphery) %>%
    left_join(physical_length_chromosomes, by = "LG")

periphery_bias %>%
  group_by(infMask) %>%
  summarise(mean_PB = round(mean(periphery_bias),2))

general_map_summary <- periphery_bias %>%
  mutate(LG = paste0("chr", LG)) %>%
  select(LG, periphery_bias, infMask) %>%
  pivot_wider(id_cols = LG, 
              values_from = periphery_bias, 
              names_from = infMask, 
              names_prefix = "periphery_bias_") %>%
  right_join(general_map_summary, by = c("LG"))


ggplot(periphery_bias, aes(LG, periphery_bias, fill = infMask)) +
  geom_col(position = "dodge") +
  ggtitle("Periphery bias")


### Is the difference in periphery bias between the sexes significant?
  # check distribution
hist(periphery_bias$periphery_bias[periphery_bias$infMask == "female map"])
hist(periphery_bias$periphery_bias[periphery_bias$infMask == "male map"])

#### paired test -> N.S.
periphery_bias_for_test <- periphery_bias %>%
  filter(LG != "chr6") %>%
  select(-c(mean_RR_periphery, mean_chr_rec_rate)) %>%
  mutate(infMask = str_replace_all(infMask, " ", "_")) %>%
  pivot_wider(names_from = infMask, values_from = periphery_bias)

# Run paired Wilcoxon test
wilcox.test(periphery_bias_for_test$female_map, periphery_bias_for_test$male_map, paired = TRUE)

#### permutation test ---> N.S
obs_diff_PB <- abs(mean(periphery_bias$periphery_bias[periphery_bias$infMask == "female map"]) - mean(periphery_bias$periphery_bias[periphery_bias$infMask == "male map"]))

# unitialize a vector 
perm_diffs_PB <- vector("numeric", 1000)

for(i in 1:1000) {
  # sample group labels
  perm_groups <- sample(periphery_bias$infMask)
  
  # calculate mean permuted difference
  perm_diffs_PB[i] <- abs(mean(periphery_bias$periphery_bias[perm_groups == "female map"]) - mean(periphery_bias$periphery_bias[perm_groups == "male map"]))
}

# calculate p-value (proportion of perm_diffs >= obs_diff)
p_value_PB <- sum(perm_diffs_PB >= obs_diff_PB) / length(perm_diffs_PB)

print(p_value_PB)

ggplot() +
  geom_histogram(aes(perm_diffs_PB), bins = 100) +
  geom_vline(xintercept = obs_diff_PB, linetype = "dashed") +
  ggtitle("periphery bias") +
  labs(y = "count", x="Difference of means") +
  theme_bw() +
  annotate("text", x = 0.15, y = 390, size = 5, label = paste("p-value =", round(p_value_PB, 3))) +
   labs(tag = "E", fontface = "bold") +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag.position = c(0.05, 0.95), 
          plot.tag = element_text(size = 20, face="bold" ))

LM_PB_full <- lm(periphery_bias ~ physical_length * infMask, data = periphery_bias)
drop1(LM_PB_full, test = "F") # drop interaction

LM_PB_1 <- lm(periphery_bias ~ physical_length + infMask, data = periphery_bias)
summary(LM_PB_1)
drop1(LM_PB_1, test = "F")

##===========================
## 18. analyse repeat classes
##===========================
repeat_classes <- read.table("~/assembly/23.EarlGrey_TE_annotation/IW24_minLen5kb_minQ10_contaminantsRemoved_polished_chromosomeScale/rhizoglyphusRobini_EarlGrey/rhizoglyphusRobini_summaryFiles/rhizoglyphusRobini.familyLevelCount.txt", h = T,  fill = T, comment.char = "")

repeat_classes <- repeat_classes %>%
  separate(name, into = c("family_name", "family"), sep = "#")

repeat_classes_summarised <- repeat_classes %>%
  group_by(family) %>%
  summarise(copy_number = sum(copy_number), 
            coverage = sum(coverage)) %>% 
  separate(family, into = c("family", "subfamily"), sep = "/")
