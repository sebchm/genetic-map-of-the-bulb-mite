map <- read.table("~/LM_new_genome/results/10.OM2_after_LA/data.Filtering_dataTol0.05_familyInfLimit4_complementedFamLim3_missingLimit0.25_RepeatsRemovedFlankingRegions200bp_minMAF0.15_lodSC21/summary_all_maps.txt", col.names = c("contig", "position_contig", "position_genetic", "LG", "infMask"))

map <- map %>%
  group_by(infMask, LG) %>%
  mutate(marker_id = row_number(), 
         infMask = factor(infMask))

ggplot(map, aes(marker_id, position_genetic)) +
  geom_point() +
  facet_grid(infMask ~ LG, scales = "free")

### map diagnostics:
  # check if first or last markers do not inflate the map length by more than 10%
# check genetic length od each LG
LG_length <-  map %>%
  group_by(LG, infMask) %>%
  summarise(gen_length = max(position_genetic))

# check first marker
marker_inflation_map_start <- map %>%
  group_by(infMask, LG) %>% # for each LG and infMask
  slice_min(n = 2, order_by = marker_id) %>% # take 2 first/last markers
  mutate(diff = position_genetic - lag(position_genetic)) %>% # check difference of genetic position between these markers
  filter(!(is.na(diff))) %>% # remove nas
  left_join(LG_length, by = c("LG", "infMask")) %>% # add total length of each linkage group
  mutate(percent_inflate = abs(diff) * 100 /gen_length) # calculate % of inflation

# check last marker
marker_inflation_map_end <- map %>%
  group_by(infMask, LG) %>%
  slice_max(n = 2, order_by = marker_id) %>%
  mutate(diff = position_genetic - lag(position_genetic)) %>%
  filter(!(is.na(diff))) %>%
  left_join(LG_length, by = c("LG", "infMask")) %>%
  mutate(percent_inflate = abs(diff) * 100 /gen_length)

# I did not remove ANY marker, as these are not indicative of map inflation- the maximum inflation is at female map (infMask23), at the end of LG8, but it's only 1.7 cM, and edge marker is from the same contig as neighbouring markers. 

