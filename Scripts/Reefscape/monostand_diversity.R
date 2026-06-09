library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
### Monostand analysis
kal <- read_delim("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/Monostand annotations Kalki 23 rename.txt", 
                  delim = ",", escape_double = FALSE, 
                  trim_ws = TRUE)
sna <- read_delim("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/Monostand annotations Snakebay 23 rename.txt", 
                  delim = ",", escape_double = FALSE, 
                  trim_ws = TRUE)
sea <- read_delim("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/Monostand annotations Seaquarium 21 rename.txt", 
                  delim = ",", escape_double = FALSE, 
                  trim_ws = TRUE)
sna$Shape <- paste0(sna$Shape, "_sna")
kal$Shape <- paste0(kal$Shape, "_kal")
sea$Shape <- paste0(sea$Shape, "_sea")
all_shapes <- bind_rows(sna, kal, sea)
all_shapes <- all_shapes[-2]
colnames(all_shapes) <- c("shape_id", "Samples")
all_shapes <- all_shapes %>%
  mutate(Size = lengths(strsplit(Samples, " ")))


monostand_df <- read.table(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/Monostand_metrics_combined.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

monostand_df$Site <- factor(monostand_df$Site)
monostand_df$Environment <- factor(monostand_df$Environment)
rugosity_lookup <- tibble(
  Site = c("Playa Kalki","Playa Kalki","Playa Kalki",
           "Snakebay","Snakebay","Snakebay",
           "Sea Aquarium","Sea Aquarium"),
  Environment = c("Slope","Crest","Terrace",
                  "Slope","Crest","Terrace",
                  "Slope","Crest"),
  rugosity = c(1.785,1.784,1.401,
               1.703,1.479,1.282,
               1.791,1.718)
)

monostand_df <- monostand_df %>%
  left_join(rugosity_lookup, by = c("Site", "Environment"))
monostand_df <- monostand_df %>%
  group_by(across(-c(Environment, rugosity))) %>%
  summarise(
    Environment = paste(sort(unique(Environment)), collapse = ". "),
    rugosity = mean(rugosity, na.rm = TRUE),
    .groups = "drop"
  )

monostands <- full_join(all_shapes, monostand_df, by = "shape_id")

#counting stand sizes
sum(monostands$area_2d > 3.14)
sum(monostands$area_2d > 10)
sum(monostands$area_2d <= 1)
sum(monostands$area_2d <= 3.14 & monostands$area_2d > 1)

monostands_large <- monostands %>%
  dplyr::filter(area_2d > 3.14)

monostands_large %>%
  summarise(
    mean_perimeter = mean(perimeter, na.rm = TRUE),
    sd_perimeter   = sd(perimeter, na.rm = TRUE),
    mean_area      = mean(area_2d, na.rm = TRUE),
    sd_area        = sd(area_2d, na.rm = TRUE),
    n_stands       = n()
  )

monostands_large %>%
  group_by(Site) %>%
  summarise(
    mean_perimeter = mean(perimeter, na.rm = TRUE),
    sd_perimeter   = sd(perimeter, na.rm = TRUE),
    mean_area      = mean(area_2d, na.rm = TRUE),
    sd_area        = sd(area_2d, na.rm = TRUE),
    n_stands       = n()
  )

kruskal.test(area_2d ~ Site, data = monostands_large)
pairwise.wilcox.test(monostands_large$area_2d,
                     monostands_large$Site,
                     p.adjust.method = "BH")
kruskal.test(perimeter ~ Site, data = monostands_large)
pairwise.wilcox.test(monostands_large$perimeter,
                     monostands_large$Site,
                     p.adjust.method = "BH")

ggplot(monostands_large, aes(x = Site, y = area_2d, label = Site)) +
  geom_boxplot() +
  geom_text(nudge_y = 1) +
  theme_classic() +
  labs(
    x = "Site",
    y = "Monospecific stand area (m²)"
  )


metadata <- read_table("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/MADR_reference_popfile_G1_full_clonal_groups.txt")
annotations <- read_csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/combined_3d_annotations_2123_aligned_renamed.coords.csv")

# Extract the ID part before "_MMIR"
metadata$ID <- sub("_SHA", "", metadata$INDV)

metadata <- merge(
  metadata,
  annotations[, c("ID", "world_z")],
  by = "ID",
  all.x = TRUE
)
monostands_long <- monostands %>%
  separate_rows(Samples, sep = " \\s*") %>% 
  rename(ID = Samples)
metadata2 <- metadata %>%
  left_join(monostands_long, by = "ID")

stand_df <- metadata2 %>%
  group_by(shape_id) %>%
  summarise(
    area_2d = first(area_2d),
    perimeter = first(perimeter),
    Depth_min = min(world_z),
    Depth_max = max(world_z),
    Depth_mean = mean(world_z),
    Site = first(Site),
    Environment = first(Environment),
    rugosity = first(rugosity)
  )


cor.test(stand_df$area_2d, stand_df$rugosity, method = "spearman")
cor.test(stand_df$area_2d, stand_df$fractal, method = "spearman")


ggplot(stand_df, aes(x = as.character(rugosity), y = area_2d, label = Site)) +
  geom_boxplot() +
  geom_text(nudge_y = 1) +
  theme_classic() +
  labs(
    x = "Structural rugosity",
    y = "Monospecific stand area (m²)"
  )


###stand comparison
metadata2 <- metadata2 %>%
  mutate(
    sample_group = case_when(
      Size == 1        ~ "isolated",
      Size <= 4        ~ "patch",
      Size > 4 & area_2d >=3.14        ~ "monostand",
      Size > 4 & area_2d <3.14        ~ "patch"
    )
  )

metadata2 <- metadata2 %>%
  mutate(
    genotype_id = ifelse(
      clone_lineage == 0,
      paste0("U_", INDV),
      as.character(clone_lineage)
    )
  )

#count of samples per stand
sum(metadata2$sample_group == "monostand", na.rm = TRUE)
metadata2 %>%
  group_by(Site) %>%
  summarise(n_monostand = sum(sample_group == "monostand", na.rm = TRUE))

metadata2 %>%
  group_by(Site) %>%
  summarise(
    n_isolated = sum(sample_group == "isolated", na.rm = TRUE),
    n_patch    = sum(sample_group == "patch", na.rm = TRUE)
  )


unique_genotypes_per_shape <- metadata2 %>%
  group_by(shape_id) %>%
  summarise(
    # Count all true unique genotypes
    n_unique = sum(clone_lineage == 0),
    
    # Count one representative per clonal lineage (clone_lineage > 0)
    n_clonal = n_distinct(clone_lineage[clone_lineage > 0]),
    
    # Total unique genotypes
    total_unique_genotypes = n_unique + n_clonal,
    # Total number of samples in this Shape
    total_samples = n(),
    
    # Unique / total ratio
    unique_ratio = total_unique_genotypes / total_samples,
    area_2d = first(area_2d),
    Site = first(Site)
  )
unique_genotypes_per_shape <- unique_genotypes_per_shape %>% filter(!is.na(shape_id))

cor.test(unique_genotypes_per_shape$area_2d, unique_genotypes_per_shape$total_unique_genotypes)
#Pearson's product-moment correlation

#data:  unique_genotypes_per_shape$area_2d and unique_genotypes_per_shape$total_unique_genotypes
#t = 23.453, df = 236, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7937984 0.8710285
#sample estimates:
#     cor 
#0.836521 


unique_genotypes_per_shape %>%
  group_by(Site) %>%
  summarize(
    correlation = cor(area_2d, total_unique_genotypes),
    n = n()
  )



ggplot(unique_genotypes_per_shape, aes(x = area_2d, y = total_unique_genotypes, color = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal(base_size = 14) +
  labs(title = "Genotypic Diversity vs. Stand Size",
       x = "Area in m2",
       y = "Total Unique Genotypes")

# average ratio across stands (multi-sample)
metadata2 %>%
  group_by(sample_group) %>%
  summarise(
    N_samples = n(),
    N_genotypes = n_distinct(genotype_id),
    NgN_ratio = N_genotypes/N_samples
  )

big_stands <- metadata2 %>%
  filter(area_2d >= 3.14 & Size > 4)
big_patches <- metadata2 %>%
  filter(area_2d < 3.14 & Size > 4)


### pie plot 1
library(dplyr)
library(ggplot2)
library(ggforce)
# Build pie geome"ggplot2"# Build pie geometry df ----------------------------
stand_comp <- big_stands %>%
  group_by(shape_id, genotype_id) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(shape_id) %>%
  mutate(
    total = sum(N),
    fraction = N / total,
    angle_end = cumsum(fraction) * 2*pi,
    angle_start = lag(angle_end, default = 0),
    radius = sqrt(total)   # pie size proportional to stand size
  )
stand_comp$clone_lineage <- as.factor(stand_comp$genotype_id)
write_csv2(stand_comp, "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_S2_monostand_genotypic_diversity/monostand_comp.csv")

# Get the largest radius across all pies
max_radius <- max(stand_comp$radius)

# Distinct colors for genotypes

library(randomcoloR)
distinct_palette <- distinctColorPalette(length(unique(stand_comp$clone_lineage)))
ggplot(stand_comp) +
  geom_point(
    aes(x = Inf, y = Inf, size = total),
    alpha = 0
  ) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = radius,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage,
      size = total      # drives size legend
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  coord_fixed(
    xlim = c(-max_radius, max_radius),
    ylim = c(-max_radius, max_radius)
  ) +
  scale_size_continuous(
    name = "Stand Size (n samples)",
    range = c(0.5, 3)
  ) +
  scale_fill_manual(
    values = distinct_palette,
    name = "Genotype (clone_lineage)"
  ) +
  theme_void() +
  theme(strip.text = element_text(size = 6))

max_fraction_per_stand <- stand_comp %>%
  group_by(shape_id) %>%
  summarise(
    max_fraction = max(fraction, na.rm = TRUE),
    .groups = "drop"
  )

dominance_of_genotypes <- metadata2 %>%
  filter(sample_group == "monostand") %>%
  group_by(shape_id) %>%
  summarise(
    total_samples = n(),
    n_genotypes   = n_distinct(genotype_id),
    dominant_genotype = genotype_id[which.max(table(genotype_id))],
    dominant_count    = max(table(genotype_id)),
    dominance_pct = dominant_count / total_samples * 100
  )

mean(dominance_of_genotypes$dominance_pct)
min(dominance_of_genotypes$dominance_pct)
max(dominance_of_genotypes$dominance_pct)

####pie plot 2

patch_comp <- big_patches %>%
  group_by(shape_id, genotype_id) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(shape_id) %>%
  mutate(
    total = sum(N),
    fraction = N / total,
    angle_end = cumsum(fraction) * 2*pi,
    angle_start = lag(angle_end, default = 0),
    radius = sqrt(total)   # pie size proportional to stand size
  )



# Get the largest radius across all pies
max_radius <- max(patch_comp$radius)

# Distinct colors for genotypes
patch_comp$clone_lineage <- as.factor(patch_comp$genotype_id)
write_csv2(patch_comp, "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_S2_monostand_genotypic_diversity/patch_comp.csv")
distinct_palette <- distinctColorPalette(length(unique(patch_comp$clone_lineage)))
ggplot(patch_comp) +
  geom_point(
    aes(x = Inf, y = Inf, size = total),
    alpha = 0
  ) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = radius,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage,
      size = total      # drives size legend
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  coord_fixed(
    xlim = c(-max_radius, max_radius),
    ylim = c(-max_radius, max_radius)
  ) +
  scale_size_continuous(
    name = "Stand Size (n samples)",
    range = c(0.5, 3)
  ) +
  scale_fill_manual(
    values = distinct_palette,
    name = "Genotype (clone_lineage)"
  ) +
  theme_void() +
  theme(strip.text = element_text(size = 6))

max_fraction_per_patch <- patch_comp %>%
  group_by(shape_id) %>%
  summarise(
    max_fraction = max(fraction, na.rm = TRUE),
    .groups = "drop"
  )


ng_n_by_stand <- metadata2 %>%
  filter(area_2d > 3.14 & Size > 4) %>%              # only larger stands
  group_by(shape_id, LOC, Size) %>%
  summarise(
    N_ramets = n(),
    Ng_genets = n_distinct(clone_lineage),
    Ng_N = Ng_genets / N_ramets,
    .groups = "drop"
  )

ng_n_summary_group <- ng_n_by_stand %>%
  group_by(clone_lineage) %>%
  summarise(
    n_stands = n(),
    mean_NgN = mean(Ng_N),
    sd_NgN = sd(Ng_N),
    se_NgN = sd_NgN / sqrt(n_stands),
    .groups = "drop"
  )

ng_n_summary_site <- ng_n_by_stand %>%
  group_by(LOC) %>%
  summarise(
    n_stands = n(),
    mean_NgN = mean(Ng_N),
    sd_NgN = sd(Ng_N),
    se_NgN = sd_NgN / sqrt(n_stands),
    .groups = "drop"
  )

kruskal.test(Ng_N ~ LOC, data = ng_n_by_stand)

pairwise.wilcox.test(
  ng_n_by_stand$Ng_N,
  ng_n_by_stand$LOC,
  p.adjust.method = "BH"
)

####Comparison of genotypic diversity
unique_genotypes_per_shape <- unique_genotypes_per_shape %>%
  mutate(
    sample_group = case_when(
      total_samples == 1        ~ "isolated",
      total_samples <= 4        ~ "patch",
      total_samples > 4         ~ "monostand"
    ),
    unique_pct = 100 * n_unique / total_samples
  )
summary_by_group <- unique_genotypes_per_shape %>%
  group_by(sample_group) %>%
  summarise(
    n_stands = n(),
    mean_unique_pct = mean(unique_pct, na.rm = TRUE),
    sd_unique_pct   = sd(unique_pct, na.rm = TRUE),
    se_unique_pct   = sd_unique_pct / sqrt(n_stands),
    mean_ngn        = mean(unique_ratio),
    sd_ngn        = sd(unique_ratio)
  )

weighted_summary <- unique_genotypes_per_shape %>%
  group_by(sample_group) %>%
  summarise(
    total_unique = sum(n_unique),
    total_samples = sum(total_samples),
    weighted_unique_pct = 100 * total_unique / total_samples
  )
kruskal.test(unique_pct ~ sample_group, data = unique_genotypes_per_shape)
pairwise.wilcox.test(
  unique_genotypes_per_shape$unique_pct,
  unique_genotypes_per_shape$sample_group,
  p.adjust.method = "BH"
)
m_binom <- glm(
  cbind(n_unique, n_clonal) ~ sample_group,
  family = binomial,
  data = unique_genotypes_per_shape
)

summary(m_binom)

emmeans(m_binom, ~ sample_group, type = "response")



#### random sampling simulation
set.seed(123)
groups <- unique(metadata2$sample_group)

group_sets <- c(
  combn(groups, 1, simplify = FALSE),
  combn(groups, 2, simplify = FALSE),
  combn(groups, 3, simplify = FALSE)
)
simulate_unique_draws <- function(df, pool_size = 100, n_draw = 20, n_iter = 100000) {
  
  unique_counts <- replicate(n_iter, {
    
    # Step 1: standardise pool size
    pool <- df %>%
      slice_sample(n = pool_size, replace = FALSE)
    
    # Step 2: draw from that pool
    draws <- pool %>%
      slice_sample(n = n_draw, replace = FALSE)
    
    length(unique(draws$genotype_id))
  })
  
  tibble(
    prob_all_unique = mean(unique_counts == n_draw),
    avg_unique_genotypes = mean(unique_counts)
  )
}

results <- lapply(group_sets, function(gset) {
  subset_metadata2 <- metadata2 %>%
    filter(!is.na(genotype_id)) %>%
    filter(!is.na(sample_group)) %>%
    filter(sample_group %in% gset)
  if(nrow(subset_metadata2) < 150) return(NULL)
  sim <- simulate_unique_draws(subset_metadata2, pool_size = 100, n_draw = 20)
  sim %>%
    mutate(sample_group = paste(gset, collapse = "+")) %>%
    select(sample_group, everything())
})

results <- bind_rows(results)
results
# Include combinations of sampling groups
df <- metadata2 %>%
  filter(!is.na(genotype_id),
         !is.na(sample_group))

# get group names
groups <- unique(df$sample_group)

# generate all combinations (1,2,3 groups)
group_sets <- c(
  combn(groups, 1, simplify = FALSE),
  combn(groups, 2, simplify = FALSE),
  combn(groups, 3, simplify = FALSE)
)

### simulation within each site
results <- lapply(group_sets, function(gset) {
  
  subset_df <- metadata2 %>%
    filter(sample_group %in% gset)
  
  site_results <- subset_df %>%
    group_by(Site) %>%
    group_modify(~ {
      
      # skip small sites
      if(nrow(.x) < 20) return(NULL)
      
      simulate_unique_draws(.x, pool_size = 20, n_draw = 20)
    })
  
  site_results %>%
    mutate(sample_group = paste(gset, collapse = "+")) %>%
    select(sample_group, Site, everything())
})

results <- bind_rows(results)
results %>%
  arrange(avg_unique_genotypes) %>%
  print(n = 21)
####Genotypes shared between stands

genotype_presence <- metadata2 %>%
  dplyr::select(genotype_id, sample_group) %>%
  distinct() %>%                
  mutate(present = 1) %>%
  pivot_wider(
    names_from = sample_group,
    values_from = present,
    values_fill = 0
  )
stand_types <- c("monostand", "patch", "isolated")

shared_pct <- expand.grid(
  from = stand_types,
  to   = stand_types,
  stringsAsFactors = FALSE
) %>%
  filter(from != to) %>%
  rowwise() %>%
  mutate(
    total_from = sum(genotype_presence[[from]] == 1),
    shared     = sum(genotype_presence[[from]] == 1 &
                       genotype_presence[[to]] == 1),
    shared_pct = 100 * shared / total_from
  ) %>%
  ungroup()



#genotypes shared between stands
genotype_stand <- metadata2 %>%
  dplyr::select(genotype_id, shape_id, sample_group) %>%
  distinct()
genotype_by_stand <- genotype_stand %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = shape_id,
    values_from = present,
    values_fill = 0
  )

genotype_stand %>%
  group_by(genotype_id) %>%
  summarise(
    n_stands = n_distinct(shape_id),
    n_groups = n_distinct(sample_group),
    stands = paste(unique(shape_id), collapse = ", "),
    groups = paste(unique(sample_group), collapse = ", ")
  ) %>%
  arrange(desc(n_stands))

genotype_stand <- genotype_stand %>%
  filter(!is.na(sample_group))  # remove rows with NA in sample_group if needed

# Summarise per genotype
genotype_summary <- genotype_stand %>%
  group_by(genotype_id) %>%
  summarise(
    n_monostand = sum(sample_group == "monostand"),
    n_patch = sum(sample_group == "patch"),
    n_isolated = sum(sample_group == "isolated"),
    n_unique_groups = n_distinct(sample_group[sample_group %in% c("monostand", "patch", "isolated")]),
    n_unique_stands = n_distinct(shape_id[!is.na(shape_id)])
  )

# Percentage of genotypes in more than 1 monostand / patch / isolated
perc_multi_monostand <- mean(genotype_summary$n_monostand > 1) * 100
perc_multi_patch <- mean(genotype_summary$n_patch > 1) * 100
perc_multi_isolated <- mean(genotype_summary$n_isolated > 1) * 100

# Percentage of genotypes found in multiple types of sample_group combinations
perc_monostand_patch <- mean(genotype_summary$n_monostand > 0 & genotype_summary$n_patch > 0) * 100
perc_monostand_isolated <- mean(genotype_summary$n_monostand > 0 & genotype_summary$n_isolated > 0) * 100
perc_patch_isolated <- mean(genotype_summary$n_patch > 0 & genotype_summary$n_isolated > 0) * 100

# Print results
perc_multi_monostand
perc_multi_patch
perc_multi_isolated
perc_monostand_patch
perc_monostand_isolated
perc_patch_isolated



##### Genetic diversity in Monospecific stands
library("hierfstat")
library("adegenet")
library("ape")
library("vegan")
library("vcfR")
library("readr")
library("pegas")
library("poppr")
library("mmod")
library("PopGenReport")
gen_div_stands <- metadata2 %>%
  # Keep only Stands with Size > 4
  filter(Size > 4 & area_2d > 3.14) %>%
  group_by(shape_id) %>%
  summarise(
    clonal_lineages = paste(sort(unique(clone_lineage)), collapse = ","),
    samples = paste(INDV, collapse = ","),
    .groups = "drop"
  )

# Preparing the data
setwd("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/F1b_overall_pop_stats/")
MADR <- read.vcfR("MADR_reference_f1_all.vcf")
MADR_maj <- read.vcfR("MADR_reference_f1_major.vcf")
MADR_min <- read.vcfR("MADR_reference_f1_minor.vcf")
popfile <- read_delim("MADR_reference_popfile_f1_clades.txt", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
colnames(popfile) <- c("INDV", "LOC", "CLADE")
popfile$Site <- substr(popfile$LOC, 1, 3)
popfile_major <- subset(popfile, popfile$CLADE != 'minor')
popfile_minor <- subset(popfile, popfile$CLADE == 'minor')
popfile$ID <- sub("_SHA$", "", popfile$INDV)


# --- 2. Merge stand info into the popfile ---
# Ensures all individuals in the VCF stay included
popfile_stand <- popfile %>%
  left_join(monostands_long, by = "ID") %>%
  mutate(
    Stand = ifelse(is.na(shape_id), "0", shape_id)   # assign 0 if not in stand
  )

# --- 3. Assign these stand values as populations in your genind object ---
MADR_genind_sites <- vcfR2genind(MADR)

# Replace pop slot with stand membership
pop(MADR_genind_sites) <- as.factor(popfile_stand$LOC)
pop_SNA <- popfile_stand[popfile_stand$LOC == "SNA", ]
gen_SNA <- MADR_genind_sites[pop(MADR_genind_sites) == "SNA", ]
pop_SNA_stand <- pop_SNA[
  match(indNames(gen_SNA), pop_SNA$INDV),
]
popfile_stand <- popfile_stand %>%
  mutate(
    Stand = if_else(!is.na(Size) & Size > 4 & area_2d > 3.14, Stand, "no_stand")
  )

pop_SNA_stand <- pop_SNA_stand %>%
  mutate(
    Stand = if_else(!is.na(Size) & Size > 4 & area_2d > 3.14, Stand, "no_stand")
  )

pop(gen_SNA) <- as.factor(pop_SNA_stand$Stand)

pop_SEA <- popfile_stand[popfile_stand$LOC == "SEA", ]
gen_SEA <- MADR_genind_sites[pop(MADR_genind_sites) == "SEA", ]
pop_SEA <- pop_SEA %>%
  mutate(
    Stand = if_else(!is.na(Size) & Size > 4 & area_2d > 3.14, Stand, "no_stand")
  )
pop(gen_SEA) <- as.factor(pop_SEA$Stand)

pop_KAL <- popfile_stand[popfile_stand$LOC == "KAL", ]
gen_KAL <- MADR_genind_sites[pop(MADR_genind_sites) == "KAL", ]
pop_KAL <- pop_KAL %>%
  mutate(
    Stand = if_else(!is.na(Size) & Size > 4 & area_2d > 3.14, Stand, "no_stand")
  )
pop(gen_KAL) <- as.factor(pop_KAL$Stand)


fst_Stands <- pairwise.WCfst(MADR_genind_sites)
nei_dist <- dist.genpop(genind2genpop(MADR_genind_sites), method = 1)
hist(nei_dist)
stand_pop <- genind2genpop(MADR_genind_sites)
stand_tab <- stand_pop@tab
centroid <- colMeans(stand_tab)

dist_to_centroid <- apply(
  stand_tab,
  1,
  function(x) dist(rbind(x, centroid))
)
dist_to_centroid


fst_Stands_SNA <- pairwise.WCfst(gen_SNA)
nei_dist_SNA <- dist.genpop(genind2genpop(gen_SNA), method = 1)
hist(nei_dist_SNA)
stand_pop_sna <- genind2genpop(gen_SNA)
stand_tab <- stand_pop_sna@tab
centroid <- colMeans(stand_tab)

dist_to_centroid <- apply(
  stand_tab,
  1,
  function(x) dist(rbind(x, centroid))
)
dist_to_centroid

fst_Stands_SEA <- pairwise.WCfst(gen_SEA)


nei_dist_SEA <- dist.genpop(genind2genpop(gen_SEA), method = 1)
hist(nei_dist_SEA)
stand_pop_sea <- genind2genpop(gen_SEA)
stand_tab <- stand_pop_sea@tab
centroid <- colMeans(stand_tab)

dist_to_centroid <- apply(
  stand_tab,
  1,
  function(x) dist(rbind(x, centroid))
)
dist_to_centroid

fst_Stands_KAL <- pairwise.WCfst(gen_KAL)
nei_dist_KAL <- dist.genpop(genind2genpop(gen_KAL), method = 1)
hist(nei_dist_KAL)

nei_df_SNA <- as.data.frame(as.table(as.matrix(nei_dist_SNA)))
colnames(nei_df_SNA) <- c("ind1", "ind2", "nei")
nei_df_SNA$type <- case_when(
  grepl("no_stand", nei_df_SNA$ind1) | grepl("no_stand", nei_df_SNA$ind2) ~ "no_stand_pairs",
  substr(nei_df_SNA$ind1,1,5) == substr(nei_df_SNA$ind2,1,5) ~ "within_stand",
  TRUE ~ "between_stand"
)
nei_df_SNA %>%
  group_by(type) %>%
  summarise(
    mean = mean(nei, na.rm = TRUE),
    sd = sd(nei, na.rm = TRUE),
    n = n()
  )

nei_df_SEA <- as.data.frame(as.table(as.matrix(nei_dist_SEA)))
colnames(nei_df_SEA) <- c("ind1", "ind2", "nei")
nei_df_SEA$type <- case_when(
  grepl("no_stand", nei_df_SEA$ind1) | grepl("no_stand", nei_df_SEA$ind2) ~ "no_stand_pairs",
  substr(nei_df_SEA$ind1,1,5) == substr(nei_df_SEA$ind2,1,5) ~ "within_stand",
  TRUE ~ "between_stand"
)
nei_df_SEA %>%
  group_by(type) %>%
  summarise(
    mean = mean(nei, na.rm = TRUE),
    sd = sd(nei, na.rm = TRUE),
    n = n()
  )
