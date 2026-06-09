library(WGCNA) # for coloring clonal groups
library(sparcl) # for ColorDendrogram
library("reshape2")
library("readr")
library(dendextend)   # colour branches
library(ggdendro)     # convert dendrogram → data frames for ggplot
library(ggplot2)
library(patchwork)    # combine the two ggplots
library(RColorBrewer)
library(ggtreeExtra)
library(ggstance)
library(ggtree)
library(ape)
library(ggnewscale)  # allows multiple fill scales if needed
library(ggtreeExtra) # for additional panels
library(dplyr)


setwd("~/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_S1_clonality_assessment/")

# reading list of bam files = order of samples in IBS matrix
bams=read.table("sample_names_major_bam.txt",header=F)

# reading IBS matrix based on SNPs with allele frequency >= 0.05:
ma = as.matrix(read.table("MADR_reference_e1_major_ANGSD_outputs.ibsMat"))
samples= bams$V1
dimnames(ma)=list(samples,samples)

# plotting hierarchical clustering tree
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.6)
abline(h=0.0156,col="red",lty=3)

# "artificial clones" (genotyping replicates)

cl1 = c("DX0054_MMIR_SNA_SHA", "DH0054_MMIR_SNA_SHA")
cl2 = c("DX0114_MMIR_KAL_SHA", "DH0114_MMIR_KAL_SHA")
cl3 = c("DX0452_MMIR_SNA", "DH0452_MMIR_SNA")
cl4 = c("DX0453_MMIR_SNA", "DH0453_MMIR_SNA")
cl5 = c("DX0454_MMIR_SNA", "DH0454_MMIR_SNA")

# color artificial clones red
art.clones=rep("black",nrow(ma))
art.clones[bams$V1 %in% c(cl1,cl2,cl3,cl4,cl5)]="red"
ColorDendrogram(hclust(as.dist(ma),"ave"), y = art.clones,branchlength=0.0022)

# sorting samples into clonal groups; singletons go into the same "color" group (for plotting) but different "cn" groups (for GLM modeling later)
cutoff <- 0.0156
abline(h = cutoff, col = "red", lty = 3)

# Cut tree by threshold
cc <- cutree(hc, h = cutoff)
tc <- table(cc)
singletons <- as.numeric(names(tc)[tc == 1])
cc[cc %in% singletons] <- singletons[1]

# Get unique group count
n_groups <- length(unique(cc))

# Generate visually distinct colors
color_palette <- colorRampPalette(brewer.pal(min(8, n_groups), "Set2"))(n_groups)
group_colors <- setNames(color_palette, sort(unique(cc)))
clones <- group_colors[as.character(cc)]
clones[names(clones) == "5"] <- "green"
# Final dendrogram with group colors
ColorDendrogram(hc, y = clones, labels = FALSE, branchlength = 0.0022)
abline(h = cutoff, col = "red", lty = 3)
# Historgram of node height to visualise the cutoff point

# Histogram of pairwise IBS distances
hist(hc$height,
     breaks = 50,
     main = "Histogram of node height IBS dendrogram major clade",
     xlab = "Node Height",
     col = "blue",
     border = "white")

# Add vertical threshold line
abline(v = cutoff, col = "red", lty = 3)

## Save clonal lineage in popfile

pop_with_lineage_maj  <- tibble(
  INDV          = names(cc),
  clone_lineage = as.integer(cc)   # keep it numeric
)

pop_with_lineage_maj <- pop_with_lineage_maj %>%
  mutate(LOC = str_extract(INDV, "(?<=MMIR_)[^_]+"),
         Depth = case_when(
           str_ends(INDV, "_SHA") ~ "reef flat",
           TRUE                 ~ "reef slope"
         )
  )

pop_with_lineage_maj$CLADE <-"Major"
pop_with_lineage_maj$clone_lineage[pop_with_lineage_maj$clone_lineage == 5] <- 0
write.table(
  pop_with_lineage_maj,
  file      = "MADR_reference_popfile_e1b_major.txt",
  quote     = FALSE,
  row.names = FALSE,
  sep       = "\t"
)

#Adding missing data
miss <- read_table("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/E1B_ANGSD/indv_miss_MADR_ANGSD_major.txt",
                   col_names = FALSE)

colnames(miss) <- c("INDV", "missing", "sites")   # <- easier name than `miss%`
miss$missing   <- as.numeric(miss$missing)        # make sure it’s numeric
miss <- miss |>
  group_by(INDV) |>
  slice(1) |>
  ungroup()

phy <- as.phylo(hc)
tip_order <- hc$labels[hc$order]
tree_reordered <- keep.tip(phy, tip_order[tip_order %in% phy$tip.label])
plot(tree_reordered, show.tip.label = F)

# Create a data frame that maps tip labels to clone colors
clone_df <- data.frame(label = names(clones), clone_color = clones)

# Build the tree and color the tips by clone_color
tree <- ggtree(tree_reordered) %<+% clone_df +
  geom_tree(aes(color = clone_color)) +
  theme_tree2() +
  theme(legend.position = "none")
tree

miss <- miss |>
  filter(INDV %in% tree$data$label)          
# Reorder 'miss' so that INDV matches the order of tip labels
miss_ordered <- miss[match(tip_order, miss$INDV), ]
tip_order <- tree$data %>%                # the data slot of the ggtree object
  filter(isTip) %>%                       # keep only the tips
  arrange(y) %>%                          # arrange by plotted y-position (top→bottom)
  pull(label)                             # grab the labels
miss_ordered$label <- miss_ordered$INDV

facet_plot(
  tree,
  panel   = "missing",
  data    = miss_ordered,
  geom    = ggstance::geom_barh,
  mapping = aes(x = missing),   
  stat    = "identity",
  width   = 0.9,
  color = "black")

ng_n_summary <- pop_with_lineage_maj %>%
  group_by(LOC) %>%
  summarise(
    N  = n(),  # total individuals
    Ng = sum(clone_lineage == 0) +           # unique samples
      n_distinct(clone_lineage[clone_lineage != 0]),  # one per clone lineage
    Ng_N = Ng / N
  )
ng_n_summary
#LOC       N    Ng  Ng_N
#<chr> <int> <int> <dbl>
# 1 KAL     227    85 0.374
# 2 SEA      89    27 0.303
# 3 SNA     300    70 0.233

ng_n_summary <- pop_with_lineage_maj %>%
  group_by(Depth) %>%
  summarise(
    N  = n(),  # total individuals
    Ng = sum(clone_lineage == 0) +           # unique samples
      n_distinct(clone_lineage[clone_lineage != 0]),  # one per clone lineage
    Ng_N = Ng / N
  )


###Minor clade
# reading list of bam files = order of samples in IBS matrix
bams=read.table("sample_names_minor_bam.txt",header=F)

# reading IBS matrix based on SNPs with allele frequency >= 0.05:
ma = as.matrix(read.table("MADR_reference_e1_minor_ANGSD_outputs.ibsMat"))
samples= bams$V1
dimnames(ma)=list(samples,samples)

# plotting hierarchical clustering tree
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.6)

# sorting samples into clonal groups; singletons go into the same "color" group (for plotting) but different "cn" groups (for GLM modeling later)
cutoff <- 0.018
abline(h = cutoff, col = "red", lty = 3)

# Cut tree by threshold
cc <- cutree(hc, h = cutoff)
tc <- table(cc)
singletons <- as.numeric(names(tc)[tc == 1])
cc[cc %in% singletons] <- singletons[1]

# Get unique group count
n_groups <- length(unique(cc))

# Generate visually distinct colors
color_palette <- colorRampPalette(brewer.pal(min(8, n_groups), "Set2"))(n_groups)
group_colors <- setNames(color_palette, sort(unique(cc)))
clones <- group_colors[as.character(cc)]
#changing the colors to match the major clade
clones[names(clones) == "2"] <- "green"

# Final deprogram with group colors
ColorDendrogram(hc, y = clones, labels = FALSE, branchlength = 0.0019)
abline(h = cutoff, col = "red", lty = 3)
# Histogram of node height to visualize the cutoff point

# Histogram of pairwise IBS distances
hist(hc$height,
     breaks = 30,
     main = "Histogram of node height IBS dendrogram major clade",
     xlab = "Node Height",
     col = "blue",
     border = "white")

# Add vertical threshold line
abline(v = cutoff, col = "red", lty = 3)

miss <- read_table("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/E1B_ANGSD/indv_miss_MADR_ANGSD_minor.txt",
                   col_names = FALSE)

colnames(miss) <- c("INDV", "missing", "sites")   # <- easier name than `miss%`
miss$missing   <- as.numeric(miss$missing)        # make sure it’s numeric
miss <- miss |>
  group_by(INDV) |>
  slice(1) |>
  ungroup()

phy <- as.phylo(hc)
tip_order <- hc$labels[hc$order]
tree_reordered <- keep.tip(phy, tip_order[tip_order %in% phy$tip.label])
plot(tree_reordered, show.tip.label = F)

tree <- ggtree(tree_reordered, aes(color = clone)) %<+% clone_df +
  geom_tree() +
  theme_tree2() +
  theme(legend.position = "none")
tree

miss <- miss |>
  filter(INDV %in% tree$data$label)          
# Reorder 'miss' so that INDV matches the order of tip labels
miss_ordered <- miss[match(tip_order, miss$INDV), ]
tip_order <- tree$data %>%                # the data slot of the ggtree object
  filter(isTip) %>%                       # keep only the tips
  arrange(y) %>%                          # arrange by plotted y-position (top→bottom)
  pull(label)                             # grab the labels
miss_ordered$label <- miss_ordered$INDV

facet_plot(
  tree,
  panel   = "missing",
  data    = miss_ordered,
  geom    = ggstance::geom_barh,
  mapping = aes(x = missing),   
  stat    = "identity",
  width   = 0.9,
  color = "black")

## Save clonal lineage in popfile

pop_with_lineage_min  <- tibble(
  INDV          = names(cc),
  clone_lineage = as.integer(cc)   # keep it numeric
)
pop_with_lineage_min <- pop_with_lineage_min %>%
  mutate(LOC = str_extract(INDV, "(?<=MMIR_)[^_]+"),
         Depth = case_when(
           str_ends(INDV, "_SHA") ~ "reef flat",
           TRUE                 ~ "reef slope"
         )
  )

pop_with_lineage_min$CLADE <-"Minor"
pop_with_lineage_min$clone_lineage[pop_with_lineage_min$clone_lineage == 2] <- 0
write.table(
  pop_with_lineage_min,
  file      = "MADR_reference_popfile_e1b_minor.txt",
  quote     = FALSE,
  row.names = FALSE,
  sep       = "\t"
)

ng_n_summary <- pop_with_lineage_min %>%
  group_by(LOC) %>%
  summarise(
    N  = n(),  # total individuals
    Ng = sum(clone_lineage == 0) +           # unique samples
      n_distinct(clone_lineage[clone_lineage != 0]),  # one per clone lineage
    Ng_N = Ng / N
  )

ng_n_summary
#  LOC       N    Ng  Ng_N
#<chr> <int> <int> <dbl>
#  1 KAL       5     4 0.8  
#  2 SEA      21     6 0.286
#  3 SNA      42    18 0.429

###combine the two clades
# Step 1: Find the max clone lineage in clade1 (excluding unique genotypes marked as 0)
max_clade1 <- max(pop_with_lineage_maj$clone_lineage[pop_with_lineage_maj$clone_lineage != 0])

# Step 2: Offset clade2 clone lineages (only non-zero ones)
pop_with_lineage_min <- pop_with_lineage_min %>%
  mutate(clone_lineage = ifelse(clone_lineage == 0, 0, clone_lineage + max_clade1))

# Step 3: Combine
combined_df <- bind_rows(pop_with_lineage_maj, pop_with_lineage_min)
combined_df <- combined_df %>%
  arrange(INDV, CLADE) %>%  
  distinct(INDV, .keep_all = TRUE)

hybrids <- c("DH0513_MMIR_SEA", "DH0514_MMIR_SEA")

combined_df <- combined_df %>%
  mutate(CLADE = if_else(INDV %in% hybrids, "Hybrid", CLADE))

ng_n_summary <- combined_df %>%
  group_by(LOC) %>%
  summarise(
    N  = n(),  # total individuals
    Ng = sum(clone_lineage == 0) +           # unique samples
      n_distinct(clone_lineage[clone_lineage != 0]),  # one per clone lineage
    Ng_N = Ng / N
  )

ng_n_summary
#LOC       N    Ng  Ng_N
#1 KAL     232    89 0.384
#2 SEA     110    33 0.3  
#3 SNA     339    87 0.257

write.table(combined_df, "MADR_reference_popfile_G1_full_clonal_groups.txt", quote = F)

annotations <- read_csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/combined_3d_annotations_2123_aligned_renamed.coords.csv")
combined_df$ID <- sub("_SHA*$", "", combined_df$INDV)
metadata <- combined_df %>%
  left_join(annotations, by = "ID")

metadata <- metadata %>% rename(Habitat = Depth)
metadata <- metadata %>% rename(Depth = world_z)

metadata <- metadata %>%
  mutate(
    Habitat = case_when(
      Depth >= -6.09 ~ "Terrace",
      Depth >= -12.48 & Depth < -6.09 ~ "Crest",
      Depth < -12.48 ~ "Slope"
    )
  )

ng_n_summary <- metadata %>%
  group_by(CLADE) %>%
  summarise(
    N  = n(),  # total individuals
    Ng = sum(clone_lineage == 0) +           # unique samples
      n_distinct(clone_lineage[clone_lineage != 0]),  # one per clone lineage
    Ng_N = Ng / N
  )

ng_n_summary


ng_n_summary <- metadata %>%
  group_by(Habitat) %>%
  summarise(
    N  = n(),  # total individuals
    Ng = sum(clone_lineage == 0) +           # unique samples
      n_distinct(clone_lineage[clone_lineage != 0]),  # one per clone lineage
    Ng_N = Ng / N
  )

ng_n_summary
# A tibble: 4 × 4
Habitat     N    Ng  Ng_N
<chr>   <int> <int> <dbl>
  1 Crest     313   110 0.351
2 Slope     233    72 0.309
3 Terrace    69    26 0.377
4 NA         66    44 0.667

#per site
ng_n_depth_summary <- metadata %>%
  group_by(LOC, Habitat) %>%
  summarise(
    N  = n(),  
    Ng = sum(clone_lineage == 0) +
      n_distinct(clone_lineage[clone_lineage != 0]),
    Ng_N = Ng / N,
    .groups = "drop"
    )
ng_n_depth_summary
# A tibble: 11 × 5
LOC   Habitat     N    Ng  Ng_N
<chr> <chr>   <int> <int> <dbl>
  1 KAL   Crest     144    62 0.431
2 KAL   Slope      33    12 0.364
3 KAL   Terrace    33    16 0.485
4 KAL   NA         22    14 0.636
5 SEA   Crest       8     4 0.5  
6 SEA   Slope      95    27 0.284
7 SEA   NA          7     6 0.857
8 SNA   Crest     161    46 0.286
9 SNA   Slope     105    33 0.314
10 SNA   Terrace    36    10 0.278
11 SNA   NA         37    25 0.676

ng_n_depth_summary <- metadata %>%
  group_by(LOC, CLADE) %>%
  summarise(
    N  = n(),  
    Ng = sum(clone_lineage == 0) +
      n_distinct(clone_lineage[clone_lineage != 0]),
    Ng_N = Ng / N,
    .groups = "drop"
  )
ng_n_depth_summary
# A tibble: 7 × 5
LOC   CLADE      N    Ng  Ng_N
<chr> <chr>  <int> <int> <dbl>
  1 KAL   Major    226    84 0.372
2 KAL   Minor      5     4 0.8  
3 SEA   Major     90    27 0.3  
4 SEA   Minor     21     6 0.286
5 SNA   Hybrid     3     2 0.667
6 SNA   Major    297    70 0.236
7 SNA   Minor     39    17 0.436

ng_n_depth_summary <- metadata %>%
  group_by(LOC, CLADE, Habitat) %>%
  summarise(
    N  = n(),  
    Ng = sum(clone_lineage == 0) +
      n_distinct(clone_lineage[clone_lineage != 0]),
    Ng_N = Ng / N,
    .groups = "drop"
  )
print(n = 25, ng_n_depth_summary)


library(tidyverse)
library(glmmTMB)

genet_data <- tribble(
  ~site, ~environment, ~lineage, ~Ng, ~N,
  "SNA", "full",    "major", 70, 300,
  "KAL", "full",    "major", 85, 227,
  "SEA", "full",    "major", 26, 87,
  "ALL", "full",    "major", 176, 614,
  
  "SNA", "full",    "minor", 17, 39,
  "KAL", "full",    "minor", 4, 5,
  "SEA", "full",    "minor", 6, 21,
  "ALL", "full",    "minor", 28, 65,
  
  "SNA", "terrace", "major", 10, 36,
  "KAL", "terrace", "major", 15, 32,
  
  "KAL", "terrace", "minor", 1, 1,
  
  "SNA", "crest",   "major", 44, 149,
  "KAL", "crest",   "major", 61, 143,
  "SEA", "crest",   "major", 2, 3,
  
  "SNA", "crest",   "minor", 2, 12,
  "KAL", "crest",   "minor", 1, 1,
  "SEA", "crest",   "minor", 2, 5,
  
  "SNA", "slope",   "major", 21, 82,
  "KAL", "slope",   "major", 10, 30,
  "SEA", "slope",   "major", 21, 79,
  
  "SNA", "slope",   "minor", 12, 23,
  "KAL", "slope",   "minor", 2, 3,
  "SEA", "slope",   "minor", 5, 14
) %>%
  filter(site != "ALL") %>%            # remove pooled rows
  mutate(ratio = Ng / N)
complexity <- tribble(
  ~site, ~environment, ~rugosity, ~fractal,
  "KAL", "slope",   1.785, 2.14,
  "KAL", "crest",   1.784, 2.09,
  "KAL", "terrace", 1.401, 2.02,
  "SNA", "slope",   1.703, 2.14,
  "SNA", "crest",   1.479, 2.09,
  "SNA", "terrace", 1.282, 2.01,
  "SEA", "slope",   1.791, 2.06,
  "SEA", "crest",   1.718, 2.04
)

df <- genet_data %>%
  left_join(complexity, by = c("site", "environment")) %>%
  drop_na(rugosity)   # removes environment/site combos without complexity

m1 <- glmmTMB(
  cbind(Ng, N - Ng) ~ rugosity + environment + lineage + (1 | site),
  family = binomial,
  data = df
)

summary(m1)

df_2 <- genet_data %>% filter(lineage != 'major')
m2 <- glmmTMB(
  cbind(Ng, N - Ng) ~ rugosity + environment + (1 | site),
  family = binomial,
  data = df
)

summary(m2)


m1 <- glmmTMB(
  cbind(Ng, N - Ng) ~ rugosity + environment + (1 | site),
  family = binomial,
  data = df
)

summary(m1) 

m1 <- glmmTMB(
  cbind(Ng, N - Ng) ~ lineage + environment + (1 | site),
  family = binomial,
  data = df
)

summary(m1) 

anova(
  glmmTMB(cbind(Ng, N-Ng) ~ environment + lineage + (1|site), family=binomial, data=df),
  glmmTMB(cbind(Ng, N-Ng) ~ environment + lineage + rugosity + (1|site), family=binomial, data=df)
)

library(emmeans)

emmeans(m1, ~ lineage, type = "response")
ggplot(df, aes(x = rugosity, y = ratio, size = N, colour = lineage)) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "glm",
    method.args = list(family = "binomial"),
    se = TRUE
  ) +
  scale_size_continuous(name = "Sample size (N)") +
  labs(
    x = "Structural rugosity",
    y = "Genet–ramet ratio (Ng / N)"
  ) +
  theme_classic()

m_fractal <- glmmTMB(
  cbind(Ng, N - Ng) ~ fractal + environment + lineage + (1 | site),
  family = binomial,
  data = df
)

summary(m_fractal)





