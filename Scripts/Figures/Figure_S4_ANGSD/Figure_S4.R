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


setwd("~/Documents/Inkfish-Phd/MADR_CH1/Figure_2_clonality_assessment/")

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
miss <- read_table("~/Documents/Inkfish-Phd/MADR_CH1/E1B_ANGSD/indv_miss_MADR_ANGSD_major.txt",
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

#LOC       N    Ng  Ng_N
#<chr> <int> <int> <dbl>
# 1 KAL     227    85 0.374
# 2 SEA      89    27 0.303
# 3 SNA     300    70 0.233

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

miss <- read_table("~/Documents/Inkfish-Phd/MADR_CH1/E1B_ANGSD/indv_miss_MADR_ANGSD_minor.txt",
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

hybrids <- c("DH0293_MMIR_SNA", "DH0360_MMIR_SNA", "DH0371_MMIR_SNA")

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

#Depth          N    Ng  Ng_N
#1 reef flat    176    61 0.347
#2 reef slope   505   165 0.327
