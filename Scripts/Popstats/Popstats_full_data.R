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
library(graph4lg)

# Preparing the data
setwd("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/F1b_overall_pop_stats/")
MADR <- read.vcfR("MADR_reference_f1_all.vcf")
MADR_maj <- read.vcfR("MADR_reference_f1_major.vcf")
MADR_min <- read.vcfR("MADR_reference_f1_minor.vcf")
popfile <- read_delim("MADR_reference_popfile_f1_clades.txt", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
colnames(popfile) <- c("INDV", "LOC", "CLADE")
popfile$depth <- ifelse(grepl("_SHA", popfile$LOC), "flat", "slope")
popfile$site <- substr(popfile$LOC, 1, 3)
popfile_major <- subset(popfile, popfile$CLADE != 'minor')
popfile_minor <- subset(popfile, popfile$CLADE == 'minor')

#popstats per site among clades
MADR_genind_sites <- vcfR2genind(MADR)
pop(MADR_genind_sites) <- as.factor(popfile$site)
gen_hier_sites <- genind2hierfstat(MADR_genind_sites)
basic_stats_sites <- basic.stats(gen_hier_sites)
basic_stats_sites
$overall
Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
0.1066 0.1308 0.1309 0.0002 0.1310 0.0002 0.0012 0.0019 0.1850 0.0003 

#within clade major
MADR_genind_sites_maj <- vcfR2genind(MADR_maj)
pop(MADR_genind_sites_maj) <- as.factor(popfile_major$site)
gen_hier_sites_major <- genind2hierfstat(MADR_genind_sites_maj)
basic_stats_sites_major <- basic.stats(gen_hier_sites_major)

basic_stats_sites_major
$overall
Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
0.1301 0.1516 0.1516 0.0000 0.1516 0.0000 0.0000 0.0000 0.1420 0.0000

#within clade minor
MADR_genind_sites_min <- vcfR2genind(MADR_min)
pop(MADR_genind_sites_min) <- as.factor(popfile_minor$site)
gen_hier_sites_minor <- genind2hierfstat(MADR_genind_sites_min)
basic_stats_sites_minor <- basic.stats(gen_hier_sites_minor)

basic_stats_sites_minor
$overall
Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
0.1706  0.2155  0.2151 -0.0004  0.2150 -0.0006 -0.0018 -0.0026  0.2084 -0.0007 


#popstats between clade
MADR_genind_clades <- vcfR2genind(MADR)
pop(MADR_genind_clades) <- as.factor(popfile$CLADE)
gen_hier_clades <- genind2hierfstat(MADR_genind_clades)
basic_stats_clades <- basic.stats(gen_hier_clades)
basic_stats_clades
$overall
Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
0.1128 0.1382 0.1465 0.0083 0.1512 0.0129 0.0566 0.0856 0.1839 0.0150 

# Calculating basic population statistics over depth
MADR_genind_depth <- vcfR2genind(MADR)
pop(MADR_genind_depth) <- as.factor(popfile$depth)
gen_hier_depth <- genind2hierfstat(MADR_genind_depth)
basic_stats_depth <- basic.stats(gen_hier_depth)
basic_stats_depth
$overall
Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
0.1072 0.1282 0.1284 0.0003 0.1287 0.0005 0.0020 0.0041 0.1635 0.0006 

wc_results <- wc(gen_hier_sites)
fst <- wc_results$FST
fis <- wc_results$FIS
fst
[1] 0.002588539
fis
[1] 0.1822069

### --- Fst & Fis between sites (ignoring clade) ---
hierf_site <- genind2hierfstat(MADR_genind_sites)

fst_by_site <- pairwise.WCfst(hierf_site)
fst_by_site
fis_by_site <- wc(hierf_site)$FIS
fis_by_site

### --- Fst & Fis between clades ---
hierf_clade <- genind2hierfstat(MADR_genind_clades)

fst_by_clade <- pairwise.WCfst(hierf_clade)
fst_by_clade
fis_by_clade <- wc(hierf_clade)$FIS
fis_by_clade

### --- Fst & Fis between sites within each clade --- NOTE use separate vcfs
unique_clades <- unique(popfile$CLADE)
fst_within_clades <- list()
fis_within_clades <- list()

for (clade in unique_clades) {
  idx <- which(popfile$CLADE == clade)
  gen_sub <- MADR_genind[idx]
  pop(gen_sub) <- popfile$site[idx]
  gen_sub_hier <- genind2hierfstat(gen_sub)
  
  fst_within_clades[[clade]] <- pairwise.WCfst(gen_sub_hier)
  fis_within_clades[[clade]] <- wc(gen_sub_hier)$FIS
}
fst_within_clades
fis_within_clades

### --- Print results ---
cat("Fst between sampling sites (ignoring clade):\n")
print(fst_by_site)
cat("\nFis per site (ignoring clade):\n")
print(fis_by_site)

cat("\nFst between clades:\n")
print(fst_by_clade)
cat("\nFis per clade:\n")
print(fis_by_clade)

cat("\nFst between sites within each clade:\n")
print(fst_within_clades)
cat("\nFis within each clade:\n")
print(fis_within_clades)

# Observed and Expected Heterozygosity
#overall
sumstats <- summary(MADR_genind)
Ho <- sumstats$Hobs
He <- sumstats$Hexp

mean_Ho <- mean(Ho, na.rm = TRUE)
mean_He <- mean(He, na.rm = TRUE)

# mean observed heterozygosity per population (clade)
mean_Ho_per_pop <- apply(basic_stats$Ho, 2, mean, na.rm = TRUE)

# mean expected heterozygosity per population (clade)
mean_He_per_pop <- apply(basic_stats$Hs, 2, mean, na.rm = TRUE)

mean_Ho_per_pop
mean_He_per_pop

jost_d_matrix <- pairwise_D(MADR_genind_clades)

# Htp (average within-pop heterozygosity)
htp <- apply(basic_stats$Ho, 2, mean)

# Print results
print(list(
  Fst = fst_overall,
  Fis = fis_overall,
  Ho = Ho,
  He = He,
  Htp = htp,
  Dst = dst_stats,
  Fstp = fstp_stats,
  Dest = dest_stats
))

strata(MADR_genind) <- data.frame(CLADE = popfile$CLADE)
amova_res <- poppr.amova(MADR_genind, ~CLADE, within = FALSE)

  
#basic stats per clade
pop_list <- seppop(MADR_genind_clades)
maj_hier <- genind2hierfstat(pop_list$major)
min_hier <- genind2hierfstat(pop_list$minor)
basic_stats_maj <- basic.stats(maj_hier)
basic_stats_min <- basic.stats(min_hier)
basic_stats_maj$overall
basic_stats_min$overall

#Amova
library(poppr)

# Assign strata (groups)
strata(MADR_genind) <- data.frame(CLADE = popfile$CLADE)

# Run AMOVA
amova_res <- poppr.amova(MADR_genind, ~CLADE)
print(amova_res)

# Test significance
set.seed(123)
amova_test <- randtest(amova_res, nrepet = 999)
print(amova_test)
plot(amova_test)


#allelic richness
gen <- vcfR2genind(
  MADR,
  ploidy = 2,
  sep = "/"
)

# Assign clade
pop(gen) <- popfile$CLADE
hf <- genind2hierfstat(gen)
ar <- allelic.richness(hf)
apply(ar$Ar, 2, mean, na.rm = TRUE)
boxplot(ar$Ar, las = 2, ylab = "Allelic richness")

###individually

MADR_maj <- read.vcfR("MADR_reference_f1_major.vcf")
MADR_min <- read.vcfR("MADR_reference_f1_minor.vcf")
popfile_major <- subset(popfile, popfile$CLADE != 'minor')
popfile_minor <- subset(popfile, popfile$CLADE == 'minor')

gen_maj <- vcfR2genind(
  MADR_maj,
  ploidy = 2,
  sep = "/"
)
pop(gen_maj) <- popfile_major$site
hf_maj <- genind2hierfstat(gen_maj)
ar_maj <- allelic.richness(hf_maj)
bs_maj <- basic.stats(hf_maj)
#$overall
#Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
#0.1301 0.1516 0.1516 0.0000 0.1516 0.0000 0.0000 0.0000 0.1420 0.0000 
AR_maj <- mean(
  apply(ar_maj$Ar, 2, mean, na.rm = TRUE),
  na.rm = TRUE
)
Hs_maj <- mean(bs_maj$Hs, na.rm = TRUE)
apply(ar_maj$Ar, 2, mean, na.rm = TRUE)
#SNA      KAL      SEA 
#1.564272 1.565355 1.603068 
ar_vals <- c(
  SNA = 1.564272,
  KAL = 1.565355,
  SEA = 1.603068
)
mean_ar <- mean(ar_vals)
sd_ar   <- sd(ar_vals)

boxplot(ar_maj$Ar, las = 2, ylab = "Allelic richness")

gen_min <- vcfR2genind(
  MADR_min,
  ploidy = 2,
  sep = "/"
)
pop(gen_min) <- popfile_minor$site
hf_min <- genind2hierfstat(gen_min)
ar_min <- allelic.richness(hf_min)
bs_min <- basic.stats(hf_min)
#$overall
#Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.1706  0.2155  0.2151 -0.0004  0.2150 -0.0006 -0.0018 -0.0026  0.2084 -0.0007 
AR_min <- mean(
  apply(ar_min$Ar, 2, mean, na.rm = TRUE),
  na.rm = TRUE
)
Hs_min <- mean(bs_min$Hs, na.rm = TRUE)
apply(ar_min$Ar, 2, mean, na.rm = TRUE)
#KAL      SNA      SEA 
#1.209830 1.214093 1.207698 
ar_vals <- c(
  SNA = 1.214093,
  KAL = 1.209830,
  SEA = 1.207698
)
mean_ar <- mean(ar_vals)
sd_ar   <- sd(ar_vals)
boxplot(ar_min$Ar, las = 2, ylab = "Allelic richness")

data.frame(
  Clade = c("major", "minor"),
  Allelic_Richness = c(AR_maj, AR_min),
  Hs = c(Hs_maj, Hs_min)
)

# compute minimum sample size (in genes) for rarefaction:
# hierfstat::allelic.richness expects 'min.n' in number of gene copies (i.e. 2 * individuals for diploids)
n_per_clade <- table(pop(MADR_genind_clades))
min_ind <- min(n_per_clade)
min_n_genes <- min_ind * 2

# compute allelic richness with rarefaction to smallest clade
ar_res <- allelic.richness(hierf_clade, min.n = min_n_genes)

# 'ar_res$Ar' usually contains per-population, per-locus richness; 'ar_res$mean.richness' or column means give summary
# Example: mean allelic richness per population across loci
mean_ar_per_pop <- colMeans(ar_res$Ar, na.rm = TRUE)
mean_ar_per_pop



### preparing for NeEstimator
MADR_maj <- read.vcfR("MADR_reference_f1_major.vcf")
fix <- as.data.frame(MADR_maj@fix)
fix$POS <- as.numeric(fix$POS)

thin_snps <- fix %>%
  group_by(CHROM) %>%
  arrange(POS) %>%
  filter(POS - lag(POS, default = -Inf) >= 5000)
vcf_ids <- MADR_maj@fix[, "ID"]
keep_rows <- which(vcf_ids %in% thin_snps$ID)
MADR_maj_thinned <- MADR_maj[keep_rows, ]

MADR_min <- read.vcfR("MADR_reference_f1_minor.vcf")
fix <- as.data.frame(MADR_min@fix)
fix$POS <- as.numeric(fix$POS)

thin_snps <- fix %>%
  group_by(CHROM) %>%
  arrange(POS) %>%
  filter(POS - lag(POS, default = -Inf) >= 5000)
vcf_ids <- MADR_min@fix[, "ID"]
keep_rows <- which(vcf_ids %in% thin_snps$ID)
MADR_min_thinned <- MADR_min[keep_rows, ]
popfile_major <- subset(popfile, popfile$CLADE != 'minor')
popfile_minor <- subset(popfile, popfile$CLADE == 'minor')
gen_maj <- vcfR2genind(MADR_maj, ploidy = 2, sep = "/")
gen_min <- vcfR2genind(MADR_min, ploidy = 2, sep = "/")
gen_maj <- gen_maj[, minorAllele(gen_maj) >= 0.05]
gen_min <- gen_min[, minorAllele(gen_min) >= 0.05]
gen_maj

pop(gen_maj) <- factor(rep("ALL", nInd(gen_maj)))
hierf_maj <- genind2hierfstat(gen_maj)
hierf_maj <- hierf_maj[hierf_maj$pop != "dumpop", ]
hierfstat::write.fstat(hierf_maj, fname = "Major_fstat.dat")

pop(gen_min) <- factor(rep("ALL", nInd(gen_min)))
hierf_min <- genind2hierfstat(gen_min)
hierf_min <- hierf_min[hierf_min$pop != "dumpop", ]
hierfstat::write.fstat(hierf_min, fname = "Minor_fstat.dat")



