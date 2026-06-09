library(reshape2)
library(dplyr)
library(ggplot2)
library(readr)
library(svglite)
setwd("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure")

# -------------------------
# Load metadata with LOC info
# -------------------------
popfile <- read_delim("MADR_reference_popfile_DAPC.txt", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = TRUE, trim_ws = TRUE)
annotations <- read_csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/combined_3d_annotations_2123_aligned.coords.csv")

# Extract the ID part before "_MMIR"
popfile$genotype <- sub("_MMIR.*", "", popfile$INDV)

popfile <- merge(
  popfile,
  annotations[, c("genotype", "world_z")],
  by = "genotype",
  all.x = TRUE
)

popfile <- popfile %>%
  rename(DEPTH = world_z)

popfile <- popfile %>%
  mutate(
    HABITAT = case_when(
      DEPTH >= -5.5 ~ "Terrace",
      DEPTH >= -11 & DEPTH < -5.5 ~ "Crest",
      DEPTH < -11 ~ "Slope"
    )
  )

# Option 1: create a new LOC_HAB column combining location and habitat
popfile$LOC <- substr(popfile$LOC, 1, 3)
popfile$LOC_HAB <- paste(popfile$LOC, popfile$HABITAT, sep = "_")

meta <- popfile[, c("INDV", "LOC")]
# -------------------------
# Load STRUCTURE (CLUMPP) output
# -------------------------
K2 <- read.csv("clumpp_K2.out.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(K2)[1] <- "INDV"
K2 <- K2[ , -2]              # remove 2nd column
colnames(K2)[2:3] <- c("pop1", "pop2")

# Melt to long format
K2_melt <- melt(K2, id.vars = "INDV", variable.name = "cluster", value.name = "q")


# -------------------------
# Merge STRUCTURE results + metadata

# -------------------------
plot_df <- left_join(K2_melt, meta, by = "INDV")

# Sort INDV by LOC (optional)
plot_df$INDV <- factor(plot_df$INDV, levels = unique(plot_df$INDV[order(plot_df$LOC)]))

# -------------------------
# Plot
# -------------------------
plot_df$LOC <- factor(
  plot_df$LOC,
  levels = c("KAL", "SNA", "SEA")
)
structure <- ggplot(plot_df, aes(x = INDV, y = q, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("firebrick", "#56B4E9")) +
  facet_grid(~ LOC, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, "lines")
  ) +
  labs(
    x = "INDV",
    y = "Q ancestry proportion",
    fill = "Cluster"
  )

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/Structure_plot.svg",
  width = 15,
  height = 6,
  bg = "transparent"
)

print(
  structure +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


K3 <- read.csv("clumpp_K3.out.csv", header = FALSE, na.strings = "", stringsAsFactors = FALSE)
colnames(K3)[1] <- "INDV"
K3 <- K3[,-2]
colnames(K3)[2] <- "pop1"
colnames(K3)[3] <- "pop2"
colnames(K3)[4] <- "pop3"
# Melt to long format
K3_melt <- melt(K3, id.vars = "INDV", variable.name = "cluster", value.name = "q")
plot_df <- left_join(K3_melt, meta, by = "INDV")

# Sort INDV by LOC (optional)
plot_df$INDV <- factor(plot_df$INDV, levels = unique(plot_df$INDV[order(plot_df$LOC)]))
plot_df$LOC <- factor(
  plot_df$LOC,
  levels = c("KAL", "SNA", "SEA")
)

# -------------------------
# Plot
# -------------------------
ggplot(plot_df, aes(x = INDV, y = q, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("firebrick", "#009E73", "#56B4E9")) +
  facet_grid(~ LOC, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, "lines")
  ) +
  labs(
    x = "INDV",
    y = "Q ancestry proportion",
    fill = "Cluster"
  )

K4 <- read.csv("clumpp_K4.out.csv", header = FALSE, na.strings = "", stringsAsFactors = FALSE)
colnames(K4)[1] <- "INDV"
K4 <- K4[,-2]
colnames(K4)[2] <- "pop4"
colnames(K4)[3] <- "pop1"
colnames(K4)[4] <- "pop2"
colnames(K4)[5] <- "pop3"
K4_melt <- melt(K4, id.vars = "INDV", variable.name = "cluster", value.name = "q")

plot_df <- left_join(K4_melt, meta, by = "INDV")

# Sort INDV by LOC (optional)
plot_df$INDV <- factor(plot_df$INDV, levels = unique(plot_df$INDV[order(plot_df$LOC)]))
plot_df$LOC <- factor(
  plot_df$LOC,
  levels = c("KAL", "SNA", "SEA")
)

# -------------------------
# Plot
# -------------------------
clade_colors <- c("hybrid" = "#BC7FB3", "major" = "#a1111a", "minor" = "#49a4e3")
ggplot(plot_df, aes(x = INDV, y = q, fill = cluster)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("#009E73", "#a1111a", "#49a4e3", "#f46036")) +
  facet_grid(~ LOC, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.2, "lines")
  ) +
  labs(
    x = "INDV",
    y = "Q ancestry proportion",
    fill = "Cluster"
  )

pie_data <- popfile %>%
  group_by(LOC, CLADE) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(LOC) %>%
  mutate(fraction = N / sum(N))

ggplot(pie_data, aes(x = "", y = fraction, fill = CLADE)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  facet_wrap(~ LOC) +
  theme_void() +
  theme(aspect.ratio = 1) +
  scale_fill_manual(values = clade_colors)

#####DEPTH ANALYSIS
clade_colors <- c("major" = "#a1111a", "minor" = "#49a4e3")
#plot depths
popfile_basic <- popfile[popfile$CLADE != "hybrid", ]
# Make violin plot
ggplot(popfile_basic, aes(x = CLADE, y = DEPTH, fill = CLADE)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  labs(title = "Depth distribution of coral samples by clade",
       y = "Depth (m)", x = "Genetic Clade") +
  theme_minimal() +
  scale_fill_manual(values = clade_colors) +
  theme(plot.title = element_text(hjust = 0.5))

# Make box plot
depth_box <- ggplot(popfile_basic, aes(x = CLADE, y = DEPTH, fill = CLADE)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  labs(title = "Depth distribution of coral samples by clade",
       y = "Depth (m)", x = "Genetic Clade") +
  theme_minimal() +
  scale_fill_manual(values = clade_colors) +
  theme(plot.title = element_text(hjust = 0.5))

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/Depth_box_plot.svg",
  width =8,
  height = 8,
  bg = "transparent"
)

print(
  depth_box +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

##### clade map (clean version)

library(tidyterra)
library(ggplot2)
library(randomcoloR)
library(terra)
library(dplyr)
library(sf)

# --- Load data ---
sna_ann <- read.csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_sna/combined_3d_annotations_2123_aligned_renamed_Snakebay_agisoft_clones.csv")
sna_ann <- sna_ann %>%
  mutate(CLADE = ifelse(CLADE == "Hybrid", "Minor", CLADE))
ortho <- rast("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sna_med_2023_ortho_small.tif")
shp   <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sna_lar_01122023_v2_monostands.gpkg")

# no CRS (local system)
terra::crs(ortho) <- NA
terra::crs(shp)   <- NA

# --- 1. Rotate raster ONCE (180° in your case) ---
ortho_rot <- terra::flip(terra::flip(ortho, "vertical"), "horizontal")

# --- 2. Convert to grayscale ---
ortho_gray <- 0.2989 * ortho_rot[[1]] +
  0.5870 * ortho_rot[[2]] +
  0.1140 * ortho_rot[[3]]

# --- 3. Crop ---
bbox <- ext(c(-40, 25, -290, -210))
ortho_gray_crop <- crop(ortho_gray, bbox)

# --- 4. Rotate points using same transformation ---
e <- ext(ortho_rot)
xc <- (e[1] + e[2]) / 2
yc <- (e[3] + e[4]) / 2

sna_ann_rot <- sna_ann %>%
  mutate(
    x = 2 * xc - agisoft_x,
    y = 2 * yc - agisoft_y
  )

# --- 5. Convert points to sf ---
pts_sf <- st_as_sf(sna_ann_rot, coords = c("x", "y"), crs = NA)

# --- 6. Rotate shapefile using SAME transform ---
shp_sf <- st_as_sf(shp)

rotate_geom <- function(geom, xc, yc) {
  coords <- st_coordinates(geom)
  
  coords[,1] <- 2 * xc - coords[,1]
  coords[,2] <- 2 * yc - coords[,2]
  
  st_polygon(list(coords))
}

shp_rot <- st_geometry(shp_sf) |>
  lapply(rotate_geom, xc = xc, yc = yc) |>
  st_sfc()

shp_sf_rot <- st_sf(shp_sf, geometry = shp_rot)

p_map <- ggplot() +
  geom_spatraster(data = ortho_gray_crop, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  
  geom_sf(
    data = shp_sf_rot,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  
  geom_sf(
    data = pts_sf,
    aes(color = factor(CLADE)),
    size = 2
  ) +
  
  coord_sf() +
  theme_minimal()

p_map
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/lineage_map_distribution_SNA.svg",
  width =8,
  height = 12,
  bg = "transparent"
)

print(
  p_map +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

#########kalki
# --- Load data ---
ortho <- rast("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_kal_med_2023_ortho_small.tif")
shp   <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_kal_lar_02122023_v2_monostands.gpkg")
terra::crs(ortho) <- NA
terra::crs(shp)   <- NA
kal_ann <- kal_ann %>%
  mutate(CLADE.y = ifelse(CLADE.y == "Hybrid", "Minor", CLADE.y))

# no CRS (local system)
terra::crs(ortho) <- NA
terra::crs(shp)   <- NA

# --- 2. Convert to grayscale ---
ortho_gray <- 0.2989 * ortho[[1]] +
  0.5870 * ortho[[2]] +
  0.1140 * ortho[[3]]

bbox <- ext(c(-30, 50, 0, 70))
ortho_gray_cropped <- crop(ortho_gray, bbox)
# --- 4. Rotate points using same transformation ---
e <- ext(ortho_gray_cropped)
xc <- (e[1] + e[2]) / 2
yc <- (e[3] + e[4]) / 2

# --- 5. Convert points to sf ---
pts_sf <- st_as_sf(kal_ann, coords = c("agisoft_x", "agisoft_y"), crs = NA)

# --- 6. Rotate shapefile using SAME transform ---
shp_sf <- st_as_sf(shp)

p_map <- ggplot() +
  geom_spatraster(data = ortho_gray_cropped, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  
  geom_sf(
    data = pts_sf,
    aes(color = factor(CLADE.y)),
    size = 2
  ) +
  
  coord_sf() +
  theme_minimal()

p_map
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/lineage_map_distribution_KAL.svg",
  width =8,
  height = 12,
  bg = "transparent"
)

print(
  p_map +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

#################Seaquarium

# --- Load data ---
sea_ann <- read.csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_sea/combined_3d_annotations_2123_aligned_renamed_Seaquarium_agisoft_clones.csv")
sea_ann <- sea_ann %>%
  mutate(CLADE = ifelse(CLADE == "Hybrid", "Minor", CLADE))

ortho <- rast("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sea_med_2021_ortho_small.tif")
shp   <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sea_med_2021_agisoft_monostands.gpkg")

# no CRS (local system)
terra::crs(ortho) <- NA
terra::crs(shp)   <- NA

# --- 2. Convert to grayscale ---
ortho_gray <- 0.2989 * ortho[[1]] +
  0.5870 * ortho[[2]] +
  0.1140 * ortho[[3]]


bbox <- ext(c(-45, 22, -20, 70))
ortho_gray_cropped <- crop(ortho_gray, bbox)

# --- 4. Rotate points using same transformation ---
e <- ext(ortho_gray_cropped)
xc <- (e[1] + e[2]) / 2
yc <- (e[3] + e[4]) / 2

# --- 5. Convert points to sf ---
pts_sf <- st_as_sf(sea_ann, coords = c("agisoft_x", "agisoft_y"), crs = NA)

# --- 6. Rotate shapefile using SAME transform ---
shp_sf <- st_as_sf(shp)

p_map <- ggplot() +
  geom_spatraster(data = ortho_gray_cropped, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  
  geom_sf(
    data = pts_sf,
    aes(color = factor(CLADE)),
    size = 2
  ) +
  
  coord_sf() +
  theme_minimal()

p_map
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/lineage_map_distribution_SEA.svg",
  width =8,
  height = 12,
  bg = "transparent"
)

print(
  p_map +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

######### Curacao map
library(sf)
library(ggplot2)
library(maptiles)
library(marmap)
library(metR)
library(terra)
library(cowplot)
library(rnaturalearth)
library(ggspatial)

# =====================================================
# Bounding box
# =====================================================

bbox <- st_bbox(c(
  xmin = -70,
  xmax = -67,
  ymin = 11,
  ymax = 13
), crs = 4326)

bbox_sf <- st_as_sfc(bbox)

# =====================================================
# Grayscale basemap
# =====================================================

base <- get_tiles(
  bbox_sf,
  provider = "CartoDB.VoyagerNoLabels",
  zoom = 11,
  crop = TRUE
)


# =====================================================
# NOAA bathymetry
# =====================================================

bathy <- getNOAA.bathy(
  lon1 = -70,
  lon2 = -67,
  lat1 = 11,
  lat2 = 13,
  resolution = 0.01
)

bathy_df <- fortify.bathy(bathy)

# =====================================================
# World inset
# =====================================================

world <- ne_countries(
  scale = "large",
  returnclass = "sf"
)

# Bounding box rectangle
bbox_rect <- st_as_sfc(bbox)

inset_map <- ggplot(world) +
  
  geom_sf(
    fill = "grey90",
    color = "grey70",
    linewidth = 0.2
  ) +
  
  geom_sf(
    data = bbox_rect,
    fill = NA,
    color = "red",
    linewidth = 0.8
  ) +
  
  coord_sf(
    xlim = c(-100, 20),
    ylim = c(-10, 50),
    expand = FALSE
  )
  
# =====================================================
# Sample sites
# =====================================================

sites <- data.frame(
  site = c("Playa Kalki", "Snakebay", "Seaquarium"),
  lon = c(
    -(69 + 9/60 + 32.69/3600),
    -(68 + 59/60 + 51.50/3600),
    -(68 + 53/60 + 54.38/3600)
  ),
  lat = c(
    12 + 22/60 + 31.46/3600,
    12 + 8/60 + 20.35/3600,
    12 + 5/60 + 3.72/3600
  )
)

sites_sf <- st_as_sf(
  sites,
  coords = c("lon", "lat"),
  crs = 4326
) 

# =====================================================
# Final plot
# =====================================================


curacao_map <- ggplot() +
  
  geom_raster(
    data = bathy_df,
    aes(x = x, y = y, fill = z),
    alpha = 0.85
  ) +
  
  scale_fill_gradientn(
    colours = c(
      "#edf2f4",
      "#d7e3ea",
      "#bfd3df",
      "#9fb9c9",
      "#6f8fa3"
    ),
    values = scales::rescale(c(
      0,
      -200,
      -1000,
      -2000,
      -4000
    )),
    guide = "none"
  ) +
  
  geom_contour(
    data = bathy_df,
    aes(x = x, y = y, z = z),
    breaks = c(-2000, -1500, -1000, -500, -200, -50),
    color = "#5c7c99",
    linewidth = 0.25,
    alpha = 0.45
  ) +
  
  #geom_contour(
  #  data = bathy_df,
  #  aes(x = x, y = y, z = z),
  #  breaks = c(0),
  #  color = "black",
  #  linewidth = 0.5,
  #  alpha = 0.8
  #) +
  
  layer_spatial(base, alpha = 0.65) +
  
  # ---------------------------------
# Sampling sites
# ---------------------------------

geom_sf(
  data = sites_sf,
  shape = 21,
  fill = "#d62828",
  color = "white",
  size = 3,
  stroke = 0.8
) +
  
  geom_sf_text(
    data = sites_sf,
    aes(label = site),
    size = 3,
    fontface = "bold",
    nudge_y = 0.05,
    nudge_x = 0.10
  ) +
  
  coord_sf(
    xlim = c(-70, -67),
    ylim = c(11, 13),
    expand = FALSE,
    label_axes = list(
      bottom = "E",
      left = "N"
    )
  ) +
  
# -------------------------
# Scale bar
# -------------------------
annotation_scale(
  location = "bl",
  width_hint = 0.25,
  text_cex = 0.8,
  line_width = 0.7
) +
  
  # -------------------------
# North arrow
# -------------------------
annotation_north_arrow(
  location = "tl",
  which_north = "true",
  pad_x = unit(0.2, "cm"),
  pad_y = unit(0.2, "cm"),
  style = north_arrow_fancy_orienteering
) +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_line(
      color = scales::alpha("grey40", 0.3),
      linewidth = 0.2
    ),
    
    axis.title = element_blank(),
    
    axis.text = element_text(
      size = 9,
      color = "black"
    ),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )


final_map <- ggdraw() +
  draw_plot(curacao_map) +
  draw_plot(
    inset_map,
    x = 0.78,
    y = 0.70,
    width = 0.25,
    height = 0.25
    ) +     theme(
      plot.background = element_rect(
        fill = scales::alpha("white", 0.85),
        color = "black"
      )
  )

final_map

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/Curacao_map.svg",
  width =12,
  height = 8,
  bg = "transparent"
)

print(
  final_map +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


#################PCA and DAPC##################################################
#############################################################
#Dependencies
library("here")
library("tidyr")
library("stringr")
library("ggplot2")
library("scales")
library("vcfR")
library("adegenet")
library("RcppCNPy")
library("readr")


#############################################################
# Set the main working directory:
setwd("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Analyses/F1b_overall_pop_stats/")



#############################################################
# Read in vcf
main <- read.vcfR("MADR_reference_f1_all.vcf") #31713 snps

# Convert to genind object
main_genind <- vcfR2genind(main)

#### Extract information from labels & put into genind object ####
names <- row.names(main_genind@tab)

metadata<- read_delim("MADR_reference_popfile_f1_clades.txt", 
                      delim = "\t", escape_double = FALSE, 
                      col_names = FALSE, trim_ws = TRUE)
colnames(metadata) <- c("INDV", "LOC", "CLADE")
strata(main_genind) <- metadata


#### PCA ####
sum(is.na(main_genind$tab))
X <- scaleGen(main_genind, NA.method = "mean")
dim(X)
#[1]   204 64244
class (X)
#[1] "matrix" "array"


# Run the PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=10)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent[1:3]
#[1] 6.73 4.47 1.60


# Convert PCA scores to a data frame
pca_coord <- as.data.frame(pca1$li)
pca_coord$names <- rownames(pca_coord)
pca_coord$clades <- with(metadata,CLADE[match(pca_coord$names,names)])
pca_coord$locations <- with(metadata,LOC[match(pca_coord$names,names)])


# Set colors
colors <- c("minor"="#0496ff","major"="#d7263d", "hybrid"="#7A3E9D")

clean_theme <-  theme(axis.title = element_text(size=20),
                      axis.line = element_line(colour = "black"),
                      axis.text = element_text(size=20),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank(),
                      legend.position="right")






# Plotting PC1 and PC2 by locations
library(ggplot2)
library(ggrepel)

pca_A1A2 <- ggplot(
  pca_coord,
  aes(
    x = Axis1,
    y = Axis2,
    color = clades
  )
) +
  
  # confidence ellipses
  stat_ellipse(
    aes(fill = clades),
    geom = "polygon",
    alpha = 0.12,
    color = NA,
    level = 0.95
  ) +
  
  # points
  geom_point(
    size = 3,
    alpha = 0.9,
    stroke = 0.4
  ) +
  
  # colors
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  
  # axes
  labs(
    x = "PC1 (6.73% variance explained)",
    y = "PC2 (4.47% variance explained)",
    color = "Clade",
    fill = "Clade"
  ) +
  
  # publication theme
  theme_classic(base_size = 12) +
  
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    
    legend.position = c(0.88, 0.88),
    legend.background = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    
    panel.border = element_blank(),
    
    plot.margin = margin(10, 10, 10, 10)
  )

pca_A1A2



# Plotting PC2 and PC3 by locations
pca_A2A3 <- ggplot(pca_coord, aes(x = Axis2, y = Axis3, label = rownames(pca_coord), color = clades)) +
  geom_point(size = 3) +
  geom_text(vjust = 0, hjust = 1.1, size = 2.5) +
  clean_theme+
  scale_color_manual(values = colors)+
  labs(x = "PC2 (4.47%)", y = "PC3 (1.60%)")


pca_A2A3

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/PCA1.svg",
  width =10,
  height = 8,
  bg = "transparent"
)

print(
  pca_A1A2 +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/PCA2.svg",
  width =10,
  height = 8,
  bg = "transparent"
)

print(
  pca_A2A3 +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()
########################################################################################################
# de novo DAPC
## Dependencies ==============================================================
suppressMessages(library("adegenet"))
suppressMessages(library("vcfR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("svglite"))
suppressMessages(library("gridExtra"))


## Main code =================================================================
## Directory
setwd("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Analyses/F1b_overall_pop_stats/")


## Main code =================================================================
## Import with vcfR and popfile (as medatada, no prior)
coral.genlight <- vcfR2genlight(main)          #643 loci with more than 2 alelles, ignored


####################################################################################################
#K=2
grp2 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 2,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#2

pop2 <- grp2$grp
max_PCAs2 <- as.integer(length(pop2) / 3) # as <= N/3 advised
dapc2 <- dapc(coral.genlight, pop2, n.pca = max_PCAs2, n.da = 10)
optimum_score2 <- optim.a.score(dapc2)              #1 PCs

dapc.pop2 <- dapc(coral.genlight, pop2, n.pca = optimum_score2$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord2 <- as.data.frame(dapc.pop2$ind.coord)
dapc_assig2 <- as.data.frame(dapc.pop2$assign)                    #Posterior assignment
dapc_coord2$indv <- rownames(dapc_coord2)
dapc_coord2 <- bind_cols(dapc_coord2, dapc_assig2)
colnames(dapc_coord2)[3] <- "assignment"                    #Hardcoded
dapc_coord2$pop <- with(metadata,CLADE[match(dapc_coord2$indv,INDV)])

#Rename assignment
dapc_coord2$assignment <- as.character(dapc_coord2$assignment)
dapc_coord2$assignment[dapc_coord2$assignment == "1"] <- "major"
dapc_coord2$assignment[dapc_coord2$assignment == "2"] <- "minor"

#Setting site order 
dapc_coord2$pop <- as.factor(dapc_coord2$pop)
levels(dapc_coord2$pop)

colors_clades <- c("hybrid" = "#BC7FB3", "major" = "firebrick", "minor" = "skyblue")
colors_depth <- c("flat" = "skyblue", "slope" = "darkblue")
#Extract variance explained by eigenvectors
percent_dapc.pop2 <- as.data.frame(dapc.pop2$eig/sum(dapc.pop2$eig)*100)



#Common across graphs
clean_theme <-  theme(axis.title = element_text(size=16),
                      axis.line = element_line(colour = "black"),
                      axis.text = element_text(size=16),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.grid = element_blank(),
                      panel.background = element_blank())

# clean_theme2 <-  theme(axis.title = element_text(size=16),
#                        axis.line = element_line(colour = "black"),
#                        axis.text = element_text(size=16),
#                        panel.grid.major = element_blank(),
#                        panel.grid.minor = element_blank(),
#                        panel.border = element_blank(),
#                        panel.grid = element_blank(),
#                        panel.background = element_blank(),
#                        legend.position ="none")
# 
# theme2 <- theme(legend.position = "none",
#                 axis.text = element_text(size = 24),
#                 axis.title = element_text(size = 24),
#                 axis.ticks = element_line(),
#                 panel.grid = element_blank(),
#                 panel.border = element_rect(size = 1, fill = NA),
#                 panel.background = element_blank(),
#                 axis.text.y = element_text(size=24),
#                 axis.text.x = element_text(size=24))


#Graph 
DAPC_K2 <- ggplot(data = dapc_coord2,aes(x=LD1, fill=pop)) + 
  geom_area(stat="bin", alpha=0.6,binwidth = 1)+ 
  scale_fill_manual("Clades", values=colors_clades)+
  clean_theme+
  labs(x = "Discriminant function 1", y = "Density")+
  theme(legend.position=c(0.9,0.9))#+
#  scale_x_continuous(limits = c(-8,4), breaks= c(-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4))

DAPC_K2

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/DAPCK2.svg",
  width =10,
  height = 8,
  bg = "transparent"
)

print(
  DAPC_K2 +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

scatter(dapc.pop2, col=colors_clades)
ggplot(dapc_coord2, aes(x = LD1, y = LD2, color = pop)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = colors_clades) +   
  labs(x = "Discriminant Function 1", 
       y = "Discriminant Function 2", 
       color = "Clade") +
  theme_minimal(base_size = 14)

grp3 <- find.clusters(coral.genlight, stat = "BIC", 
                      n.iter = 100000, n.pca = 120, 
                      choose.n.clust = 3,
                      n.start = 1000, set.seed(2))

#Choose the number of clusters (>=2 - assess BIC): 
#3

pop3 <- grp3$grp

max_PCAs3 <- as.integer(length(pop3) / 3) # as <= N/3 advised
dapc3 <- dapc(coral.genlight, pop3, n.pca = max_PCAs3, n.da = 10)
optimum_score3 <- optim.a.score(dapc3)              #6 PCs
dapc.pop3 <- dapc(coral.genlight, pop3, n.pca = optimum_score3$best, n.da = 10)

#Extract coordinates and add metadata
dapc_coord3 <- as.data.frame(dapc.pop3$ind.coord)
dapc_assig3 <- as.data.frame(dapc.pop3$assign)                    #Posterior assignment
dapc_coord3$indv <- rownames(dapc_coord3)
dapc_coord3 <- bind_cols(dapc_coord3, dapc_assig3)
colnames(dapc_coord3)[4] <- "assignment"                    #Hardcoded
dapc_coord3$pop <- with(metadata,CLADE[match(dapc_coord3$indv,INDV)])
dapc_coord3$depth <- with(metadata,depth[match(dapc_coord3$indv,INDV)])

#Rename assignment
dapc_coord3$assignment <- as.character(dapc_coord3$assignment)
dapc_coord3$assignment[dapc_coord3$assignment == "1"] <- "Cluster 1"
dapc_coord3$assignment[dapc_coord3$assignment == "2"] <- "Cluster 2"
dapc_coord3$assignment[dapc_coord3$assignment == "3"] <- "Cluster 3"

#Setting site order 
dapc_coord3$pop <- as.factor(dapc_coord3$pop)
levels(dapc_coord3$pop)

#Graph 
DAPC_K3 <- ggplot(data = dapc_coord3,aes(x=LD1, fill=pop)) + 
  geom_area(stat="bin", alpha=0.6,binwidth = 1)+ 
  scale_fill_manual("K = 3", values=colors_clades)+
  clean_theme+
  labs(x = "Discriminant function 1", y = "Density")+
  theme(legend.position=c(0.9,0.9))#+
#  scale_x_continuous(limits = c(-8,4), breaks= c(-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4))

DAPC_K3

scatter(dapc.pop3, col=colors_clades)


DAPC_K3_plot <- ggplot(dapc_coord3, aes(x = LD1, y = LD2, color = pop, shape = assignment)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = colors_clades) +   
  labs(x = "Discriminant Function 1", 
       y = "Discriminant Function 2", 
       color = "Depth (m)", 
       shape = "Cluster") +
  theme_minimal(base_size = 14)
DAPC_K3_plot

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/DAPCK3.svg",
  width =10,
  height = 8,
  bg = "transparent"
)

print(
  DAPC_K3_plot +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_2_overall_structure/DAPCK3_scatter.svg",
  width =10,
  height = 8,
  bg = "transparent"
)
  par(
    bg = NA,
    bty = "n",
    mar = c(5, 5, 1, 1)
  )
  scatter(
    dapc.pop4,
    col = colors_clades)
dev.off()
