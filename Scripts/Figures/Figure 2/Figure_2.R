library(readr)
library(dplyr)
library(FNN)
library(purrr)
library(ggplot2)

sna_ann <- sna_ann <- read_csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_sna/combined_3d_annotations_2123_aligned_renamed_Snakebay_agisoft.csv")
sea_ann <- read_csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_sea/combined_3d_annotations_2123_aligned_renamed_Seaquarium_agisoft.csv")
kal_ann <- read_csv("~/Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_kal/combined_3d_annotations_2123_aligned_renamed_Kalki_agisoft.csv")

metadata <- read_table("~/Documents/Inkfish-Phd/MADR_CH1/Analyses/G1_reefscape_genomics/MADR_reference_popfile_G1_full_clonal_groups.txt")
metadata$ID <- sub("_SHA$", "", metadata$INDV)


sna_ann <- sna_ann %>%
  left_join(
    metadata %>% select(ID, clone_lineage, CLADE),
    by = c("genotype" = "ID")
  )

sea_ann <- sea_ann %>%
  left_join(
    metadata %>% select(ID, clone_lineage, CLADE),
    by = c("genotype" = "ID")
  )

kal_ann <- kal_ann %>%
  left_join(
    metadata %>% select(ID, clone_lineage, CLADE),
    by = c("genotype" = "ID")
  )

write.csv(sna_ann, "Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_sna/combined_3d_annotations_2123_aligned_renamed_Snakebay_agisoft_clones.csv", row.names = FALSE)
write.csv(sea_ann, "Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_sea/combined_3d_annotations_2123_aligned_renamed_Seaquarium_agisoft_clones.csv", row.names = FALSE)
write.csv(kal_ann, "Documents/Inkfish-Phd/MADR_CH1/Models/annotations/3D annotations/cur_kal/combined_3d_annotations_2123_aligned_renamed_Kalki_agisoft_clones.csv", row.names = FALSE)


## ORthomosaics
library(tidyterra)
library(ggplot2)
library(randomcoloR)
library(terra)
library(sf)
library(ggrepel)
library(svglite)

distinct_palette <- distinctColorPalette(length(unique(sna_ann$clone_lineage)))
clone_levels <- sort(unique(sna_ann$clone_lineage))
pal <- distinctColorPalette(length(clone_levels))
names(pal) <- clone_levels

pal["0"] <- "black"

clone_scale <- scale_color_manual(
  values = pal,
  name = "Genotype (clone_lineage)"
)

clone_fill <- scale_fill_manual(
  values = pal,
  name = "Genotype (clone_lineage)"
)

ortho <- rast("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sna_med_2023_ortho_small.tif")
shp   <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sna_lar_01122023_v2_monostands.gpkg")
terra::crs(ortho) <- NA
terra::crs(shp)   <- NA

ortho_rotated <- terra::flip(terra::flip(ortho, "vertical"), "horizontal")
ortho_gray <- 0.2989 * ortho_rotated[[1]] +
  0.5870 * ortho_rotated[[2]] +
  0.1140 * ortho_rotated[[3]]
bbox <- ext(c(-55, 45, -300, -210))
ortho_gray_cropped <- crop(ortho_gray, bbox)

ggplot() +
  geom_spatraster(data = ortho_gray_cropped, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none")

e <- ext(ortho_rotated)
xc <- (e[1] + e[2]) / 2  # center x
yc <- (e[3] + e[4]) / 2  # center y

sna_ann_rotated <- sna_ann %>%
  mutate(
    agisoft_x = 2*xc - agisoft_x,
    agisoft_y = 2*yc - agisoft_y
  )

coords <- sna_ann_rotated %>%
  select(agisoft_x, agisoft_y) %>%
  as.matrix()

pts_sf <- st_as_sf(sna_ann_rotated, coords = c("agisoft_x", "agisoft_y"), crs = NA)

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

###Connector lines
k <- 1  # number of neighbours
nn <- get.knn(coords, k = k)

edges <- data.frame(
  from = rep(1:nrow(coords), k),
  to   = as.vector(nn$nn.index)
)

library(igraph)
edges_df <- sna_ann_rotated %>%
  filter(clone_lineage != 0) %>%
  group_split(clone_lineage) %>%
  purrr::map_dfr(function(df) {
    
    coords <- as.matrix(df[, c("agisoft_x","agisoft_y")])
    n <- nrow(coords)
    
    if(n < 2) return(NULL)
    
    # Distance matrix
    d <- as.matrix(dist(coords))
    
    g <- graph_from_adjacency_matrix(
      d,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    
    mst_g <- mst(g, weights = E(g)$weight)
    
    edgelist <- as_edgelist(mst_g)
    edgelist <- matrix(as.numeric(edgelist), ncol = 2)  # 🔑 FIX
    
    data.frame(
      x = coords[edgelist[,1],1],
      y = coords[edgelist[,1],2],
      xend = coords[edgelist[,2],1],
      yend = coords[edgelist[,2],2],
      clone_lineage = df$clone_lineage[1]
    )
  })

### adding cropped model outline
crop_box <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sna_med_dec_2023_cropped_model.shp")
terra::crs(crop_box) <- NA
crop_box_sf <- st_as_sf(crop_box)

crop_box_rot <- st_geometry(crop_box_sf) |>
  lapply(rotate_geom, xc = xc, yc = yc) |>
  st_sfc()

crop_box_sf_rot <- st_sf(crop_box_sf, geometry = crop_box_rot)



p_map <- ggplot() +
  geom_spatraster(data = ortho_gray_cropped, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  
  geom_segment(
    data = edges_df,
    aes(x = x, y = y, xend = xend, yend = yend,
        color = factor(clone_lineage)),
    alpha = 0.5,
    linewidth = 0.7
  ) +
  
  geom_sf(
    data = shp_sf_rot,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = sna_ann_rotated,
    aes(x = agisoft_x, y = agisoft_y, color = factor(clone_lineage)),
    size = 3
  ) +
  
  geom_sf(
    data = crop_box_sf_rot,
    fill = NA,
    color = "orange",
    linewidth = 0.5,
    linetype = "solid"
  ) +
  
  clone_scale +
  coord_sf() +
  labs(x = "X", y = "Y") +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
p_map

####exporting seperatly for inkscape alignment
p_raster <- ggplot() +
  geom_spatraster(
    data = ortho_gray_cropped,
    maxcell = 2e6
  ) +
  scale_fill_gradient(
    low = "black",
    high = "white",
    guide = "none"
  ) +
  coord_sf() +
  theme_void()

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/raster_layer.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  p_raster +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()
p_vector <- ggplot() +
  
  geom_sf(
    data = shp_sf_rot,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = edges_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = factor(clone_lineage)
    ),
    alpha = 0.5
  ) +
  
  geom_point(
    data = sna_ann_rotated,
    aes(
      x = agisoft_x,
      y = agisoft_y,
      color = factor(clone_lineage)
    ),
    size = 2
  ) +
  
  geom_sf(
    data = crop_box_sf_rot,
    fill = NA,
    color = "orange",
    linewidth = 1
  ) +
  
  coord_sf() +
  theme_void()

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/overlay_vectors.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  p_vector +
    clone_scale +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


####pie charts#####
library(ggforce)
shp_sf_rot <- shp_sf_rot %>% mutate(shape_id = paste0(NAME, "_sna"))
shp_sf_rot <- st_zm(shp_sf_rot, drop = TRUE, what = "ZM")
centroids <- shp_sf_rot %>%
  st_centroid() %>%
  mutate(
    x0 = st_coordinates(.)[,1],
    y0 = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(shape_id, x0, y0)

stand_comp <- read_delim("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_S2_monostand_genotypic_diversity/monostand_comp.csv", delim = ";")
stand_comp_sna <- stand_comp %>%
  left_join(centroids, by = "shape_id") %>%
  mutate(across(
    c(angle_start, angle_end, radius, total),
    ~ as.numeric(gsub(",", ".", .))
  )) %>%
  filter(!is.na(x0))

scale_factor <- 1
max_radius <- max(stand_comp_sna$radius, na.rm = TRUE)
stand_comp_sna <- stand_comp_sna %>%
  filter(!is.na(x0)) %>%
  arrange(total) %>%
  mutate(
    x_pos = seq_along(total)
  )

labels_df <- stand_comp_sna %>%
  group_by(shape_id) %>%
  summarise(
    total = first(total),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", total))

label_size = 3

pie_vector <- ggplot(stand_comp_sna) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = 1,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  geom_text(
    data = labels_df,
    aes(x = 0, y = -max_radius * 1.1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_fixed(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  clone_fill +
  theme_void() +
  theme(
    strip.text = element_text(size = label_size * 2.8, family = "sans"),
    legend.position = "none"
  )
pie_vector
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Snakebay/pie_vectors.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  pie_vector +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

#########Pie charts based on size:
pie_vector <- ggplot(stand_comp_sea) +
  
  geom_point(
    aes(x = Inf, y = Inf, size = total),
    alpha = 0
  ) +
  
  geom_arc_bar(
    aes(
      x0 = x_pos, y0 = 0,
      r0 = 0,
      r = radius,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  geom_text(
    data = labels_df,
    aes(x = 0, y = -max_radius * 1.1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_fixed(
    xlim = c(-max_radius, max_radius),
    ylim = c(-max_radius * 1.2, max_radius)
  ) +
  scale_size_continuous(
    name = "Stand Size (n samples)",
    range = c(0.5, 3)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  clone_fill +
  theme_void() +
  theme(
    strip.text = element_text(size = label_size * 2.8, family = "sans"),
    legend.position = "none"
  )
pie_vector

####Patches
patch_comp <- read_delim("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_S2_monostand_genotypic_diversity/patch_comp.csv", delim = ";")
patch_comp_sna <- patch_comp %>%
  left_join(centroids, by = "shape_id") %>%
  mutate(across(
    c(angle_start, angle_end, radius, total),
    ~ as.numeric(gsub(",", ".", .))
  ))

scale_factor <- 1
max_radius <- max(patch_comp_sna$radius, na.rm = TRUE)
patch_comp_sna <- patch_comp_sna %>%
  filter(!is.na(x0)) %>%
  arrange(total) %>%
  mutate(
    x_pos = seq_along(total)
  )

labels_df <- patch_comp_sna %>%
  group_by(shape_id) %>%
  summarise(
    total = first(total),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", total))

label_size = 3

pie_vector <- ggplot(patch_comp_sna) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = 1,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  geom_text(
    data = labels_df,
    aes(x = 0, y = -max_radius * 1.1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_fixed(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  clone_fill +
  theme_void() +
  theme(
    strip.text = element_text(size = label_size * 2.8, family = "sans"),
    legend.position = "none"
  )
pie_vector
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Snakebay/pie_vectors_patches.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  pie_vector +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


###overall cover
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

total_survey_area <- 2475
total_cropped_area <- 450

area_by_class_sna <- subset(area_by_class, Site == "Snakebay")
area_by_class_sna_cropped <- subset(area_by_class_cropped, Site == "Snakebay")

fill_cover_area <- scale_fill_manual(
  values = c(
    "stand" = "#C96B00",
    "patch" = "#E9A03B",
    "individual" = "#F4D06F",
    "other" = "darkgrey"
  )
)

make_pie_df <- function(df, total_area) {
  
  df %>%
    select(size_class, total_area_m2) %>%
    rename(area = total_area_m2) %>%
    
    bind_rows(
      tibble(
        size_class = "other",
        area = total_area - sum(df$total_area_m2)
      )
    ) %>%
    
    mutate(
      fraction = area / total_area,
      ymax = cumsum(fraction),
      ymin = lag(ymax, default = 0),
      label = paste0(round(fraction * 100, 1), "%")
    )
}

pie_df_total <- make_pie_df(area_by_class_sna, total_survey_area)
pie_df_crop  <- make_pie_df(area_by_class_sna_cropped, total_cropped_area)

make_pie_plot <- function(df, title_text) {
  
  ggplot(df) +
    
    geom_rect(
      aes(
        xmin = 0,
        xmax = 1,
        ymin = ymin,
        ymax = ymax,
        fill = size_class
      ),
      color = "white"
    ) +
    
    coord_polar(theta = "y") +
    
    geom_text(
      aes(
        x = 0.5,
        y = (ymin + ymax) / 2,
        label = ifelse(fraction > 0.05, label, "")
      ),
      size = 4
    ) +
    
    fill_cover_area +
    
    labs(title = title_text) +
    
    theme_void() +
    
    theme(
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        size = 14
      )
    )
}

cover_pie_total <- make_pie_plot(
  pie_df_total,
  "Total Survey Area"
)

cover_pie_crop <- make_pie_plot(
  pie_df_crop,
  "Cropped Area"
)

bar_df <- area_by_class_sna %>%
  mutate(
    fraction = total_area_m2 / sum(total_area_m2),
    label = paste0(round(fraction * 100, 1), "%")
  )

cover_bar <- ggplot(bar_df) +
  
  geom_bar(
    aes(
      x = "Cover Composition",
      y = fraction,
      fill = size_class
    ),
    stat = "identity",
    width = 0.6
  ) +
  
  geom_text(
    aes(
      x = "Cover Composition",
      y = cumsum(fraction) - fraction / 2,
      label = ifelse(fraction > 0.03, label, "")
    ),
    size = 4
  ) +
  
  coord_flip() +
  
  fill_cover_area +
  
  scale_y_continuous(
    labels = scales::percent
  ) +
  
  labs(
    title = "Relative Cover Composition",
    x = NULL,
    y = "Percent of Coral Cover"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )

final_cover_plot <-
  (cover_pie_total + cover_pie_crop) /
  cover_bar

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/cover_vectors.svg",
  width = 10,
  height = 10,
  bg = "transparent"
)

print(
  final_cover_plot &
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

#######################################################################################################################
###KALKI#####
distinct_palette <- distinctColorPalette(length(unique(kal_ann$clone_lineage)))
clone_levels <- sort(unique(kal_ann$clone_lineage))
pal <- distinctColorPalette(length(clone_levels))
names(pal) <- clone_levels

pal["0"] <- "black"

clone_scale <- scale_color_manual(
  values = pal,
  name = "Genotype (clone_lineage)"
)
clone_fill <- scale_fill_manual(
  values = pal,
  name = "Genotype (clone_lineage)"
)

ortho <- rast("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_kal_med_2023_ortho_small.tif")
shp   <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_kal_lar_02122023_v2_monostands.gpkg")
terra::crs(ortho) <- NA
terra::crs(shp)   <- NA
ortho_gray <- 0.2989 * ortho[[1]] +
  0.5870 * ortho[[2]] +
  0.1140 * ortho[[3]]
bbox <- ext(c(-30, 50, 0, 70))
ortho_gray_cropped <- crop(ortho_gray, bbox)

ggplot() +
  geom_spatraster(data = ortho_gray, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none")

e <- ext(ortho_gray)
xc <- (e[1] + e[2]) / 2  # center x
yc <- (e[3] + e[4]) / 2  # center y


coords <- kal_ann %>%
  select(agisoft_x, agisoft_y) %>%
  as.matrix()

shp_sf <- st_as_sf(shp)

###Connector lines
k <- 1  # number of neighbours
nn <- get.knn(coords, k = k)

edges <- data.frame(
  from = rep(1:nrow(coords), k),
  to   = as.vector(nn$nn.index)
)

edges_df <- kal_ann %>%
  filter(clone_lineage != 0) %>%
  group_split(clone_lineage) %>%
  purrr::map_dfr(function(df) {
    
    coords <- as.matrix(df[, c("agisoft_x", "agisoft_y")])
    
    if(nrow(coords) < 2) return(NULL)
    
    n_unique <- nrow(unique(as.data.frame(coords)))
    
    if(n_unique < 2) return(NULL)
    
    d <- as.matrix(dist(coords))
    
    g <- graph_from_adjacency_matrix(
      d,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    
    mst_g <- mst(g, weights = E(g)$weight)
    
    edgelist <- as_edgelist(mst_g)
    
    # 🚨 SECOND SAFETY CHECK (this fixes your crash)
    if(length(edgelist) == 0) return(NULL)
    
    edgelist <- matrix(as.numeric(edgelist), ncol = 2)
    
    data.frame(
      x = coords[edgelist[,1],1],
      y = coords[edgelist[,1],2],
      xend = coords[edgelist[,2],1],
      yend = coords[edgelist[,2],2],
      clone_lineage = df$clone_lineage[1]
    )
  })

### adding cropped model outline
crop_box <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_kal_lar_02122023_v2_cropped_model.shp")
terra::crs(crop_box) <- NA
crop_box_sf <- st_as_sf(crop_box)

p_map <- ggplot() +
  geom_spatraster(data = ortho_gray, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  
  geom_segment(
    data = edges_df,
    aes(x = x, y = y, xend = xend, yend = yend,
        color = factor(clone_lineage)),
    alpha = 0.5,
    linewidth = 0.7
  ) +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = kal_ann,
    aes(x = agisoft_x, y = agisoft_y, color = factor(clone_lineage)),
    size = 3
  ) +
  
  geom_sf(
    data = crop_box_sf,
    fill = NA,
    color = "orange",
    linewidth = 0.5,
    linetype = "solid"
  ) +
  
  clone_scale +
  coord_sf() +
  labs(x = "X", y = "Y") +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
p_map

####exporting seperatly for inkscape alignment
p_raster <- ggplot() +
  geom_spatraster(
    data = ortho_gray,
    maxcell = 2e6
  ) +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = edges_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = factor(clone_lineage)
    ),
    alpha = 0.5
  ) +
  
  geom_point(
    data = kal_ann,
    aes(
      x = agisoft_x,
      y = agisoft_y,
      color = factor(clone_lineage)
    ),
    size = 2
  ) +
  scale_fill_gradient(
    low = "black",
    high = "white",
    guide = "none"
  ) +
  coord_sf() +
  theme_void()

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Kalki/raster_layer.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  p_raster +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

p_vector <- ggplot() +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = edges_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = factor(clone_lineage)
    ),
    alpha = 0.5
  ) +
  
  geom_point(
    data = kal_ann,
    aes(
      x = agisoft_x,
      y = agisoft_y,
      color = factor(clone_lineage)
    ),
    size = 2
  ) +
  
  coord_sf() +
  theme_void()

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Kalki/overlay_vectors.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  p_vector +
    clone_scale +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


####pie charts#####
shp_sf <- shp_sf %>% mutate(shape_id = paste0(NAME, "_kal"))
shp_sf <- st_zm(shp_sf, drop = TRUE, what = "ZM")
centroids <- shp_sf %>%
  st_centroid() %>%
  mutate(
    x0 = st_coordinates(.)[,1],
    y0 = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(shape_id, x0, y0)

stand_comp <- read_delim("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_S2_monostand_genotypic_diversity/monostand_comp.csv", delim = ";")
stand_comp_kal <- stand_comp %>%
  left_join(centroids, by = "shape_id") %>%
  mutate(across(
    c(angle_start, angle_end, radius, total),
    ~ as.numeric(gsub(",", ".", .))
  ))

scale_factor <- 1
max_radius <- max(stand_comp_kal$radius, na.rm = TRUE)
stand_comp_kal <- stand_comp_kal %>%
  filter(!is.na(x0)) %>%
  arrange(total) %>%
  mutate(
    x_pos = seq_along(total)
  )

labels_df <- stand_comp_kal %>%
  group_by(shape_id) %>%
  summarise(
    total = first(total),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", total))

label_size = 3

pie_vector <- ggplot(stand_comp_kal) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = 1,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  geom_text(
    data = labels_df,
    aes(x = 0, y = -max_radius * 1.1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_fixed(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  clone_fill +
  theme_void() +
  theme(
    strip.text = element_text(size = label_size * 2.8, family = "sans"),
    legend.position = "none"
  )
pie_vector
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Kalki/pie_vectors.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  pie_vector +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

#######Patches
patch_comp_kal <- patch_comp %>%
  left_join(centroids, by = "shape_id") %>%
  mutate(across(
    c(angle_start, angle_end, radius, total),
    ~ as.numeric(gsub(",", ".", .))
  ))

scale_factor <- 1
max_radius <- max(patch_comp_kal$radius, na.rm = TRUE)
patch_comp_kal <- patch_comp_kal %>%
  filter(!is.na(x0)) %>%
  arrange(total) %>%
  mutate(
    x_pos = seq_along(total)
  )

labels_df <- patch_comp_kal %>%
  group_by(shape_id) %>%
  summarise(
    total = first(total),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", total))

label_size = 3

pie_vector <- ggplot(patch_comp_kal) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = 1,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  geom_text(
    data = labels_df,
    aes(x = 0, y = -max_radius * 1.1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_fixed(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  clone_fill +
  theme_void() +
  theme(
    strip.text = element_text(size = label_size * 2.8, family = "sans"),
    legend.position = "none"
  )
pie_vector
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Kalki/pie_vectors_patches.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  pie_vector +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


###overall cover

total_survey_area <- 2300
total_cropped_area <- 450

area_by_class_kal <- subset(area_by_class, Site == "Playa Kalki")
area_by_class_kal_cropped <- subset(area_by_class_cropped, Site == "Playa Kalki")

fill_cover_area <- scale_fill_manual(
  values = c(
    "stand" = "#C96B00",
    "patch" = "#E9A03B",
    "individual" = "#F4D06F",
    "other" = "darkgrey"
  )
)

pie_df_total <- make_pie_df(area_by_class_kal, total_survey_area)
pie_df_crop  <- make_pie_df(area_by_class_kal_cropped, total_cropped_area)

make_pie_plot <- function(df, title_text) {
  
  ggplot(df) +
    
    geom_rect(
      aes(
        xmin = 0,
        xmax = 1,
        ymin = ymin,
        ymax = ymax,
        fill = size_class
      ),
      color = "white"
    ) +
    
    coord_polar(theta = "y") +
    
    geom_text(
      aes(
        x = 0.5,
        y = (ymin + ymax) / 2,
        label = ifelse(fraction > 0.05, label, "")
      ),
      size = 4
    ) +
    
    fill_cover_area +
    
    labs(title = title_text) +
    
    theme_void() +
    
    theme(
      legend.position = "none",
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        size = 14
      )
    )
}

cover_pie_total <- make_pie_plot(
  pie_df_total,
  "Total Survey Area"
)

cover_pie_crop <- make_pie_plot(
  pie_df_crop,
  "Cropped Area"
)

bar_df <- area_by_class_kal %>%
  mutate(
    fraction = total_area_m2 / sum(total_area_m2),
    label = paste0(round(fraction * 100, 1), "%")
  )

cover_bar <- ggplot(bar_df) +
  
  geom_bar(
    aes(
      x = "Cover Composition",
      y = fraction,
      fill = size_class
    ),
    stat = "identity",
    width = 0.6
  ) +
  
  geom_text(
    aes(
      x = "Cover Composition",
      y = cumsum(fraction) - fraction / 2,
      label = ifelse(fraction > 0.03, label, "")
    ),
    size = 4
  ) +
  
  coord_flip() +
  
  fill_cover_area +
  
  scale_y_continuous(
    labels = scales::percent
  ) +
  
  labs(
    title = "Relative Cover Composition",
    x = NULL,
    y = "Percent of Coral Cover"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )

final_cover_plot <-
  (cover_pie_total + cover_pie_crop) /
  cover_bar

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Kalki/cover_vectors.svg",
  width = 10,
  height = 10,
  bg = "transparent"
)

print(
  final_cover_plot &
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

#######################################################################################################
###Sea Aquarium#####
# plotting hierarchical clustering tree
distinct_palette <- distinctColorPalette(length(unique(sea_ann$clone_lineage)))
clone_levels <- sort(unique(sea_ann$clone_lineage))
pal <- distinctColorPalette(length(clone_levels))
names(pal) <- clone_levels

pal["0"] <- "black"

clone_scale <- scale_color_manual(
  values = pal,
  name = "Genotype (clone_lineage)"
)
clone_fill <- scale_fill_manual(
  values = pal,
  name = "Genotype (clone_lineage)"
)


ortho <- rast("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sea_med_2021_ortho_small.tif")
shp   <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sea_med_2021_agisoft_monostands.gpkg")
terra::crs(ortho) <- NA
terra::crs(shp)   <- NA
ortho_gray <- 0.2989 * ortho[[1]] +
  0.5870 * ortho[[2]] +
  0.1140 * ortho[[3]]
bbox <- ext(c(-45, 22, -20, 70))
ortho_gray_cropped <- crop(ortho_gray, bbox)

e <- ext(ortho_gray)
xc <- (e[1] + e[2]) / 2  # center x
yc <- (e[3] + e[4]) / 2  # center y

coords <- sea_ann%>%
  select(agisoft_x, agisoft_y) %>%
  as.matrix()

shp_sf <- st_as_sf(shp)

k <- 1  # number of neighbours
nn <- get.knn(coords, k = k)

edges <- data.frame(
  from = rep(1:nrow(coords), k),
  to   = as.vector(nn$nn.index)
)

edges_df <- sea_ann %>%
  filter(clone_lineage != 0) %>%
  group_split(clone_lineage) %>%
  purrr::map_dfr(function(df) {
    
    coords <- as.matrix(df[, c("agisoft_x","agisoft_y")])
    n <- nrow(coords)
    
    if(n < 2) return(NULL)
    
    d <- as.matrix(dist(coords))
    
    g <- graph_from_adjacency_matrix(
      d,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )
    
    mst_g <- mst(g, weights = E(g)$weight)
    
    edgelist <- as_edgelist(mst_g)
    
    if(length(edgelist) == 0) return(NULL)
    
    edgelist <- matrix(as.numeric(edgelist), ncol = 2)
    
    data.frame(
      x = coords[edgelist[,1],1],
      y = coords[edgelist[,1],2],
      xend = coords[edgelist[,2],1],
      yend = coords[edgelist[,2],2],
      clone_lineage = df$clone_lineage[1]
    )
  })

### adding cropped model outline
crop_box <- terra::vect("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Models/orthos/cur_sea_lar_2021_v2_cropped_model.shp")
terra::crs(crop_box) <- NA
crop_box_sf <- st_as_sf(crop_box)

p_map <- ggplot() +
  geom_spatraster(data = ortho_gray, maxcell = 2e6) +
  scale_fill_gradient(low = "black", high = "white", guide = "none") +
  
  geom_segment(
    data = edges_df,
    aes(x = x, y = y, xend = xend, yend = yend,
        color = factor(clone_lineage)),
    alpha = 0.5,
    linewidth = 0.7
  ) +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = sea_ann,
    aes(x = agisoft_x, y = agisoft_y, color = factor(clone_lineage)),
    size = 3
  ) +
  
  geom_sf(
    data = crop_box_sf,
    fill = NA,
    color = "orange",
    linewidth = 0.5,
    linetype = "solid"
  ) +
  
  clone_scale +
  coord_sf() +
  labs(x = "X", y = "Y") +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))
p_map

####exporting seperatly for inkscape alignment
p_raster <- ggplot() +
  geom_spatraster(
    data = ortho_gray,
    maxcell = 2e6
  ) +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = edges_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = factor(clone_lineage)
    ),
    alpha = 0.5
  ) +
  
  geom_point(
    data = sea_ann,
    aes(
      x = agisoft_x,
      y = agisoft_y,
      color = factor(clone_lineage)
    ),
    size = 2
  ) +
  scale_fill_gradient(
    low = "black",
    high = "white",
    guide = "none"
  ) +
  coord_sf() +
  theme_void()

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Seaquarium/raster_layer.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  p_raster +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

p_vector <- ggplot() +
  
  geom_sf(
    data = shp_sf,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_segment(
    data = edges_df,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = factor(clone_lineage)
    ),
    alpha = 0.5
  ) +
  
  geom_point(
    data = sea_ann,
    aes(
      x = agisoft_x,
      y = agisoft_y,
      color = factor(clone_lineage)
    ),
    size = 2
  ) +
  
  coord_sf() +
  theme_void()

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Seaquarium/overlay_vectors.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  p_vector +
    clone_scale +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

####pie charts#####
shp_sf <- shp_sf %>% mutate(shape_id = paste0(NAME, "_sea"))
shp_sf <- st_zm(shp_sf, drop = TRUE, what = "ZM")
centroids <- shp_sf %>%
  st_centroid() %>%
  mutate(
    x0 = st_coordinates(.)[,1],
    y0 = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  select(shape_id, x0, y0)

stand_comp <- read_delim("/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_S2_monostand_genotypic_diversity/monostand_comp.csv", delim = ";")
stand_comp_sea <- stand_comp %>%
  left_join(centroids, by = "shape_id") %>%
  mutate(across(
    c(angle_start, angle_end, radius, total),
    ~ as.numeric(gsub(",", ".", .))
  ))

scale_factor <- 1
max_radius <- max(stand_comp_sea$radius, na.rm = TRUE)
stand_comp_sea <- stand_comp_sea %>%
  filter(!is.na(x0)) %>%
  arrange(total) %>%
  mutate(
    x_pos = seq_along(total)
  )

labels_df <- stand_comp_sea %>%
  group_by(shape_id) %>%
  summarise(
    total = first(total),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", total))

label_size = 3

pie_vector <- ggplot(stand_comp_sea) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = 1,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  geom_text(
    data = labels_df,
    aes(x = 0, y = -max_radius * 1.1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_fixed(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  clone_fill +
  theme_void() +
  theme(
    strip.text = element_text(size = label_size * 2.8, family = "sans"),
    legend.position = "none"
  )
pie_vector
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Seaquarium/pie_vectors.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  pie_vector +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

#####Patches
patch_comp_sea <- patch_comp %>%
  left_join(centroids, by = "shape_id") %>%
  mutate(across(
    c(angle_start, angle_end, radius, total),
    ~ as.numeric(gsub(",", ".", .))
  ))

scale_factor <- 1
max_radius <- max(patch_comp_sea$radius, na.rm = TRUE)
patch_comp_sea <- patch_comp_sea %>%
  filter(!is.na(x0)) %>%
  arrange(total) %>%
  mutate(
    x_pos = seq_along(total)
  )

labels_df <- patch_comp_sea %>%
  group_by(shape_id) %>%
  summarise(
    total = first(total),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", total))

label_size = 3

pie_vector <- ggplot(patch_comp_sea) +
  geom_arc_bar(
    aes(
      x0 = 0, y0 = 0,
      r0 = 0,
      r = 1,
      start = angle_start,
      end = angle_end,
      fill = clone_lineage
    ),
    color = "white",
    linewidth = 0.3
  ) +
  facet_wrap(~ shape_id) +
  geom_text(
    data = labels_df,
    aes(x = 0, y = -max_radius * 1.1, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  coord_fixed(
    xlim = c(-1.1, 1.1),
    ylim = c(-1.1, 1.1)
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  clone_fill +
  theme_void() +
  theme(
    strip.text = element_text(size = label_size * 2.8, family = "sans"),
    legend.position = "none"
  )
pie_vector
svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Seaquarium/pie_vectors_patches.svg",
  width = 8,
  height = 8,
  bg = "transparent"
)

print(
  pie_vector +
    theme(
      panel.background = element_blank(),
      plot.background  = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()


###overall cover

total_survey_area <- 560
total_cropped_area <- 450

area_by_class_sea <- subset(area_by_class, Site == "Sea Aquarium")
area_by_class_sea_cropped <- subset(area_by_class_cropped, Site == "Sea Aquarium")

fill_cover_area <- scale_fill_manual(
  values = c(
    "stand" = "#C96B00",
    "patch" = "#E9A03B",
    "individual" = "#F4D06F",
    "other" = "darkgrey"
  )
)

pie_df_total <- make_pie_df(area_by_class_sea, total_survey_area)
pie_df_crop  <- make_pie_df(area_by_class_sea_cropped, total_cropped_area)

cover_pie_total <- make_pie_plot(
  pie_df_total,
  "Total Survey Area"
)

cover_pie_crop <- make_pie_plot(
  pie_df_crop,
  "Cropped Area"
)

bar_df <- area_by_class_sea %>%
  mutate(
    fraction = total_area_m2 / sum(total_area_m2),
    label = paste0(round(fraction * 100, 1), "%")
  )

cover_bar <- ggplot(bar_df) +
  
  geom_bar(
    aes(
      x = "Cover Composition",
      y = fraction,
      fill = size_class
    ),
    stat = "identity",
    width = 0.6
  ) +
  
  geom_text(
    aes(
      x = "Cover Composition",
      y = cumsum(fraction) - fraction / 2,
      label = ifelse(fraction > 0.03, label, "")
    ),
    size = 4
  ) +
  
  coord_flip() +
  
  fill_cover_area +
  
  scale_y_continuous(
    labels = scales::percent
  ) +
  
  labs(
    title = "Relative Cover Composition",
    x = NULL,
    y = "Percent of Coral Cover"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      size = 14
    )
  )

final_cover_plot <-
  (cover_pie_total + cover_pie_crop) /
  cover_bar

svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Seaquarium/cover_vectors.svg",
  width = 10,
  height = 10,
  bg = "transparent"
)

print(
  final_cover_plot &
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()

####################################################################################
Creating the map:
  library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)

library(geodata)
library(sf)

curacao <- gadm(country = "CUW", level = 0, path = tempdir())
curacao <- st_as_sf(curacao)
curacao <- st_make_valid(curacao)
geom_sf(data = curacao, fill = "gray95", color = "gray40", size = 0.4)
# Location color mapping
location_colors <- c("SNA" = "#2e294e", "SEA" = "#1b998b", "KAL" = "#9b2226")

# Sample site data with corresponding location codes
samples <- data.frame(
  SiteName = c("Snakebay", "Sea Aquarium", "Playa Kalki"),
  LocationCode = c("SNA", "SEA", "KAL"),
  Longitude = c(-68.99730274399919, -68.89843572928432, -69.15876448304232),  # negative = west
  Latitude = c(12.139165840463166, 12.084367178152423, 12.375245235353116)    # positive = north
)

# Convert to spatial points
samples_sf <- st_as_sf(samples, coords = c("Longitude", "Latitude"), crs = 4326)


# Plot
map <- ggplot() +
  geom_sf(data = curacao, fill = "gray95", color = "gray40", size = 0.4) +
  geom_sf(data = samples_sf, aes(fill = LocationCode), 
          shape = 21, size = 5, color = "white", stroke = 1) +
  geom_text(
    data = samples, 
    aes(x = Longitude, y = Latitude, label = SiteName), 
    size = 3.5,               # slightly smaller text
    fontface = "plain",       # not bold, lighter style
    hjust = -0.1,             # push text just right of the point
    nudge_x = 0.02            # adds spacing to the right
  ) +
  coord_sf(xlim = c(-69.4, -68.6), ylim = c(11.9, 12.45), expand = FALSE) +
  scale_fill_manual(values = location_colors, name = "Sample site") +
  annotation_scale(location = "bl", width_hint = 0.3, line_width = 0.6) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_fancy_orienteering, 
                         height = unit(1, "cm"), width = unit(1, "cm")) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#e6f2fa", color = NA),  # Softer blue ocean
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.key = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 10)
  ) +
  labs(title = "Sample Sites on Curaçao", 
       x = "Longitude", 
       y = "Latitude")
map 



svglite::svglite(
  "/Users/dvan216/Documents/Inkfish-Phd/MADR_CH1/Figures/Figure_3_genotype_spread/Curacao_map.svg",
  width = 10,
  height = 10,
  bg = "transparent"
)

print(
  map &
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
)

dev.off()
