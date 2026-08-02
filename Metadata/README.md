# Data Availability

This directory contains the metadata associated with the data presented in **van Hulten et al. (in preparation)**.

## Metadata Mastersheet

The `Metadata_mastersheet.csv` file contains metadata for every sample included in this study. The variables are described below.

| Variable | Description |
|----------|-------------|
| `targetid` | Unique sequencing identifier assigned by DArT for each sequencing run. Where a sample was sequenced multiple times, reads from sequential runs were concatenated for downstream analyses. |
| `ordernumber` | DArT batch identifier corresponding to the sequencing submission (representing samples submitted in 2022 and 2024, respectively). |
| `genus` | Scientific name of the genus. |
| `species` | Scientific name of the species. |
| `sample_id` | Unique identifier assigned to each sample. Standard notation follows `DHxxxx_MMIR_{SITE}`, with technical replicates denoted as `DXxxxx_MMIR_{SITE}`. |
| `site` | Sampling site of origin. |
| `annotation` | Identifier assigned to the sample annotation within the Viscore 3D models. These identifiers are unique within each site but not across sites. |
| `orig_x`, `orig_y`, `orig_z` | Original (untransformed) sample coordinates as annotated in the Viscore `.ply` model. |
| `world_x`, `world_y`, `world_z` | Transformed sample coordinates aligned across both survey years. The `z`-axis corresponds to real-world depth, while the `x`-axis is oriented parallel to the shoreline. |
| `agisoft_x`, `agisoft_y`, `agisoft_z` | Original sample coordinates in Agisoft Metashape Professional, used for alignment with polygon annotations on the orthomosaics. |
| `clone_group` | Clonal group assignment based on ANGSD analyses. A value of `0` indicates a unique genotype that was not assigned to a clonal group. |
| `lineage` | Genetic lineage assignment based on concordant results from K-clustering, PCA, and DAPC analyses. |
| `shape_id` | Identifier linking each sample to a polygon delineated on the orthomosaics, representing aggregations of *Madracis auretenra* branches. |
| `lat` | Latitude of the sampling site (decimal degrees). |
| `lon` | Longitude of the sampling site (decimal degrees). |

## Monostand Metrics

The monostand metrics files contain demographic measurements for all *Madracis auretenra* colonies included in this study. Separate files are provided for:

- the complete dataset, and
- the standardized survey areas used for comparative analyses.