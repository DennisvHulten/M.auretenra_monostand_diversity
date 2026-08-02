# Data Availability

This directory contains the processed genomic datasets used for the population genomic analyses presented in **van Hulten et al. (in preparation)**.

The datasets are provided in two file formats:

- **VCF (`.vcf.gz`)** – Filtered SNP genotype data in Variant Call Format, used for population genomic analyses including population structure, genetic diversity, and lineage assignment.
- **IBS (`.ibs`)** – Pairwise Identity-by-State matrix generated using **ANGSD**, used for clone identification and the assignment of samples to clonal groups.

For both file types, the same naming convention is used:

| File suffix | Description |
|-------------|-------------|
| `*_all` | Complete dataset containing all samples included in the study. |
| `*_major` | Subset containing only samples assigned to the major genetic lineage. |
| `*_minor` | Subset containing only samples assigned to the minor genetic lineage. |

The accompanying `Metadata` directory contains detailed information for each sample, including sample identifiers, geographic coordinates, spatial positions within the photogrammetric models, lineage assignments, and clonal group assignments.