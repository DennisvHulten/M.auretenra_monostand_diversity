# *Madracis auretenra* RADseq/DArTseq Population Genomics Analysis

Analysis notebook documenting the genomic workflow used to accompany **van Hulten (in prep.)** — a population genomics study of the scleractinian coral *Madracis auretenra* using DArTseq (reduced representation SNP data).

This README walks through every stage of the analysis, from raw read processing through to reefscape-level genetic diversity assessment, along with the exact commands, scripts, and parameter files used at each step.

## Overview

Data were collected across two expeditions to Curaçao: an initial expedition in 2021 (493 samples) and a follow-up in 2023 (184 samples), for a total of **681 individuals**. All samples were shipped to Diversity Arrays Technology (Canberra, Australia) for DArTseq sequencing. Voucher/backup samples are held at the California Academy of Sciences, San Francisco, USA.

## Table of Contents

- [A1a. Data Preparation](#a1a-data-preparation--concatenating-raw-fastq-files)
- [A1b. Raw Read Quality Control](#a1b-raw-read-quality-control-fastqc-multiqc-trim-galore)
- [B1. Reference-Based Assembly with ipyrad](#b1-reference-based-assembly-with-ipyrad)
- [C1. VCF Post-Processing](#c1-vcf-post-processing)
  - [C1a. Cleaning Up Sample Names](#c1a-cleaning-up-sample-names-in-the-vcf)
  - [C1b. VCF Quality Control](#c1b-vcf-quality-control)
- [C2. Removing Symbiont Contamination](#c2-removing-symbiont-contamination)
- [C3. Screening for Other Contamination](#c3-screening-for-other-contamination-blast-against-ncbi-nt)
- [D1. Initial Genetic Assessment](#d1-initial-genetic-assessment)
  - [D1a. Removing Underperforming Individuals & SNPs](#d1a-removing-underperforming-individuals--snps)
  - [D1b. Building Initial Neighbour-Joining Trees](#d1b-building-initial-neighbour-joining-trees)
  - [D1c. Genetic Structure, Clone Detection & STRUCTURE Analysis](#d1c-genetic-structure-clone-detection--structure-analysis)
- [E1. Clone Detection](#e1-clone-detection)
  - [E1a. Initial Genetic Similarity Comparison](#e1a-initial-genetic-similarity-comparison)
  - [E1b. ANGSD IBS-Based Clonality Assessment](#e1b-angsd-ibs-based-clonality-assessment)
- [F1. Lineage Comparison](#f1-lineage-comparison)
  - [F1a. Private and Fixed Alleles](#f1a-private-and-fixed-alleles)
  - [F1b. Population Statistical Analysis](#f1b-population-statistical-analysis-pca--dapc)
  - [F1c. Genetic Distance Between Lineages](#f1c-genetic-distance-between-lineages-hierfstat)
- [F2. Outlier Detection with PCAdapt & RDA](#f2-outlier-detection-with-pcadapt--rda)
- [G1. Reefscape Assessment](#g1-reefscape-assessment)
  - [G1a. Demographic Analysis of Stands](#g1a-demographic-analysis-of-stands)
  - [G1b. Genotypic Diversity of Stands](#g1b-genotypic-diversity-of-stands)
  - [G1c. Genetic Diversity of Stands](#g1c-genetic-diversity-of-stands)

## Software & Versions

| Tool | Version |
|---|---|
| FastQC | v0.12.1 |
| BWA | 0.7.18-r1243-dirty |
| MultiQC | 1.27 |
| ipyrad | 0.9.102 |
| VCFtools | 0.1.16 |
| Trim Galore | 0.6.10-1 |
| BLASTN | 2.16.0+ |
| STRUCTURE | 2.3.4 |
| adegenet (R) | 2.1.10 |
| poppr (R) | 2.9.6 |
| ANGSD | htslib 1.17 (build Aug 27 2023) |
| NewHybrids | 2.0+ Developmental (Jul/Aug 2007) |

## Scripts

Custom scripts referenced throughout this notebook (e.g. `vcf_missing_data.py`, `pyrad2fasta.py`, `vcf_gdmatrix.py`, `gdmatrix2tree.py`) are available at [pimbongaerts/radseq](https://github.com/pimbongaerts/radseq).

## Images

Figures referenced below are stored in the [`Images/`](./Images) folder of this repository.

---

## A1a. Data Preparation — Concatenating raw FASTQ files

**Concatenate fastq files for each individual.**

```bash
#!/bin/bash
# Define input CSV file
csv_file="MADR_targets_22_24.csv"
FASTQ_path="/home/deepcat/rawseq/MADR_dart/MADR_full"
FASTQ_out="/home/deepcat/rawseq/MADR_dart/MADR_merge"
# Loop through unique genotypes in the CSV file
awk -F',' 'NR>1 {print $5}' "$csv_file" | sort | uniq |
while read genotype; do
if [ -z "$genotype" ]; then
continue # Skip empty genotype values
fi
# Get all targetid values for the current genotype
target_ids=$(awk -F',' -v g="$genotype" 'NR>1 && $5 == g {print
$1}' "$csv_file")
# Create an array of FASTQ files corresponding to these target_ids
fastq_files=() # Initialize an empty array
for target_id in $target_ids; do
file_path="${FASTQ_path}/${target_id}.FASTQ.gz"
if [ -f "$file_path" ]; then
fastq_files+=("$file_path") # Add the file path to the array
fi
done
# Skip if no FASTQ files are found
if [ ${#fastq_files[@]} -eq 0 ]; then
echo "No FASTQ files found for genotype $genotype"
continue
fi
# Combine the FASTQ files for this genotype
output_file="${FASTQ_out}/${genotype}.FASTQ.gz"
echo "Combining files for genotype $genotype into $output_file"
cat "${fastq_files[@]}" > "$output_file"
done
```

**Total samples**

```bash
ls /home/deepcat/rawseq/MADR_dart/MADR_merge/ | wc -l
#681 total samples
```

## A1b. Raw Read Quality Control (FastQC, MultiQC, Trim Galore)

```bash
# Run FastQC on all FASTQ files in the directory
% for file in $FASTQ_PATH/*.FASTQ.gz; do fastqc -o
~/Projects/MADR/A1a_QC/fastqc
$file; done
# Collate with MultiQC
$ multiqc fastqc/
# found 681 reports
$ zip multiqc*
```

![](Images/image22.png)

![](Images/image5.png)

![](Images/image57.png)

Most of my sequences looks like this, they are supposed to run parallel... most sequences look similar in the position of the peaks but differ on peak intensity, this could be biologically relevant but could also be contamination?

![](Images/image53.png)

![](Images/image39.png)

And warnings for sequence length dist....

![](Images/image47.png)

![](Images/image9.png)

![](Images/image33.png)

**Trim Galore trimming step**

```bash
trim_galore --cores 8 -q 20 --length 30 --phred33 -o fastq_trimmed/
$FASTQ_PATH/*.FASTQ.gz
```

## B1. Reference-Based Assembly with ipyrad

Assembly was run against a reference genome for *Madracis auretenra* (NCBI BioProject [PRJEB76086](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB76086/), Sanger Institute).

```bash
% ipyrad -n MADR_reference
% nano params-MADR_reference
```

Key parameters set in the resulting `params-MADR_reference.txt` file:

```text
------- ipyrad params file (v.0.9.102) -------------------------------------
MADR_REF                                                      ## [0] [assembly_name]: Assembly name. Used to name output directories for assembly steps
/home/dvanhulten/Projects/MADR/B2_reference_assembly          ## [1] [project_dir]: Project dir (made in curdir if not present)
                                                                ## [2] [raw_fastq_path]: Location of raw non-demultiplexed fastq files
                                                                ## [3] [barcodes_path]: Location of barcodes file
/home/dvanhulten/Projects/MADR/A1b_qc/fastq_trimmed/*.fq.gz    ## [4] [sorted_fastq_path]: Location of demultiplexed/sorted fastq files
reference                                                      ## [5] [assembly_method]: Assembly method (denovo, reference)
~/Projects/MADR/B2_reference_assembly/madr_auretenra_ref.fasta ## [6] [reference_sequence]: Location of reference sequence file
ddrad                                                          ## [7] [datatype]: Datatype (see docs): rad, gbs, ddrad, etc.
TGCAG,                                                         ## [8] [restriction_overhang]: Restriction overhang (cut1,) or (cut1, cut2)
5                                                               ## [9] [max_low_qual_bases]: Max low quality base calls (Q<20) in a read
33                                                              ## [10] [phred_Qscore_offset]: phred Q score offset (33 is default and very standard)
6                                                               ## [11] [mindepth_statistical]: Min depth for statistical base calling
6                                                               ## [12] [mindepth_majrule]: Min depth for majority-rule base calling
10000                                                          ## [13] [maxdepth]: Max cluster depth within samples
0.85                                                           ## [14] [clust_threshold]: Clustering threshold for de novo assembly
0                                                               ## [15] [max_barcode_mismatch]: Max number of allowable mismatches in barcodes
2                                                               ## [16] [filter_adapters]: Filter for adapters/primers (1 or 2=stricter)
35                                                              ## [17] [filter_min_trim_len]: Min length of reads after adapter trim
2                                                               ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences
0.05                                                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus
0.05                                                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus
4                                                               ## [21] [min_samples_locus]: Min # samples per locus for output
0.2                                                             ## [22] [max_SNPs_locus]: Max # SNPs per locus
8                                                               ## [23] [max_Indels_locus]: Max # of indels per locus
0.5                                                             ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus
0, 0, 0, 0                                                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)
0, 0, 0, 0                                                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)
p, s, l, v                                                     ## [27] [output_formats]: Output formats (see docs)
                                                                ## [28] [pop_assign_file]: Path to population assignment file
                                                                ## [29] [reference_as_filter]: Reads mapped to this reference are removed in step 3
```

Run the assembly across all 7 ipyrad steps using an ipcluster of 30 CPUs:

```bash
% ipcluster start -n 30 --daemonize --profile="MADR_ipyrad_reference"; sleep 60
% ipyrad -p params-MADR_reference.txt -s1234567 -c 30 --ipcluster MADR_ipyrad_reference
```

## C1. VCF Post-Processing

### C1a. Cleaning Up Sample Names in the VCF

Remove .FASTQ.gz_trimmed from the sample names

```bash
awk 'BEGIN{OFS="\t"} /^#CHROM/ {for (i=10; i<=NF; i++)
sub(/\.FASTQ.gz_trimmed$/, "", $i)}1' MADR_reference.vcf >
```

MADR_reference_C1a.vcf

**###**

### C1b. VCF Quality Control

```bash
% wc -l MADR_reference_C1a.vcf | awk '{print $1-12}'
# 428.665 #number of SNPs
#vcf_missing_data.py from https://github.com/pimbongaerts/radseq
/home/pbongaerts/Github/radseq/vcf_missing_data.py
```

MADR_reference_C1a.vcf > MADR_reference_C1b_qc.txt

```bash
$ sort -gk 5 MADR_reference_C1b_qc.txt
$ awk '
BEGIN { less_10 = 0; less_30 = 0; less_1000_geno = 0; total_percent = 0; count = 0 }
NR > 1 {                                  # Skip header
    if ($5 < 10)   less_10++;             # Check percentage genotyped
    if ($5 < 30)   less_30++;
    if ($3 < 1000) less_1000_geno++;      # Check GENO column
    total_percent += $5;                  # Add % genotyped to total
    count++;
}
END {
    print "Samples with <10% genotyped: "   less_10;
    print "Samples with <30% genotyped: "   less_30;
    print "Samples with <1000 SNPs (GENO): " less_1000_geno;
    print "Average % genotyped: " total_percent / count;
}

' MADR_reference_C1b_qc.txt

```bash
#total 681 samples
#Samples with <10% genotyped: 19
#Samples with <30% genotyped: 293
#Samples with <1000 SNPs (GENO): 3
#Average % genotyped: 32.9087
```

## C2. Removing Symbiont Contamination

```bash
#pyrad2fasta.py from
[https://github.com/pimbongaerts/radseq](https://github.com/pimbongaerts/radseq)
%
/home/pbongaerts/Github/radseq/[pyrad2fasta.py](http://pyrad2fasta.py)
MADR_reference.loci > MADR_reference.fasta
% grep ">" MADR_reference.fasta | wc -l
# 53.725 number of loci
#Remove hyphens for blastn
% sed 's/-//g' MADR_reference.fasta > MADR_reference_cleaned.fasta
#align with bwa
% bwa mem -t 32 -M /home/deepcat/genomes/breviolum/breviolum
```

MADR_reference_cleaned.fasta > MADR_reference_to_brev.sam

```bash
#Exploring sam output:
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | wc -l
#Number of reads: 53733 (more than in the fasta?)
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | awk '$2
== 0' | wc -l
#Mapped reads forward: 96
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | awk '$2
== 16' | wc -l
#Mapped reads forward: 103
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | awk '$2
== 4' | wc -l
#Unmapped reads: 53526
#Extract a list of succesfully mapped loci with mapping quality of >=20
% /home/pbongaerts/Github/radseq/mapping_get_bwa_matches.py
bwa_outputs/MADR_denovo_to_brev.sam >
bwa_outputs/MADR_denovo_to_brev_q20.txt
% wc -l bwa_outputs/MADR_denovo_to_brev_q20.txt
#69 loci to remove
% /home/pbongaerts/Github/radseq/mapping_get_bwa_matches.py
bwa_outputs/MADR_reference_to_brev.sam >
bwa_outputs/MADR_reference_to_brev_q20.txt
% wc -l bwa_outputs/MADR_reference_to_brev_q20.txt
#15 loci to remove
#Repeat for other genomes
GENOMES=("symbiodinium" "cladocopium" "durusdinium")
GENOME_DIR=/home/deepcat/genomes/
INPUT_FASTA=MADR_denovo_cleaned.fasta
OUTPUT_DIR=bwa_outputs
for GENOME in "${GENOMES[@]}"; do
GENOME_SHORT=${GENOME:0:4}
bwa mem -t 32 -M "${GENOME_DIR}/${GENOME}/${GENOME}"
"$INPUT_FASTA" >
"${OUTPUT_DIR}/MADR_denovo_to_${GENOME_SHORT}.sam"
/home/pbongaerts/Github/radseq/mapping_get_bwa_matches.py
"${OUTPUT_DIR}/MADR_denovo_to_${GENOME_SHORT}.sam" >
"${OUTPUT_DIR}/MADR_denovo_to_${GENOME_SHORT}_q20.txt"
echo "Processing for ${GENOME} completed."
done
cat bwa_outputs/MADR_denovo_*.txt | cut -f1 | sort | uniq >
bwa_outputs/MADR_denovo_sym_loci_to_remove.txt
wc -l bwa_outputs/MADR_denovo_sym_loci_to_remove.txt
#237 loci to remove
#repeat for reference
wc -l bwa_outputs/reference/MADR_reference_sym_loci_to_remove.txt
#55 loci to remove
```

Compare results with blastn method

```bash
**#1. run blastn against the databases**
$ INPUT_FASTA=MADR_reference_cleaned.fasta
$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue
```

bitscore"

```bash
GENOMES=("symbiodinium" "cladocopium" "durusdinium" "breviolum")
$ for GENOME in "${GENOMES[@]}"; do
GENOME_SHORT=${GENOME:0:4}
blastn -query $INPUT_FASTA\
-db "${GENOME_DIR}/${GENOME}/${GENOME}" \
```

-task blastn \ -outfmt "$OUTPUT_FORMAT" 

```bash
-out "**blastn_outputs/**MADR_reference_to_${GENOME_SHORT}.txt"
done
```

**#2. extract matches with E-value lower than 10^-15^**

```bash
$ MAX_E_VALUE="0.000000000000001"
for GENOME in "${GENOMES[@]}"; do
GENOME_SHORT=${GENOME:0:4}
python3 /home/pbongaerts/Github/radseq/mapping_get_blastn_matches.py
"**blastn_outputs/**MADR_reference_to_${GENOME_SHORT}.txt"
$MAX_E_VALUE
done
#Symbiodinium
Matches: 12 | Min.length: 78.0 bp | Min. nident: 66.0 bp | Min.
```

pident: 73.729 %

```bash
#Cladocopium
Matches: 10 | Min.length: 72.0 bp | Min. nident: 62.0 bp | Min.
```

pident: 78.095 %

```bash
#Durusdinium
Matches: 2 | Min.length: 77.0 bp | Min. nident: 73.0 bp | Min.
```

pident: 86.869 %

```bash
#Breviolum
Matches: 15 | Min.length: 65.0 bp | Min. nident: 60.0 bp | Min.
```

pident: 73.171 % Total of 39 loci that match symbiont genome

## C3. Screening for Other Contamination (BLAST against NCBI nt)

```bash
#Start a screen
$ screen -S MADR_denovo_blastn
blastn -query
~/Projects/MADR/C2_remove_symbiont_contamination/reference/MADR_reference_cleaned.fasta
```

-db ../blast_db_12_26_2023/nt -task blastn -evalue 0.0001 -max_target_seqs 10 -outfmt "7 qseqid sseqid length nident pident evalue bitscore staxids stitle" -out

```bash
blastn_outputs/MADR_reference_all.txt -num_threads 32
```

**Repeat for denovo**

Looking at the results

```bash
grep "Query" blastn_outputs/MADR_reference_all.txt | wc -l
#53725 Loci
$ mapping_identify_blast_matches_modAH.py
blastn_outputs/MADR_reference_all.txt 0.0001
$ cut -f 8 blastn_outputs/MADR_reference_all_match0.0001.txt | sort |
uniq -c > MADR_reference_other_loci_hit_summary.txt
#Reference results
30 Annelida
6 Apicomplexa
177 Arthropoda
4 Ascomycota
3 Bacillota
2 Bacteroidota
2 Brachiopoda
4 Bryozoa
237 Chordata
21233 Cnidaria
30 Echinodermata
74 Mollusca
1 Mucoromycota
2 Nematoda
5 NOT_FOUND
1 Peploviricota
1 Planctomycetota
5 Platyhelminthes
3 Porifera
1 Priapulida
5 Pseudomonadota
30 Streptophyta
1 Uroviricota
# 2. Create a list of loci that match phyla outside Cnidaria and
NOT_FOUND
$ grep -v -E "Cnidaria|NOT_FOUND"
blastn_outputs/MADR_reference_all_match0.0001.txt >
```

MADR_reference_other_loci_to_remove.txt

```bash
$ wc -l MADR_reference_other_loci_to_remove.tx
#674 other loci to remove
Remove all contaminants:
# 1. Merge the symbiotns and other contaminants
$ cat MADR_reference_sym_loci_to_remove.txt
MADR_reference_other_loci_to_remove.txt | sort | uniq >
```

MADR_reference_all_loci_to_remove_temp.txt

```bash
$ wc -l MADR_reference_all_loci_to_remove_temp.txt
#729 loci to remove
# for reference selection we need the name of the chrom
awk '{split($2, acc, "[|.]"); print "ENA|" acc[4] "|"
acc[4] ".1";}' MADR_reference_all_loci_to_remove_temp.txt >
```

MADR_reference_all_loci_to_remove.txt

```bash
# 3. Remove from vcf file
$ /home/pbongaerts/Github/radseq/vcf_remove_chrom.py
```

MADR_reference_c3_temp.vcf MADR_reference_all_loci_to_remove.txt > MADR_reference_c3.vcf

## D1. Initial Genetic Assessment

### D1a. Removing Underperforming Individuals & SNPs

Start with qc control of snps and the vcf

```bash
vcf_pos_count_MODref.py reference/MADR_reference_c3.vcf >
ref_pos_count.txt
awk '{count[$2]++} END {for (value in count) print value ": "
count[value]}' ref_pos_count.txt | sort -n\
```

4: 1 5: 82.245 6: 1.855 7: 176 8: 10 9: 1

```bash
# %missing data
% /home/pbongaerts/Github/radseq/vcf_missing_data.py
MADR_reference_C3.vcf > MADR_reference_preformance_C3.txt
% sort -gk 5 -o MADR_reference_preformance_C3.txt
```

MADR_reference_preformance_C3.txt

```bash
#INDIVIDUAL MISS GENO TOTAL % GENOTYPED
#DH0293_MMIR_SNA 418633 0 418633 0.0
#DH0360_MMIR_SNA 418633 0 418633 0.0
#DH0371_MMIR_SNA 418614 19 418633 0.0
#DH0373_MMIR_SNA 400754 17879 418633 4.27
#DH0308_MMIR_SNA 392283 26350 418633 6.29
#DH0277_MMIR_SNA 388740 29893 418633 7.14
#DH0627_MMIR_KAL 388546 30087 418633 7.19
#DH0552_MMIR_SEA 388343 30290 418633 7.24
#DH0553_MMIR_SEA 388260 30373 418633 7.26
#Samples with <10% genotyped: 19
#Samples with <30% genotyped: 293
#Samples with <1000 SNPs (GENO): 3
#Average % genotyped: 32.9087
#less than 1000 snps
```

INDIVIDUAL MISS GENO TOTAL % GENOTYPED

```bash
DH0293_MMIR_SNA 428665 0 428665 0.0
DH0360_MMIR_SNA 428665 0 428665 0.0
DH0371_MMIR_SNA 428649 16 428665 0.0
DH0371_MMIR_SNA_x 428661 4 428665 0.0
# save all names with <1000 SNPS genotyped to remove file
$ awk '$3 < 1000 {print $1}' MADR_reference_preformance_c3.txt >
```

indv_to_remove_reference_c3.txt

```bash
3 individuals to remove
```

**remove underperforming snps and individuals**

```bash
$ vcftools --vcf MADR_reference_c3.vcf --remove
```

indv_to_remove_reference_C3.txt --max-missing 0.8 --mac 1 --minDP 8 --recode --stdout > MADR_reference_d1_strict.vcf

```bash
#After filtering, kept 678 out of 681 Individuals
#After filtering, kept 34.466 out of a possible 418633 Sites
$ /home/pbongaerts/Github/radseq/vcf_missing_data.py
MADR_reference_d1_temp.vcf > MADR_reference_performance_d1_temp.txt
$ sort -gk 5 -o MADR_reference_performance_d1_temp.txt
```

MADR_reference_performance_d1_temp.txt INDIVIDUAL MISS GENO TOTAL % GENOTYPED

```bash
DH0373_MMIR_SNA 110633 9964 120597 8.26
DH0308_MMIR_SNA 103181 17416 120597 14.44
DH0552_MMIR_SEA 101742 18855 120597 15.63
DH0311_MMIR_SNA 99322 21275 120597 17.64
DH0553_MMIR_SEA 98731 21866 120597 18.13
popfile_from_vcf_MODdennis.py MADR_reference_d1.vcf 13 >
MADR_reference_popfile_d1.txt
#Build the GD matrix
/home/pbongaerts/Github/radseq/vcf_gdmatrix.py MADR_reference_d1.vcf
```

MADR_reference_popfile_d1.txt > MADR_reference_gd_d1.txt

### D1b. Building Initial Neighbour-Joining Trees

Dart tree

```bash
/home/pbongaerts/Github/radseq/gdmatrix2tree.py
```

MADR_reference_gd_d1.txt MADR_reference_d1.tre

![](Images/image63.png)

QC of individuals using the tree

```bash
#repeat for all vcf variants
vcftools --vcf MADR_reference_d1.vcf --het
vcftools --vcf MADR_reference_d1.vcf --missing-indv
vcftools --vcf MADR_reference_d1.vcf --freq
vcftools --vcf MADR_reference_d1.vcf --depth
```

Adding depth to the tree:

```bash
% cat cur_kal_med.co2.csv cur_sna_med_full.co2.csv
cur_sea_med_full.co2.csv > cur_annotations_full.csv
% awk -F',' 'NR==FNR {map[$4] = $2; next} {print $0,
(map[$6] ? map[$6] : "NA")}' OFS=','
```

MADR_sample_record.csv cur_annotations_full.csv > cur_annotations_temp.csv For figures, see Tree_figures_d1.R

![](Images/image72.png)

![](Images/image45.png)

### D1c. Genetic Structure, Clone Detection & STRUCTURE Analysis

Initial clone detection

```bash
$ popfile_from_vcf_MODdennis.py MADR_reference_d1.vcf 13 >
MADR_reference_popfile_d1.txt
# 1. Calculate allelic similarity
$ vcf_clone_detect_npMOD.py -v MADR_reference_d1.vcf -p
```

MADR_reference_popfile_d1.txt -o MADR_reference_d1_clones.txt Check highest scores for replicates:

```bash
awk '
{
match($0, /DH([0-9]+)/, dh);
match($0, /DX([0-9]+)/, dx);
if (dh[1] == dx[1]) print $0;
```

}' clone_detect.out

```bash
#99.81 0.0 [SNA] DH0453_MMIR_SNA vs DX0453_MMIR_SNA 92485.0/92660
114223 99034
#99.79 0.0 [SNA] DH0452_MMIR_SNA vs DX0452_MMIR_SNA 92828.0/93022
112484 101135
#99.69 0.0 [SNA_SHA] DH0054_MMIR_SNA_SHA vs DX0054_MMIR_SNA_SHA
101450.0/101767 117210 105154
#99.67 0.0 [KAL_SHA] DH0114_MMIR_KAL_SHA vs DX0114_MMIR_KAL_SHA
```

101121.5/101460 117087 104970

```bash
#99.12 0.0 [SNA] DH0454_MMIR_SNA vs DX0454_MMIR_SNA 57379.0/57889
117749 60737
```

Based on these results we can put a threshold around 98% similarity for our initial clone detection

```bash
python3 /home/pbongaerts/Github/radseq/detect_clones_vcf.py -v
```

MADR_reference_d1.vcf -p MADR_reference_popfile_d1.txt -t 98 > MADR_reference_d1_cloneout_98.txt Use the output to select clonal lineages (between ###4 and ###5)

```bash
sed -n '/###4/,/###5/p' MADR_reference_d1_cloneout_98.txt >
MADR_reference_d1_clonal_groups.txt
$ python3 ~/scripts/add_clonal_group_to_popfile.py
```

MADR_reference_popfile_d1.txt MADR_reference_d1_clonal_groups.txt MADR_reference_d1_popfile_temp.txt Visualise clonal groups in the tree using R

![](Images/image67.png)

![](Images/image69.png)

Initial results look good enough to remove these individuals from the data and run structure. Get the recommended clones to remove based on %missing from the cloneout

```bash
file.
sed -n '/###5/,$p' MADR_reference_d1_cloneout_98.txt >
```

MADR_reference_d1_clones_to_remove.txt

**Remove clones from vcf**

```bash
vcftools --vcf MADR_reference_d1.vcf --remove
```

MADR_reference_d1_clones_to_remove.txt --recode --out MADR_reference_d1_no_clones

```bash
#kept 309 individuals -reference
```

Run through NJ tree building steps again for the data without clones (See D1b) Run Structure

```bash
$ scp /home/deepcat/STRUCTURE_params/extraparams .
$ scp /home/deepcat/STRUCTURE_params/mainparams .
$ vcftools --vcf ../MADR_reference_d1_no_clones.vcf
#After filtering, kept 309 out of 309 Individuals
#After filtering, kept 121930 out of a possible 121930 Sites
$ nano mainparams
```

The `mainparams` file for STRUCTURE was edited as follows (comments explain each field; see the STRUCTURE documentation for full details):

```text
Basic Program Parameters
#define MAXPOPS 7      // (int) number of populations assumed
#define BURNIN 100000    // (int) length of burnin period
#define NUMREPS 50000     // (int) number of MCMC reps after burnin

Input/Output files
#define INFILE infile     // (str) name of input data file
#define OUTFILE outfile   // (str) name of output data file

Data file format
#define NUMINDS 296     // (int) number of diploid individuals in data file
#define NUMLOCI 34466    // (int) number of loci in data file
#define PLOIDY 2          // (int) ploidy of data
#define MISSING -9        // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line
#define LABEL 1           // (B) Input file contains individual labels
#define POPDATA 1         // (B) Input file contains a population identifier
#define POPFLAG 0         // (B) Input file contains a flag saying whether to use popinfo when USEPOPINFO==1
#define LOCDATA 0         // (B) Input file contains a location identifier
#define PHENOTYPE 0       // (B) Input file contains phenotype information
#define EXTRACOLS 0       // (int) Number of additional columns of data before the genotype data start
#define MARKERNAMES 0     // (B) data file contains row of marker names
#define RECESSIVEALLELES 0 // (B) data file contains dominant markers (eg AFLPs) and a row to indicate which alleles are recessive
#define MAPDISTANCES 0    // (B) data file contains row of map distances between loci

Advanced data file options
#define PHASED 0          // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO 0       // (B) the data for each individual contains a line indicating phase (linkage model)
#define MARKOVPHASE 0     // (B) the phase info follows a Markov model
#define NOTAMBIGUOUS -999 // (int) for use in some analyses of polyploid data
```

```bash
#set up interactive shell with 8 cpu's and 2gb ram
$ ~/radseq/structure_mp.py MADR_reference_d1_strict_no_clones.vcf
MADR_reference_d1_strict_popfile_no_clones.txt 4 10 32
K = 2: MedMeaK 2.0 MaxMeaK 2 MedMedK 2.0 MaxMedK 2
```

MADR_reference_d1_no_clones.vcf

```bash
K = 3: MedMeaK 1.0 MaxMeaK 2 MedMedK 1.0 MaxMedK 2
MADR_reference_d1_no_clones.vcf
K = 4: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1
```

MADR_reference_d1_no_clones.vcf

```bash
K = 5: MedMeaK 0.0 MaxMeaK 1 MedMedK 0.0 MaxMedK 1
MADR_reference_d1_no_clones.vcf
K = 6: MedMeaK 0.0 MaxMeaK 1 MedMedK 0.0 MaxMedK 1
```

MADR_reference_d1_no_clones.vcf

```bash
K = 7: MedMeaK 0.0 MaxMeaK 1 MedMedK 0.0 MaxMedK 1
```

MADR_reference_d1_no_clones.vcf For strict dataset: Subsample SNPs (one random SNP per locus)... [4957 SNPs/loci] Outputting 10 STRUCTURE files...10 reps DONE Executing 32 parallel STRUCTURE runs for K = 2 ...10 reps DONE Executing 32 parallel STRUCTURE runs for K = 3 ...10 reps DONE Executing 32 parallel STRUCTURE runs for K = 4 ...10 reps DONE Running CLUMPP on replicates for K = 2 ... Running CLUMPP on replicates for K = 3 ... Running CLUMPP on replicates for K = 4 ...

```bash
K = 2: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1
MADR_reference_d1_strict_no_clones.vcf
K = 3: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1
```

MADR_reference_d1_strict_no_clones.vcf

```bash
K = 4: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1
```

MADR_reference_d1_strict_no_clones.vcf Running snapclust and structure on d1_no_clones snapclust.R

```bash
#warning
> MADR_snap_k2 <- snapclust(MADR_genind, 2)
Large dataset syndrome:
for 308 individuals, differences in log-likelihoods exceed computer
```

precision; group membership probabilities are approximated (only trust clear-cut values)

![](Images/image61.png)

![](Images/image24.png)

![](Images/image32.png)

![](Images/image27.png)

Check to see if we get the same assignment of individuals from the denovo and dart vcfs.

```bash
# use pop assignment script to assign populations based on structure output
python3 ~/scripts/MADR_compare_strucout.py MADR_reference_d1_strict_no_clones_1739223814/*.out.csv
```

Compared assignments with lineages based on the NJ tree using iTOL, for both the reference and denovo assemblies: found no inconsistencies in samples belonging to the minor or major lineage.

```bash
%awk -F ': ' '{gsub(/\(.*\)/, "", $2); gsub(/ /, "", $2); print $2}' \
    MADR_reference_d1_strict_clonal_groups.txt > MADR_reference_d1_strict_clonal_groups_clean.txt
```

```bash
% python3 map_clones_to_popfile.py
```

MADR_reference_d1_strict_clonal_groups_clean.txt MADR_reference_d1_strict_no_clones_popfile_strucK2.txt MADR_reference_d1_strict_popfile_strucK2_with_clones.txt

```bash
%awk '{if ($4 == "Major_lineage") print $1 >
"MADR_reference_d1_major_lineage.txt"; else if ($4 ==
"Minor_lineage") print $1 >
```

"MADR_reference_d1_minor_lineage.txt"}' MADR_reference_d1_strict_popfile_strucK2_with_clones.txt

```bash
613 MADR_reference_d1_major_lineage.txt
68 MADR_reference_d1_minor_lineage.txt
263 MADR_reference_d1_major_no_clones
33 MADR_reference_d1_minor_no_clones
```

**#**

## E1. Clone Detection

### E1a. Initial Genetic Similarity Comparison

Remove monomorphic sites and low genotyped sites and create vcf files

```bash
for each lineage
$ vcftools --vcf MADR_reference_C3.vcf --remove
MADR_reference_d1_major_lineage.txt --recode --out MADR_reference_e1_minor_full_temp
#After filtering, kept 68 out of 681 Individuals
#Outputting VCF file...
#After filtering, kept 418633 out of a possible 418633 Sites
$ vcftools --vcf MADR_reference_e1_minor_full.vcf --remove
```

indv_to_remove_reference_c3.txt --recode --out MADR_reference_e1_minor_full_temp.vcf

```bash
#After filtering, kept 65 out of 68 Individuals
$ vcftools --vcf MADR_reference_C3.vcf --remove
MADR_reference_d1_minor_lineage.txt --recode --out MADR_reference_e1_major_full
#After filtering, kept 616 out of 681 Individuals
#Outputting VCF file...
#After filtering, kept 418633 out of a possible 418633 Sites
$ vcftools --vcf MADR_reference_e1_major_full.vcf --remove
```

indv_to_remove_reference_c3.txt --recode --out MADR_reference_e1_major_full_temp.vcf

```bash
#After filtering, kept 613 out of 616 Individuals
vcftools --vcf MADR_reference_d1_strict_minor.vcf --mac 1
```

--max-missing 0.8 --minDP 8 --recode --out MADR_reference_e1_minor After filtering, kept 65 out of 65 Individuals After filtering, kept 12199 out of a possible 34466 Sites

```bash
vcftools --vcf MADR_reference_d1_strict_major.vcf --mac 1
```

--max-missing 0.8 --minDP 8 --recode --out MADR_reference_e1_major After filtering, kept 613 out of 613 Individuals After filtering, kept 28409 out of a possible 34466 Sites

```bash
vcftools --vcf MADR_reference_e1_minor_full.vcf --mac 1 --max-missing
```

0.8 --minDP 8 --recode --out MADR_reference_e1_minor_full_filtered After filtering, kept 65 out of 65 Individuals After filtering, kept 14752 out of a possible 418633 Sites

```bash
vcftools --vcf MADR_reference_e1_major_full.vcf --mac 1 --max-missing
```

0.8 --minDP 8 --recode --out MADR_reference_e1_major_full_filtered After filtering, kept 613 out of 613 Individuals After filtering, kept 30334 out of a possible 418633 Sites

```bash
$ python3 popfile_from_vcf_MODdennis.py MADR_reference_e1_minor.vcf 13
> MADR_reference_popfile_e1_minor.txt
$ python3 popfile_from_vcf_MODdennis.py MADR_reference_e1_major.vcf 13
```

> MADR_reference_popfile_e1_major.txt

```bash
# 1. Calculate allelic similarity
$ screen -S clonality
python3 /home/pbongaerts/Github/radseq/detect_clones_vcf.py -v
```

MADR_reference_e1_minor_full.vcf -p MADR_reference_popfile_e1_minor.txt -o MADR_reference_e1_minor_ful_genetic_sim.txt > MADR_reference_e1_minor_full_clonedet_out2.txt

```bash
python3 /home/pbongaerts/Github/radseq/detect_clones_vcf.py -v
```

MADR_reference_e1_major_full.vcf -p MADR_reference_popfile_e1_major.txt -o MADR_reference_e1_major_full_genetic_sim.txt > MADR_reference_e1_major_full_clonedet_out2.txt Building trees for each lineage

```bash
$ /home/pbongaerts/Github/radseq/vcf_gdmatrix.py
```

MADR_reference_e1_minor_full.vcf MADR_reference_popfile_e1_minor.txt > MADR_reference_e1_minor_full_gdmat.txt

```bash
$ /home/pbongaerts/Github/radseq/gdmatrix2tree.py
MADR_reference_e1_minor_full_gdmat.txt MADR_reference_e1_minor_full.tre
$ /home/pbongaerts/Github/radseq/vcf_gdmatrix.py
```

MADR_reference_e1_major_full.vcf MADR_reference_popfile_e1_major.txt > MADR_reference_e1_major_full_gdmat.txt

```bash
$ /home/pbongaerts/Github/radseq/gdmatrix2tree.py
```

MADR_reference_e1_major_full_gdmat.txt MADR_reference_e1_major_full.tre Clone detection using clone_detect.R for Major lineage with the general

```bash
filtered data
```

![](Images/image18.png)

Most appropriate threshold for Major lineages seems to be the Average threshold dictated by poppr pictured in Blue in the graph above. Important notice: the purple line indicates the genetic similarity of the replicate pair with the highest

![](Images/image30.png)

Visualise the genetic similarity between pairs in a histogram to determine the appropriate threshold. Average similarity and average threshold respectively for each measures seem to lie around the same point, just after the last big peak and might be useful as a threshold. However Using the lowest similarity between known replicates could be an easy, more lenient way to select clones as this value lies much lower in both cases. This would lead to a threshold of 99% similarity and a distance threshold of 0.002748542 for the 28409 site dataset. Comparing clone calls between the two methods:

```bash
python3 /home/pbongaerts/Github/radseq/detect_clones_vcf.py -v
```

MADR_reference_e1_major.vcf -p MADR_reference_popfile_e1_major.txt -o MADR_reference_e1_major_genetic_sim2_t99.txt -t 99 > MADR_reference_e1_major_clonedet_out2_t99.txt Results based on the lineage specific filtered data:

![](Images/image36.png)

![](Images/image13.png)

```bash
#minor lineage
```

![](Images/image49.png)

![](Images/image56.png)

Given the results we have decided to go with the lowest genetic similarity of the replicates and set the threshold for the definition of clones at 99% similarity

```bash
clonal_pairs="DH0452_MMIR_SNA DX0452_MMIR_SNA
DH0453_MMIR_SNA DX0453_MMIR_SNA
DH0454_MMIR_SNA DX0454_MMIR_SNA
DH0054_MMIR_SNA_SHA DX0054_MMIR_SNA_SHA
DH0114_MMIR_KAL_SHA DX0114_MMIR_KAL_SHA"
while IFS=$'\n' read -r pair; do
# Split the pair into two variables
read -r ind1 ind2 <<< "$pair"
awk -F',' -v i1="$ind1" -v i2="$ind2" '
($1 == i1 && $2 == i2) || ($1 == i2 && $2 == i1) {
print $0
}
```

' MADR_reference_e1_major_full_genetic_sim.txt

```bash
done <<< "$clonal_pairs"
# ind1,ind2,ind1_snps,ind2_snps,both_snps,match,match_perc,pop
#DH0452_MMIR_SNA,DX0452_MMIR_SNA,29965,29507,29138,29116.000000,99.920000,SNA
#DH0453_MMIR_SNA,DX0453_MMIR_SNA,30056,29200,28922,28907.500000,99.950000,SNA
#DH0454_MMIR_SNA,DX0454_MMIR_SNA,30084,22045,21795,21579.500000,99.010000,SNA
#DH0054_MMIR_SNA_SHA,DX0054_MMIR_SNA_SHA,30083,29411,29160,29116.000000,99.850000,SNA_SHA
#DH0114_MMIR_KAL_SHA,DX0114_MMIR_KAL_SHA,30013,29246,28925,28865.500000,99.790000,KAL_SHA
```

![](Images/image25.png)

With the initial threshold of 99.01% similarity we can find this interesting laderised pattern in the NJ tree. Upon further inspection combining the missing data we can see there is a clear overlap between these "difficult cases" and missing data. This is likely caused by the small denominator effect where due to the lower shared sites each "error" carries more weight compared to complete data. To get a better understanding of the sequencing error rate of our data we will look into the samples sequenced multiple times by DArT. First we make a list of all samples that were sequenced multiple times.

```bash
input_file="$HOME/Projects/MADR/A1a_data_prep/MADR_targets_22_24.csv"
output_file="MADR_targets_multi_run.csv"
awk -F, '
BEGIN {
```

OFS = ","

```bash
}
NR == 1 {
```

print $0 # Print header as is next

```bash
}
{
count[$5]++
if (count[$5] > 1) {
$5 = $5 "_x" count[$5] # Correct way to append "_x{number}"
}
print $0
}' "$input_file" > "$output_file"
```

Renaming the FASTQ_files

```bash
#!/bin/bash
# Define input CSV file
csv_file="MADR_targets_multi_run.csv"
FASTQ_path="/home/deepcat/rawseq/MADR_dart/MADR_full"
FASTQ_out="FASTQ_repeats"
# Loop through unique genotypes in the CSV file
awk -F',' 'NR>1 {print $5}' "$csv_file" | sort | uniq |
while read genotype; do
if [ -z "$genotype" ]; then
continue # Skip empty genotype values
fi
# Get all targetid values for the current genotype
target_ids=$(awk -F',' -v g="$genotype" 'NR>1 && $5 == g {print
$1}' "$csv_file")
# Create an array of FASTQ files corresponding to these target_ids
fastq_files=() # Initialize an empty array
for target_id in $target_ids; do
file_path="${FASTQ_path}/${target_id}.FASTQ.gz"
if [ -f "$file_path" ]; then
fastq_files+=("$file_path") # Add the file path to the array
fi
done
# Skip if no FASTQ files are found
if [ ${#fastq_files[@]} -eq 0 ]; then
echo "No FASTQ files found for genotype $genotype"
continue
fi
# Combine the FASTQ files for this genotype
output_file="${FASTQ_out}/${genotype}.FASTQ.gz"
echo "Combining files for genotype $genotype into $output_file"
cat "${fastq_files[@]}" > "$output_file"
done
```

QC, lineage split and Ipyrad run

```bash
trim_galore --cores 8 -q 20 --length 30 --phred33 -o fastq_trimmed/
FASTQ_PATH/*.FASTQ.gz
```

Filter the data in a few ways to get close to the original dataset.

```bash
vcftools --vcf MADR_reference_e1_replicates.vcf --mac 1
--max-missing 0.5 --recode --out MADR_reference_e1_replicates_loose_filter
vcftools --vcf MADR_reference_e1_replicates.vcf --mac 1 --max-missing
```

0.8 --minDP 8 --recode --out MADR_reference_e1_replicates_strict_filter

**Take out underperforming individuals**

```bash
$ /home/pbongaerts/Github/radseq/vcf_missing_data.py
MADR_reference_e1_replicates.vcf > MADR_reference_e1_replicates_qc.txt
$ sort -gk 5 -o MADR_reference_e1_replicates_qc.txt
```

MADR_reference_e1_replicates_qc.txt

```bash
$ awk '
BEGIN { less_10 = 0; less_30 = 0; less_1000_geno = 0; total_percent = 0; count = 0 }
NR > 1 {                                  # Skip header
    if ($5 < 10)   less_10++;             # Check percentage genotyped
    if ($5 < 30)   less_30++;
    if ($3 < 1000) less_1000_geno++;      # Check GENO column
    total_percent += $5;                  # Add % genotyped to total
    count++;
}
END {
    print "Samples with <10% genotyped: "   less_10;
    print "Samples with <30% genotyped: "   less_30;
    print "Samples with <1000 SNPs (GENO): " less_1000_geno;
    print "Average % genotyped: " total_percent / count;
}

' MADR_reference_e1_replicates_qc.txt

```bash
#Samples with <10% genotyped: 22
#Samples with <30% genotyped: 291
#Samples with <1000 SNPs (GENO): 4
#Average % genotyped: 32.3641
$ awk '$3 < 1000 {print $1}' MADR_reference_e1_replicates_qc.txt
```

> indv_to_remove_reference_e1_replicates.txt

```bash
vcftools --vcf MADR_reference_e1_replicates.vcf --remove
```

indv_to_remove_reference_e1_replicates.txt --recode --out MADR_reference_e1_replicates_temp.vcf Run similarity analisis on both

```bash
python3 /home/pbongaerts/Github/radseq/vcf_clone_detect.py -v
```

MADR_reference_e1_replicates.vcf -p MADR_reference_popfile_e1_replicates.txt -o MADR_reference_e1_replicates_genetic_sim.txt > MADR_reference_e1_replicates_clonedet_out.txt

```bash
python3/home/pbongaerts/Github/radseq/vcf_clone_detect.py -v
```

MADR_reference_e1_replicates_strict.vcf -p MADR_reference_popfile_e1_replicates.txt -o MADR_reference_e1_replicates_genetic_sim_strict.txt > MADR_reference_e1_replicates_clonedet_strict_out.txt

```bash
python3 /home/pbongaerts/Github/radseq/vcf_clone_detect.py -v
```

MADR_reference_e1_replicates_loose.vcf -p MADR_reference_popfile_e1_replicates.txt -o MADR_reference_e1_replicates_genetic_sim_loose.txt > MADR_reference_e1_replicates_clonedet_loose_out.txt Compare replicates

```bash
$ awk -F',' '
$1 ~ /_x[0-9]/ && $2 == substr($1, 1, length($1)-3) ||
$2 ~ /_x[0-9]/ && $1 == substr($2, 1, length($2)-3) {
print $0
}
```

' MADR_reference_e1_replicates_genetic_sim_strict.txt > MADR_reference_e1_replicates_strict_genetic_sim_double_sequenced_pairs.txt

```bash
$ tail -n 15
MADR_reference_e1_replicates_strict_genetic_sim_double_sequenced_pairs.txt
# ind1,ind2,ind1_snps,ind2_snps,both_snps,match,match_perc,pop
#DH0527_MMIR_SEA,DH0527_MMIR_SEA_x2,27258,30108,25239,25022.500000,99.140000,SEA-SEA_x2
#DH0525_MMIR_SEA,DH0525_MMIR_SEA_x2,25887,29945,24708,24408.000000,98.790000,SEA-SEA_x2
#DH0327_MMIR_SNA,DH0327_MMIR_SNA_x2,13945,15663,9134,9018.500000,98.740000,SNA-SNA_x2
#DH0575_MMIR_SEA,DH0575_MMIR_SEA_x2,26476,29328,24093,23759.000000,98.610000,SEA-SEA_x2
#DH0370_MMIR_SNA,DH0370_MMIR_SNA_x2,28781,19552,17822,17461.500000,97.980000,SNA-SNA_x2
#DH0369_MMIR_SNA,DH0369_MMIR_SNA_x2,26606,25991,21512,21037.000000,97.790000,SNA-SNA_x2
#DH0423_MMIR_SNA,DH0423_MMIR_SNA_x2,28805,18250,16814,16293.500000,96.900000,SNA-SNA_x2
#DH0326_MMIR_SNA,DH0326_MMIR_SNA_x2,19391,13805,8653,8380.500000,96.850000,SNA-SNA_x2
#DH0421_MMIR_SNA,DH0421_MMIR_SNA_x2,18437,25031,14747,14259.000000,96.690000,SNA-SNA_x2
#DH0325_MMIR_SNA,DH0325_MMIR_SNA_x2,22475,12923,9562,9185.500000,96.060000,SNA-SNA_x2
#DH0258_MMIR_SNA,DH0258_MMIR_SNA_x2,23252,16331,12020,11536.000000,95.970000,SNA-SNA_x2
#DH0211_MMIR_SNA,DH0211_MMIR_SNA_x2,17748,6137,3661,3510.500000,95.890000,SNA-SNA_x2
#DH0212_MMIR_SNA,DH0212_MMIR_SNA_x2,19160,14632,9069,8693.500000,95.860000,SNA-SNA_x2
#DH0419_MMIR_SNA,DH0419_MMIR_SNA_x2,15384,12487,6471,6195.500000,95.740000,SNA-SNA_x2
#DH0372_MMIR_SNA,DH0372_MMIR_SNA_x2,21448,12851,8720,8324.000000,95.460000,SNA-SNA_x2
```

The lowest replicate that shares at least 10.000 sites is

```bash
#DH0258_MMIR_SNA,DH0258_MMIR_SNA_x2,23252,16331,12020,11536.000000,95.970000,SNA-SNA_x2
```

If we use this as our threshold:

![](Images/image2.png)

![](Images/image71.png)

### E1b. ANGSD IBS-Based Clonality Assessment

Clone detection based on ANGSD IBS distance matrix /ccg/bin/angsd -b bam_list.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 475 -setMinDepth 8 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_ANGSD_outputs

**-> Output filenames:**

->"MADR_reference_e1_ANGSD_outputs.arg" ->"MADR_reference_e1_ANGSD_outputs.ibs.gz" ->"MADR_reference_e1_ANGSD_outputs.ibsMat" -> Thu Mar 13 15:01:00 2025 -> Arguments and parameters for all analysis are located in .arg file -> Total number of sites analyzed: 43341802 -> Number of sites retained after filtering: 2164473 [ALL done] cpu-time used = 6760.01 sec [ALL done] walltime used = 6843.00 sec Create bam order sample name list:

```bash
awk -F'/' '{split($NF, a, ".FASTQ"); print a[1]}' bam_list.txt > sample_names_bam_order.txt
```

Initial IBS results:

![](Images/image58.png)

```bash
/ccg/bin/angsd -b bam_list.txt -GL 1 -uniqueOnly 1 -remove_bads 1
```

-minMapQ 20 -minQ 20 -minInd 543 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_ANGSD_outputs

```bash
-> Output filenames:
->"MADR_reference_e1_ANGSD_outputs.arg"
->"MADR_reference_e1_ANGSD_outputs.ibs.gz"
->"MADR_reference_e1_ANGSD_outputs.ibsMat"
-> Thu Mar 13 20:56:28 2025
```

-> Arguments and parameters for all analysis are located in .arg

```bash
file
-> Total number of sites analyzed: 43341802
-> Number of sites retained after filtering: 1561519
[ALL done] cpu-time used = 6030.48 sec
[ALL done] walltime used = 6121.00 sec
```

Splitting the lineages:

```bash
BAM_LIST="bam_list.txt" # File containing BAM file paths
REMOVE_LIST="$HOME/Projects/MADR/D1c_genetic_structure/reference/strict/MADR_reference_d1_minor_lineage.txt"
# File containing sample names to remove
OUTPUT_FILE="bam_list_major_lineage.txt" # Output file
# Create a grep pattern from REMOVE_LIST
PATTERN=$(awk '{print $1}' "$REMOVE_LIST" | paste -sd "|" -)
# Filter BAM list, keeping only lines that do NOT match the pattern
grep -Ev "$PATTERN" "$BAM_LIST" > "$OUTPUT_FILE"
echo "Filtered BAM list saved to $OUTPUT_FILE"
#repeat for minor lineage
/ccg/bin/angsd -b bam_list_major_lineage.txt -GL 1 -uniqueOnly 1
```

-remove_bads 1 -minMapQ 20 -minQ 20 -minInd 494 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_major_ANGSD_outputs

```bash
-> Output filenames:
->"MADR_reference_e1_major_ANGSD_outputs.arg"
->"MADR_reference_e1_major_ANGSD_outputs.ibs.gz"
->"MADR_reference_e1_major_ANGSD_outputs.ibsMat"
-> Thu Mar 13 21:01:07 2025
```

-> Arguments and parameters for all analysis are located in .arg

```bash
file
-> Total number of sites analyzed: 40679307
-> Number of sites retained after filtering: 1619345
[ALL done] cpu-time used = 5531.98 sec
[ALL done] walltime used = 5592.00 sec
```

/ccg/bin/angsd -b bam_list_minor_lineage.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 54 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_minor_ANGSD_outputs

```bash
-> Output filenames:
->"MADR_reference_e1_minor_ANGSD_outputs.arg"
->"MADR_reference_e1_minor_ANGSD_outputs.ibs.gz"
->"MADR_reference_e1_minor_ANGSD_outputs.ibsMat"
-> Thu Mar 13 19:32:58 2025
```

-> Arguments and parameters for all analysis are located in .arg

```bash
file
-> Total number of sites analyzed: 17133579
-> Number of sites retained after filtering: 1079849
[ALL done] cpu-time used = 345.06 sec
[ALL done] walltime used = 345.00 sec
#adding options to infer missingness running the same filters but not
```

calculating ibs to save on time /ccg/bin/angsd -b bam_list_major_lineage.txt -GL 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 2 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 494 -setMinDepth 20 -P 4 -out geno_text_output_major

```bash
gzip -d geno_text_output_major.geno.gz
```

/ccg/bin/angsd -b bam_list_minor_lineage.txt -GL 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 2 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 54 -setMinDepth 20 -P 4 -out geno_text_output_minor

```bash
gzip -d geno_text_output_minor.geno.gz
```

Devided into clonal groups using Figure_S1_clonality.R Initial IBS results minor lineage: setting the threshold according to the split in node heights

![](Images/image42.png)

![](Images/image55.png)

Mapping the IBS clonal groups to the NJ tree (keeping colors consistent)

![](Images/image17.png)

![](Images/image68.png)

![](Images/image38.png)

![](Images/image8.png)

![](Images/image35.png)

![](Images/image1.png)

Using replicates to set a threshold:

```bash
/ccg/bin/angsd -b bam_list_replicates.txt -GL 1 -uniqueOnly 1
```

-remove_bads 1 -minMapQ 20 -minQ 20 -minInd 543 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_replicates_ANGSD_outputs

```bash
-> Output filenames:
->"MADR_reference_e1_replicates_ANGSD_outputs.arg"
->"MADR_reference_e1_replicates_ANGSD_outputs.ibs.gz"
->"MADR_reference_e1_replicates_ANGSD_outputs.ibsMat"
-> Sun Mar 16 16:29:51 2025
```

-> Arguments and parameters for all analysis are located in .arg

```bash
file
-> Total number of sites analyzed: 43342207
-> Number of sites retained after filtering: 1504952
[ALL done] cpu-time used = 7326.38 sec
[ALL done] walltime used = 7372.00 sec
```

**#Setting threshold to 0.016 based on the lowest scoring replicate**

cluster ibs.

![](Images/image34.png)

Find the threshold for each individually:

```bash
/ccg/bin/angsd -b bam_list_replicates_minor_lineage.txt -GL 1
```

-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 66 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_replicates_minor_ANGSD_outputs

```bash
-> Output filenames:
->"MADR_reference_e1_replicates_minor_ANGSD_outputs.arg"
->"MADR_reference_e1_replicates_minor_ANGSD_outputs.ibs.gz"
->"MADR_reference_e1_replicates_minor_ANGSD_outputs.ibsMat"
-> Mon Mar 17 14:05:00 2025
```

-> Arguments and parameters for all analysis are located in .arg

```bash
file
-> Total number of sites analyzed: 17129578
-> Number of sites retained after filtering: 1270478
[ALL done] cpu-time used = 522.16 sec
[ALL done] walltime used = 523.00 sec
```

/ccg/bin/angsd -b bam_list_replicates_major_lineage.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 563 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_replicates_major_ANGSD_outputs

```bash
-> Output filenames:
->"MADR_reference_e1_replicates_major_ANGSD_outputs.arg"
->"MADR_reference_e1_replicates_major_ANGSD_outputs.ibs.gz"
->"MADR_reference_e1_replicates_major_ANGSD_outputs.ibsMat"
-> Mon Mar 17 15:44:52 2025
```

-> Arguments and parameters for all analysis are located in .arg

```bash
file
-> Total number of sites analyzed: 40677711
-> Number of sites retained after filtering: 1621212
[ALL done] cpu-time used = 6412.29 sec
[ALL done] walltime used = 6420.00 sec
```

USE missing % from VCF to select clones to keep:

```bash
samples_to_keep <- clone_group_df %>%
# Keep all samples with clone_group == 0
filter(clone_group == 0) %>%
# Combine with the samples from other groups with the lowest F_miss
```

bind_rows( clone_group_df %>%

```bash
filter(clone_group != 0) %>%
```

group_by(clone_group) %>%

```bash
slice_min(F_miss, n = 1, with_ties = FALSE) %>% # Select the lowest
```

F_miss without ties ungroup() )

```bash
#177 samples left for Major lineage
#28 samples left for Minor lineage
```

## F1. Lineage Comparison

Create new datasets without clonal samples:

```bash
vcftools --vcf MADR_reference_C3.vcf --keep
```

MADR_reference_e1_major_samples_to_keep.txt --recode --out MADR_reference_f1_major_full

```bash
VCFtools - 0.1.16
Keeping individuals in 'keep' list
After filtering, kept 177 out of 681 Individuals
Outputting VCF file...
After filtering, kept 418633 out of a possible 418633 Sites
Run Time = 196.00 seconds
vcftools --vcf MADR_reference_C3.vcf --keep
```

MADR_reference_e1_minor_samples_to_keep.txt --recode --out MADR_reference_f1_minor_full

```bash
Keeping individuals in 'keep' list
After filtering, kept 28 out of 681 Individuals
Outputting VCF file...
After filtering, kept 418633 out of a possible 418633 Sites
vcftools --vcf MADR_reference_C3.vcf --keep
```

MADR_reference_e1_all_samples_to_keep.txt --recode --out MADR_reference_f1_full

```bash
Keeping individuals in 'keep' list
After filtering, kept 205 out of 681 Individuals
Outputting VCF file...
After filtering, kept 418633 out of a possible 418633 Sites
Run Time = 205.00 seconds
#filter data
vcftools --vcf MADR_reference_f1_major_full.vcf --mac 1 --max-missing
```

0.8 --minDP 8 --recode --out MADR_reference_f1_major After filtering, kept 177 out of 177 Individuals Outputting VCF file... After filtering, kept 27,201 out of a possible 418,633 Sites

```bash
vcftools --vcf MADR_reference_f1_minor_full.vcf --mac 1 --max-missing
```

0.8 --minDP 8 --recode --out MADR_reference_f1_minor After filtering, kept 28 out of 28 Individuals Outputting VCF file... After filtering, kept 13,498 out of a possible 418,633 Sites

```bash
vcftools --vcf MADR_reference_f1_full.vcf --mac 1 --max-missing 0.8
```

--minDP 8 --recode --out MADR_reference_f1_full After filtering, kept 205 out of 205 Individuals Outputting VCF file... After filtering, kept 31,713 out of a possible 418,633 Sites

```bash
vcftools --vcf MADR_reference_f1_full_allele_analysis.vcf --mac 1
```

--max-missing 0.05 --minDP 8 --recode --out MADR_reference_f1_full_loose_filter After filtering, kept 205 out of 205 Individuals Outputting VCF file... After filtering, kept 169,029 out of a possible 418,633 Sites QC

```bash
/home/pbongaerts/Github/radseq/vcf_missing_data.py
~/Projects/MADR/F1_populations_statistics/MADR_reference_f1_major.vcf
> MADR_reference_f1_major_qc.txt
/home/pbongaerts/Github/radseq/vcf_missing_data.py
~/Projects/MADR/F1_populations_statistics/MADR_reference_f1_minor.vcf
```

> MADR_reference_f1_minor_qc.txt

```bash
/home/pbongaerts/Github/radseq/vcf_missing_data.py
~/Projects/MADR/F1_populations_statistics/MADR_reference_f1_all.vcf >
MADR_reference_f1_all_qc.txt
$ sort -gk 5 -o MADR_reference_f1_major_qc.txt MADR_reference_f1_major_qc.txt
$ awk '
BEGIN { less_10 = 0; less_30 = 0; less_1000_geno = 0; total_percent = 0; count = 0 }
NR > 1 {                                  # Skip header
    if ($5 < 10)   less_10++;             # Check percentage genotyped
    if ($5 < 30)   less_30++;
    if ($3 < 1000) less_1000_geno++;      # Check GENO column
    total_percent += $5;                  # Add % genotyped to total
    count++;
}
END {
    print "Samples with <10% genotyped: "   less_10;
    print "Samples with <30% genotyped: "   less_30;
    print "Samples with <1000 SNPs (GENO): " less_1000_geno;
    print "Average % genotyped: " total_percent / count;
}
' MADR_reference_f1_major_qc.txt
# repeated for the minor lineage and the full (all) dataset
```

| Dataset | <10% genotyped | <30% genotyped | <1000 SNPs (GENO) | Average % genotyped |
|---|---|---|---|---|
| Major lineage | 0 | 2 | 0 | 90.18% |
| Minor lineage | 1 | 1 | 1 | 87.12% |
| All (1 individual removed) | 1 | 3 | 1 | 89.42% |

individual to remove: DH0293_MMIR_SNA Building NJ trees

```bash
$ python3 popfile_from_vcf_MODdennis.py MADR_reference_f1_major.vcf 13
> MADR_reference_popfile_f1_major.txt
$ /home/pbongaerts/Github/radseq/vcf_gdmatrix.py
```

MADR_reference_f1_major.vcf MADR_reference_popfile_f1_major.txt > MADR_reference_f1_major_gd_mat.txt

```bash
$ /home/pbongaerts/Github/radseq/gdmatrix2tree.py
MADR_reference_f1_major_gd_mat.txt MADR_reference_f1_major.tre
#repeat for minor and all
```

### F1a. Private and Fixed Alleles

Hybrids not separated

```bash
# For major lineage
vcftools --vcf MADR_reference_f1_full_loose_filter.recode.vcf --keep
```

MADR_reference_e1_major_samples_to_keep.txt --freq --out allele_frequencies/major_lineage_allele_freq

```bash
# For minor lineage
vcftools --vcf MADR_reference_f1_full_loose_filter.recode.vcf --keep
```

MADR_reference_e1_minor_samples_to_keep.txt --freq --out allele_frequencies/minor_lineage_allele_freq

```bash
Finding private and fixed alleles using {compare_allele_frequencies.py}
compare_allele_frequencies.py freq_list.txt 0.96 0.01
MADR_reference_f1_loose
#freq_list.txt = txt file containing .frq file path + num_indv
#0.96 % of individuals genotyped for that site
#0.01 leniancy for genotyping error
wc -l major/minor_lineage_allele_freq.frq
```

169030

```bash
grep "major" MADR_reference_f1_loose_unique_fixed_alleles.txt | wc -l
```

2533 (1.5% of total sites)

```bash
grep "minor" MADR_reference_f1_loose_unique_fixed_alleles.txt | wc -l
```

866 (0.5% of total sites)

```bash
wc -l MADR_reference_f1_loose_fixed_alleles.txt
```

3810 (- 3399 = 411 shared fixed alleles) (2.25% of total sites)

```bash
grep "major" MADR_reference_f1_loose_private_alleles.txt | wc -l
```

49047 (29.0 % of total sites)

```bash
grep "minor" MADR_reference_f1_loose_private_alleles.txt | wc -l
```

31014 (18.3 % of total sites)

```bash
wc -l MADR_reference_f1_loose_private_alleles.txt
```

78855 (46.6% of total sites) (49047+31014)-78855=1206 sites with alternating private alleles

```bash
grep "minor" MADR_reference_f1_loose_private_sites.txt | wc -l
```

22 (0.013% of total sites)

```bash
grep "major" MADR_reference_f1_loose_private_sites.txt | wc -l
```

8970 (5.30% of total sites)

**\**

Hybrids separated (See F1B)

```bash
# For major lineage
vcftools --vcf
```

../hybrids_not_separated/MADR_reference_f1_full_loose_filter.recode.vcf --keep DAPC_major_samples_to_keep.txt --freq --out major_lineage_allele_freq kept 176 samples

```bash
# For minor lineage
vcftools --vcf
```

../hybrids_not_separated/MADR_reference_f1_full_loose_filter.recode.vcf --keep DAPC_minor_samples_to_keep.txt --freq --out minor_lineage_allele_freq kept 23 samples

```bash
# For hybrids
vcftools --vcf
```

../hybrids_not_separated/MADR_reference_f1_full_loose_filter.recode.vcf --keep DAPC_hybrid_samples_to_keep.txt --freq --out hybrids_allele_freq kept 5 samples

```bash
Finding private and fixed alleles using {compare_allele_frequencies.py}
~/scripts/compare_allele_frequencies.py freq_list_maj_min.txt 0.96
0.01 MADR_reference_f1_min_maj
#freq_list.txt = txt file containing .frq file path + num_indv
#0.96 % of individuals genotyped for that site
#0.01 leniancy for genotyping error
wc -l major/minor_lineage_allele_freq.frq
```

169030

```bash
grep "major" MADR_reference_f1_min_maj_unique_fixed_alleles.txt | wc
```

-l 1913 (xx% of total sites)

```bash
grep "minor" MADR_reference_f1_min_maj_unique_fixed_alleles.txt | wc
```

-l 2871 (xx% of total sites)

```bash
wc -l MADR_reference_f1_min_maj_fixed_alleles.txt
```

5839 (- 4784 = 1055 shared fixed alleles) (xx% of total sites)

```bash
wc -l MADR_reference_f1_min_maj_unique_fixed_alleles.txt
```

4785 (- 4784 -1) = 0 alternatively fixed alleles) (xx% of total sites)

```bash
grep "major" MADR_reference_f1_min_maj_private_alleles.txt | wc -l
```

49619 (xx % of total sites)

```bash
grep "minor" MADR_reference_f1_min_maj_private_alleles.txt | wc -l
```

29761 (xx % of total sites)

```bash
wc -l MADR_reference_f1_min_maj_private_alleles.txt
```

78136 (xx% of total sites) (49047+31014)-78855=1206 sites with alternating private alleles

```bash
grep "minor" MADR_reference_f1_min_maj_private_sites.txt | wc -l
```

31 (xx% of total sites)

```bash
grep "major" MADR_reference_f1_min_maj_private_sites.txt | wc -l
```

9343 (xx% of total sites)

```bash
#loose filter for fixed alleles:
~/scripts/compare_allele_frequencies.py freq_list_maj_min.txt 0.6 0.05
MADR_reference_f1_min_maj_loose
grep "major" MADR_reference_f1_min_maj_loose_unique_fixed_alleles.txt
```

| wc -l 15275 (xx% of total sites)

```bash
grep "minor" MADR_reference_f1_min_maj_loose_unique_fixed_alleles.txt
```

| wc -l 9321 (xx% of total sites)

```bash
wc -l MADR_reference_f1_min_maj_loose_unique_fixed_alleles.txt
```

24593 (- 24596) = -3 alternatively fixed alleles) (xx% of total sites)

```bash
grep "major" MADR_reference_f1_min_maj_loose_private_alleles.txt | wc
```

-l 29108 (xx % of total sites)

```bash
grep "minor" MADR_reference_f1_min_maj_loose_private_alleles.txt | wc
```

-l 30193 (xx % of total sites)

```bash
wc -l MADR_reference_f1_min_maj_loose_private_alleles.txt
```

57857 (xx% of total sites) (49047+31014)-78855=1206 sites with alternating private alleles

```bash
grep "minor" MADR_reference_f1_min_maj_loose_private_sites.txt | wc
```

-l 31 (xx% of total sites)

```bash
grep "major" MADR_reference_f1_min_maj_loose_private_sites.txt | wc
```

-l 9343 (xx% of total sites) Finding most divergent sites:

```bash
~/scripts/compare_allele_frequencies.py freq_list_maj_min.txt 0.96
0.01 MADR_reference_f1_min_maj 20000
~/Projects/MADR/F1_populations_statistics/F1a_Private_and_fixed_alleles/reference/hybrids_separated$
awk '$4 > 0.5' MADR_reference_f1_min_maj_most_divergent_loci.txt |
wc -l
```

4962

```bash
awk '$4 > 0.6' MADR_reference_f1_min_maj_most_divergent_loci.txt |
wc -l
3730
awk '$4 > 0.6' MADR_reference_f1_min_maj_most_divergent_loci.txt >
```

MADR_reference_f1_most_divergent_loci.txt

### F1b. Population Statistical Analysis (PCA / DAPC)

Adding pop to popfile

```bash
awk 'NR==FNR {samples[$1]="major"; next} {print $0, ($1 in
samples) ? samples[$1] : "minor"}'
~/Projects/MADR/E1_clone_detection/MADR_reference_e1_major_samples_to_keep.txt
MADR_reference_popfile_f1_all.txt > MADR_reference_popfile_f1_lineages.txt Running popstats_MADR.R
#PCA with NA = 'mean'\
```

![](Images/image4.png)

![](Images/image66.png)

![](Images/image29.png)

DAPC

![](Images/image19.png)

![](Images/image26.png)

![](Images/image60.png)

Following the DAPC posterior assignment:

![](Images/image16.png)

![](Images/image48.png)

![](Images/image15.png)

Re running structure:

```bash
$ vcftools --vcf ../MADR_reference_f1_all.vcf
#After filtering, kept 204 out of 204 Individuals
#After filtering, kept 31713 out of a possible 31713 Sites
$ nano mainparams
```

The `mainparams` file for STRUCTURE was edited as follows (comments explain each field; see the STRUCTURE documentation for full details):

```text
Basic Program Parameters
#define MAXPOPS 4      // (int) number of populations assumed
#define BURNIN 100000    // (int) length of burnin period
#define NUMREPS 50000     // (int) number of MCMC reps after burnin

Input/Output files
#define INFILE infile     // (str) name of input data file
#define OUTFILE outfile   // (str) name of output data file

Data file format
#define NUMINDS 204     // (int) number of diploid individuals in data file
#define NUMLOCI 31713    // (int) number of loci in data file
#define PLOIDY 2          // (int) ploidy of data
#define MISSING -9        // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line
#define LABEL 1           // (B) Input file contains individual labels
#define POPDATA 1         // (B) Input file contains a population identifier
#define POPFLAG 0         // (B) Input file contains a flag saying whether to use popinfo when USEPOPINFO==1
#define LOCDATA 0         // (B) Input file contains a location identifier
#define PHENOTYPE 0       // (B) Input file contains phenotype information
#define EXTRACOLS 0       // (int) Number of additional columns of data before the genotype data start
#define MARKERNAMES 0     // (B) data file contains row of marker names
#define RECESSIVEALLELES 0 // (B) data file contains dominant markers (eg AFLPs) and a row to indicate which alleles are recessive
#define MAPDISTANCES 0    // (B) data file contains row of map distances between loci

Advanced data file options
#define PHASED 0          // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO 0       // (B) the data for each individual contains a line indicating phase (linkage model)
#define MARKOVPHASE 0     // (B) the phase info follows a Markov model
#define NOTAMBIGUOUS -999 // (int) for use in some analyses of polyploid data
```

```bash
#set up interactive shell with 8 cpu's and 2gb ram
$ ~/scripts/structure_mp_ref_update.py ../MADR_reference_f1_all.vcf
../MADR_reference_popfile_f1_all.txt 4 10 32
Initialise indivs and pops for ../MADR_reference_f1_all.vcf...
Subsample SNPs (one random SNP per locus)... [5492 SNPs/loci]
Outputting 10 STRUCTURE files...10 reps DONE
Executing 32 parallel STRUCTURE runs for K = 2 ...10 reps DONE
Executing 32 parallel STRUCTURE runs for K = 3 ...10 reps DONE
Executing 32 parallel STRUCTURE runs for K = 4 ...10 reps DONE
Running CLUMPP on replicates for K = 2 ...
Running CLUMPP on replicates for K = 3 ...
Running CLUMPP on replicates for K = 4 ...
```

K = 2: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1

```bash
../MADR_reference_f1_all.vcf
```

K = 3: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1

```bash
../MADR_reference_f1_all.vcf
```

K = 4: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1 ../MADR_reference_f1_all.vcf

![](Images/image23.png)

![](Images/image64.png)

Get max likelihood for k struct:

```bash
output_file="Evanos_K_values.txt"
for file in *.log; do
if grep -q "Estimated Ln Prob of Data" "$file"; then
# Extract the value
lnpd=$(grep "Estimated Ln Prob of Data" "$file" | awk -F'= '
'{print $2}')
echo "$file,$lnpd" >> "$output_file"
else
echo "$file,N/A" >> "$output_file"
fi
done
for suffix in "_2.log" "_3.log" "_4.log"; do
echo -n "Average for $suffix: "
awk -F',' -v sfx="$suffix" '
$1 ~ sfx {sum += $2; count++}
END {if (count > 0) print sum / count; else print "No entries"}
```

' "$output_file"

```bash
done
```

Average for _2.log: -410874 Average for _3.log: -410190 Average for _4.log: -421119 Test to see true hybrid assignment Running PCA/DAPC with simulated data using adegenet hybridize: To do:

![](Images/image40.png)

![](Images/image59.png)

![](Images/image62.png)

MAX PCAs =2 (more similar to results without added

![](Images/image52.png)

![](Images/image54.png)

![](Images/image10.png)

"Suspected" hybrids in green

![](Images/image31.png)

![](Images/image28.png)

![](Images/image3.png)

Running newhybrids with simulated and actual data:

```bash
Program New Hybrids completed at Thu Mar 27 16:49:44 2025
Data = MADR_reference_f1_nhyb.txt, PiPrior = JEFFREYS, ThetaPrior =
JEFFREYS, Seeds = 27 13, NumBurnIn = 50000, NumRepsAfterBurnIn = 100000
```

SeedsAtEnding = 1940924971 560041219 Program New Hybrids completed at Fri Mar 28 09:05:38 2025

```bash
Data = SIM_hybrids_nhyb.txt, PiPrior = JEFFREYS, ThetaPrior = JEFFREYS,
Seeds = 27 13, NumBurnIn = 50000, NumRepsAfterBurnIn = 100000
```

SeedsAtEnding = 594100171 307743162

![](Images/image51.png)

![](Images/image6.png)

![](Images/image43.png)

Nj tree It still seems like there is a fair bit of noise which could be due to the selection of samples to hybridize.

![](Images/image37.png)

HC clust

![](Images/image70.png)

NJ tree

![](Images/image14.png)

Full SIM data with only pure hybrids

![](Images/image12.png)

![](Images/image50.png)

![](Images/image73.png)

![](Images/image44.png)

![](Images/image46.png)

### F1c. Genetic Distance Between Lineages (hierfstat)

Genetic distance between lineages was calculated using `hierfstat` in the `Popstants_full_data.R` script.

![](Images/image41.png)

| Stat | Meaning | Interpretation |
|------|---------|-----------------|
| **Ho** | Observed heterozygosity | Mean proportion of individuals that are *heterozygous* at loci across the populations. (Here: **0.0998**) |
| **Hs** | Expected heterozygosity within populations | Expected heterozygosity assuming Hardy-Weinberg within each *population*. (Here: **0.1383**) |
| **Ht** | Total expected heterozygosity | Expected heterozygosity if you treated *all populations together* as one big group. (Here: **0.1499**) |
| **Dst** | Genetic divergence between populations | How much of the total genetic diversity is **due to differences between populations** (Ht − Hs). (Here: **0.0115**) |
| **Htp** | Ht corrected for sampling size | A corrected version of Ht to account for unequal or small sample sizes. (Here: **0.1557**) |
| **Dstp** | Dst corrected for sampling size | Same correction for Dst. (Here: **0.0174**) |
| **Fst** | Unbiased Weir & Cockerham global Fst | Proportion of genetic variation **due to differences between populations**. (Here: **0.0770**, or 7.7%) |
| **Fstp** | Fst corrected for sampling size | A corrected Fst. (Here: **0.1117**) |
| **Fis** | Inbreeding coefficient within populations | Measures *departure from Hardy-Weinberg* inside populations: **positive = excess of homozygosity**. (Here: **0.2786**) |
| **Dest** | Jost's D, an alternative to Fst | A measure of *actual allele differentiation* between pops. (Here: **0.0202**) |

## F2. Outlier Detection with PCAdapt & RDA

Outlier loci were identified using the `pcadapt` R package, run from `MADR_pcadapt.R`.

![](Images/image11.png)

![](Images/image20.png)

![](Images/image65.png)

![](Images/image74.png)

RDA

![](Images/image7.png)

![](Images/image21.png)

## G1. Reefscape Assessment

### G1a. Demographic Analysis of Stands

```bash
monstand_metrics.r
```

### G1b. Genotypic Diversity of Stands

monostand_diversity.R

### G1c. Genetic Diversity of Stands

**To do:** genotypes per monostand; distribution of genotypes within and among stands.
