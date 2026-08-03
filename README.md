Notebook for initial analysis of Madracis auretenra (DarT)

Data was collected during two expeditions to Curacao, first in 2021 when 493 samples were sequenced, and second in 2023 when 184 samples were collected and sequenced. 

Log into one of the genomics servers: e.g. Alice  
Use rawseq data from deepcat

Set up environment:  
Version information:  
FastQC v0.12.1  
Bwa(0.7.18-r1243-dirty)  
multiqc, version 1.27  
ipyrad(0.9.102)  
vcftools(0.1.16)  
trimgalore(0.6.10-1)  
blastn(2.16.0+)  
structure(2.3.4)  
adegenet (2.1.10)  
poppr(2.9.6)  
 angsd version:  (htslib: 1.17) build(Aug 27 2023 15:34:05)  
newhybrids \-- VERSION: 2.0+ Developmental.  July/August 2007

```py
#set up miniconda env for needed programs:
%  ~/miniconda3/bin/mamba create --name ipyrad_workshop python=2.7
 % source ~/miniconda3/etc/profile.d/conda.sh          
% conda activate ipyrad_workshop
% source ~/miniconda3/etc/profile.d/conda.sh        
%  ~/miniconda3/bin/mamba install -c bioconda multiqc  
%  ~/miniconda3/bin/mamba install -c bioconda fastqc   
```

A0 Setting up working env  
First copy all FASTQ files into one dir so we can work with the whole data instead of ‘22 or ‘24 data alone:

```py
#create new folder in rawseq/MADR_dart/
mkdir MADR_merge
% cp MADR_dart24/RAW/*.FASTQ.gz MADR_merge/
% cp MADR_dart22/RAW/*.FASTQ.gz MADR_merge/

```

# 

# A1a\_data\_prep

Some individuals have been sequenced multiple times first we contatenate the fastq files so that there is only 1 fastq per individual. Run code below to combine fastq files per genotype

```py
#!/bin/bash

# Define input CSV file
csv_file="MADR_targets_22_24.csv"
FASTQ_path="/home/deepcat/rawseq/MADR_dart/MADR_full"
FASTQ_out="/home/deepcat/rawseq/MADR_dart/MADR_merge"

# Loop through unique genotypes in the CSV file
awk -F',' 'NR>1 {print $5}' "$csv_file" | sort | uniq | while read genotype; do
    if [ -z "$genotype" ]; then
        continue # Skip empty genotype values
    fi

    # Get all targetid values for the current genotype
    target_ids=$(awk -F',' -v g="$genotype" 'NR>1 && $5 == g {print $1}' "$csv_file")

    # Create an array of FASTQ files corresponding to these target_ids
    fastq_files=()  # Initialize an empty array
    for target_id in $target_ids; do
        file_path="${FASTQ_path}/${target_id}.FASTQ.gz"
        if [ -f "$file_path" ]; then
            fastq_files+=("$file_path")  # Add the file path to the array
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

```py
ls /home/deepcat/rawseq/MADR_dart/MADR_merge/ | wc -l
#681 total samples

```

# A1b\_qc

```py
% screen -S dennis_ipyradworkshop zsh
% conda activate ipyrad_workshop
# environmental variables:RAW_SEQ_PATH="/home/deepcat/rawseq/MADR_dart"
FASTQ_PATH="/home/deepcat/rawseq/MADR_dart/MADR_merge"
% cd Projects/MADR/A1a_qc/
```

```py
# Run FastQC on all FASTQ files in the directory
% for file in $FASTQ_PATH/*.FASTQ.gz; do fastqc -o ~/Projects/MADR/A1a_QC/fastqc 
 $file; done
# temporarily close screen with CTRL-A-D

# Collate with MultiQC
$ multiqc fastqc/
# found 791 reports
$ zip multiqc*
```

```py
# download locally for viewing:  scp dvanhulten@alice:/home/dvanhulten/Projects/MADR/Ipyrad_trial/A1_fastqc/multiqc_data.zip [location_local]
```

![][image1]

## ![][image2]

## 

![][image3]

Most of my sequences looks like this, they are supposed to run parallel… most sequences look similar in the position of the peaks but differ on peak intensity, this could be biologically relevant but could also be contamination?

![][image4]

## ![][image5]

And warnings for sequence length dist….  
![][image6]  
![][image7]

![][image8]

```py

%  ~/miniconda3/bin/mamba create --name ipyrad_env python=3.9
% conda install ipyrad -c conda-forge -c bioconda
% ~/miniconda3/bin/mamba install -c bioconda trim-galore
```

Use trimgalore to remove sequences with quality lower than 20 and remove reads shorter than 30 bp after initial filtering. We might do more stringent filtering after.

```py
trim_galore --cores 8 -q 20 --length 30 --phred33 -o fastq_trimmed/ $FASTQ_PATH/*.FASTQ.gz
```

Now do another qc check

```py
multiqc trimmed_qc
```

`#file saved @ /Users/dvan216/Documents/Inkfish-Phd/Project_Mir/MADR_dart/multiqc_report_trimmed`

# B1\_ipyrad\_denovo\_assembly

`% cd ~/Projects/MADR/B1_ipyrad_denovo_assembly/ipyrad_run_merged`  
`% ipyrad -n MADR_denovo`  
 `# New file 'params-MADR.txt' created in /home/dvanhulten/Projects/MADR/B1_ipyrad_denovo_assembly/ipyrad_run_merged`  
`% nano params-MADR_denovo.txt` 

`------- ipyrad params file (v.0.9.102)------------------------------------------`  
`MADR_denovo                           ## [0] [assembly_name]: Assembly name. Used to n>`  
`/home/dvanhulten/Projects/MADR/Ipyrad_trial/ipyrad_run1 ## [1] [project_dir]: P>`  
                               `## [2] [raw_fastq_path]: Location of raw non-dem>`  
                               `## [3] [barcodes_path]: Location of barcodes file`  
`/home/dvanhulten/Projects/MADR/A1b_qc/fastq_trimmed/*.fq.gz             >`  
`denovo                         ## [5] [assembly_method]: Assembly method (denov>`  
                               `## [6] [reference_sequence]: Location of referen>`  
`ddrad                            ## [7] [datatype]: Datatype (see docs): rad, g>`  
`TGCAG,                         ## [8] [restriction_overhang]: Restriction overh>`  
`5                              ## [9] [max_low_qual_bases]: Max low quality bas>`  
`33                             ## [10] [phred_Qscore_offset]: phred Q score off>`  
`6                              ## [11] [mindepth_statistical]: Min depth for st>`  
`6                              ## [12] [mindepth_majrule]: Min depth for majori>`  
`10000                          ## [13] [maxdepth]: Max cluster depth within sam>`  
`0.85                           ## [14] [clust_threshold]: Clustering threshold >`  
`0                              ## [15] [max_barcode_mismatch]: Max number of al>`  
`2                              ## [16] [filter_adapters]: Filter for adapters/p>`  
`35                             ## [17] [filter_min_trim_len]: Min length of rea>`  
`2                              ## [18] [max_alleles_consens]: Max alleles per site in consensus sequences`  
`0.05                           ## [19] [max_Ns_consens]: Max N's (uncalled bases) in consensus`  
`0.05                           ## [20] [max_Hs_consens]: Max Hs (heterozygotes) in consensus`  
`4                              ## [21] [min_samples_locus]: Min # samples per locus for output`  
`0.2                            ## [22] [max_SNPs_locus]: Max # SNPs per locus`  
`8                              ## [23] [max_Indels_locus]: Max # of indels per locus`  
`0.5                            ## [24] [max_shared_Hs_locus]: Max # heterozygous sites per locus`  
`0, 0, 0, 0                     ## [25] [trim_reads]: Trim raw read edges (R1>, <R1, R2>, <R2) (see docs)`  
`0, 0, 0, 0                     ## [26] [trim_loci]: Trim locus edges (see docs) (R1>, <R1, R2>, <R2)`  
`p, s, l, v                        	 ## [27] [output_formats]: Output formats (see docs)`  
                               `## [28] [pop_assign_file]: Path to population assignment file`  
                               `## [29] [reference_as_filter]: Reads mapped to this reference are removed in` 

Run ipyrad

`### IMPORTANT!!!! run in a screen to ensure it keeps running!!!!`  
`% screen -S dennis_ipyrad`  
`% ipcluster start -n 30 --daemonize --profile="MADR_ipyrad_dennis_denovo"; sleep 60`  
`% ipyrad -p params-MADR_denovo.txt  -s1234567 -c 30 --ipcluster MADR_ipyrad_dennis_denovo`  
`(ipyrad_env) dvanhulten@alice:~/Projects/MADR/B1a_ipyrad$ ipyrad -p params-MADR.txt -s1234567 -c 30 --ipcluster MADR_ipyrad_dennis`  
 `-------------------------------------------------------------`  
  `ipyrad [v.0.9.102]`  
  `Interactive assembly and analysis of RAD-seq data`  
 `-------------------------------------------------------------`  
  `Parallel connection | alice: 30 cores`  
  `Step 1: Loading sorted fastq data to Samples`  
  `[####################] 100% 0:01:47 | loading reads`  
  `791 fastq files loaded to 791 Samples.`  
  `Step 2: Filtering and trimming reads`  
  `[####################] 100% 0:19:15 | processing reads`  
  `Step 3: Clustering/Mapping reads within samples`  
  `[####################] 100% 0:06:41 | dereplicating`  
  `[####################] 100% 1:33:35 | clustering/mapping`  
  `[####################] 100% 0:00:03 | building clusters`  
  `[####################] 100% 0:00:01 | chunking clusters`  
  `[####################] 100% 5:20:58 | aligning clusters`  
  `[####################] 100% 0:04:20 | concat clusters`  
  `[####################] 100% 0:06:46 | calc cluster stats`  
  `Step 4: Joint estimation of error rate and heterozygosity`  
  `[####################] 100% 0:12:38 | inferring [H, E]`  
  `Step 5: Consensus base/allele calling`  
  `Mean error  [0.00356 sd=0.00066]`  
  `Mean hetero [0.01732 sd=0.00319]`  
  `[####################] 100% 0:05:04 | calculating depths`  
  `[####################] 100% 0:04:57 | chunking clusters`  
  `[####################] 100% 2:48:48 | consens calling`  
  `[####################] 100% 0:10:30 | indexing alleles`  
  `Step 6: Clustering/Mapping across samples`  
  `[####################] 100% 0:03:51 | concatenating inputs`  
  `[####################] 100% 1:26:41 | clustering across`  
  `[####################] 100% 0:02:47 | building clusters`  
  `[####################] 100% 2:02:44 | aligning clusters`  
  `Step 7: Filtering and formatting output files`  
  `[####################] 100% 0:01:29 | applying filters`  
  `[####################] 100% 1:15:05 | building arrays`  
  `[####################] 100% 0:58:42 | writing conversions`  
  `[############        ]  60% 1:56:51 | indexing vcf depths`  

# B2\_ipyrad\_reference\_assembly

```py
%mkdir B2_ipyrad_reference_assembly
%conda activate ipyrad_env
%ipyrad -n MADR_reference
%nano params-MADR_reference
#NOTE, had to rename madr_auretenra_ref.fasta.gz to remove gz 
```

\------- ipyrad params file (v.0.9.102)------------------------------------------  
MADR\_REF                       \#\# \[0\] \[assembly\_name\]: Assembly name. Used to name output directories for assembly steps  
/home/dvanhulten/Projects/MADR/B2\_reference\_assembly \#\# \[1\] \[project\_dir\]: Project dir (made in curdir if not present)  
                               \#\# \[2\] \[raw\_fastq\_path\]: Location of raw non-demultiplexed fastq files  
                               \#\# \[3\] \[barcodes\_path\]: Location of barcodes file  
/home/dvanhulten/Projects/MADR/A1b\_qc/fastq\_trimmed/\*.fq.gz    
                            \#\# \[4\] \[sorted\_fastq\_path\]: Location of demultiplexed/sorted fastq files  
reference                         \#\# \[5\] \[assembly\_method\]: Assembly method (denovo, reference)  
\~/Projects/MADR/B2\_reference\_assembly/madr\_auretenra\_ref.fasta                              \#\# \[6\] \[reference\_sequence\]: Location of reference sequence file  
ddrad                            \#\# \[7\] \[datatype\]: Datatype (see docs): rad, gbs, ddrad, etc.  
TGCAG,                         \#\# \[8\] \[restriction\_overhang\]: Restriction overhang (cut1,) or (cut1, cut2)  
5                              \#\# \[9\] \[max\_low\_qual\_bases\]: Max low quality base calls (Q\<20) in a read  
33                             \#\# \[10\] \[phred\_Qscore\_offset\]: phred Q score offset (33 is default and very standard)  
6                              \#\# \[11\] \[mindepth\_statistical\]: Min depth for statistical base calling  
6                              \#\# \[12\] \[mindepth\_majrule\]: Min depth for majority-rule base calling  
10000                          \#\# \[13\] \[maxdepth\]: Max cluster depth within samples  
0.85                           \#\# \[14\] \[clust\_threshold\]: Clustering threshold for de novo assembly  
0                              \#\# \[15\] \[max\_barcode\_mismatch\]: Max number of allowable mismatches in barcodes  
2                              \#\# \[16\] \[filter\_adapters\]: Filter for adapters/primers (1 or 2=stricter)  
35                             \#\# \[17\] \[filter\_min\_trim\_len\]: Min length of reads after adapter trim  
2                              \#\# \[18\] \[max\_alleles\_consens\]: Max alleles per site in consensus sequences  
0.05                           \#\# \[19\] \[max\_Ns\_consens\]: Max N's (uncalled bases) in consensus  
0.05                           \#\# \[20\] \[max\_Hs\_consens\]: Max Hs (heterozygotes) in consensus  
4                              \#\# \[21\] \[min\_samples\_locus\]: Min \# samples per locus for output  
0.2                            \#\# \[22\] \[max\_SNPs\_locus\]: Max \# SNPs per locus  
8                              \#\# \[23\] \[max\_Indels\_locus\]: Max \# of indels per locus  
0.5                            \#\# \[24\] \[max\_shared\_Hs\_locus\]: Max \# heterozygous sites per locus  
0, 0, 0, 0                     \#\# \[25\] \[trim\_reads\]: Trim raw read edges (R1\>, \<R1, R2\>, \<R2) (see docs)  
0, 0, 0, 0                     \#\# \[26\] \[trim\_loci\]: Trim locus edges (see docs) (R1\>, \<R1, R2\>, \<R2)  
p, s, l, v                        \#\# \[27\] \[output\_formats\]: Output formats (see docs)  
                               \#\# \[28\] \[pop\_assign\_file\]: Path to population assignment file  
                               \#\# \[29\] \[reference\_as\_filter\]: Reads mapped to this reference are removed in step 3

`### IMPORTANT!!!! run in a screen to ensure it keeps running!!!!`

```py
% screen -S dennis_ipyrad
% ipcluster start -n 30 --daemonize --profile="MADR_ipyrad_dennis_reference"; sleep 60
% ipyrad -p params-MADR_reference.txt -s1234567 -c 30 --ipcluster MADR_ipyrad_dennis_reference

```

# C1\_vcf\_analysis

C1a\_vcf\_edit  
Remove .FASTQ.gz\_trimmed from the sample names

```py
awk 'BEGIN{OFS="\t"} /^#CHROM/ {for (i=10; i<=NF; i++) sub(/\.FASTQ.gz_trimmed$/, "", $i)}1' ~/Projects/MADR/B1_ipyrad_denovo_assembly/ipyrad_run_merged/MADR_denovo_outfiles/MADR_denovo.vcf > ~/Projects/MADR/C1a_vcf_edit/MADR_denovo_C1a.vcf
```

C1b\_vcf\_QC

```py
% wc -l ~/Projects/MADR/C1a_vcf_edit/MADR_denovo_C1a.vcf | awk '{print $1-12}'
# 892.416 #number of SNPs 
% wc -l ~/Projects/MADR/C1a_vcf_comparison/MADR_reference_C1a.vcf | awk '{print $1-12}'
# 428.665 #number of SNPs

#Assess sample performance through percentage of missing data
$ /home/pbongaerts/Github/radseq/vcf_missing_data.py ~/Projects/MADR/C1a_vcf_edit/MADR_denovo_C1a.vcf > MADR_denovo_C1b_qc.txt
$ sort -gk 5 -o MADR_denovo_C1b_qc.txt MADR_denovo_C1b_qc.txt
$ awk '
BEGIN {
    less_10 = 0;
    less_30 = 0;
    less_1000_geno = 0;
    total_percent = 0;
    count = 0;
}
NR > 1 {  # Skip header
    # Check percentage genotyped
    if ($5 < 10) less_10++;
    if ($5 < 30) less_30++;
    # Check GENO column
    if ($3 < 1000) less_1000_geno++;
    # Add % genotyped to total and increment sample count
    total_percent += $5;
    count++;
}
END {
    print "Samples with <10% genotyped: " less_10;
    print "Samples with <30% genotyped: " less_30;
    print "Samples with <1000 SNPs (GENO): " less_1000_geno;
    print "Average % genotyped: " total_percent / count;
}
' MADR_denovo_C1b_qc.txt
#total 681 samples
#Samples with <10% genotyped: 134
#Samples with <30% genotyped: 588
#Samples with <1000 SNPs (GENO): 3 
#Average % genotyped: 16.0669


$ /home/pbongaerts/Github/radseq/vcf_missing_data.py MADR_reference_C1a.vcf > MADR_reference_C1b_qc.txt
$ sort -gk 5 MADR_reference_C1b_qc.txt    
#run same awk as before:
#total 681 samples
#Samples with <10% genotyped: 19  
#Samples with <30% genotyped: 293  
#Samples with <1000 SNPs (GENO): 3 
#Average % genotyped: 32.9087 

```

# C2\_removing\_symbiont\_contamination

check for symbiont genome contamination  
Extract a single reference sequence (using the first sample) for each RAD locus:

```py
% /home/pbongaerts/Github/radseq/pyrad2fasta.py ~/Projects/MADR/B1_ipyrad_denovo_assembly/ipyrad_run_merged/MADR_denovo_outfiles/MADR_denovo.loci > MADR_denovo.fasta
% /home/pbongaerts/Github/radseq/pyrad2fasta.py ~/Projects/MADR/B2_reference_assembly/MADR_reference_outfiles/MADR_reference.loci > MADR_reference.fasta
% grep ">" MADR_denovo.fasta | wc -l
# 153.721 #number of loci
% grep ">" MADR_reference.fasta | wc -l
# 53.725 number of loci 

#Remove hyphens for blastn
% sed 's/-//g' MADR_denovo.fasta > MADR_denovo_cleaned.fasta
% sed 's/-//g' MADR_reference.fasta > MADR_reference_cleaned.fasta
```

```py
#align with bwa
% bwa mem -t 32 -M  /home/deepcat/genomes/breviolum/breviolum MADR_denovo_cleaned.fasta > bwa_outputs/MADR_denovo_to_brev.sam

#Exploring sam output:
$ grep -v '^@' bwa_outputs/MADR_denovo_to_brev.sam | wc -l
#Number of reads: 153730 (10 more than in the fasta?)
$ grep -v '^@' bwa_outputs/MADR_denovo_to_brev.sam | awk '$2 == 0' | wc -l
#Mapped reads forward: 362
$ grep -v '^@' bwa_outputs/MADR_denovo_to_brev.sam | awk '$2 == 16' | wc -l
#Mapped reads forward: 354
$ grep -v '^@' bwa_outputs/MADR_denovo_to_brev.sam | awk '$2 == 4' | wc -l
#Unmapped reads: 153005


% bwa mem -t 32 -M  /home/deepcat/genomes/breviolum/breviolum MADR_reference_cleaned.fasta > bwa_outputs/MADR_reference_to_brev.sam

#Exploring sam output:
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | wc -l
#Number of reads: 53733 (more than in the fasta?)
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | awk '$2 == 0' | wc -l
#Mapped reads forward: 96
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | awk '$2 == 16' | wc -l
#Mapped reads forward: 103
$ grep -v '^@' bwa_outputs/MADR_reference_to_brev.sam | awk '$2 == 4' | wc -l
#Unmapped reads: 53526 

###NOTE since we used reference for the ipyrad assembly all sequences that do not match the reference should have been discarded this means that the reference sequence must contain some symbiont sequences?

#Extract a list of succesfully mapped loci with mapping quality of >=20
% /home/pbongaerts/Github/radseq/mapping_get_bwa_matches.py bwa_outputs/MADR_denovo_to_brev.sam > bwa_outputs/MADR_denovo_to_brev_q20.txt
% wc -l bwa_outputs/MADR_denovo_to_brev_q20.txt
#69 loci to remove

% /home/pbongaerts/Github/radseq/mapping_get_bwa_matches.py bwa_outputs/MADR_reference_to_brev.sam > bwa_outputs/MADR_reference_to_brev_q20.txt
% wc -l bwa_outputs/MADR_reference_to_brev_q20.txt
#15 loci to remove

#Repeat for other genomes
GENOMES=("symbiodinium" "cladocopium" "durusdinium")
GENOME_DIR=/home/deepcat/genomes/
INPUT_FASTA=MADR_denovo_cleaned.fasta
OUTPUT_DIR=bwa_outputs

for GENOME in "${GENOMES[@]}"; do
    GENOME_SHORT=${GENOME:0:4}
    
    bwa mem -t 32 -M "${GENOME_DIR}/${GENOME}/${GENOME}" "$INPUT_FASTA" > "${OUTPUT_DIR}/MADR_denovo_to_${GENOME_SHORT}.sam"
    
    /home/pbongaerts/Github/radseq/mapping_get_bwa_matches.py "${OUTPUT_DIR}/MADR_denovo_to_${GENOME_SHORT}.sam" > "${OUTPUT_DIR}/MADR_denovo_to_${GENOME_SHORT}_q20.txt"
    
    echo "Processing for ${GENOME} completed."
done
cat bwa_outputs/MADR_denovo_*.txt | cut -f1 | sort | uniq > bwa_outputs/MADR_denovo_sym_loci_to_remove.txt
wc -l bwa_outputs/MADR_denovo_sym_loci_to_remove.txt
#237 loci to remove

#repeat for reference
wc -l bwa_outputs/reference/MADR_reference_sym_loci_to_remove.txt
#55 loci to remove
```

Compare results with blastn method

```py
#1. run blastn against the databases
$ mkdir blastn_outputs
$ INPUT_FASTA=MADR_denovo_cleaned.fasta
$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue bitscore"
GENOMES=("symbiodinium" "cladocopium" "durusdinium" "breviolum")
$ for GENOME in "${GENOMES[@]}"; do
GENOME_SHORT=${GENOME:0:4}
blastn -query $INPUT_FASTA\
       -db "${GENOME_DIR}/${GENOME}/${GENOME}" \
       -task blastn \
       -outfmt "$OUTPUT_FORMAT" \
       -out "blastn_outputs/MADR_denovo_to_${GENOME_SHORT}.txt"
done

#2. extract matches with E-value lower than 10-15
$ MAX_E_VALUE="0.000000000000001"
for GENOME in "${GENOMES[@]}"; do
GENOME_SHORT=${GENOME:0:4}
python3 /home/pbongaerts/Github/radseq/mapping_get_blastn_matches.py "blastn_outputs/MADR_denovo_to_${GENOME_SHORT}.txt" $MAX_E_VALUE
done
#Symbiodinium
Matches: 24 | Min.length: 62.0 bp | Min. nident: 56.0 bp | Min. pident: 76.471 
#Cladocopium
Matches: 27 | Min.length: 57.0 bp | Min. nident: 54.0 bp | Min. pident: 76.991 
#Durusdinium
Matches: 8 | Min.length: 63.0 bp | Min. nident: 56.0 bp | Min. pident: 77.193 
#Breviolum
Matches: 29 | Min.length: 55.0 bp | Min. nident: 51.0 bp | Min. pident: 75.862
Total of 88 loci that match symbiont genome
```

```py
#1. run blastn against the databases
$ INPUT_FASTA=MADR_reference_cleaned.fasta
$ OUTPUT_FORMAT="7 qseqid sseqid length nident pident evalue bitscore"
GENOMES=("symbiodinium" "cladocopium" "durusdinium" "breviolum")
$ for GENOME in "${GENOMES[@]}"; do
GENOME_SHORT=${GENOME:0:4}
blastn -query $INPUT_FASTA\
       -db "${GENOME_DIR}/${GENOME}/${GENOME}" \
       -task blastn \
       -outfmt "$OUTPUT_FORMAT" \
       -out "blastn_outputs/MADR_reference_to_${GENOME_SHORT}.txt"
done

#2. extract matches with E-value lower than 10-15
$ MAX_E_VALUE="0.000000000000001"
for GENOME in "${GENOMES[@]}"; do
GENOME_SHORT=${GENOME:0:4}
python3 /home/pbongaerts/Github/radseq/mapping_get_blastn_matches.py "blastn_outputs/MADR_reference_to_${GENOME_SHORT}.txt" $MAX_E_VALUE
done
#Symbiodinium 
Matches: 12 | Min.length: 78.0 bp | Min. nident: 66.0 bp | Min. pident: 73.729 %
#Cladocopium  
Matches: 10 | Min.length: 72.0 bp | Min. nident: 62.0 bp | Min. pident: 78.095 %
#Durusdinium   
Matches: 2 | Min.length: 77.0 bp | Min. nident: 73.0 bp | Min. pident: 86.869 %
#Breviolum  
Matches: 15 | Min.length: 65.0 bp | Min. nident: 60.0 bp | Min. pident: 73.171 %
Total of 39 loci that match symbiont genome
```

# C3 Checking for other potential contaminations:

To allow for more queries on NCBI blastn we need an API key  
Key for NCBI dennis: 

```py
c2242b822cd5dc7f5c3f7002d57388cae608
```

```py
#Start a screen
$ screen -S MADR_denovo_blastn

blastn -query ~/Projects/MADR/C2_remove_symbiont_contamination/reference/MADR_reference_cleaned.fasta -db ../blast_db_12_26_2023/nt -task blastn -evalue 0.0001 -max_target_seqs 10 -outfmt "7 qseqid sseqid length nident pident evalue bitscore staxids stitle" -out blastn_outputs/MADR_reference_all.txt -num_threads 32
```

Repeat for denovo  
Looking at the results

```py
grep "Query" blastn_outputs/MADR_denovo_all.txt | wc -l
#153721 Loci
grep "Query" blastn_outputs/MADR_reference_all.txt | wc -l
#53725 Loci

$ python3 ../mapping_identify_blast_matches_modAH.py blastn_outputs/MADR_denovo_all.txt 0.0001 dvanhulten@calacademy.org c2242b822cd5dc7f5c3f7002d57388cae608

$ cut -f 8 blastn_outputs/MADR_denovo_all_match0.0001.txt | sort | uniq -c > MADR_denovo_other_loci_hit_summary.txt
#denovo results
    2 Acidobacteriota
    107 Actinomycetota
    123 Annelida
     55 Apicomplexa
    672 Arthropoda
     46 Ascomycota
     10 Bacillota
     68 Bacteroidota
      9 Basidiomycota
      1 Bdellovibrionota
     11 Bryozoa
      2 Campylobacterota
      1 Chlamydiota
      1 Chloroflexota
    129 Chlorophyta
  28345 Chordata
      2 Chytridiomycota
  31767 Cnidaria
     14 Cyanobacteriota
     39 Echinodermata
     19 ERROR
      6 Euglenozoa
      1 Kitrinoviricota
      3 Methanobacteriota
    248 Mollusca
      1 Myxococcota
      5 Nematoda
     16 NOT_FOUND
      3 Nucleocytoviricota
      1 Phixviricota
      3 Planctomycetota
      9 Platyhelminthes
     33 Porifera
      1 Priapulida
   2826 Pseudomonadota
      1 Rhodophyta
      1 Rhodothermota
      1 Rotifera
      1 Spirochaetota
    157 Streptophyta
      6 Thermodesulfobacteriota
      1 Uroviricota
      1 Verrucomicrobiota
      4 Xenacoelomorpha



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

# 2. Create a list of loci that match phyla outside Cnidaria and NOT_FOUND
$ grep -v -E "Cnidaria|NOT_FOUND" blastn_outputs/MADR_denovo_all_match0.0001.txt > MADR_denovo_other_loci_to_remove.txt
$ wc -l MADR_denovo_other_loci_to_remove.txt
# 32969 loci to remove

```

Repeat for reference:

```py
#674 other loci to remove
```

Remove all contaminants:

```py
# 1. Merge the symbiotns and other contaminants
$ cat ~/Projects/MADR/C2_remove_symbiont_contamination/denovo/bwa_outputs/MADR_denovo_sym_loci_to_remove.txt ~/Projects/MADR/C3_checking_for_other_contamination/denovo/MADR_denovo_other_loci_to_remove.txt | sort | uniq > MADR_denovo_all_loci_to_remove_temp.txt
$ wc -l MADR_denovo_all_loci_to_remove_temp.txt
#33206 loci to remove denovo

#Repeat for reference
#729 loci to remove 

# 2. Adding the prefix RAD for the denovo assembly
$ awk '{print "RAD_" $1}' MADR_denovo_all_loci_to_remove_temp.txt > MADR_denovo_all_loci_to_remove.txt
# for reference selection we need the name of the chrom 
awk '{split($2, acc, "[|.]"); print "ENA|" acc[4] "|" acc[4] ".1";}' MADR_reference_all_loci_to_remove_temp.txt > MADR_reference_all_loci_to_remove.txt

# 3. Remove from vcf file
$ cp ~/Projects/MADR/C1a_vcf_edit/MADR_denovo_C1a.vcf MADR_denovo_c3_temp.vcf

$ /home/pbongaerts/Github/radseq/vcf_remove_chrom.py MADR_denovo_c3_temp.vcf MADR_denovo_all_loci_to_remove.txt > MADR_denovo_c3.vcf


$ vcftools --vcf MADR_denovo_c3.vcf --out output_prefix --site-depth
#After filtering, kept 681 out of 681 Individuals 
#After filtering, kept 682.854 out of a possible 682.854 Sites 
#209.562 sites removed after filtering
```

```py
#Compare with original DART vcf
vcftools --vcf /home/deepcat/rawseq/MADR_dart/MADR_dart24/RAW/Report_DMadr24-8994_1_moreOrders_SNP_vcf_2.csv 
 --out output_prefix --site-depth
#After filtering, kept 681 out of 681 Individuals  
#After filtering, kept 43.903 out of a possible 43.903 Sites
```

# D1\_initial\_NJ\_tree

### D1a\_removal\_of\_underperforming\_ind\_and\_snps

Start with qc control of snps and the vcf

```py
~/radseq/vcf_pos_count.py denovo/MADR_denovo_c3.vcf
1-5	30758	*********************************************************************
6-10	31422	**********************************************************************
11-15	32028	************************************************************************
16-20	32453	*************************************************************************
21-25	32640	*************************************************************************
26-30	32785	*************************************************************************
31-35	32734	*************************************************************************
36-40	32334	************************************************************************
41-45	31400	**********************************************************************
46-50	29754	*******************************************************************
51-55	29221	*****************************************************************
56-60	28724	****************************************************************
61-65	27533	**************************************************************
66-70	26976	************************************************************
71-75	26352	***********************************************************
76-80	25500	*********************************************************
81-85	25362	*********************************************************
86-90	25468	*********************************************************
91-95	23749	*****************************************************
96-100	23249	****************************************************
101-105	22765	***************************************************
106-110	22806	***************************************************
111-115	23724	*****************************************************
116-120	25002	********************************************************
121-125	10175	**********************
126-130	1034	**
131-135	220	
136-140	163	
141-145	159	
146-150	123	
151-155	103	
156-160	102	
161-165	65	
166-170	74	
171-175	37	
176-180	12	
181-185	12	
186-190	2	
191-191	1	

vcf_pos_count_MODref.py reference/MADR_reference_c3.vcf > ref_pos_count.txt

awk '{count[$2]++} END {for (value in count) print value ": " count[value]}' ref_pos_count.txt | sort -n4: 1
5: 82.245
6: 1.855
7: 176
8: 10
9: 1

```

```py
# %missing data
$ ~/radseq/vcf_missing_data.py MADR_denovo_c3.vcf > MADR_denovo_preformance_c3.txt

$ sort -gk 5 -o MADR_denovo_preformance_c3.txt MADR_denovo_preformance_c3.txt
#INDIVIDUAL      MISS    GENO    TOTAL   % GENOTYPED
#DH0360_MMIR_SNA 682852  2       682854  0.0
#DH0293_MMIR_SNA 682557  297     682854  0.04
#DH0371_MMIR_SNA 682448  406     682854  0.06
#DH0373_MMIR_SNA 663268  19586   682854  2.87
#DH0308_MMIR_SNA 655393  27461   682854  4.02
#DH0552_MMIR_SEA 652152  30702   682854  4.5
#DH0277_MMIR_SNA 651800  31054   682854  4.55
#DH0623_MMIR_KAL 651303  31551   682854  4.62
#DH0627_MMIR_KAL 651214  31640   682854  4.63

$ awk '
BEGIN {
    less_10 = 0;
    less_30 = 0;
    less_1000_geno = 0;
    total_percent = 0;
    count = 0;
}
NR > 1 {  # Skip header
    # Check percentage genotyped
    if ($5 < 10) less_10++;
    if ($5 < 30) less_30++;
    # Check GENO column
    if ($3 < 1000) less_1000_geno++;
    # Add % genotyped to total and increment sample count
    total_percent += $5;
    count++;
}
END {
    print "Samples with <10% genotyped: " less_10;
    print "Samples with <30% genotyped: " less_30;
    print "Samples with <1000 SNPs (GENO): " less_1000_geno;
    print "Average % genotyped: " total_percent / count;
}
' MADR_denovo_preformance_c3.txt
#Samples with <10% genotyped: 82
#Samples with <30% genotyped: 587
#Samples with <1000 SNPs (GENO): 3
#Average % genotyped: 19.4148

# save all names with <1000 SNPS genotyped to remove file
$ awk '$3 < 1000 {print $1}' MADR_denovo_preformance_c3.txt > indv_to_remove_denovo_c3.txt
3 individuals to remove

#repeat for reference
% ~/radseq/vcf_missing_data.py MADR_reference_C3.vcf > MADR_reference_preformance_C3.txt
% sort -gk 5 -o MADR_reference_preformance_C3.txt MADR_reference_preformance_C3.txt
#INDIVIDUAL	MISS	GENO	TOTAL	% GENOTYPED
#DH0293_MMIR_SNA	418633	0	418633	0.0
#DH0360_MMIR_SNA	418633	0	418633	0.0
#DH0371_MMIR_SNA	418614	19	418633	0.0
#DH0373_MMIR_SNA	400754	17879	418633	4.27
#DH0308_MMIR_SNA	392283	26350	418633	6.29
#DH0277_MMIR_SNA	388740	29893	418633	7.14
#DH0627_MMIR_KAL	388546	30087	418633	7.19
#DH0552_MMIR_SEA	388343	30290	418633	7.24
#DH0553_MMIR_SEA	388260	30373	418633	7.26

#repeat sort +awk command above:

#Samples with <10% genotyped: 19
#Samples with <30% genotyped: 293
#Samples with <1000 SNPs (GENO): 3 
#Average % genotyped: 32.9087 
#less than 1000 snps
INDIVIDUAL      MISS    GENO    TOTAL   % GENOTYPED
DH0293_MMIR_SNA 428665  0       428665  0.0
DH0360_MMIR_SNA 428665  0       428665  0.0
DH0371_MMIR_SNA 428649  16      428665  0.0
DH0371_MMIR_SNA_x       428661  4       428665  0.0

# save all names with <1000 SNPS genotyped to remove file
$ awk '$3 < 1000 {print $1}' MADR_reference_preformance_c3.txt > indv_to_remove_reference_c3.txt
3 individuals to remove
```

remove underperforming snps and individuals

```py

$ vcftools --vcf MADR_denovo_c3.vcf --remove indv_to_remove_denovo_c3.txt --max-missing 0.5 --recode --stdout > MADR_denovo_d1_temp.vcf
#After filtering, kept 678 out of 681 Individuals
#After filtering, kept 96.825 out of a possible 682.854 Sites
 
$ vcftools --vcf MADR_reference_c3.vcf --remove indv_to_remove_reference_C3.txt --max-missing 0.5 --recode --stdout > MADR_reference_d1_temp.vcf
#After filtering, kept 678 out of 791 Individuals
#After filtering, kept 120.597 out of a possible 418.633 Sites

#testing dart file
$ vcftools --vcf MADR_dart_d1.vcf 
#After filtering, kept 661 out of 661 Individuals
#After filtering, kept 43.903 out of a possible 43.903 Sites
```

After meeting with Ale: might have too many sites left, will drag out the STRUCTURE run. Can do more aggressive filtering since we have enough sites that should help pull appart the two clades

```py
$ vcftools --vcf MADR_denovo_c3.vcf --remove indv_to_remove_denovo_c3.txt --max-missing 0.8 --mac 1 --minDP 8 --recode --stdout > MADR_denovo_d1_strict.vcf
#After filtering, kept 678 out of 681 Individuals
#After filtering, kept 25.995 out of a possible 682854 Sites
 
$ vcftools --vcf MADR_reference_c3.vcf --remove indv_to_remove_reference_C3.txt --max-missing 0.8 --mac 1 --minDP 8 --recode --stdout > MADR_reference_d1_strict.vcf
#After filtering, kept 678 out of 681 Individuals
#After filtering, kept 34.466 out of a possible 418633 Sites

#testing dart file
$ vcftools --vcf MADR_dart_d1.vcf 
#After filtering, kept 661 out of 661 Individuals
#After filtering, kept 43.903 out of a possible 43.903 Sites
```

Create a popfile using modded script to incorporate both sna and sna\_sha for the shallow samples

```py
###repeat 
$ ~/radseq/vcf_missing_data.py MADR_denovo_d1_temp.vcf > MADR_denovo_performance_d1_temp.txt
$ sort -gk 5 -o MADR_denovo_performance_d1_temp.txt MADR_denovo_performance_d1_temp.txt
INDIVIDUAL      MISS    GENO    TOTAL   % GENOTYPED
DH0373_MMIR_SNA 88288   8537    96825   8.82
DH0308_MMIR_SNA 81893   14932   96825   15.42
DH0552_MMIR_SEA 81186   15639   96825   16.15
DH0311_MMIR_SNA 78893   17932   96825   18.52
DH0420_MMIR_SNA 78651   18174   96825   18.77


$ ~/radseq/vcf_missing_data.py MADR_reference_d1_temp.vcf > MADR_reference_performance_d1_temp.txt
$ sort -gk 5 -o MADR_reference_performance_d1_temp.txt MADR_reference_performance_d1_temp.txt
INDIVIDUAL      MISS    GENO    TOTAL   % GENOTYPED
DH0373_MMIR_SNA 110633  9964    120597  8.26
DH0308_MMIR_SNA 103181  17416   120597  14.44
DH0552_MMIR_SEA 101742  18855   120597  15.63
DH0311_MMIR_SNA 99322   21275   120597  17.64
DH0553_MMIR_SEA 98731   21866   120597  18.13


#test dart vcf as well
$ ~/radseq/vcf_missing_data.py dart/MADR_dart_d1.vcf > dart/MADR_dart_performance_d1.txt
$ sort -gk 5 -o dart/MADR_dart_performance_d1.txt dart/MADR_dart_performance_d1.txt
#Samples with <10% genotyped: 0
#Samples with <30% genotyped: 30
#Samples with <1000 SNPs (GENO): 0
#Average % genotyped: 61.9826


#rename MADR_[method]_d1_temp.vcf to MADR_[method]_d1.vcf 
```

 Try the gd matrix again, at this point we will also include the Dart vcf file

```py
% python3 ~/scripts/popfile_from_vcf_MODdennis.py MADR_denovo_d1.vcf 13 > MADR_denovo_popfile_d1.txt

% python3 ~/scripts/popfile_from_vcf_MODdennis.py MADR_reference_d1.vcf 13 > MADR_reference_popfile_d1.txt

% python3 ~/scripts/popfile_from_vcf_MODdennis.py dart/MADR_dart_d1.vcf 13 > dart/MADR_dart_popfile_d1.txt

#Build the GD matrix
~/radseq/vcf_gdmatrix.py MADR_denovo_d1.vcf MADR_denovo_popfile_d1.txt > MADR_denovo_gd_d1.txt 

~/radseq/vcf_gdmatrix.py MADR_reference_d1.vcf MADR_reference_popfile_d1.txt > MADR_reference_gd_d1.txt

~/radseq/vcf_gdmatrix.py dart/MADR_dart_d1.vcf dart/MADR_dart_popfile_d1.txt > dart/MADR_dart_gd_d1.txt


```

### D1b\_building\_initial\_nj\_trees

Dart tree

```py
$ ~/radseq/gdmatrix2tree.py MADR_denovo_gd_d1.txt MADR_denovo_d1.tre

$ ~/radseq/gdmatrix2tree.py MADR_reference_gd_d1.txt MADR_reference_d1.tre

$ ~/radseq/gdmatrix2tree.py ../../D1a_removal_of_underperforming_ind_and_snps/dart/MADR_dart_gd_d1.txt MADR_dart_d1.tre
```

![][image9]  
\#61 out of 661 individuals in the minor/blue clade   
\#2 individuals seem to fall between the two clades   
\#600 out of 661 individuals fall in the major or red clade 

QC of individuals using the tree

```py
#repeat for all vcf variants
vcftools --vcf MADR_reference_d1.vcf --het
vcftools --vcf MADR_reference_d1.vcf --missing-indv
vcftools --vcf MADR_reference_d1.vcf --freq
vcftools --vcf MADR_reference_d1.vcf --depth
```

Adding depth to the tree:

```py
% cat cur_kal_med.co2.csv cur_sna_med_full.co2.csv cur_sea_med_full.co2.csv > cur_annotations_full.csv

% awk -F',' 'NR==FNR {map[$4] = $2; next} {print $0, (map[$6] ? map[$6] : "NA")}' OFS=',' 
MADR_sample_record.csv cur_annotations_full.csv > cur_annotations_temp.csv


```

For figures, see Tree\_figures\_d1.R

![][image10]

Denovo tree

![][image11]

![][image12]

Reference tree  
Interesting observation: the “Hybrid” samples are always 2 samples but in the dart and denovo tree are grouped with the Major clade, however, in the reference tree these 2 samples are grouped with the Minor clade.![][image13]  
![][image14]

### D1c\_genetic\_structure

Initial clone detection  
\#repeate methods below for both Denovo and Dart data

```py
$ python3 ~/scripts/popfile_from_vcf_MODdennis.py MADR_reference_d1.vcf 13 > MADR_reference_popfile_d1.txt

# 1. Calculate allelic similarity 
$ screen -S clonality
$ ~/scripts/vcf_clone_detect_npMOD.py -v MADR_reference_d1.vcf -p MADR_reference_popfile_d1.txt -o MADR_reference_d1_clones.txt

Error: clone matrix not properly sorted
```

Check highest scores for replicates:

```py
awk '
{
    match($0, /DH([0-9]+)/, dh);
    match($0, /DX([0-9]+)/, dx);
    if (dh[1] == dx[1]) print $0;
}' clone_detect.out 
#99.81   0.0     [SNA]   DH0453_MMIR_SNA vs DX0453_MMIR_SNA      92485.0/92660   114223  99034
#99.79   0.0     [SNA]   DH0452_MMIR_SNA vs DX0452_MMIR_SNA      92828.0/93022   112484  101135
#99.69   0.0     [SNA_SHA]       DH0054_MMIR_SNA_SHA vs DX0054_MMIR_SNA_SHA      101450.0/101767 117210  105154
#99.67   0.0     [KAL_SHA]       DH0114_MMIR_KAL_SHA vs DX0114_MMIR_KAL_SHA      101121.5/101460 117087  104970
#99.12   0.0     [SNA]   DH0454_MMIR_SNA vs DX0454_MMIR_SNA      57379.0/57889   117749  60737

```

Based on these results we can put a threshold around 98% similarity for our initial clone detection

```py
python3 ~/scripts/detect_clones_vcf.py -v MADR_reference_d1.vcf -p MADR_reference_popfile_d1.txt -t 98 > MADR_reference_d1_cloneout_98.txt
```

Use the output to select clonal lineages (between \#\#\#4 and \#\#\#5)

```py
sed -n '/###4/,/###5/p' MADR_reference_d1_cloneout_98.txt > MADR_reference_d1_clonal_groups.txt
```

Check if any clonal groups belong to different populations

```py
$ python3 ~/scripts/get_populations_clonal_groups.py MADR_reference_d1_clonal_groups.txt 
```

`#reference results`  
`###4 - Clonal groups (threshold: 98.0)`  
`KAL_SHA-SNA_SHA: DH0002_MMIR_SNA_SHA, DH0087_MMIR_KAL_SHA (99.82 %)`  
`KAL_SHA-SNA_SHA: DH0005_MMIR_SNA_SHA, DH0086_MMIR_KAL_SHA, DH0074_MMIR_SNA_SHA (99.79-99.8 %)`  
`SNA_SHA_x-KAL: DH0002_MMIR_SNA_SHA_x, DH0756_MMIR_KAL (99.77 %)`  
`KAL_SHA-SNA_SHA: DH0074_MMIR_SNA_SHA, DH0015_MMIR_SNA_SHA, DH0086_MMIR_KAL_SHA, DH0005_MMIR_SNA_SHA, DH0004_MMIR_SNA_SHA, DH0007_MMIR_SNA_SHA, DH0010_MMIR_SNA_SHA, DH0003_MMIR_SNA_SHA, DH0073_MMIR_SNA_SHA, DH0012_MMIR_SNA_SHA, DH0006_MMIR_SNA_SHA, DH0011_MMIR_SNA_SHA, DH0085_MMIR_SNA_SHA (99.88-99.95 %)`  
`KAL-SEA: DH0582_MMIR_SEA, DH0591_MMIR_SEA, DH0585_MMIR_SEA, DH0586_MMIR_SEA, DH0584_MMIR_SEA, DH0595_MMIR_KAL, DH0590_MMIR_SEA, DH0541_MMIR_SEA (99.58-99.92 %)`

\#denovo results  
`###4 - Clonal groups (threshold: 98.0)`  
`SNA_SHA-KAL_SHA-KAL_SHA_x-SNA_SHA_x: DH0085_MMIR_SNA_SHA, DH0015_MMIR_SNA_SHA, DH0005_MMIR_SNA_SHA, DH0006_MMIR_SNA_SHA, DH0086_MMIR_KAL_SHA_x, DH0012_MMIR_SNA_SHA_x, DH0086_MMIR_KAL_SHA, DH0074_MMIR_SNA_SHA, DH0003_MMIR_SNA_SHA_x, DH0007_MMIR_SNA_SHA, DH0074_MMIR_SNA_SHA_x (99.61-99.81 %)`  
`SNA_SHA-KAL_SHA-KAL-SNA_SHA_x: DH0002_MMIR_SNA_SHA_x, DH0087_MMIR_KAL_SHA, DH0002_MMIR_SNA_SHA, DH0756_MMIR_KAL (99.74-99.81 %)`  
`SEA-KAL: DH0541_MMIR_SEA, DH0595_MMIR_KAL, DH0591_MMIR_SEA, DH0586_MMIR_SEA, DH0585_MMIR_SEA, DH0582_MMIR_SEA, DH0584_MMIR_SEA, DH0590_MMIR_SEA (99.58-99.93 %)`  
`SNA_SHA-KAL_SHA: DH0086_MMIR_KAL_SHA, DH0006_MMIR_SNA_SHA, DH0074_MMIR_SNA_SHA, DH0085_MMIR_SNA_SHA, DH0015_MMIR_SNA_SHA (99.89-99.91 %)`  
`SNA_SHA-KAL-KAL_SHA: DH0002_MMIR_SNA_SHA, DH0753_MMIR_KAL, DH0087_MMIR_KAL_SHA (99.88-99.91 %)`

Two clonal groups shared between locations, DH0086 and DH0087\_KAL\_SHA seem to be the culprits (see results E1)

\#dart results  
`KAL_SHA-SNA_SHA: DH0010_MMIR_SNA_SHA, DH0085_MMIR_SNA_SHA, DH0086_MMIR_KAL_SHA, DH0073_MMIR_SNA_SHA, DH0004_MMIR_SNA_SHA, DH0015_MMIR_SNA_SHA, DH0006_MMIR_SNA_SHA, DH0012_MMIR_SNA_SHA, DH0005_MMIR_SNA_SHA, DH0003_MMIR_SNA_SHA, DH0011_MMIR_SNA_SHA, DH0007_MMIR_SNA_SHA, DH0074_MMIR_SNA_SHA (98.17-99.52 %)`  
Again DH0086\_KAL\_SHA in the data

```py
$ python3 ~/scripts/add_clonal_group_to_popfile.py MADR_reference_popfile_d1.txt MADR_reference_d1_clonal_groups.txt MADR_reference_d1_popfile_temp.txt
```

Visualise clonal groups in the tree using R

Initial results look good enough to remove these individuals from the data and run structure. 

Get the recommended clones to remove based on %missing from the cloneout file.

```py
sed -n '/###5/,$p' MADR_reference_d1_cloneout_98.txt > MADR_reference_d1_clones_to_remove.txt
```

 Remove clones from vcf

```py
vcftools --vcf MADR_reference_d1.vcf --remove MADR_reference_d1_clones_to_remove.txt --recode --out MADR_reference_d1_no_clones
#kept 302 individuals -denovo
#kept 302 individuals -dart
#kept 300 individuals -reference
```

Run through NJ tree building steps again for the data without clones (See D1b)

Run Structure

```py
$ scp /home/deepcat/STRUCTURE_params/extraparams .
$ scp /home/deepcat/STRUCTURE_params/mainparams .

$ vcftools --vcf ../MADR_reference_d1_no_clones.vcf
#After filtering, kept 309 out of 309 Individuals
#After filtering, kept 121930 out of a possible 121930 Sites 

$ nano mainparams
KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.
"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!)
Basic Program Parameters
#define MAXPOPS    7      // (int) number of populations assumed
#define BURNIN    100000   // (int) length of burnin period
#define NUMREPS   50000   // (int) number of MCMC reps after burnin
Input/Output files
#define INFILE   infile   // (str) name of input data file
#define OUTFILE  outfile  //(str) name of output data file
Data file format
#define NUMINDS    296    // (int) number of diploid individuals in data file
#define NUMLOCI    34466    // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line
#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says 
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier
#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data 
                             before the genotype data start.
#define MARKERNAMES      0  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances 
                            // between loci
Advanced data file options
#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data

#set up interactive shell with 8 cpu's and 2gb ram

$ ~/radseq/structure_mp.py MADR_reference_d1_strict_no_clones.vcf MADR_reference_d1_strict_popfile_no_clones.txt 4 10 32

K = 2: MedMeaK 2.0 MaxMeaK 2 MedMedK 2.0 MaxMedK 2      MADR_reference_d1_no_clones.vcf
K = 3: MedMeaK 1.0 MaxMeaK 2 MedMedK 1.0 MaxMedK 2      MADR_reference_d1_no_clones.vcf
K = 4: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      MADR_reference_d1_no_clones.vcf
K = 5: MedMeaK 0.0 MaxMeaK 1 MedMedK 0.0 MaxMedK 1      MADR_reference_d1_no_clones.vcf
K = 6: MedMeaK 0.0 MaxMeaK 1 MedMedK 0.0 MaxMedK 1      MADR_reference_d1_no_clones.vcf
K = 7: MedMeaK 0.0 MaxMeaK 1 MedMedK 0.0 MaxMedK 1      MADR_reference_d1_no_clones.vcf
```

Strange structure results for ref assembly; could it be due to missing data/lenient filtering?  
The single line in the major clade seems to be a 3 sample clonal pair…

\#\#denovo  
After filtering, kept 300 out of 300 Individuals  
After filtering, kept 98416 out of a possible 98416 Sites  
Initialise indivs and pops for MADR\_denovo\_d1\_strict\_no\_clones.vcf...  
Subsample SNPs (one random SNP per locus)... \[4416 SNPs/loci\]  
Outputting 10 STRUCTURE files...10 reps DONE  
Executing 32 parallel STRUCTURE runs for K \= 2 ...10 reps DONE  
Executing 32 parallel STRUCTURE runs for K \= 3 ...10 reps DONE  
Executing 32 parallel STRUCTURE runs for K \= 4 ...10 reps DONE  
Running CLUMPP on replicates for K \= 2 ...  
Running CLUMPP on replicates for K \= 3 ...  
Running CLUMPP on replicates for K \= 4 ...  
K \= 2: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      MADR\_denovo\_d1\_strict\_no\_clones.vcf  
K \= 3: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      MADR\_denovo\_d1\_strict\_no\_clones.vcf  
K \= 4: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      MADR\_denovo\_d1\_strict\_no\_clones.vcf

\#\#\# repeat for dart to see if problem with program or data…  
After filtering, kept 302 out of 302 Individuals  
After filtering, kept 43903 out of a possible 43903 Sites

For strict dataset:

Subsample SNPs (one random SNP per locus)... \[4957 SNPs/loci\]  
Outputting 10 STRUCTURE files...10 reps DONE  
Executing 32 parallel STRUCTURE runs for K \= 2 ...10 reps DONE  
Executing 32 parallel STRUCTURE runs for K \= 3 ...10 reps DONE  
Executing 32 parallel STRUCTURE runs for K \= 4 ...10 reps DONE  
Running CLUMPP on replicates for K \= 2 ...  
Running CLUMPP on replicates for K \= 3 ...  
Running CLUMPP on replicates for K \= 4 ...  
K \= 2: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      MADR\_reference\_d1\_strict\_no\_clones.vcf  
K \= 3: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      MADR\_reference\_d1\_strict\_no\_clones.vcf  
K \= 4: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      MADR\_reference\_d1\_strict\_no\_clones.vcf

Running snapclust and structure on d1\_no\_clones snapclust.R

```py
#warning
> MADR_snap_k2 <- snapclust(MADR_genind, 2)
Large dataset syndrome:
 for 308 individuals, differences in log-likelihoods exceed computer precision;
 group membership probabilities are approximated
 (only trust clear-cut values)
```

![][image15]  
![][image16]

Check to see if we get the same assignment of individuals from the denovo and dart vcfs.

```py
#use pop assignment script to assign populations based on structure output
python3 ~/scripts/MADR_compare_strucout.py MADR_reference_d1_strict_no_clones_1739223814/*.out.csv

Compare with clades based on NJ tree using ITOL
### found no inconsistencies in samples beloning to minor or major clade assignments line up with the NJ tree for reference assembly and for Denovo assembly 
```

```py
%awk -F ': ' '{gsub(/\(.*\)/, "", $2); gsub(/ /, "", $2); print $2}' MADR_reference_d1_strict_clonal_groups.txt > MADR_reference_d1_strict_clonal_groups_clean.txt

% python3 map_clones_to_popfile.py MADR_reference_d1_strict_clonal_groups_clean.txt MADR_reference_d1_strict_no_clones_popfile_strucK2.txt MADR_reference_d1_strict_popfile_strucK2_with_clones.txt%awk '{if ($4 == "Major_Clade") print $1 > "MADR_reference_d1_major_clade.txt"; else if ($4 == "Minor_Clade") print $1 > "MADR_reference_d1_minor_clade.txt"}' MADR_reference_d1_strict_popfile_strucK2_with_clones.txt613 MADR_reference_d1_major_clade.txt
68 MADR_reference_d1_minor_clade.txt
263 MADR_reference_d1_major_no_clones
33 MADR_reference_d1_minor_no_clones

613 MADR_denovo_d1_major_clade.txt
65 MADR_denovo_d1_minor_clade.txt


```

Split the vcf into two based on the population assignment

```py
#remove all minor clade samples
vcftools --vcf MADR_reference_d1_strict.vcf --remove MADR_reference_d1_minor_clade.txt --recode --out MADR_reference_d1_strict_major
#After  filtering kept 701 out of 786 individuals

#remove all major clade samples
vcftools --vcf MADR_reference_d1_strict.vcf --remove MADR_reference_d1_major_clade.txt --recode --out MADR_reference_d1_strict_minor
#After  filtering kept 83 out of 786 individuals
```

# 

# E1\_Clone\_detection

### E1a\_defining\_clonal\_threshold

First we will have a look at the data using R and poppr version 2.9.6 to get an understanding of clonality in the dataset. See Clone\_detection.R script.

For clone detection it is important to filter so that there are no monomorphic sites (sites that share the same state over all samples) as these sites will always count toward genetic similarity, and make sure we have a low % of missing data. Because the minor clade is so much smaller we should go back to the original data before filtering, pull the clades appart and then filter again.

```py
$ vcftools --vcf ~/Projects/MADR/C3_checking_for_other_contamination/reference/MADR_reference_C3.vcf --remove ~/Projects/MADR/D1c_genetic_structure/reference/strict/MADR_reference_d1_major_clade.txt --recode --out MADR_reference_e1_minor_full_temp
#After filtering, kept 68 out of 681 Individuals
#Outputting VCF file... 
#After filtering, kept 418633 out of a possible 418633 Sites 
#Still have to remove underperforming individuals from D1 
$ vcftools --vcf MADR_reference_e1_minor_full.vcf --remove  ~/Projects/MADR/D1a_removal_of_underperforming_ind_and_snps/reference/indv_to_remove_reference_c3.txt --recode --out MADR_reference_e1_minor_full_temp.vcf 
#After filtering, kept 65 out of 68 Individuals
$ vcftools --vcf ~/Projects/MADR/C3_checking_for_other_contamination/reference/MADR_reference_C3.vcf --remove ~/Projects/MADR/D1c_genetic_structure/reference/strict/MADR_reference_d1_minor_clade.txt --recode --out MADR_reference_e1_major_full
#After filtering, kept 616 out of 681 Individuals
#Outputting VCF file... 
#After filtering, kept 418633 out of a possible 418633 Sites
$ vcftools --vcf MADR_reference_e1_major_full.vcf --remove  ~/Projects/MADR/D1a_removal_of_underperforming_ind_and_snps/reference/indv_to_remove_reference_c3.txt --recode --out MADR_reference_e1_major_full_temp.vcf
#After filtering, kept 613 out of 616 Individuals

```

```py
vcftools --vcf MADR_reference_d1_strict_minor.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_e1_minor
After filtering, kept 65 out of 65 Individuals
After filtering, kept 12199 out of a possible 34466 Sites

vcftools --vcf MADR_reference_d1_strict_major.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_e1_major
After filtering, kept 613 out of 613 Individuals
After filtering, kept 28409 out of a possible 34466 Sites
##############################################################################
vcftools --vcf MADR_reference_e1_minor_full.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_e1_minor_full_filtered
After filtering, kept 65 out of 65 Individuals
After filtering, kept 14752 out of a possible 418633 Sites 

vcftools --vcf MADR_reference_e1_major_full.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_e1_major_full_filtered
After filtering, kept 613 out of 613 Individuals
After filtering, kept 30334 out of a possible 418633 Sites

```

```py
$ python3 ~/scripts/popfile_from_vcf_MODdennis.py MADR_reference_e1_minor.vcf 13 > MADR_reference_popfile_e1_minor.txt

$ python3 ~/scripts/popfile_from_vcf_MODdennis.py MADR_reference_e1_major.vcf 13 > MADR_reference_popfile_e1_major.txt

# 1. Calculate allelic similarity 
$ screen -S clonality

python3 ~/scripts/detect_clones_vcf.py -v MADR_reference_e1_minor_full.vcf -p MADR_reference_popfile_e1_minor.txt -o MADR_reference_e1_minor_ful_genetic_sim.txt > MADR_reference_e1_minor_full_clonedet_out2.txt

python3 ~/scripts/detect_clones_vcf.py -v MADR_reference_e1_major_full.vcf -p MADR_reference_popfile_e1_major.txt -o MADR_reference_e1_major_full_genetic_sim.txt > MADR_reference_e1_major_full_clonedet_out2.txt
```

Building trees for each clade

```py
$ ~/radseq/vcf_gdmatrix.py MADR_reference_e1_minor_full.vcf MADR_reference_popfile_e1_minor.txt > MADR_reference_e1_minor_full_gdmat.txt
$ ~/radseq/gdmatrix2tree.py MADR_reference_e1_minor_full_gdmat.txt MADR_reference_e1_minor_full.tre

$ ~/radseq/vcf_gdmatrix.py MADR_reference_e1_major_full.vcf MADR_reference_popfile_e1_major.txt > MADR_reference_e1_major_full_gdmat.txt
$ ~/radseq/gdmatrix2tree.py MADR_reference_e1_major_full_gdmat.txt MADR_reference_e1_major_full.tre
```

Clone detection using clone\_detect.R for Major clade with the general filtered data  
![][image17]  
Most appropriate threshold for Major clades seems to be the Average threshold dictated by poppr pictured in Blue in the graph above. Important notice: the purple line indicates the genetic similarity of the replicate pair with the highest distance…..  
Visualise the genetic similarity between pairs in a histogram to determine the appropriate threshold. Average similarity and average threshold respectively for each measures seem to lie around the same point, just after the last big peak and might be useful as a threshold. However Using the lowest similarity between known replicates could be an easy, more lenient way to select clones as this value lies much lower in both cases. This would lead to a threshold of 99% similarity and a distance threshold of 0.002748542 for the 28409 site dataset. 

Comparing clone calls between the two methods:

```py
 python3 ~/scripts/detect_clones_vcf.py -v MADR_reference_e1_major.vcf -p MADR_reference_popfile_e1_major.txt -o MADR_reference_e1_major_genetic_sim2_t99.txt -t 99 > MADR_reference_e1_major_clonedet_out2_t99.txt

```

Results based on the clade specific filtered data:  
![][image18]  
![][image19]  
   
\#minor clade

![][image20]  
![][image21]

Given the results we have decided to go with the lowest genetic similarity of the replicates and set the threshold for the definition of clones at 99% similarity

```py
clonal_pairs="DH0452_MMIR_SNA DX0452_MMIR_SNA
DH0453_MMIR_SNA DX0453_MMIR_SNA
DH0454_MMIR_SNA DX0454_MMIR_SNA
DH0054_MMIR_SNA_SHA DX0054_MMIR_SNA_SHA
DH0114_MMIR_KAL_SHA DX0114_MMIR_KAL_SHA"

# Ensure the while loop reads each line correctly
while IFS=$'\n' read -r pair; do
    # Split the pair into two variables
    read -r ind1 ind2 <<< "$pair"

    # Search for the matching pair in the file
    awk -F',' -v i1="$ind1" -v i2="$ind2" '
        ($1 == i1 && $2 == i2) || ($1 == i2 && $2 == i1) {
            print $0
        }
    ' MADR_reference_e1_major_full_genetic_sim.txt
done <<< "$clonal_pairs"

# ind1,ind2,ind1_snps,ind2_snps,both_snps,match,match_perc,pop
#DH0452_MMIR_SNA,DX0452_MMIR_SNA,29965,29507,29138,29116.000000,99.920000,SNA
#DH0453_MMIR_SNA,DX0453_MMIR_SNA,30056,29200,28922,28907.500000,99.950000,SNA 
#DH0454_MMIR_SNA,DX0454_MMIR_SNA,30084,22045,21795,21579.500000,99.010000,SNA 
#DH0054_MMIR_SNA_SHA,DX0054_MMIR_SNA_SHA,30083,29411,29160,29116.000000,99.850000,SNA_SHA 
#DH0114_MMIR_KAL_SHA,DX0114_MMIR_KAL_SHA,30013,29246,28925,28865.500000,99.790000,KAL_SHA


```

With the initial threshold of 99.01% similarity we can find this interesting laderised pattern in the NJ tree. Upon further inspection combining the missing data we can see there is a clear overlap between these “difficult cases” and missing data. This is likely caused by the small denominator effect where due to the lower shared sites each “error” carries more weight compared to complete data.

To get a better understanding of the sequencing error rate of our data we will look into the samples sequenced multiple times by DArT.

First we make a list of all samples that were sequenced multiple times.

```py
input_file="$HOME/Projects/MADR/A1a_data_prep/MADR_targets_22_24.csv"
output_file="MADR_targets_multi_run.csv"

awk -F, '
BEGIN {
    OFS = ","
}
NR == 1 {
    print $0  # Print header as is
    next
}
{
    count[$5]++
    if (count[$5] > 1) {
        $5 = $5 "_x" count[$5]  # Correct way to append "_x{number}"
    }
    print $0
}' "$input_file" > "$output_file"


```

Renaming the FASTQ\_files

```py
#!/bin/bash

# Define input CSV file
csv_file="MADR_targets_multi_run.csv"
FASTQ_path="/home/deepcat/rawseq/MADR_dart/MADR_full"
FASTQ_out="FASTQ_repeats"

# Loop through unique genotypes in the CSV file
awk -F',' 'NR>1 {print $5}' "$csv_file" | sort | uniq | while read genotype; do
    if [ -z "$genotype" ]; then
        continue # Skip empty genotype values
    fi

    # Get all targetid values for the current genotype
    target_ids=$(awk -F',' -v g="$genotype" 'NR>1 && $5 == g {print $1}' "$csv_file")

    # Create an array of FASTQ files corresponding to these target_ids
    fastq_files=()  # Initialize an empty array
    for target_id in $target_ids; do
        file_path="${FASTQ_path}/${target_id}.FASTQ.gz"
        if [ -f "$file_path" ]; then
            fastq_files+=("$file_path")  # Add the file path to the array
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

QC, clade split and Ipyrad run

```py
trim_galore --cores 8 -q 20 --length 30 --phred33 -o fastq_trimmed/ FASTQ_PATH/*.FASTQ.gz
```

Filter the data in a few ways to get close to the original dataset.

```py
vcftools --vcf MADR_reference_e1_replicates.vcf --mac 1 --max-missing 0.5 --recode --out MADR_reference_e1_replicates_loose_filter
vcftools --vcf MADR_reference_e1_replicates.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_e1_replicates_strict_filter

```

Take out underperforming individuals

```py
$ /home/pbongaerts/Github/radseq/vcf_missing_data.py MADR_reference_e1_replicates.vcf > MADR_reference_e1_replicates_qc.txt
$ sort -gk 5 -o  MADR_reference_e1_replicates_qc.txt MADR_reference_e1_replicates_qc.txt
$ awk '
BEGIN {
    less_10 = 0;
    less_30 = 0;
    less_1000_geno = 0;
    total_percent = 0;
    count = 0;
}
NR > 1 {  # Skip header
    # Check percentage genotyped
    if ($5 < 10) less_10++;
    if ($5 < 30) less_30++;
    # Check GENO column
    if ($3 < 1000) less_1000_geno++;
    # Add % genotyped to total and increment sample count
    total_percent += $5;
    count++;
}
END {
    print "Samples with <10% genotyped: " less_10;
    print "Samples with <30% genotyped: " less_30;
    print "Samples with <1000 SNPs (GENO): " less_1000_geno;
    print "Average % genotyped: " total_percent / count;
}
'  MADR_reference_e1_replicates_qc.txt
#Samples with <10% genotyped: 22 
#Samples with <30% genotyped: 291 
#Samples with <1000 SNPs (GENO): 4
#Average % genotyped: 32.3641

$ awk '$3 < 1000 {print $1}' MADR_reference_e1_replicates_qc.txt > indv_to_remove_reference_e1_replicates.txt

vcftools --vcf MADR_reference_e1_replicates.vcf --remove indv_to_remove_reference_e1_replicates.txt --recode --out MADR_reference_e1_replicates_temp.vcf

```

Run similarity analisis on both

```py
 python3 ~/scripts/vcf_clone_detect.py -v MADR_reference_e1_replicates.vcf -p MADR_reference_popfile_e1_replicates.txt -o MADR_reference_e1_replicates_genetic_sim.txt > MADR_reference_e1_replicates_clonedet_out.txt

 python3 ~/scripts/vcf_clone_detect.py -v MADR_reference_e1_replicates_strict.vcf -p MADR_reference_popfile_e1_replicates.txt -o MADR_reference_e1_replicates_genetic_sim_strict.txt > MADR_reference_e1_replicates_clonedet_strict_out.txt

 python3 ~/scripts/vcf_clone_detect.py -v MADR_reference_e1_replicates_loose.vcf -p MADR_reference_popfile_e1_replicates.txt -o MADR_reference_e1_replicates_genetic_sim_loose.txt > MADR_reference_e1_replicates_clonedet_loose_out.txt
```

Compare replicates

```py
$ awk -F',' '                                                        
    $1 ~ /_x[0-9]/ && $2 == substr($1, 1, length($1)-3) ||
    $2 ~ /_x[0-9]/ && $1 == substr($2, 1, length($2)-3) {
        print $0
    }
' MADR_reference_e1_replicates_genetic_sim_strict.txt > MADR_reference_e1_replicates_strict_genetic_sim_double_sequenced_pairs.txt

$ tail -n 15 MADR_reference_e1_replicates_strict_genetic_sim_double_sequenced_pairs.txt
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
`#DH0258_MMIR_SNA,DH0258_MMIR_SNA_x2,23252,16331,12020,11536.000000,95.970000,SNA-SNA_x2`  
If we use this as our threshold:

![][image22]

![][image23]

### E1b\_ANGSD IBS clonality assessment

```py

/ccg/bin/angsd -b bam_list.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 475 -setMinDepth 8 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_ANGSD_outputs
```

        \-\> Output filenames:  
                \-\>"MADR\_reference\_e1\_ANGSD\_outputs.arg"  
                \-\>"MADR\_reference\_e1\_ANGSD\_outputs.ibs.gz"  
                \-\>"MADR\_reference\_e1\_ANGSD\_outputs.ibsMat"  
        \-\> Thu Mar 13 15:01:00 2025  
        \-\> Arguments and parameters for all analysis are located in .arg file  
        \-\> Total number of sites analyzed: 43341802  
        \-\> Number of sites retained after filtering: 2164473  
        \[ALL done\] cpu-time used \=  6760.01 sec  
        \[ALL done\] walltime used \=  6843.00 sec  
Create bam order sample name list:

```py
awk -F'/' '{split($NF, a, ".FASTQ"); print a[1]}' bam_list.txt > sample_names_bam_order.txt
```

Initial IBS results:  
![][image24]

```py
/ccg/bin/angsd -b bam_list.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 543 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_ANGSD_outputs

 -> Output filenames:
                ->"MADR_reference_e1_ANGSD_outputs.arg"
                ->"MADR_reference_e1_ANGSD_outputs.ibs.gz"
                ->"MADR_reference_e1_ANGSD_outputs.ibsMat"
        -> Thu Mar 13 20:56:28 2025
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 43341802
        -> Number of sites retained after filtering: 1561519
        [ALL done] cpu-time used =  6030.48 sec
        [ALL done] walltime used =  6121.00 sec

```

Splitting the clades:

```py
BAM_LIST="bam_list.txt"          # File containing BAM file paths
REMOVE_LIST="$HOME/Projects/MADR/D1c_genetic_structure/reference/strict/MADR_reference_d1_minor_clade.txt"    # File containing sample names to remove
OUTPUT_FILE="bam_list_major_clade.txt"  # Output file
# Create a grep pattern from REMOVE_LIST
PATTERN=$(awk '{print $1}' "$REMOVE_LIST" | paste -sd "|" -)
# Filter BAM list, keeping only lines that do NOT match the pattern
grep -Ev "$PATTERN" "$BAM_LIST" > "$OUTPUT_FILE"
echo "Filtered BAM list saved to $OUTPUT_FILE"

#repeat for minor clade
```

```py
/ccg/bin/angsd -b bam_list_major_clade.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 494 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_major_ANGSD_outputs
     -> Output filenames:
                ->"MADR_reference_e1_major_ANGSD_outputs.arg"
                ->"MADR_reference_e1_major_ANGSD_outputs.ibs.gz"
                ->"MADR_reference_e1_major_ANGSD_outputs.ibsMat"
-> Thu Mar 13 21:01:07 2025
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 40679307
        -> Number of sites retained after filtering: 1619345
        [ALL done] cpu-time used =  5531.98 sec
        [ALL done] walltime used =  5592.00 sec


/ccg/bin/angsd -b bam_list_minor_clade.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 54 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_minor_ANGSD_outputs
        -> Output filenames:
                ->"MADR_reference_e1_minor_ANGSD_outputs.arg"
                ->"MADR_reference_e1_minor_ANGSD_outputs.ibs.gz"
                ->"MADR_reference_e1_minor_ANGSD_outputs.ibsMat"
        -> Thu Mar 13 19:32:58 2025
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 17133579
        -> Number of sites retained after filtering: 1079849 
        [ALL done] cpu-time used =  345.06 sec
        [ALL done] walltime used =  345.00 sec
```

```py
#adding options to infer missingness running the same filters but not calculating ibs to save on time
/ccg/bin/angsd -b bam_list_major_clade.txt -GL 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 2 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 494 -setMinDepth 20 -P 4 -out geno_text_output_major
gzip -d geno_text_output_major.geno.gz


/ccg/bin/angsd -b bam_list_minor_clade.txt -GL 1 -doMajorMinor 1 -doMaf 1 -doPost 1 -doGeno 2 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 54 -setMinDepth 20 -P 4 -out geno_text_output_minor
gzip -d geno_text_output_minor.geno.gz


```

Devided into clonal groups using Figure\_S1\_clonality.R

Initial IBS results minor clade: setting the threshold according to the split in node heights

![][image25]

![][image26]

Mapping the IBS clonal groups to the NJ tree (keeping colors consistent)  
![][image27]

![][image28]

![][image29]![][image30]  
Using replicates to set a threshold:

```py
/ccg/bin/angsd -b bam_list_replicates.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 543 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_replicates_ANGSD_outputs
-> Output filenames:
                ->"MADR_reference_e1_replicates_ANGSD_outputs.arg"
                ->"MADR_reference_e1_replicates_ANGSD_outputs.ibs.gz"
                ->"MADR_reference_e1_replicates_ANGSD_outputs.ibsMat"
        -> Sun Mar 16 16:29:51 2025
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 43342207
        -> Number of sites retained after filtering: 1504952
        [ALL done] cpu-time used =  7326.38 sec
        [ALL done] walltime used =  7372.00 sec
```

\#Setting threshold to 0.016 based on the lowest scoring replicate cluster ibs. 

![][image31]

Find the threshold for each individually:

```py
/ccg/bin/angsd -b bam_list_replicates_minor_clade.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 66 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_replicates_minor_ANGSD_outputs
 -> Output filenames:
                ->"MADR_reference_e1_replicates_minor_ANGSD_outputs.arg"
                ->"MADR_reference_e1_replicates_minor_ANGSD_outputs.ibs.gz"
                ->"MADR_reference_e1_replicates_minor_ANGSD_outputs.ibsMat"
        -> Mon Mar 17 14:05:00 2025
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 17129578
        -> Number of sites retained after filtering: 1270478
        [ALL done] cpu-time used =  522.16 sec
        [ALL done] walltime used =  523.00 sec

/ccg/bin/angsd -b bam_list_replicates_major_clade.txt -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -minInd 563 -setMinDepth 20 -doCounts 1 -doMajorMinor 1 -doIBS 1 -makeMatrix 1 -P 1 -out MADR_reference_e1_replicates_major_ANGSD_outputs
-> Output filenames:
                ->"MADR_reference_e1_replicates_major_ANGSD_outputs.arg"
                ->"MADR_reference_e1_replicates_major_ANGSD_outputs.ibs.gz"
                ->"MADR_reference_e1_replicates_major_ANGSD_outputs.ibsMat"
        -> Mon Mar 17 15:44:52 2025
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 40677711
        -> Number of sites retained after filtering: 1621212 
        [ALL done] cpu-time used =  6412.29 sec
        [ALL done] walltime used =  6420.00 sec
```

USE missing % from VCF to select clones to keep:

```py
samples_to_keep <- clone_group_df %>%
  # Keep all samples with clone_group == 0
  filter(clone_group == 0) %>%
  # Combine with the samples from other groups with the lowest F_miss
  bind_rows(
    clone_group_df %>%
      filter(clone_group != 0) %>%
      group_by(clone_group) %>%
      slice_min(F_miss, n = 1, with_ties = FALSE) %>% # Select the lowest F_miss without ties
      ungroup()
  )
#177 samples left for Major clade
#28 samples left for Minor clade
```

# F1\_Clade\_comparison

Create new datasets without clonal samples:

```py
vcftools --vcf ~/Projects/MADR/C3_checking_for_other_contamination/reference/MADR_reference_C3.vcf --keep ~/Projects/MADR/E1_clone_detection/MADR_reference_e1_major_samples_to_keep.txt  --recode --out MADR_reference_f1_major_full
VCFtools - 0.1.16
Keeping individuals in 'keep' list
After filtering, kept 177 out of 681 Individuals
Outputting VCF file...
After filtering, kept 418633 out of a possible 418633 Sites
Run Time = 196.00 seconds

vcftools --vcf ~/Projects/MADR/C3_checking_for_other_contamination/reference/MADR_reference_C3.vcf --keep ~/Projects/MADR/E1_clone_detection/MADR_reference_e1_minor_samples_to_keep.txt  --recode --out MADR_reference_f1_minor_full
Keeping individuals in 'keep' list
After filtering, kept 28 out of 681 Individuals
Outputting VCF file...
After filtering, kept 418633 out of a possible 418633 Sites

vcftools --vcf ~/Projects/MADR/C3_checking_for_other_contamination/reference/MADR_reference_C3.vcf --keep ~/Projects/MADR/E1_clone_detection/MADR_reference_e1_all_samples_to_keep.txt  --recode --out MADR_reference_f1_full
Keeping individuals in 'keep' list
After filtering, kept 205 out of 681 Individuals
Outputting VCF file...
After filtering, kept 418633 out of a possible 418633 Sites
Run Time = 205.00 seconds
```

```py
#filters
vcftools --vcf MADR_reference_f1_major_full.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_f1_major
After filtering, kept 177 out of 177 Individuals
Outputting VCF file...
After filtering, kept 27,201 out of a possible 418,633 Sites

vcftools --vcf MADR_reference_f1_minor_full.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_f1_minor
After filtering, kept 28 out of 28 Individuals
Outputting VCF file...
After filtering, kept 13,498 out of a possible 418,633 Sites


vcftools --vcf MADR_reference_f1_full.vcf --mac 1 --max-missing 0.8 --minDP 8 --recode --out MADR_reference_f1_full
After filtering, kept 205 out of 205 Individuals
Outputting VCF file...
After filtering, kept 31,713 out of a possible 418,633 Sites

vcftools --vcf MADR_reference_f1_full_allele_analysis.vcf --mac 1 --max-missing 0.05 --minDP 8 --recode --out MADR_reference_f1_full_loose_filter
After filtering, kept 205 out of 205 Individuals
Outputting VCF file...
After filtering, kept 169,029 out of a possible 418,633 Sites
 
```

QC

```py
/home/pbongaerts/Github/radseq/vcf_missing_data.py ~/Projects/MADR/F1_populations_statistics/MADR_reference_f1_major.vcf > MADR_reference_f1_major_qc.txt

/home/pbongaerts/Github/radseq/vcf_missing_data.py ~/Projects/MADR/F1_populations_statistics/MADR_reference_f1_minor.vcf > MADR_reference_f1_minor_qc.txt

/home/pbongaerts/Github/radseq/vcf_missing_data.py ~/Projects/MADR/F1_populations_statistics/MADR_reference_f1_all.vcf > MADR_reference_f1_all_qc.txt

$ sort -gk 5 -o MADR_reference_f1_major_qc.txt MADR_reference_f1_major_qc.txt
$ awk '
BEGIN {
    less_10 = 0;
    less_30 = 0;
    less_1000_geno = 0;
    total_percent = 0;
    count = 0;
}
NR > 1 {  # Skip header
    # Check percentage genotyped
    if ($5 < 10) less_10++;
    if ($5 < 30) less_30++;
    # Check GENO column
    if ($3 < 1000) less_1000_geno++;
    # Add % genotyped to total and increment sample count
    total_percent += $5;
    count++;
}
END {
    print "Samples with <10% genotyped: " less_10;
    print "Samples with <30% genotyped: " less_30;
    print "Samples with <1000 SNPs (GENO): " less_1000_geno;
    print "Average % genotyped: " total_percent / count;
}
'  MADR_reference_f1_major_qc.txt
Samples with <10% genotyped: 0
Samples with <30% genotyped: 2
Samples with <1000 SNPs (GENO): 0
Average % genotyped: 90.1801

#repeat:' MADR_reference_f1_minor_qc.txt
Samples with <10% genotyped: 1
Samples with <30% genotyped: 1
Samples with <1000 SNPs (GENO): 1
Average % genotyped: 87.1186
#Remove 1 indv

' MADR_reference_f1_all_qc.txt
Samples with <10% genotyped: 1
Samples with <30% genotyped: 3
Samples with <1000 SNPs (GENO): 1
Average % genotyped: 89.4151

individual to remove: DH0293_MMIR_SNA

```

Building NJ trees

```py
$ python3 ~/scripts/popfile_from_vcf_MODdennis.py MADR_reference_f1_major.vcf 13 > MADR_reference_popfile_f1_major.txt

$ ~/radseq/vcf_gdmatrix.py MADR_reference_f1_major.vcf MADR_reference_popfile_f1_major.txt > MADR_reference_f1_major_gd_mat.txt

$ ~/radseq/gdmatrix2tree.py MADR_reference_f1_major_gd_mat.txt MADR_reference_f1_major.tre

#repeat for minor and all
```

### F1a\_Private and fixed alleles

Hybrids not separated

```py
# For Clade 1
vcftools --vcf MADR_reference_f1_full_loose_filter.recode.vcf --keep ../E1_clone_detection/MADR_reference_e1_major_samples_to_keep.txt --freq --out allele_frequencies/major_clade_allele_freq

# For Clade 2
vcftools --vcf MADR_reference_f1_full_loose_filter.recode.vcf --keep ../E1_clone_detection/MADR_reference_e1_minor_samples_to_keep.txt --freq --out allele_frequencies/minor_clade_allele_freq

```

Finding private and fixed alleles using {compare\_allele\_frequencies.py}

```py
~/scripts/compare_allele_frequencies.py freq_list.txt 0.96 0.01 MADR_reference_f1_loose
#freq_list.txt = txt file containing .frq file path + num_indv
#0.96 % of individuals genotyped for that site
#0.01 leniancy for genotyping error  

wc -l major/minor_clade_allele_freq.frq 
169030

grep "major" MADR_reference_f1_loose_unique_fixed_alleles.txt | wc -l
2533 (1.5% of total sites)
grep "minor" MADR_reference_f1_loose_unique_fixed_alleles.txt | wc -l
866 (0.5% of total sites)
wc -l MADR_reference_f1_loose_fixed_alleles.txt 
3810 (- 3399 = 411 shared fixed alleles) (2.25% of total sites)
grep "major" MADR_reference_f1_loose_private_alleles.txt | wc -l
49047 (29.0 % of total sites)
grep "minor" MADR_reference_f1_loose_private_alleles.txt | wc -l
31014 (18.3 % of total sites)
wc -l MADR_reference_f1_loose_private_alleles.txt 
78855 (46.6% of total sites)
(49047+31014)-78855=1206 sites with alternating private alleles
grep "minor" MADR_reference_f1_loose_private_sites.txt | wc -l
22 (0.013% of total sites)
grep "major" MADR_reference_f1_loose_private_sites.txt | wc -l
8970 (5.30% of total sites)
```

Private alleles with Denovo assembly:

```py
vcftools --vcf ~/Projects/MADR/C3_checking_for_other_contamination/denovo/MADR_denovo_c3.vcf --keep ~/Projects/MADR/E1_clone_detection/MADR_reference_e1_all_samples_to_keep.txt  --recode --out MADR_denovo_f1_full
Keeping individuals in 'keep' list
After filtering, kept 204 out of 681 Individuals
Outputting VCF file...
After filtering, kept 682854 out of a possible 682854 Sites

vcftools --vcf MADR_denovo_f1_full.vcf --mac 1 --max-missing 0.05 --minDP 8 --recode --out MADR_denovo_f1_full_loose_filter
After filtering, kept 204 out of 204 Individuals
Outputting VCF file...
After filtering, kept 204,408 out of a possible 682,854 Sites
Run Time = 99.00 seconds


# For Clade 1
vcftools --vcf MADR_denovo_f1_full_loose_filter.recode.vcf --keep ../../../E1_clone_detection/MADR_reference_e1_major_samples_to_keep.txt --freq --out major_clade_allele_freq_denovo


# For Clade 2
vcftools --vcf MADR_denovo_f1_full_loose_filter.recode.vcf --keep ../../../E1_clone_detection/MADR_reference_e1_minor_samples_to_keep.txt --freq --out minor_clade_allele_freq_denovo

~/scripts/compare_allele_frequencies.py freq_list.txt 0.96 0.01 MADR_denovo_f1_loose
#freq_list.txt = txt file containing .frq file path + num_indv
#0.96 % of individuals genotyped for that site
#0.01 leniancy for genotyping error  

wc -l major/minor_clade_allele_freq.frq 
204409

grep "major" MADR_denovo_f1_loose_unique_fixed_alleles.txt | wc -l
1457 (0.7% of total sites)

grep "minor" MADR_denovo_f1_loose_unique_fixed_alleles.txt | wc -l
714 (0.3% of total sites)
wc -l MADR_denovo_f1_loose_fixed_alleles.txt 
2381 (- 2171 = 210 shared fixed alleles) (0.1% of total sites)
grep "major" MADR_denovo_f1_loose_private_alleles.txt | wc -l
77602 (37.96 % of total sites)
grep "minor" MADR_denovo_f1_loose_private_alleles.txt | wc -l
34223 (16.7 % of total sites)
wc -l MADR_denovo_f1_loose_private_alleles.txt 
109680 (53.6% of total sites)
(77602+34223)-109680=2145 sites with alternating private alleles
grep "minor" MADR_denovo_f1_loose_private_sites.txt | wc -l
37 (0.018% of total sites)
grep "major" MADR_denovo_f1_loose_private_sites.txt | wc -l
13959 (6.82% of total sites)
```

Hybrids separated

```
# For Clade 1
vcftools --vcf ../hybrids_not_separated/MADR_reference_f1_full_loose_filter.recode.vcf --keep DAPC_major_samples_to_keep.txt --freq --out major_clade_allele_freq
kept 176 samples
# For Clade 2
vcftools --vcf ../hybrids_not_separated/MADR_reference_f1_full_loose_filter.recode.vcf --keep DAPC_minor_samples_to_keep.txt --freq --out minor_clade_allele_freq
kept 23 samples
# For hybrids
vcftools --vcf ../hybrids_not_separated/MADR_reference_f1_full_loose_filter.recode.vcf --keep DAPC_hybrid_samples_to_keep.txt --freq --out hybrids_allele_freq
kept 5 samples
```

Finding private and fixed alleles using {compare\_allele\_frequencies.py}

```py
~/scripts/compare_allele_frequencies.py freq_list_maj_min.txt 0.96 0.01 MADR_reference_f1_min_maj 

#freq_list.txt = txt file containing .frq file path + num_indv
#0.96 % of individuals genotyped for that site
#0.01 leniancy for genotyping error  

wc -l major/minor_clade_allele_freq.frq 
169030

grep "major" MADR_reference_f1_min_maj_unique_fixed_alleles.txt | wc -l
1913 (xx% of total sites)
grep "minor" MADR_reference_f1_min_maj_unique_fixed_alleles.txt | wc -l
2871 (xx% of total sites)
wc -l MADR_reference_f1_min_maj_fixed_alleles.txt 
5839 (- 4784 = 1055 shared fixed alleles) (xx% of total sites)
wc -l MADR_reference_f1_min_maj_unique_fixed_alleles.txt 
4785 (- 4784 -1) = 0 alternatively fixed alleles) (xx% of total sites)
grep "major" MADR_reference_f1_min_maj_private_alleles.txt | wc -l
49619 (xx % of total sites)
grep "minor" MADR_reference_f1_min_maj_private_alleles.txt | wc -l
29761 (xx % of total sites)
wc -l MADR_reference_f1_min_maj_private_alleles.txt 
78136 (xx% of total sites)
(49047+31014)-78855=1206 sites with alternating private alleles
grep "minor" MADR_reference_f1_min_maj_private_sites.txt | wc -l
31 (xx% of total sites)
grep "major" MADR_reference_f1_min_maj_private_sites.txt | wc -l
9343 (xx% of total sites)
```

```py
#loose filter for fixed alleles:
 ~/scripts/compare_allele_frequencies.py freq_list_maj_min.txt 0.6 0.05 MADR_reference_f1_min_maj_loose 

grep "major" MADR_reference_f1_min_maj_loose_unique_fixed_alleles.txt | wc -l
15275 (xx% of total sites)
grep "minor" MADR_reference_f1_min_maj_loose_unique_fixed_alleles.txt | wc -l
9321 (xx% of total sites)
wc -l MADR_reference_f1_min_maj_loose_unique_fixed_alleles.txt 
24593 (- 24596) = -3 alternatively fixed alleles) (xx% of total sites)
grep "major" MADR_reference_f1_min_maj_loose_private_alleles.txt | wc -l
29108 (xx % of total sites)
grep "minor" MADR_reference_f1_min_maj_loose_private_alleles.txt | wc -l
30193 (xx % of total sites)
wc -l MADR_reference_f1_min_maj_loose_private_alleles.txt 
57857 (xx% of total sites)
(49047+31014)-78855=1206 sites with alternating private alleles
grep "minor" MADR_reference_f1_min_maj_loose_private_sites.txt | wc -l
31 (xx% of total sites)
grep "major" MADR_reference_f1_min_maj_loose_private_sites.txt | wc -l
9343 (xx% of total sites)

```

Finding most divergent sites:

```py
 ~/scripts/compare_allele_frequencies.py freq_list_maj_min.txt 0.96 0.01 MADR_reference_f1_min_maj 20000 

~/Projects/MADR/F1_populations_statistics/F1a_Private_and_fixed_alleles/reference/hybrids_separated$ awk '$4 > 0.5' MADR_reference_f1_min_maj_most_divergent_loci.txt | wc -l
4962
awk '$4 > 0.6' MADR_reference_f1_min_maj_most_divergent_loci.txt | wc -l
3730
awk '$4 > 0.6' MADR_reference_f1_min_maj_most_divergent_loci.txt > MADR_reference_f1_most_divergent_loci.txt
```

###  F1b\_Hybrid\_analysis

Adding pop to popfile

```py
awk 'NR==FNR {samples[$1]="major"; next} {print $0, ($1 in samples) ? samples[$1] : "minor"}' ~/Projects/MADR/E1_clone_detection/MADR_reference_e1_major_samples_to_keep.txt MADR_reference_popfile_f1_all.txt 
 > MADR_reference_popfile_f1_clades.txt
```

Running popstats\_MADR.R  
\#PCA with NA \= ‘mean’  
![][image32]%%  
![][image33]

![][image34]

DAPC  
![][image35]

![][image36]

Following the DAPC posterior assignment:

![][image37]

![][image38]  
![][image39]  
Re running structure:

```py
$ scp /home/deepcat/STRUCTURE_params/extraparams .
$ scp /home/deepcat/STRUCTURE_params/mainparams .

$ vcftools --vcf ../MADR_reference_f1_all.vcf
#After filtering, kept 204 out of 204 Individuals
#After filtering, kept 31713 out of a possible 31713 Sites


$ nano mainparams
KEY PARAMETERS FOR THE PROGRAM structure.  YOU WILL NEED TO SET THESE
IN ORDER TO RUN THE PROGRAM.  VARIOUS OPTIONS CAN BE ADJUSTED IN THE
FILE extraparams.
"(int)" means that this takes an integer value.
"(B)"   means that this variable is Boolean
        (ie insert 1 for True, and 0 for False)
"(str)" means that this is a string (but not enclosed in quotes!)
Basic Program Parameters
#define MAXPOPS    4      // (int) number of populations assumed
#define BURNIN    100000   // (int) length of burnin period
#define NUMREPS   50000   // (int) number of MCMC reps after burnin
Input/Output files
#define INFILE   infile   // (str) name of input data file
#define OUTFILE  outfile  //(str) name of output data file
Data file format
#define NUMINDS    204    // (int) number of diploid individuals in data file
#define NUMLOCI    31713    // (int) number of loci in data file
#define PLOIDY       2    // (int) ploidy of data
#define MISSING     -9    // (int) value given to missing genotype data
#define ONEROWPERIND 0    // (B) store data for individuals in a single line
#define LABEL     1     // (B) Input file contains individual labels
#define POPDATA   1     // (B) Input file contains a population identifier
#define POPFLAG   0     // (B) Input file contains a flag which says 
                              whether to use popinfo when USEPOPINFO==1
#define LOCDATA   0     // (B) Input file contains a location identifier
#define PHENOTYPE 0     // (B) Input file contains phenotype information
#define EXTRACOLS 0     // (int) Number of additional columns of data 
                             before the genotype data start.
#define MARKERNAMES      0  // (B) data file contains row of marker names
#define RECESSIVEALLELES 0  // (B) data file contains dominant markers (eg AFLPs)
                            // and a row to indicate which alleles are recessive
#define MAPDISTANCES     0  // (B) data file contains row of map distances 
                            // between loci
Advanced data file options
#define PHASED           0 // (B) Data are in correct phase (relevant for linkage model only)
#define PHASEINFO        0 // (B) the data for each individual contains a line
                                  indicating phase (linkage model)
#define MARKOVPHASE      0 // (B) the phase info follows a Markov model.
#define NOTAMBIGUOUS  -999 // (int) for use in some analyses of polyploid data

#set up interactive shell with 8 cpu's and 2gb ram

$ ~/scripts/structure_mp_ref_update.py ../MADR_reference_f1_all.vcf ../MADR_reference_popfile_f1_all.txt 4 10 32
Initialise indivs and pops for ../MADR_reference_f1_all.vcf...
Subsample SNPs (one random SNP per locus)... [5492 SNPs/loci]
Outputting 10 STRUCTURE files...10 reps DONE
Executing 32 parallel STRUCTURE runs for K = 2 ...10 reps DONE
Executing 32 parallel STRUCTURE runs for K = 3 ...10 reps DONE
Executing 32 parallel STRUCTURE runs for K = 4 ...10 reps DONE
Running CLUMPP on replicates for K = 2 ...
Running CLUMPP on replicates for K = 3 ...
Running CLUMPP on replicates for K = 4 ...
K = 2: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      ../MADR_reference_f1_all.vcf
K = 3: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      ../MADR_reference_f1_all.vcf
K = 4: MedMeaK 1.0 MaxMeaK 1 MedMedK 1.0 MaxMedK 1      ../MADR_reference_f1_all.vcf
```

![][image40]

Get max likelihood for k struct:

```py
output_file="Evanos_K_values.txt"
for file in *.log; do
    if grep -q "Estimated Ln Prob of Data" "$file"; then
        # Extract the value
        lnpd=$(grep "Estimated Ln Prob of Data" "$file" | awk -F'= ' '{print $2}')
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
    ' "$output_file"
done
Average for _2.log: -410874
Average for _3.log: -410190
Average for _4.log: -421119

```

Test to see true hybrid assignment

Running PCA/DAPC with simulated data using adegenet hybridize:

To do:  
![][image41]

MAX PCAs= 106

MAX PCAs \=2 (more similar to results without added hybrids)

“Suspected” hybrids in green PCA

![][image42]

Running newhybrids with simulated and actual data:

```py
Program New Hybrids completed at Thu Mar 27 16:49:44 2025
Data = MADR_reference_f1_nhyb.txt, PiPrior = JEFFREYS, ThetaPrior = JEFFREYS, Seeds = 27 13, NumBurnIn = 50000, NumRepsAfterBurnIn = 100000
SeedsAtEnding = 1940924971 560041219

Program New Hybrids completed at Fri Mar 28 09:05:38 2025
 
Data = SIM_hybrids_nhyb.txt, PiPrior = JEFFREYS, ThetaPrior = JEFFREYS, Seeds = 27 13, NumBurnIn = 50000, NumRepsAfterBurnIn = 100000
SeedsAtEnding = 594100171 307743162
```

![][image43]

![][image44]

![][image45]  
Nj tree

It still seems like there is a fair bit of noise which could be due to the selection of samples to hybridize. 

![][image46]  
HC clust  
![][image47]  
NJ tree  
![][image48]

Full SIM data with only pure hybrids

![][image49]  
![][image50]

![][image51]

![][image52]

### ![][image53]

### 

### 

### 

### 

### 

### 

### F1c\_genetic\_distance
Genetic distances and population statistics were calculated using the Popstats_full_data.R script.

# G1\_reefscape\_assessment
Reefscape genomics assessment was completed using the monostand_metrics.R, the monostand_diversity.R, and Figure_2.R scripts.

