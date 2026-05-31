# Metagenomic Processing Workflow

This document describes the core command-line workflow used for metagenomic read processing, host read depletion, taxonomic profiling, metagenomic assembly, ORF prediction, non-redundant gene catalog construction, MAG recovery, and read-based gene abundance quantification.

The workflow is associated with the manuscript:

**Ventilation-shaped farm environments link infectious coryza dissemination, infraorbital sinus microbiome collapse and mobile resistome accumulation in laying hens**

This document is intended to provide representative commands and computational parameters for reproducibility.

---

# 1. Input files and directory structure

Raw paired-end metagenomic reads were assumed to be named as follows:

```bash
raw/${SAMPLE}_R1.fastq.gz
raw/${SAMPLE}_R2.fastq.gz

Recommended directory structure:

raw/
qc/
sam/
bams/
clean_nohost/
taxonomy/
assembly/
genes/
NR_catalog/
binning/
MAGs/
abundance/
databases/

The variable ${SAMPLE} represents the sample ID.

2. Software and databases

The major tools used in this workflow include:

Trimmomatic
fastp
Bowtie2
samtools
Kraken2
Bracken
DIAMOND BLASTX
MEGAHIT
seqkit
MetaGeneMark
CD-HIT-EST
BWA
SAMtools
MetaBAT2
MaxBin2
CONCOCT
DAS Tool
CheckM
dRep
GTDB-Tk
Kraken2/Bracken-compatible taxonomic databases
Metaxa2
CoverM
FastANI
inStrain
IQ-TREE
iTOL


Required databases or indices include:

databases/human_chicken_composite
databases/MetaGeneMark_meta_3.mod
GTDB-Tk database
database/ResFinder
database/Comprehensive Antibiotic Resistance Database
database/MobileGeneticElementDatabase
database/TnCentral, ISfinder, INTEGRALL



3. Quality Control & Host Decontamination
Raw reads are processed using `fastp` for quality filtering (Phred < 20, length $\ge$ 50 bp). Host-derived contamination (Human GRCh38 and Chicken GRCg7b) is systematically depleted by mapping reads using `Bowtie2` (v2.5.1).

```bash
# 3.1 Quality Trimming
fastp -i raw/${SAMPLE}_R1.fastq.gz -I raw/${SAMPLE}_R2.fastq.gz \
      -o qc/${SAMPLE}_clean_R1.fastq.gz -O qc/${SAMPLE}_clean_R2.fastq.gz \
      -q 20 --length_required 50 -w 8

# 3.2 Host Genome Depletion via Bowtie2 (Against Human + Chicken composite index)
bowtie2 -x databases/human_chicken_composite \
        -1 qc/${SAMPLE}_clean_R1.fastq.gz -2 qc/${SAMPLE}_clean_R2.fastq.gz \
        -S sam/${SAMPLE}_host_mapped.sam --threads 16 --sensitive

# 3.3 Extract paired unmapped reads (Host-depleted data)
samtools view -@ 8 -b -f 12 -F 256 sam/${SAMPLE}_host_mapped.sam -o bams/${SAMPLE}_nohost.bam
samtools sort -@ 8 -n bams/${SAMPLE}_nohost.bam -o bams/${SAMPLE}_nohost_sorted.bam
samtools fastq -@ 8 bams/${SAMPLE}_nohost_sorted.bam \
               -1 clean_nohost/${SAMPLE}_nohost_1.fastq.gz \
               -2 clean_nohost/${SAMPLE}_nohost_2.fastq.gz





4. Taxonomic Classification and Profiling
Bacterial community composition is profiled directly from host-depleted reads using Kraken2 and refined with Bracken for species-level relative abundance.

# 4.1 Kraken2 Taxonomic Assignment
kraken2 --db databases/k2_pluspf_db --threads 16 --gzip-compressed \
        --paired clean_nohost/${SAMPLE}_nohost_1.fastq.gz clean_nohost/${SAMPLE}_nohost_2.fastq.gz \
        --report taxonomy/${SAMPLE}.k2report --output taxonomy/${SAMPLE}.kraken

# 4.2 Bracken Species-Level Estimation
bracken -d databases/k2_pluspf_db -i taxonomy/${SAMPLE}.k2report \
        -o taxonomy/${SAMPLE}.bracken -r 150 -l S -t 10




5. Assembly, ORF Prediction & Gene Catalog Construction
Host-depleted reads are assembled de novo using MEGAHIT. Open Reading Frames (ORFs) are predicted via MetaGeneMark and clustered via CD-HIT-EST to build the Non-Redundant (NR) gene catalog.

# 5.1 De novo Assembly
megahit -1 clean_nohost/${SAMPLE}_nohost_1.fastq.gz -2 clean_nohost/${SAMPLE}_nohost_2.fastq.gz \
        --num-cpu-threads 16 --out-dir assembly/${SAMPLE}_out/ --out-prefix ${SAMPLE}

# 5.2 Retain contigs longer than 300 bp
seqkit seq -m 300 assembly/${SAMPLE}_out/${SAMPLE}.final.contigs.fa > assembly/${SAMPLE}.ge300.fa

# 5.3 MetaGeneMark Gene Prediction
gmhmmp -m databases/MetaGeneMark_meta_3.mod -f G \
       -o genes/${SAMPLE}.gff -d genes/${SAMPLE}.fna -a genes/${SAMPLE}.faa \
       assembly/${SAMPLE}.ge300.fa

# 5.4 Construct Non-Redundant Catalog (95% Identity, 90% Coverage)
cd-hit-est -i all_combined_orfs.fna -o NR_catalog/nr_gene_catalog.fna \
           -c 0.95 -aS 0.90 -M 0 -T 24 -g 1




6. Ensemble Metagenomic Binning and SGB Analysis
Contigs $\ge$ 1,500 bp are binned independently using MetaBAT2, MaxBin2, and CONCOCT. Predictions are integrated via DAS Tool, dereplicated with dRep, and classified with GTDB-Tk.

# 6.1 Generate Contig Coverage Profiles via Bowtie2
seqkit seq -m 1500 assembly/${SAMPLE}.ge300.fa > binning/${SAMPLE}.ge1.5k.fa
bowtie2-build binning/${SAMPLE}.ge1.5k.fa binning/${SAMPLE}.ge1.5k
bowtie2 -x binning/${SAMPLE}.ge1.5k -1 clean_nohost/${SAMPLE}_nohost_1.fastq.gz -2 clean_nohost/${SAMPLE}_nohost_2.fastq.gz | samtools view -@ 8 -bS - | samtools sort -@ 8 - > binning/${SAMPLE}_sorted.bam

# 6.2 Multi-engine binning (MetaBAT2, MaxBin2, CONCOCT)

metabat2 -i contigs.ge1.5k.fa -a depth.txt -o metabat_out/bin
run_MaxBin.pl -contig contigs.ge1.5k.fa -abund maxbin_abund.txt -out maxbin_out/bin
concoct --composition_file contigs.ge1.5k.fa --coverage_file coverage.tsv -b concoct_out/


# 6.3 DAS Tool Ensemble Integration
DAS_Tool -i MetaBAT2.txt,MaxBin2.txt,CONCOCT.txt -l metabat,maxbin,concoct \
         -c binning/${SAMPLE}.ge1.5k.fa -o binning/${SAMPLE}_DASTool/ --write_bins

# 6.4 MAG Quality Check & dRep Dereplication (95% ANI)
checkm lineage_wf binning/${SAMPLE}_DASTool/_bins/ checkm_out/
drep dereplicate drep_out/ -g binning/*_DASTool/_bins/*.fa -sa 0.95 --comp 50 --con 10




7. Advanced Resistome and Mobilome Quantification
ARG profiling utilizes both CARD and ResFinder (unified under CARD nomenclature via customized BLASTN harmonization). MGEs are identified via DIAMOND. Total abundances are normalized to 16S rRNA gene counts parsed by Metaxa2.

# 7.1 Map host-depleted reads back to the NR catalog for quantification
bwa index NR_catalog/nr_gene_catalog.fna
bwa mem -t 16 NR_catalog/nr_gene_catalog.fna clean_nohost/${SAMPLE}_nohost_1.fastq.gz clean_nohost/${SAMPLE}_nohost_2.fastq.gz | samtools view -@ 8 -bS - | samtools sort -@ 8 - > abundance/${SAMPLE}_mapped.bam
samtools idxstats abundance/${SAMPLE}_mapped.bam > abundance/${SAMPLE}_counts.txt

# 7.2 Extract 16S rRNA Reads via Metaxa2 for Abundance Normalization
metaxa2 -i clean_nohost/${SAMPLE}_nohost_1.fastq.gz -o abundance/${SAMPLE}_16s --mode g -g SSU --threads 8


