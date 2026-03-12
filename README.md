This is a sample MitoSeq data analysis pipeline for a multiplexed run. Conda environments for longread_umi and samtools are required. 
This is an archive of data analysis pipelines used in MitoSeq. Copyright © 2026 The Trustees of Columbia University in the City of New York. All Rights Reserved.

MitoSeq Long-Read UMI Analysis Pipeline

Pipeline for analyzing MitoSeq long-read sequencing data with UMIs, including read filtering, alignment, polishing, and variant calling.

This workflow requires the longread_umi framework:

https://github.com/SorenKarst/longread_umi

The pipeline was used in the associated manuscript and tested on macOS systems.

1. Conda Environment
1.1 Intel Mac Environment (Used in Manuscript)

Environment name: longread_umi_Intel

name: longread_umi_Intel
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.6
  - seqtk=1.3
  - parallel=20191122
  - racon=1.4.10
  - minimap2=2.17
  - medaka=0.11.5
  - gawk=4.1.3
  - cutadapt=2.7
  - filtlong=0.2.0
  - bwa=0.7.17
  - samtools=1.9
  - bcftools=1.9
  - git
platform: osx-64

1.2 Apple Silicon Mac Environment (Alternative)

If the Intel version is unavailable when building the Conda environment, the following Apple Silicon (ARM) environment was tested and showed slightly improved performance.

Environment name: longread_umi

name: longread_umi
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.9
  - seqtk=1.5
  - parallel=20220722
  - racon=1.5.0
  - minimap2=2.24
  - medaka=2.1.1
  - gawk=5.3.1
  - cutadapt=5.1
  - filtlong=0.3.1
  - bwa=0.7.17
  - samtools=1.9
  - bcftools=1.19
  - git
platform: osx-arm64

2. Additional Required Software

The pipeline requires additional tools outside the main Conda environment.

Base Environment
samtools 1.21
Additional Conda Environment

Environment name: samtools_env_Intel

name: samtools_env_Intel
dependencies:
  - samtools=1.2
  - bcftools=1.15.1
  - cutadapt=3.5
  - gawk=5.1.0
  - parallel=20210822
  - minimap2=2.24
platform: osx-64

Additional Software

The pipeline also requires:

USEARCH 11

3. Operating Systems Tested

The pipeline was tested on:

macOS 15.7.1
macOS 10.14.6

4. Quick Start

Before running the pipeline, update the following variables in the scripts:

$PATH_to_basecalled_files
$PATH_to_barcode
$PATH_to_reference
$PATH_to_script_files
$Run_name

These variables must be updated in the following scripts:

part1
part2
part3

5. Running the Pipeline

From a shell terminal, run:

time sh Sample_run.part1.60K.sh

This command will execute:

part1
part2
part3

sequentially.

6. Output

The pipeline runs in three steps.

Part 1

Demultiplexing using cutadapt

UMI sequence recovery using longread_umi

Part 2

UMI clustering

Variant detection

Part 3

Consensus sequence generation and alignment to the mitochondrial reference (chrM)

Output Directories

Each demultiplexed sample is stored in a separate folder named:

$Run_name.101.dup.V01.cut1.21.m150
$Run_name.101.dup.V02.cut1.21.m150
$Run_name.101.dup.V03.cut1.21.m150
Coverage-Based Results

Each sample folder contains results generated with different read coverage thresholds:

Result	Description
≥4 reads	Consensus from clusters with ≥4 reads
≥5 reads (Ave-1)	Consensus from clusters with ≥5 reads
≥6 reads (Ave-2)	Consensus from clusters with ≥6 reads
Average coverage (Ave-cal)	Consensus based on calculated average coverage

Results using Average coverage (Ave-cal) were used in the manuscript.

Key Output Files

Aligned consensus sequences:

$Run_name.101.dup.V02.cut1.21.m150.60K.Ave-cal.con.chrM.sort.bam

Variant detection results:

cluster_94_60K_4up/Ave-cal.60.vcf.count.1.input

Notes

Ensure the correct Conda environment is activated before running the pipeline.

Verify that all required software versions are installed.

Check that all paths in the scripts are correctly set.


Reference
1. Karst SM, et al. High-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. Nat Methods 18, 165–169 (2021)
2. Wei S, et al. Development of a Clinically Applicable High-Resolution Assay for Sperm Mosaicism. J Mol Diagn 27, 525–537 (2025)
