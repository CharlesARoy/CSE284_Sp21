# CSE284 Spring 2021 — Term Project

## Investigating Structural Variant Detection in Simple Tandem Repeat Regions Using HG002 (NA24385) Nanopore Data

### Overview

The goal of this project is to test a variety of tools and tool parameters to optimize structural variant (SV) detection in simple tandem repeat (STR) regions of the human genome using Nanopore (ONT) data. This problem stems from the fact that even state-of-the-art aligners, such as Minimap2, often introduce gaps and insertions somewhat randomly in repetitive regions, resulting in patterns such as those seen in the following image. Poorly localized indels in reads can reduce the accuracy of SV calls and generally make genotyping more challenging. In trio data this issue can also result in inflated estimates of _de novo_ mutation (e.g. when SVs are called in children but are missed in parents). 


![alt text](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/misaligned_deletions.jpg?raw=true)


### Methods

#### Sources of Data

The human reference genome used in this project was obtained from the BROAD GATK Resource Bundle which is hosted [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0) (Homo_sapiens_assembly38.fasta).

The data used in this project was ultra-long GridION Nanopore data from HG002 (NA24385) that was sequenced by [Genome in a Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle). The raw fast5 files were obtained from [this Amazon bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=) but are no longer hosted there. Fastq files base called with an older version of Guppy than used in this project can be obtained from the [Human Pangenomics HG002 data freeze](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0) (data stored [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/nanopore/)).

The truth set for HG002 (NA was obtained with the following command and was then lifted over to hg38:
```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
```

A bed file of Simple Tandem Repeats was obtained from the UCSC Table Browser by seleting the following options:
* Group: Repeats
* Track: Simple Repeats
* Region: Genome
* Output Format: BED


#### Pipeline Commands

A broad overview of the pipeline is displayed below.

The raw fast5 data from HG002 was basecalled with [Guppy version 4.0.11](https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_4.0.11_linux64.tar.gz) using the following command:
```bash
guppy_basecaller --recursive --input_path /path/to/fast5/files --save_path /output/path --config dna_r9.4.1_450bps_hac.cfg --verbose_logs --device "cuda:all"
```

After basecalling the data, it was mapped to the HG38 human reference genome with Minimap2 and then sorted, indexed, and downsampled with [Samtools](http://www.htslib.org/). The data was downsampled from 50X coverage to 10X coverage to speed up testing and to better approximate sample coverage that might be more realistic in clinical testing. The downsampled data was then converted back to fastq format so that the reads could be aligned with other aligners. 
```bash
# Align
minimap2 -x map-ont -t 10 -a --MD -Y reference_file.fasta unmapped_reads.fastq -o aligned_reads.bam

# Sort and Index
samtools sort aligned_reads.bam -o aligned_reads_sorted.bam
samtools index aligned_reads_sorted.bam

# Downsample data
samtools view -s 0.2 --threads 2 -b aligned_reads_sorted.bam > aligned_reads_sorted_downsampled.bam

# Extract chromosome 21
samtools view aligned_reads_sorted_downsampled.bam "Chr21" > aligned_reads_sorted_downsampled_chr21.bam

# Convert downsampled bam back to fastq
samtools fastq aligned_reads_sorted_downsampled_chr21.bam > unmapped_reads_downsampled_chr21.fastq
```

Following this, the data was aligned with multiple tools. [Minimap2](https://github.com/lh3/minimap2) was used because it is currently the most popular long read aligner and is recommended by both Nanopore and PacBio scientists. [Minimap2 v2.20 was released during this project](https://github.com/lh3/minimap2/blob/master/NEWS.md) and is supposed to produce better alignment in highly repetitive regions, so it was tested in addition to version 2.18. [Winnowmap](https://github.com/marbl/Winnowmap) was used because it is optimized for mapping ONT and PacBio reads to repetitive reference sequences. Finally, because Winnowmap is an extension of Minimap2, [NGMLR](https://github.com/philres/ngmlr) was used as an alternative.

The alignments were performed with the following commands and were then sorted and indexed with Samtools and converted from sam to bam (not shown):
```bash
# Minimap2 v2.18 and 2.20
minimap2 -x map-ont -t 10 -a --MD -Y reference_file.fasta unmapped_reads.fastq -o aligned_reads.bam

# NGMLR
ngmlr -x ont -t 10 -r reference_file.fasta -q unmapped_reads.fastq -o aligned_reads.sam

# Winnowmap
winnowmap -x map-ont --MD -Y -a reference_file.fasta unmapped_reads.fastq -o aligned_reads.sam
```

Following alignment, for each aligner, structural variants were called with [NanoVar](https://github.com/benoukraflab/NanoVar), [Sniffles](https://github.com/fritzsedlazeck/Sniffles), and [SVIM](https://github.com/eldariont/svim) using the following commands:
```bash
# NanoVar
nanovar -x ont -f hg38 -c 1 -t 2 aligned_reads_sorted.bam Homo_sapiens_assembly38.fasta ./nanovar_work_dir

# Sniffles
sniffles -m aligned_reads_sorted.bam -v aligned_reads_sorted.vcf -t 2 -s 1 -n -1 --genotype --cluster

# SVIM
svim alignment --interspersed_duplications_as_insertions --tandem_duplications_as_insertions --max_sv_size 10000000 ./svim_work_dir aligned_reads_sorted.bam Homo_sapiens_assembly38.fasta
```

The VCF files and the truth set were processed similarly using [Snakemake](https://snakemake.readthedocs.io/en/stable/) to ensure consistency. The Snakemake file for the truth set can be seen [here](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/snakefile_preparing_truth_set) and one of the Snakemake files for processing the SV calls can be seen [here](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/snakefile_post-processing_sniffles-svs). The main processing steps included:
* Filtering SVs that were smaller than 50 bp. This was done because SVs are usually defined as 50 bp or larger but some SV callers have defaults that are smaller than that.
* Filtering SVs that are larger than 10 Mb. This was done to avoid including massive inversions that sometimes get called and slow down later processing steps.
* Convert the VCF files to BED format, keeping the Chr, Start, End, SV length, SV type, Read support, and Genotype.
* Retain SVs that intersect by 50% or more with simple tandem repeat regions.
* Split the BED file into a file of deletions and a file of insertions and duplications. The insertions and duplications were kept together because different SV callers can refer to insertions and duplications inconsistently and might e.g. call an insertion a duplication.
* Pad the insertion start/end coordinates so that they span up to 100 bp. This was done because insertions and duplications are usually just given a 1 bp location but this is impractical when intersecting calls with the truth set to identify true positives and false positives. 

After this, the final step was to intersect the SV callsets of each type (deletions vs insertions + duplications) with the calls of the same type in the truth set. To be considered in the truth set, a given call had to have 50% reciprocal overlap with a call in the truth set which is achieved with the following command. Fifty percent reciprocal overlap was chosen rather than say 100% reciprocal overlap to account for the fact that different aligners and SV callers can cause the coordinates of SVs to shift somewhat.

```bash
bedtools intersect -f 0.5 -r -wa -wb -a sv_calls.bed -b truthset.bed
```


![alt text](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/Flowchart.jpeg?raw=true)


### Results

The primary metrics used to assess performance were sensitivity and false discovery rate (FDR). Consistent with expectations, across all combinations of aligner and SV caller, relative to non-repetitive regions, in STR regions the sensitivity was found to be 8% lower and the FDR was found to be 15% higher. Sensitivity and FDR were both slightly better for deletions compared to insertions. It was also found that when calls were required to have two or more supporting reads, sensitivity dropped by 4% and FDR decreased by 35% in non-repetitive regions and 9% in STR regions.

As indicated in **Table 1**, for both deletions and insertions it was found that Winnowmap had the highest sensitivity and NGMLR had the lowest FDR. Among SV callers, it was found that SVIM typically had the highest sensitivity and the highest FDR. For deletions and insertions, the highest sensitivity was achieved with Winnowmap and SVIM. For deletions the lowest FDR was achieved with NGMLR and NanoVar and for insertions it was achieved with NGMLR and Sniffles.

**Table 1.** Sensitivity and False Discovery Rate of Various Aligners and SV Callers
|                            | Minimap2 v2.18 | Minimap2 v2.20 | NGM-LR | Winnowmap |
|----------------------------|:--------------:|:--------------:|:------:|:---------:|
|  NanoVar; DEL; Sensitivity |      0.25      |      0.25      |  0.30  |    0.25   |
| Sniffles; DEL; Sensitivity |      0.33      |      0.34      |  0.30  |    0.36   |
|   SVIM; DEL; Sensitivity   |      0.40      |      0.41      |  0.38  |    0.44   |
|      NanoVar; DEL; FDR     |      0.42      |      0.44      |  0.35  |    0.42   |
|     Sniffles; DEL; FDR     |      0.56      |      0.56      |  0.40  |    0.54   |
|       SVIM; DEL; FDR       |      0.73      |      0.67      |  0.46  |    0.62   |
|  NanoVar; INS; Sensitivity |      0.29      |      0.30      |  0.37  |    0.26   |
| Sniffles; INS; Sensitivity |      0.35      |      0.35      |  0.32  |    0.33   |
|   SVIM; INS; Sensitivity   |      0.52      |      0.52      |  0.45  |    0.56   |
|      NanoVar; INS; FDR     |      0.70      |      0.69      |  0.61  |    0.73   |
|     Sniffles; INS; FDR     |      0.66      |      0.65      |  0.50  |    0.62   |
|       SVIM; INS; FDR       |      0.70      |      0.68      |  0.61  |    0.71   |

Comparisons are restricted to calls with read support of two or greater that intersect with simple tandem repeats.
