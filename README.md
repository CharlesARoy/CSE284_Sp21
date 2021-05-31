# CSE284 Spring 2021 — Term Project

## Investigating Structural Variant Detection in Simple Tandem Repeat Regions Using HG002 (NA24385) Nanopore Data

The goal of this project is to test a variety of tools and tool parameters to optimize structural variant (SV) detection in simple tandem repeat (STR) regions of the human genome using Nanopore (ONT) data. This problem stems from the fact that even state-of-the-art aligners, such as Minimap2, often introduce gaps and insertions somewhat randomly in repetitive regions, resulting in patterns such as those seen in the following image. Poorly localized indels in reads can reduce the accuracy of SV calls and generally make genotyping more challenging. In trio data this issue can also result in inflated estimates of _de novo_ mutation (e.g. when SVs are called in children but are missed in parents). 


![alt text](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/misaligned_deletions.jpg?raw=true)


#### Methods

The data used in this project was Ultra-long GridION data sequenced by Genome in a Bottle (GIAB). The raw fast5 files were obtained from [this Amazon bucket](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=) but are no longer hosted there. Currently, fastq files can be obtained from the human pangenomics Amazon bucket found [here](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/hpp_HG002_NA24385_son_v1/nanopore/).

The truth set for this sample was obtained with the following command and was then lifted over to hg38:
```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
```


The raw fast5 data was basecalled with Guppy version 4.0.11 using the following command: <br>
```bash
guppy_basecaller --recursive --input_path /path/to/fast5/files --save_path /output/path --config dna_r9.4.1_450bps_hac.cfg --verbose_logs --device "cuda:all"
```

After basecalling the data, it was mapped to the HG38 human reference genome with Minimap2 and then sorted, indexed, and downsampled with samtools. The data was downsampled from 50X coverage to 10X coverage to speed up testing and to better approximate sample coverage that might be more realistic in clinical testing. The downsampled data was then converted back to fastq format so that the reads could be aligned with other aligners. 
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



![alt text](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/Flowchart.jpeg?raw=true)


#### Results


|                 | Minimap2 v2.18 | Minimap2 v2.20 | NGM-LR | Winnowmap |
|-----------------|:--------------:|:--------------:|:------:|:---------:|
| Sniffles; Sens. |      0.33      |      0.34      |  0.30  |    0.36   |
|   SVIM; Sens.   |      0.40      |      0.41      |  0.38  |    0.44   |
|   NanoVar; FDR  |      0.42      |      0.44      |  0.35  |    0.42   |
|  Sniffles; FDR  |      0.56      |      0.56      |  0.40  |    0.54   |
|    SVIM; FDR    |      0.73      |      0.67      |  0.46  |    0.62   |
|  NanoVar; Sens. |      0.29      |      0.30      |  0.37  |    0.26   |
| Sniffles; Sens. |      0.35      |      0.35      |  0.32  |    0.33   |
|   SVIM; Sens.   |      0.52      |      0.52      |  0.45  |    0.56   |
|   NanoVar; FDR  |      0.70      |      0.69      |  0.61  |    0.73   |
