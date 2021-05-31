# CSE284 Spring 2021 — Term Project

## Investigating Structural Variant Detection in Simple Tandem Repeat Regions Using HG002 (NA24385) Nanopore Data

The goal of this project is to test a variety of tools and tool parameters to optimize structural variant (SV) detection in simple tandem repeat (STR) regions of the human genome using Nanopore (ONT) data. As seen in the following image, this problem arises from the fact that even state-of-the-art aligners, such as Minimap2, often introduce gaps and insertions somewhat randomly in repetitive regions. Poorly localized indels in reads can reduce the accuracy of SV calls and generally makes genotyping more challenging. In trio data this issue can also result in inflated estimates of _de novo_ mutation (e.g. when SVs are called in children but are missed in parents). 



![alt text](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/misaligned_deletions.jpg?raw=true)


#### Methods


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
