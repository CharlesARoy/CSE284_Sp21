# CSE284 Spring 2021 — Term Project

## Investigating Structural Variant Detection in Simple Tandem Repeat Regions with Long Read Data

The goal of this project is to test a variety of tools and tool parameters to optimize structural variant (SV) detection in simple tandem repeat (STR) regions of the human genome using Nanopore (ONT) data. As seen in the following image, this problem arises from the fact that even state-of-the-art aligners, such as Minimap2, can introduce gaps and insertions somewhat randomly in repetitive regions, making it challenging for SV detection algorithms to properly aggregate supporting reads when calling SVs. 


![alt text](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/misaligned_deletions.jpg?raw=true)


#### Methods


![alt text](https://github.com/CharlesARoy/CSE284_Sp21/blob/main/Flowchart.png?raw=true)
