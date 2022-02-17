## APSE - C++ implementation
<img src="https://user-images.githubusercontent.com/37526521/154415881-271124f1-f514-4815-a5df-dccb5faffcba.png" width=530 height=250>

This repository provides the C++ implementation of the Assessment Program for Systematic Error causing phylogenetic incongruence of gene markers.
>**Motivation**: Software program developed to date to evaluate and improve phylogenetic reliability did not actively include analysis of the quality of molecular biological data in the early stage of analyze. As the amount of data is rapidly increasing, it remains difficult to evaluate the reliability of the phylogenetic reconstruction using only the existing analysis pipeline. The initial dataset has the potential to contain noise to some extent, and the quality of the molecular biological dataset should be evaluated in consideration of potential problems known as systematic biases that cause phylogenetic inaccuracies. Ultimately, through the information dericed from this program, a phylogenetic tree with improved reliability that is closer to reality is generated from molecular data.
>
>**Objective**: This software program estimates the reliability of the phylogenetic tree with specific values using the phylogenetic tree and input dataset (multi-aligned file format: .fas, .aln) to estimate the evolutionary relationship between individuals. **The three sub-goals set for this are**: The first goal is to develop a component that parses the sequence alignment file of multiple data into an easy-to-compute structure for performing phylogenetic analysis. The second goal is to implement the calculation and estimation of potential biases with statistical properties to infer the accuracy of phylogenetic trees using characters (base sequence) extracted from the multi-aligned file, and to provide the corresponding output values independently. The third goal is to develop a single module system by integrating data parsing logic and data calculation logic.

## Dependencies
+ C++ 11
+ C++ STL
+ Clang 12

## Downloading datasets
The datasets used in this study belong to Terebelliformia (Annelida), Daphniid (Arthoropoda) and Mammalia, and taxa with misplacement problem within each clade were included in the analysis.

1. TriAA, TriTer hypothesis - EF1α (1163 positions, 6 species), 18S rDNA (1897 positions, 7 species), mtDNA (16580 positions, 7 species), and 28S rDNA (4387 positions, 7 species)

2. Daphniid (Crustacea) - 16S rDNA (502 positions, 9 species), 28S rDNA (4702 positions, 10 species)

3. Glire hypothesis - A2AB (2391 positions, 30 species), IRBP (9719 positions, 27 species), vWF (8663 positions, 30 species) and a concatenation of three marker mitochondrial genes (2789 positions, 26 species), which are 12s rRNA, tRNA valine (MT-TV), and 16s rRNA

The datasets were collected from the National Center for Biotechnology Information (NCBI, https://www.ncbi.nlm.nih.gov). Corresponding fasta-formatted data was collected using a ncbi-acc-download version 0.2.6 Python tool.
<pre><code>ncbi-acc-download --format fasta “accession number”</code></pre>

## Results
<img width="467" alt="result" src="https://user-images.githubusercontent.com/37526521/154423671-d31ba575-efc9-417d-943d-b9d0770fe640.png">

## Citation
If you refer to this work (including dataset) useful for your work, please cite electronic resource in [Seoul National University library](https://library.snu.ac.kr/).

