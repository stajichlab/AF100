# Introduction

This repository contains the code used to conduct variant and comparative genomic analysis associated with the manuscript 
"_Aspergillus fumigatus_ In-Host HOG pathway mutation for Cystic Fibrosis Lung Microenvironment Persistence". Brandon Ross, Lotus A. Lofgren, Alix Ashare,  Jason E. Stajich, Robert A. Cramer. _Submitted_

# Usage

1. `sbatch pipeline/00_init.sh` - to download and build local databases for analysis
2. `sbatch -a 1-121 pipeline/01_aln_arrayjob.sh` - run the bwa alignment
