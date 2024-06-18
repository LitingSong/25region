#!/bin/bash


#cat /sc/arion/projects/roussp01b/resources/databases/PEC_RNAseq_REFERENCE_HG38_GENCODE30/gencode.v30.annotation.gtf | 
#  gtf2bed | grep -w gene > gencodeGenes.bed


/sc/arion/projects/roussp01b/resources/software/bedtools/bin/sortBed -i /hpc/users/songl05/PF_25BR/data/gwas_sig.avinput  > /hpc/users/songl05/PF_25BR/data/gwas_sig.avinput.sort.bed
#sortBed -i gencodeGenes.bed > gencodeGenes.sort.bed

/hpc/users/songl05/software/bedops/bin/closest-features --closest --dist  /hpc/users/songl05/PF_25BR/data/gwas_sig.avinput.sort.bed  /sc/arion/projects/CommonMind/roussp01a/INGELHEIM/BI_25/data/gencodeGenes.sort.bed > /hpc/users/songl05/PF_25BR/data/gwas_sig.answer.bed


