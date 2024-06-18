#!/bin/bash
gwas_snp=/hpc/users/songl05/PF_25BR/data/gwas_sig.avinput
/sc/arion/projects/roussp01b/resources/software/annovar/table_annovar.pl $gwas_snp /sc/arion/projects/roussp01b/resources/software/annovar/humandb/ -buildver hg38 -out /hpc/users/songl05/PF_25BR/data/gwas_sig.anno   -protocol refGene -operation gx -nastring . 
