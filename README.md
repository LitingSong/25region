# Codes for BI project (Fig 4-6, Fig S13-16,S19-25)

# 1. statictics of brain specificity of abc links (Fig 4a-c & Fig S13)

**ABC_link_stat.R**

source data: /sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata (60M)


# 2. validatation (Fig S14 & S15)

**vali_eqtl_non5.R** (for 5 and non5 links, seperately) barplot to compare 5 and non5 (Fig S15)

**valid_ABC_combined_eid.R** combined validation (Fig S14)

source data: abc link: /sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata
biccn: /sc/arion/projects/roussp01a/liting/Pf_25/biccn*pos.bedpe (19M)
gtex eqtl: /sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc/fine/GTex/GTEx_vcf.*_eqtls.txt.gz


# 3. ldsc codes (Fig 4d-e & Fig S16)

## 3.1 to get ldsc coef and p value

**LDSC_enhancer_1k.r** (in enhancers regulating 5', non5' and nonregulatory enhancers, respectively) (Fig S16)

**LDSC_promoter_1k.r** (in promoter, 5' promoter, and non5' promoter) (Fig S16)

source data: /sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata
source data: /sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/peaks/brain_Rg,'/',brain_Rg,'.bed

## 3.2 enhancers and promoters heatmap (Fig 4d-e )
LDSC_enrich_heatmap_ep.R   
source data:
Isoforms=c('5','n5','nonreg')
Enhancer:/sc/arion/projects/roussp01a/liting/ldsc_1k_,isoform,/meta-files/aggregatedPartInfo.tsv
promoter: /sc/arion/projects/roussp01a/liting/ldsc_1k_Promoter_',isoform,/meta-files/aggregatedPartInfo.tsv

# 4.Finemap (Fig 5 and Fig 6)
## step 1: gwas to enhancer-target: 

**4.1.1. gwas_2target_scz_finemap.R** （finemap snp for scz）
ABC_summary.Rdata
finemap: /sc/arion/projects/roussp01a/jaro/data/dual_assay/data/scz_finemap.tsv
closet gene: /sc/arion/projects/CommonMind/roussp01a/INGELHEIM/BI_25/data/answer.sort.bed

**4.1.2 gwas_2target_BIP_Finemap.R** （finemap snp for bip）
ABC_summary.Rdata
finemap: daner_bip_pgc3_only_finemap_only_finemap_all.txt.gz
closet gene: /sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/data/gwas_finemap_bip.answer.bed

**4.1.3 get_sig_gwas.R & gwas_2target_sigwas_trait.R & annovar.sh & find.closest.sh** (significant snp for all psychiatric disorders; both get target and statitics analysis, Fig S19-25)
gwas summary: /sc/arion/projects/roussp01b/resources/databases/gwas/',tt,'/',meta_gwas[tt,'gwasFilename']

## step 2: snp target gene heatmap  

### 4.2.1 scz snp target gene heatmap 

**scz_finemap_heatmap.R** (local) finemap snp, abc max (Fig 5A-C )
/sc/arion/projects/roussp01a/liting/Pf_25/output/SCZ_finemap_final_link.RData
/sc/arion/projects/roussp01a/jaro/data/dual_assay/data/scz_finemap.tsv
/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/POPS/files/final_pops/added_ridge_sz3.preds

### 4.2.2 bip snp target gene heatmap 

**bip_finemap_heatmap.R** (local) finemap snp, abc max (Fig 6A-C )
ABC_summary.Rdata
finemap: daner_bip_pgc3_only_finemap_only_finemap_all.txt.gz
closet gene: /sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/data/gwas_finemap_bip.answer.bed
pops: /sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/data/added_ridge_bip2.preds


## step 3: UCSC

**4.3.1 UCSC_link_scz_finemap.R** (prepare input file for ucsc)  finemap snp, abc max

**4.3.2 UCSC_link_bip_finemap.R**

get_bw.sh (now using pengfei's aws ucsc) 




