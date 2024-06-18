# Codes for BI project (Fig 4-6, Fig S13-16,S19-25)

codes: NC_codes
source data: NC_data

# 1. statictics of brain specificity of abc links (Fig 4a-c & Fig S13)

**ABC_link_stat.R**

source data: ABC_summary.Rdata 


# 2. validatation (Fig S14 & S15)

**vali_eqtl_non5.R** (for 5 and non5 links, seperately) barplot to compare 5 and non5 (Fig S15)

**valid_ABC_combined_eid.R** combined validation (Fig S14)

source data: ABC_summary.Rdata; biccn ** pos.bedpe; GTEx_vcf. ** _eqtls.txt.gz


# 3. ldsc codes (Fig 4d-e & Fig S16)

## 3.1 to get ldsc coef and p value

**LDSC_enhancer_1k.r** (in enhancers regulating 5', non5' and nonregulatory enhancers, respectively) (Fig S16)

**LDSC_promoter_1k.r** (in promoter, 5' promoter, and non5' promoter) (Fig S16)

source data: ABC_summary.Rdata; 'brain_Rg'.bed

## 3.2 enhancers and promoters heatmap (Fig 4d-e )
LDSC_enrich_heatmap_ep.R   
source data is the output of step 3.1.

# 4.Finemap (Fig 5 and Fig 6)
## step 1: gwas to enhancer-target: 

**4.1.1. gwas_2target_scz_finemap.R** （finemap snp for scz）
source data: ABC_summary.Rdata; finemap: scz_finemap.tsv; closet gene: answer.sort.bed

**4.1.2 gwas_2target_BIP_Finemap.R** （finemap snp for bip）
source data: ABC_summary.Rdata; finemap: daner_bip_pgc3_only_finemap_only_finemap_all.txt.gz; closet gene: gwas_finemap_bip.answer.bed

**4.1.3 get_sig_gwas.R & gwas_2target_sigwas_trait.R & annovar.sh & find.closest.sh** (significant snp for all psychiatric disorders; both get target and statitics analysis, Fig S19-25)
gwas summary: data can be download based on metadata.tsv

## step 2: snp target gene heatmap  

### 4.2.1 scz snp target gene heatmap 

**scz_finemap_heatmap.R** (local) finemap snp, abc max (Fig 5A-C )
source data: SCZ_finemap_final_link.RData; scz_finemap.tsv; added_ridge_sz3.preds

### 4.2.2 bip snp target gene heatmap 

**bip_finemap_heatmap.R** (local) finemap snp, abc max (Fig 6A-C )
source data: ABC_summary.Rdata; finemap: daner_bip_pgc3_only_finemap_only_finemap_all.txt.gz; closet gene: gwas_finemap_bip.answer.bed; pops: added_ridge_bip2.preds


## step 3: UCSC

**4.3.1 UCSC_link_scz_finemap.R** (prepare input file for ucsc)  finemap snp, abc max

**4.3.2 UCSC_link_bip_finemap.R**





