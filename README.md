# Codes for BI project

# 1. statictics of brain specificity of abc links (Fig 4a-c & Fig S13)

**ABC_link_stat.R**


# 2. validatation (Fig S14 & S15)

**vali_eqtl_non5.R** (for 5 and non5 links, seperately) barplot to compare 5 and non5 (Fig S15)

**valid_ABC_combined_eid.R** combined validation (Fig S14)


# 3. ldsc codes (Fig 4d-e & Fig S16)

## 3.1 to get ldsc coef and p value

## nonregulatory/active: active but not in abc link

**LDSC_enhancer_1k.r** (in enhancers regulating 5', non5' and nonregulatory enhancers, respectively) (Fig S16)

**LDSC_promoter_1k.r** (in promoter, 5' promoter, and non5' promoter) (Fig S16)


## 3.2 enhancers and promoters heatmap (Fig 4d-e )
LDSC_enrich_heatmap_ep.R   

# 4.Finemap (Fig 5 and Fig 6)
## step 1: gwas to enhancer-target: 

**4.1.1. gwas_2target_scz_finemap.R** （finemap snp for scz）

**4.1.2 gwas_2target_BIP_Finemap.R** （finemap snp for bip）

**4.1.3 get_sig_gwas.R & gwas_2target_sigwas_trait.R** (significant snp for all psychiatric disorders; both get target and statitics analysis, Fig S19-25)


## step 2: snp target gene heatmap  

### 4.2.1 scz snp target gene heatmap 

**scz_finemap_heatmap.R** (local) finemap snp, abc max (Fig 5A-C )

### 4.2.2 bip snp target gene heatmap 

**bip_finemap_heatmap.R** (local) finemap snp, abc max (Fig 6A-C )


## step 3: UCSC

**4.3.1 UCSC_link_scz_finemap.R** (prepare input file for ucsc)  finemap snp, abc max

**4.3.2 UCSC_link_bip_finemap.R**

get_bw.sh (now using pengfei's aws ucsc) 




