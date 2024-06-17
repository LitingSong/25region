# /sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/codes

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





# ucsc 
https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A72233434%2D72458433&hgsid=1941877242_0HuisATst6fZ53elf5YB9PXWAjiZ

# lead snp
track name="Lead SNP" itemRgb="On"

# ld snp
track name="LD SNP" itemRgb="On"

# candidate snp
track name="Candidate SNP" itemRgb="On"

# isoform
track name="Promoter Isoforms" itemRgb="On" useScore=1000

# inter
track type=interact name="ABC links" interactDirectional=false maxHeightPixels=200:100:50 visibility=full

# gwas
track type=bedGraph name="GWAS" visibility=full color=200,100,0 priority=20


# https://usegalaxy.org



# /hpc/users/songl05/PF_25BR/codes
# /sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/codes

# 1. statictics of brain specificity of abc links (Fig 4a-c & Fig S13)
 **ABC_link_stat.R (all)**


# 2. validatation (Fig S14 & S15)
## peak width=1000
vali_eqtl_non5.R (for 5 and non5 links, seperately) barplot to compare 5 and non5
combined validation (/hpc/users/songl05/PF_25BR/codes/valid_ABC_combined_eid.R) 


# 3. ldsc codes (Fig 4d-e )
## to get ldsc coef and p value
## nonregulatory/active: active but not in abc link
/sc/arion/projects/roussp01a/liting/Pf_25/LDSC_enhancer_1k.r (in enhancers regulating 5' non5' and nonregulatory enhancers, respectively)
/sc/arion/projects/roussp01a/liting/Pf_25/LDSC_promoter_1k.r (in promoter, 5' promoter, and non5' promoter)

## enhancers and promoters heatmap (Fig 4d-e )
LDSC_enrich_heatmap_ep.R   

# 4.Finemap (Fig 5 and Fig 6)
# step 1: gwas to enhancer-target: 
gwas_2PromoterI.R 

gwas2_target1k.R（index and ld）

**1.gwas_2target_scz_finemap.R** （finemap snp for scz）

**1.gwas_2target_BIP_Finemap.R**

**0.get_sig_gwas.R -> 1.gwas_2target_sigwas_trait.R** (significant snp for all psychiatric disorders; both get target and statitics analysis)

# gwas to promoter: 

gwas_2target_1k.R


# step 2: snp target gene heatmap  

## 2.1 scz snp target gene heatmap 

scz_combined_heatmap.R (local) index and ld

scz_ABCMAX_heatmap.R (local) index and ld, abc max

**scz_finemap_heatmap.R** (local) finemap snp, abc max


## 2.2 bip snp target gene heatmap 

**bip_finemap_heatmap.R**


# step 3: UCSC

scz_links.R (prepare input file for ucsc) index and ld

trait_links.R (index and ld ) for bd or scz

**3.1 UCSC_link_scz_finemap.R** (prepare input file for ucsc)  finemap snp, abc max

**3.2 UCSC_link_bip_finemap.R**

get_bw.sh (now using pengfei's aws ucsc) 





# ucsc 
https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A72233434%2D72458433&hgsid=1941877242_0HuisATst6fZ53elf5YB9PXWAjiZ

# lead snp
track name="Lead SNP" itemRgb="On"

# ld snp
track name="LD SNP" itemRgb="On"

# candidate snp
track name="Candidate SNP" itemRgb="On"

# isoform
track name="Promoter Isoforms" itemRgb="On" useScore=1000

# inter
track type=interact name="ABC links" interactDirectional=false maxHeightPixels=200:100:50 visibility=full

# gwas
track type=bedGraph name="GWAS" visibility=full color=200,100,0 priority=20


# https://usegalaxy.org

