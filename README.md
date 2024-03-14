chr18:79,843,656-79,942,420
chr7:72,236,306-72,464,119
chr1:2,415,168-2,481,835
chr17:40,046,026-40,068,248

https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr18%3A79843656%2D79942420&hgsid=2011572880_SS6Gi2CPTKwkh8Hk1NM7wFP25K9g

 https://genome.ucsc.edu/s/litingsong/hg38

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


# gtf 

# highlight 
#eca9cb  # i
#a9c1ec  # i target
#ecdda9  # candidate snp

# https://usegalaxy.org


track type=bigWig name="BasGan" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135f764a0e1e67e5e1a/display?to_ext=bigwig color=51,160,43
track type=bigWig name="NEX" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135ac2a4e78c725fbf4/display?to_ext=bigwig color=226,33,29
track type=bigWig name="Limbic" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135ac2a4e78c725fbf4/display?to_ext=bigwig color=251,154,153
track type=bigWig name="ARC" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721355bc30d570c5d537f/display?to_ext=bigwig color=178,89,41
track type=bigWig name="HAB" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135ac2a4e78c725fbf4/display?to_ext=bigwig color=251,191,111
track type=bigWig name="MDT" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135ac2a4e78c725fbf4/display?to_ext=bigwig color=251,126,2
track type=bigWig name="DRN" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721355bc30d570c5d537f/display?to_ext=bigwig color=202,177,214
track type=bigWig name="RMTG" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135a9cc8a99037e79fe/display?to_ext=bigwig color=127,115,173
track type=bigWig name="VTA" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135a9cc8a99037e79fe/display?to_ext=bigwig color=106,61,154

# /hpc/users/songl05/PF_25BR/codes

# statictics of brain specificity of abc links 
ABC_link_stat_uniquegene.R (unique ID); ABC_link_stat.R(all)


# validatation:
## peak width=1000
vali_eqtl_non5.R (for 5 and non5 links, seperately) barplot to compare 5 and non5
combined validation (/hpc/users/songl05/PF_25BR/codes/valid_ABC_combined_eid.R) 


# ldsc codes
## to get ldsc coef and p value
## nonregulatory/active: active but not in abc link
/sc/arion/projects/roussp01a/liting/Pf_25/LDSC_enhancer_1k.r (in enhancers regulating 5' non5' and nonregulatory enhancers, respectively)
/sc/arion/projects/roussp01a/liting/Pf_25/LDSC_promoter_1k.r (in promoter, 5' promoter, and non5' promoter)

### enhancers and promoters
LDSC_enrich_heatmap_ep.R   (final)

# gwas to enhancer-target: 
gwas_2PromoterI.R 
gwas2_target1k.R（index and ld）
*1.gwas_2target_scz_finemap.R* （finemap snp for scz）
*supp1.gwas_2target_sigwas_trait.R* (significant snp for all psychiatric disorders)


# gwas to promoter: 
gwas_2target_1k.R

# 

# scz links (ucsc)
scz_links.R (prepare input file for ucsc) index and ld
trait_links.R (index and ld ) for bd or scz
*2.UCSC_link_scz_finemap.R* (prepare input file for ucsc)  finemap snp, abc max
get_bw.sh (now using pengfei's aws ucsc) 


# scz snp target gene heatmap 
scz_combined_heatmap.R (local) index and ld
scz_ABCMAX_heatmap.R (local) index and ld, abc max
*3.scz_finemap_heatmap.R* (local) finemap snp, abc max



venn_target.R (overlap between fine regions, supple)
pops_trait.R: distribution of pops score across different diasese



track type=bigWig name="Basal glia TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135783630cff4a01d44/display?to_ext=bigwig color=51,160,43
track type=bigWig name="NEX TCF4" https://usegalaxy.org/api/datasets/f9cad7b01a472135b73edd62d47b2837/display?to_ext=bigwig" bigDataUrl= color=226,33,29
track type=bigWig name="Limbic TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135b6222e9e83db7c80/display?to_ext=bigwig color=251,154,153
track type=bigWig name="ARC TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721352ad231f9708d7a09/display?to_ext=bigwig color=249,252,154
track type=bigWig name="HAB TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135616cd5986d1ef2ed/display?to_ext=bigwig color=251,191,111
track type=bigWig name="MDT TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135179c777a85eae64e/display?to_ext=bigwig color=251,126,2
track type=bigWig name="DRN TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135287f608d821c1da4/display?to_ext=bigwig color=202,177,214
track type=bigWig name="RMTG TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721352ada6669a67677bd/display?to_ext=bigwig color=127,115,173
track type=bigWig name="VTA TCF4" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135fda6a8fcc6552d7e/display?to_ext=bigwig color=106,61,154

track type=bigWig name="Basal glia CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721351085f9dd199f68ea/display?to_ext=bigwig color=51,160,43
track type=bigWig name="NEX CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135d6bb9abf21d6cae3/display?to_ext=bigwig color=226,33,29
track type=bigWig name="Limbic CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135a62e94bec89bf162/display?to_ext=bigwig color=251,154,153
track type=bigWig name="ARC CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135e122579bd09b90ff/display?to_ext=bigwig color=249,252,154
track type=bigWig name="HAB CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721354d0e22bc369e5b81/display?to_ext=bigwig color=251,191,111
track type=bigWig name="MDT CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135800861a19f1623b1/display?to_ext=bigwig color=251,126,2
track type=bigWig name="DRN CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a47213529c5ec3d24fb33db/display?to_ext=bigwig color=202,177,214
track type=bigWig name="RMTG CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721350d0d13636607f458/display?to_ext=bigwig color=127,115,173
track type=bigWig name="VTA CACNB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135cffb1e782f654a0b/display?to_ext=bigwig color=106,61,154

track type=bigWig name="Basal glia ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721356fe240b4dda5f11f/display?to_ext=bigwig color=51,160,43
track type=bigWig name="NEX ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721359498e89b21d28bfe/display?to_ext=bigwig color=226,33,29
track type=bigWig name="Limbic ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721356c2eedcdbd5f11ca/display?to_ext=bigwig color=251,154,153
track type=bigWig name="ARC ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a47213571db2946c29f0479/display?to_ext=bigwig color=249,252,154
track type=bigWig name="HAB ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a4721352d2a9eb991ab4355/display?to_ext=bigwig color=251,191,111
track type=bigWig name="MDT ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135efee467986e4fcd6/display?to_ext=bigwig color=251,126,2
track type=bigWig name="DRN ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a47213516ddc329d46f644e/display?to_ext=bigwig color=202,177,214
track type=bigWig name="RMTG ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135e32242c65d3eaeac/display?to_ext=bigwig color=127,115,173
track type=bigWig name="VTA ZEB2" bigDataUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135ba438eec16503766/display?to_ext=bigwig color=106,61,154






