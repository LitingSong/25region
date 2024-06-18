# get significant gwas snp 
library('SNPlocs.Hsapiens.dbSNP144.GRCh37')
library('BSgenome')




setwd('~/PF_25BR/data/')
meta_gwas <- read.delim('metadata.tsv')
rownames(meta_gwas) <- meta_gwas$gwasAcronym

mytraits=list(`Neuropsychiatric and behavioral`=c("sz3",'asd', "bip2",  "mdd_wray", #"bdAndSczShared",
                                                  "adhd", "drinking", "smoking",
                                                  "insomn2", "neu2"),
              Neurological=c('alz2','als','ms','pd_without_23andMe','stroke'),
              `Non-brain related`=c('bmi','cad','cd','height','ibd','ra','dm2','uc'))

psy_t <- mytraits$`Neuropsychiatric and behavioral`


gsums_sig <- c()
for (tt in c((psy_t))){
  gwas_sum <- paste0('/sc/arion/projects/roussp01b/resources/databases/gwas/',tt,'/',meta_gwas[tt,'gwasFilename'])
  
  if (tt!="insomn2"){gsum <- read.delim(gzfile(gwas_sum))}
  if (tt=="insomn2"){gsum <- read.delim(gzfile(gwas_sum),sep=' ')}
  
  if(!tt %in%c('drinking','smoking','neu2','adhd','mdd_wray')){gsum <- gsum[,c('CHR','SNP','BP','P','A1','A2')]}
  if(tt %in%c('drinking','smoking')){gsum <- gsum[,c('CHR','MarkerName','POS','Pval','A1','A2')]}
  if(tt %in%c('neu2')){gsum <- gsum[,c('CHR','SNP','POS','P','A1','A2')]}
  if(tt %in%c('mdd_wray')){gsum1 <- gsum[,c('MARKERNAME','P','A1','A2')] %>%subset(P < 5e-8);
  gsum2 <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37,gsum1$MARKERNAME,ifnotfound="drop");
  gsum2 <- as.data.frame(gsum2);
  gsum <- merge(gsum1,gsum2, by.x='MARKERNAME',by.y='RefSNP_id')[,c('seqnames','MARKERNAME','pos','P','A1','A2')]
  
  }
  colnames(gsum) <- c('CHR','SNP','BP','P','A1','A2')
  
  
  gsum_sig <- subset(gsum,P < 5e-8)
  print(c(tt, nrow(gsum_sig)))
  gsum_sig$trait <- tt
  gsums_sig <- rbind(gsums_sig,gsum_sig)
}




save(gsums_sig, file='gsums_sig.RData')
