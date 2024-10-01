

load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/rnaseq/test/neuron_gene_norm_none_BIC_4_in_0.05_CPM_1_in_0.2/files/rmhind_residualized_DxBrainRegion_EffectKept.Rdata')

#Promoter isoform 

load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/rnaseq/test/neuron_proactive_promoter_norm_none_BIC_4_in_0.05_CPM_3_in_0.4/files/rev1_qc.Rdata')
#load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/exp_promoters_CPM.Rdata')
load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/rnaseq/promoter/proactive_readcount.Rdata')
# expression matrix

allInfo$Brain_region_H1 = 'ForeBr'
allInfo[allInfo$Brain_region %in% c("DS","NAC","GP"),"Brain_region_H1"] = "BasGan"
allInfo[allInfo$Brain_region %in% c("RMTG","VTA","DRN"),"Brain_region_H1"] = "MidDien"
allInfo[allInfo$Brain_region %in% c("HAB", "MDT", "ARC"), "Brain_region_H1"] = "MidDien"
allInfo$Brain_region_H1b = 'NEX'
allInfo$Brain_region_H1b[allInfo$Brain_region_H1=='BasGan'] = 'BasGan'
allInfo[allInfo$Brain_region %in% c("AMY",'HIPP'),"Brain_region_H1b"] = "Limbic"
allInfo[allInfo$Brain_region %in% c("RMTG","VTA","DRN"),"Brain_region_H1b"] = "MidBr"
allInfo[allInfo$Brain_region %in% c("HAB", "MDT", "ARC"), "Brain_region_H1b"] = "Dien"
load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')

load('/sc/arion/projects/roussp01a/liting/Pf_25/output/SCZ_finemap_final_link.RData')
load('/sc/arion/projects/roussp01a/liting/Pf_25/output/bip_finemap_final_link.RData')

library(limma)
voomobj = voom(readcount, design=NULL)
CPM <- voomobj$E

plot_isoexp <- function(G,trait){
  if(trait=='bip'){final_link <- bip_final_link}
  if(trait=='SCZ'){final_link <- SCZ_final_link}
  
  ENSG <- unique(subset(final_link ,TargetGene_symbol==G)[,'gene'])
  ENSG <- unique(subset(final_link ,TargetGene_symbol==G)[,'gene'])
  print(ENSG)
  exp_m <- melt(CPM[names(subset(promoter_anno, gene_id==ENSG)),])
  exp_m$order <- as.data.frame(promoter_anno)[as.character(exp_m$Var1),'order']
  exp_m$Brain_region <- allInfo[exp_m$Var2,"Brain_region_H1b"]
  exp_m$Brain_region[exp_m$Brain_region%in%c("Dien","MidBr")] <- 
    allInfo[exp_m$Var2[exp_m$Brain_region%in%c("Dien","MidBr")],'Brain_region']
  exp_m$cell <- allInfo[exp_m$Var2,"cell"]
  
  exp_m_n <- subset(exp_m, cell!='glia')
  
  library(ggpubr)
  
  
  p1 <- ggplot(subset(exp_m_n),# & order%in%c(1:2)),
         aes(x=Brain_region,y=value, fill=as.character(order)))+
    geom_boxplot(width=0.5)+#geom_violin()+
    guides(fill=guide_legend(title='Promoter isoform'))+
    theme_classic() +ylab('log2(CPM)') + xlab('Brain region')+theme(legend.position = 'top')
  ggsave(p1,file = paste0("/hpc/users/songl05/PF_25BR/figures/supp/sp4_",G,".pdf"), width = 4, height = 4)
  
  
}

plot_isoexp('WDR82',trait='bip')
plot_isoexp('FURIN',trait='bip')
plot_isoexp('ABHD2',trait='SCZ')
plot_isoexp('CALN1',trait='SCZ')

