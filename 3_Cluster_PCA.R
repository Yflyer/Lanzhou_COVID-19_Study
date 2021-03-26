source('0_function.R')

#################### save dir #######################
folder = paste0('3_Cluster')
dir.create(folder)
# This code is to investigate the cluster of OTUs or Genes
otu_pack = function(project,rank,label,note){
  load(paste0('overlap_',rank,'_',note))
  overlap = get(paste0('overlap_',rank,'_',note))
  overlap_index = c(tax_table(project)[,rank] %in% overlap) %>% rownames(tax_table(project))[.]
  OTU = otu_table(project) %>% .[overlap_index,]
  taxa_list = tax_table(project)[,rank][overlap_index]
  dt = group_sum(OTU,taxa_list,margin = 2)
  colnames(dt)=paste(label,colnames(OTU))
  #rownames(dt)= unique(taxa_list) # since group sum will use order from unique(), so they should be the same
  dt
}
######################################
pcoa_table <- function(note,rank){
  label = 'RNA:'
  load(paste0('project_metatran_',note))
  project = get(paste0('project_',prefix,'_',note))
  tran_sample = sample_names(project)
  OTU_tran = otu_pack(project,rank,label,note)

  label = 'DNA:'
  load(paste0('project_metageno_',note))
  project = get(paste0('project_',prefix,'_',note))
  project = subset_samples(project, sample_names(project) %in% tran_sample)
  OTU_geno = otu_pack(project,rank,label,note)
  OTU_geno=OTU_geno[,paste(label,tran_sample)]
  
  ### must check
  all(rownames(OTU_tran) == rownames(OTU_geno))
  all(rownames(OTU_tran) %in% rownames(OTU_geno))
  all(rownames(OTU_geno) %in% rownames(OTU_tran))
  ### notice: when pack up the overlap otu, bacterial dataset has more rows than the overlap, that's because bacterial dataset glom at genus level has replicate of genus so that more otu can match to genus within overlap. THIS OCASSION ONLY OCCURED IN BACTERIAL DATASET
  OTU_geno=OTU_geno[rownames(OTU_geno) %>% order(.),]
  OTU_tran=OTU_tran[rownames(OTU_tran) %>% order(.),]
  ### must check
  all(rownames(OTU_tran) == rownames(OTU_geno))
  OTU = cbind(OTU_geno,OTU_tran)
  #########
  dis = vegdist(t(OTU),method = 'bray')
  pca=prcomp(dis, center = TRUE, scale = TRUE)
  dt_pcoa = pca$rotation[,1:5] %>% as.data.frame(.)
  
  ###########
  prefix = 'overlap'
  dt_pcoa$Data=c(rep('DNA seq',10),rep('RNA seq',10))
  dt_pcoa$Sample=rep(tran_sample,2)
  dt_pcoa$Taxonomic_level = rank
  dt_pcoa$Taxa = note
  dt_pcoa
}

##########################################
tran_sample =c("A1","A21","A3","A47","A64","A68","A70","A72","A73")
########################################
prefix = 'overlap'
rank = 'genus'
note = 'Bacteria'
dt_pcoa = pcoa_table(note,rank)
assign(paste0('dt_pcoa_',rank,'_',note),dt_pcoa)

rank = 'phylum'
note = 'Bacteria'
dt_pcoa = pcoa_table(note,rank)
assign(paste0('dt_pcoa_',rank,'_',note),dt_pcoa)
############## Fungi #####
rank = 'genus'
note = 'Fungi'
dt_pcoa = pcoa_table(note,rank)
assign(paste0('dt_pcoa_',rank,'_',note),dt_pcoa)

rank = 'phylum'
note = 'Fungi'
dt_pcoa = pcoa_table(note,rank)
assign(paste0('dt_pcoa_',rank,'_',note),dt_pcoa)

############## Viruses #####
rank = 'genus'
note = 'Viruses'
dt_pcoa = pcoa_table(note,rank)
assign(paste0('dt_pcoa_',rank,'_',note),dt_pcoa)

rank = 'phylum'
note = 'Viruses'
dt_pcoa = pcoa_table(note,rank)
assign(paste0('dt_pcoa_',rank,'_',note),dt_pcoa)

###############################################################################
dt_pcoa = Reduce(rbind,list(dt_pcoa_genus_Bacteria,dt_pcoa_genus_Fungi,dt_pcoa_genus_Viruses,dt_pcoa_phylum_Bacteria,dt_pcoa_phylum_Fungi,dt_pcoa_phylum_Viruses))
dt_pcoa$Study = factor(dt_pcoa$Study,levels = c('Latitude (North America)','Altitude (SNNR)','Altitude (HKV)'))
dt_pcoa$Gradient = as.factor(dt_pcoa$Gradient)
###############################################################################

#### plot
PTs <- theme(panel.background=element_blank(),
             text = element_text(size=10),
             axis.title=element_blank(),
             #legend.position="none",
             panel.border = element_rect(colour = "grey", fill=NA, size=1))
cycle = stat_ellipse(geom = "polygon",aes(fill= Data),level = 0.9,alpha=0.1)
round = stat_ellipse(aes(color= Data),level = 0.9)
ssm <- scale_shape_manual(values=c(16,21))

ggplot(dt_pcoa,aes(PC1, PC2)) + 
  geom_point(aes(shape = Data, color = Data),size=3,alpha=0.5)+ssm+
  PTs+facet_grid(Taxonomic_level~Taxa,scales = 'free') +cycle+round
ggsave(filename = paste0(folder,'/Cluster_PCOA_total.jpg'),dpi = 900,width=18,height=13,units='cm')

plot = dt_pcoa %>% group_by(Taxonomic_level,Taxa) %>% 
  do(gg ={ ggplot(.,aes(PC1, PC2)) + 
      geom_point(aes(shape = Data, color = Data),size=3,alpha=0.5)+
      PTs+ssm+round+cycle#+facet_grid(~Study,scales = 'free_x')
    ggsave(filename = paste0(folder,'/Cluster_PCOA_',unique(.$Taxa),'_',unique(.$Taxonomic_level),'.jpg'),dpi = 900,width=8,height=7,units='cm')
  }) 

