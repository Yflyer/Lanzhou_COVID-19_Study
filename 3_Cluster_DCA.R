source('0_function.R')

#################### save dir #######################
folder = paste0('3_Cluster_fun')
dir.create(folder)

########################################
otu_pack = function(project,rank,label,tran_sample){
  ########################################
  project = subset_samples(project, sample_names(project) %in% tran_sample)
  overlap_index = tax_table(project)[,rank] %in% overlap
  ########################################
  dt=otu_table(project) %>% as.data.frame(.) %>% .[overlap_index,]
  colnames(dt)=paste(label,colnames(dt))
  ########################################
  dt = group_sum(dt,overlap,margin = 2)
  dt=dt[rownames(dt) %>% order(.),]
  dt
}
### must check
dca_pack=function(OTU_tran,OTU_geno){
  all(rownames(OTU_tran) == rownames(OTU_geno))
  OTU = cbind(OTU_geno,OTU_tran)
  
  ########################################
  dis = vegdist(t(OTU),method = 'bray')
  dca <- decorana(dis)
  dt_dca = dca$rproj %>% as.data.frame(.)
  
  ########################################
  dt_dca$Data=c(rep('DNA seq',ncol(OTU_geno)),rep('RNA seq',ncol(OTU_tran)))
  dt_dca$Sample=rep(tran_sample,2)
  dt_dca$Taxonomic_level = rank
  dt_dca$Taxa = note
  dt_dca
}
######################################

######################################
###################### GO ############
rank = 'GO_hit'
note = 'GO'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_GO')
label = 'RNA:'
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label,tran_sample)
load('project_metageno_GO')
label = 'DNA:'
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label,tran_sample)

dt_dca = dca_pack(OTU_tran,OTU_geno)
assign(paste0('dt_dca_',note),dt_dca)
######################################
#################### KEGG ############
rank = 'KEGG_hit'
note = 'KEGG'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_KEGG')
label = 'RNA:'
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label,tran_sample)
load('project_metageno_KEGG')
label = 'DNA:'
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label,tran_sample)

dt_dca = dca_pack(OTU_tran,OTU_geno)
assign(paste0('dt_dca_',note),dt_dca)
######################################
################## eggCOG ############
rank = 'eggNOG_hit'
note = 'eggNOG'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_eggNOG')
label = 'RNA:'
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label,tran_sample)
load('project_metageno_eggNOG')
label = 'DNA:'
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label,tran_sample)

dt_dca = dca_pack(OTU_tran,OTU_geno)
assign(paste0('dt_dca_',note),dt_dca)

###############################################################################
dt_dca = Reduce(rbind,list(dt_dca_eggNOG,dt_dca_GO,dt_dca_KEGG))
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

ggplot(dt_dca,aes(DCA1, DCA2)) + 
  geom_point(aes(shape = Data, color = Data),size=3,alpha=0.5)+ssm+
  PTs+facet_grid(~Taxa,scales = 'free') +cycle+round
ggsave(filename = paste0(folder,'/Cluster_DCA_total.jpg'),dpi = 900,width=18,height=13,units='cm')

plot = dt_dca %>% group_by(Taxa) %>% 
  do(gg ={ ggplot(.,aes(DCA1, DCA2)) + 
      geom_point(aes(shape = Data, color = Data),size=3,alpha=0.5)+
      PTs+ssm+round+cycle#+facet_grid(~Study,scales = 'free_x')
    ggsave(filename = paste0(folder,'/Cluster_DCA_',unique(.$Taxa),'.jpg'),dpi = 900,width=8,height=7,units='cm')
  }) 

