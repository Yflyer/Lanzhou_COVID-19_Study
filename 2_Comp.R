source('0_function.R')

#################### save dir #######################
folder = paste0('2_Comp')
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

########################################
merge_pack <- function(note,rank,tran_sample=tran_sample){
  load(paste0('project_metatran_',note))
  project = get(paste0('project_',prefix,'_',note))
  OTU_tran = otu_pack(project,rank,'RNA:',note)
  OTU_tran = OTU_tran[,paste('RNA:',tran_sample)]

  load(paste0('project_metageno_',note))
  project = get(paste0('project_',prefix,'_',note))
  OTU_geno = otu_pack(project,rank,'DNA:',note)
  OTU_geno=OTU_geno[,paste('DNA:',tran_sample)]
  
  ### must check in debug
  all(rownames(OTU_tran) == rownames(OTU_geno))
  all(rownames(OTU_tran) %in% rownames(OTU_geno))
  all(rownames(OTU_geno) %in% rownames(OTU_tran))
  ### notice: when pack up the overlap otu, bacterial dataset has more rows than the overlap, that's because bacterial dataset glom at genus level has replicate of genus so that more otu can match to genus within overlap. THIS OCASSION ONLY OCCURED IN BACTERIAL DATASET
  OTU_geno=OTU_geno[rownames(OTU_geno) %>% order(.),]
  OTU_tran=OTU_tran[rownames(OTU_tran) %>% order(.),]
  ### must check
  all(rownames(OTU_tran) == rownames(OTU_geno))
  OTU = cbind(OTU_geno,OTU_tran)
  OTU
}

table_pack <- function(OTU,tran_sample){
  ##############################
  if (nrow(OTU)>10) {
    abun_10_taxa = rowMeans(OTU) %>% sort(.,decreasing = TRUE) %>% .[1:10]
    OTU = OTU[rownames(OTU)%in% names(abun_10_taxa),]
  }
  Patient = rep(tran_sample,2)
  ###############################
  plot_data = data.frame()
  for (i in 1:ncol(OTU)) {
    OTU_abundance = OTU[,i] %>% as.numeric()
    Ln_OTU_abundance = (OTU[,i]+1) %>% log(.)
    Sample = rep(colnames(OTU)[i],nrow(OTU)) 
    #OTU_abundance = apply(OTU,1,function(x) sum(x[site %in% levels(site)[i]])/summary(site)[i] )
    Patient_ID = rep(Patient[i],nrow(OTU))
    Taxa=rownames(OTU)
    OOO = data.frame(OTU_abundance,Ln_OTU_abundance,Sample,Patient_ID,Taxa)
    plot_data = rbind(plot_data,OOO)
  }
  plot_data$Data=c(rep('DNA seq',nrow(OTU)*length(tran_sample)),rep('RNA seq',nrow(OTU)*length(tran_sample)))
  plot_data
}

plot_pack <- function(plot_data){
  comp = ggplot(data = plot_data,aes(x = Patient_ID, y = Ln_OTU_abundance,group = factor(Taxa)))
  barset = geom_bar(aes(fill = factor(Taxa)),stat='identity',width=1)
  area =  geom_area(aes(fill = factor(Taxa)),stat = "identity",alpha = 0.3)
  #line = geom_line(aes(colour = factor(Taxa)), position = "stack")
  comp+PTs +barset+area +sfm +facet_grid(~Data,scales = "free", space = "free")#+line
}

comp_pack <- function(note,rank,tran_sample){
  OTU = merge_pack(note,rank,tran_sample)
  plot_data = table_pack(OTU,tran_sample)
  plot_pack(plot_data)
}
########################################
PTs = theme(#axis.text.x=element_blank(),
  #axis.title=element_blank(),
  legend.position="right",
  legend.title = element_blank(),
  text = element_text(size=10),
  panel.spacing = unit(0.2, "lines"))
scm = scale_color_manual(values = c('orangered','#43a2ca','#018571','#3C5488','#f4a582','slategray','#80cdc1','#ca0020','#b2abd2','#5e3c99'))
sfm = scale_fill_manual(values = c('orangered','#43a2ca','#018571','#3C5488','#f4a582','slategray','#80cdc1','#ca0020','#b2abd2','#5e3c99'))

########################################
tran_sample =c("A1","A101","A21","A3","A47","A64","A68","A70","A72","A73")
########################################
rank = 'phylum'
note = 'Bacteria'
comp_pack(note,rank,tran_sample)
ggsave(filename = paste0(folder,'/Comp_',note,'_',rank,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank = 'genus'
note = 'Bacteria'
comp_pack(note,rank,tran_sample)
ggsave(filename = paste0(folder,'/Comp_',note,'_',rank,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank = 'phylum'
note = 'Fungi'
comp_pack(note,rank,tran_sample)
ggsave(filename = paste0(folder,'/Comp_',note,'_',rank,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank = 'genus'
note = 'Fungi'
comp_pack(note,rank,tran_sample)
ggsave(filename = paste0(folder,'/Comp_',note,'_',rank,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank = 'phylum'
note = 'Viruses'
comp_pack(note,rank,tran_sample)
ggsave(filename = paste0(folder,'/Comp_',note,'_',rank,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank = 'genus'
note = 'Viruses'
comp_pack(note,rank,tran_sample)
ggsave(filename = paste0(folder,'/Comp_',note,'_',rank,'.jpg'),dpi=900,width=18,height=19,units='cm')

