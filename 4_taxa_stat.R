source('0_function.R')

#################### save dir #######################
folder = paste0('2_Comp')
dir.create(folder)
# This code is to investigate the cluster of OTUs or Genes
otu_pack = function(project,rank,label,tran_sample,taxa_list){
  ########################################
  project = subset_samples(project, sample_names(project) %in% tran_sample)
  overlap_index = tax_table(project)[,rank] %in% overlap
  tax_table(project) = tax_table(project)[overlap_index,]
  ########################################
  OTU=otu_table(project)
  TAX=tax_table(project)
  
  dt = group_sum(OTU,TAX[,taxa_list],margin = 2)
  colnames(dt)=paste(label,colnames(dt))
  dt=dt[rownames(dt) %>% order(.),]
  dt
}

table_pack <- function(OTU,tran_sample){
  ##############################
  #if (nrow(OTU)>10) {
  #  abun_10_taxa = rowMeans(OTU) %>% sort(.,decreasing = TRUE) %>% .[1:10]
  #  OTU = OTU[rownames(OTU)%in% names(abun_10_taxa),]
  #}
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

########################################
ttest_pack <- function(plot_data){
  t.result = data.frame()
  Taxa = plot_data$Taxa
  Taxa_info = t_value = mean_difference = p_value = c()
  for (i in 1:nlevels(Taxa)) {
    test = plot_data[Taxa==levels(Taxa)[i],]
    ttest = t.test(test$Ln_OTU_abundance[1:9],test$Ln_OTU_abundance[10:18], paired = T)
    if (ttest$p.value < 0.05) {
      Taxa_info = c(Taxa_info,levels(Taxa)[i])
      t_value = c(t_value,ttest$statistic)
      mean_difference = c(mean_difference,ttest$estimate)
      p_value = c(p_value,ttest$p.value) # maybe wrong in another set
    }
  }
  t.result = data.frame(Taxa_info,t_value,mean_difference,p_value)
  t.result
}
plot_pack = function(t.result,OTU,plot_data){
  #### get the 10 abundant taxa #####
  Taxa_info = t.result[['Taxa_info']]
  OTU = OTU[rownames(OTU) %in% Taxa_info,]
  if (nrow(OTU)>10) {
    Taxa_info = rowMeans(OTU) %>% sort(.,decreasing = TRUE) %>% names(.) %>% .[1:10]
  }
  ###################################
  plot_data = filter(plot_data, Taxa %in% Taxa_info) %>% droplevels(.)
  ###################################
  reads_plot <- ggplot(aes(x=Taxa,y=Ln_OTU_abundance),data = plot_data)+
    geom_boxplot(aes(color=Data),lwd=0.8, 
                 outlier.shape = 16,
                 outlier.size = 1,
                 outlier.alpha = 0.5)+
    geom_jitter(aes(color=Data),shape=16,alpha=0.5,position=position_jitter(0.1),size=1)
  plot = reads_plot+PTs
  plot+ coord_flip()+ scale_color_npg()+ scale_fill_npg()
}
##########################################
PTs <- theme(#text = element_text(size=12),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  legend.position = "bottom",
  panel.border = element_rect(colour = "grey", fill=NA, size=1))
########################################
########################################
tran_sample =c("A1",'A21',"A3","A47","A64","A68","A70","A72","A73")
########################################
rank = taxa_list = 'phylum'
note = 'Bacteria'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_Bacteria')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label='RNA:',tran_sample,taxa_list)

load('project_metageno_Bacteria')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label='DNA:',tran_sample,taxa_list)

### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)
###################################
plot_data = table_pack(OTU,tran_sample) %>% filter(.,Taxa != 'None') %>% droplevels(.)
t.result = ttest_pack(plot_data)
write.csv(file = paste0(folder,'/Ttest_',taxa_list,'.csv'),t.result)
###################################
plot_pack(t.result,OTU,plot_data)
ggsave(filename = paste0(folder,'/Significant_pathway_boxplot_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank =taxa_list = 'genus'
note = 'Bacteria'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_Bacteria')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label='RNA:',tran_sample,taxa_list)

load('project_metageno_Bacteria')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label='DNA:',tran_sample,taxa_list)

### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)
###################################
plot_data = table_pack(OTU,tran_sample) %>% filter(.,Taxa != 'None') %>% droplevels(.)
t.result = ttest_pack(plot_data)
write.csv(file = paste0(folder,'/Ttest_',taxa_list,'.csv'),t.result)
###################################
plot_pack(t.result,OTU,plot_data)
ggsave(filename = paste0(folder,'/Significant_pathway_boxplot_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')


########################################
rank = taxa_list = 'phylum'
note = 'Fungi'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_Fungi')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label='RNA:',tran_sample,taxa_list)

load('project_metageno_Fungi')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label='DNA:',tran_sample,taxa_list)

### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)
###################################
plot_data = table_pack(OTU,tran_sample) %>% filter(.,Taxa != 'None') %>% droplevels(.)
t.result = ttest_pack(plot_data)
write.csv(file = paste0(folder,'/Ttest_',taxa_list,'.csv'),t.result)
###################################
plot_pack(t.result,OTU,plot_data)
ggsave(filename = paste0(folder,'/Significant_pathway_boxplot_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank =taxa_list = 'genus'
note = 'Fungi'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_Fungi')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label='RNA:',tran_sample,taxa_list)

load('project_metageno_Fungi')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label='DNA:',tran_sample,taxa_list)

### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)
###################################
plot_data = table_pack(OTU,tran_sample) %>% filter(.,Taxa != 'None') %>% droplevels(.)
t.result = ttest_pack(plot_data)
write.csv(file = paste0(folder,'/Ttest_',taxa_list,'.csv'),t.result)
###################################
plot_pack(t.result,OTU,plot_data)
ggsave(filename = paste0(folder,'/Significant_pathway_boxplot_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank = taxa_list = 'phylum'
note = 'Viruses'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_Viruses')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label='RNA:',tran_sample,taxa_list)

load('project_metageno_Viruses')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label='DNA:',tran_sample,taxa_list)

### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)
###################################
plot_data = table_pack(OTU,tran_sample) %>% filter(.,Taxa != 'None') %>% droplevels(.)
t.result = ttest_pack(plot_data)
write.csv(file = paste0(folder,'/Ttest_',taxa_list,'.csv'),t.result)
###################################
plot_pack(t.result,OTU,plot_data)
ggsave(filename = paste0(folder,'/Significant_pathway_boxplot_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')

########################################
rank =taxa_list = 'genus'
note = 'Viruses'
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_Viruses')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label='RNA:',tran_sample,taxa_list)

load('project_metageno_Viruses')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label='DNA:',tran_sample,taxa_list)

### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)
###################################
plot_data = table_pack(OTU,tran_sample) %>% filter(.,Taxa != 'None') %>% droplevels(.)
t.result = ttest_pack(plot_data)
write.csv(file = paste0(folder,'/Ttest_',taxa_list,'.csv'),t.result)
###################################
plot_pack(t.result,OTU,plot_data)
ggsave(filename = paste0(folder,'/Significant_pathway_boxplot_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')
