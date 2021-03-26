source('0_function.R')
library(ggvenn)
library(ggplot2)

########### mkdir #############
folder = '1_Venn_function'
dir.create(folder)
####################################

venn_pack = function(project,rank,label,group){
  ##################################
  project = subset_samples(project,Patient_Group=='COVID-19')
  ##################################
  OTU=otu_table(project)
  TAX=tax_table(project)
  ENV=sample_data(project)
  rank_index=!duplicated(TAX[,rank])
  dt = OTU[rank_index,] %>% as.data.frame(.)
  rownames(dt)=TAX[rank_index,rank]
  ##################################
  group = paste(label,ENV[[group]])
  venn.dt = group_count(dt,group) %>% as.data.frame(.)
  venn.dt$taxa = rownames(venn.dt)
  venn.dt
}
vennplot_pack = function(venn.dt){
  venn.dt[is.na(venn.dt)]=0
  a = rownames(venn.dt)[venn.dt$`DNA: COVID-19`>0]
  b = rownames(venn.dt)[venn.dt$`RNA: COVID-19`>0]
  venn.plot = list(`DNA: COVID-19`=a,`RNA: COVID-19`=b)
}
overlap_pack = function(venn.dt){
  venn.dt[is.na(venn.dt)]=0
  a = rownames(venn.dt)[venn.dt$`DNA: COVID-19`>0]
  b = rownames(venn.dt)[venn.dt$`RNA: COVID-19`>0]
  overlap = intersect(a,b)
}
###############################
rank = 'KEGG_hit'
label = 'RNA:'
group = 'Patient_Group'
load('project_metatran_KEGG')
project = get(paste0('project_',prefix,'_',note))
venn.dt.tran = venn_pack(project,rank,label,group)
tran_sample = sample_names(project)

####### match #################
label = 'DNA:'
load('project_metageno_KEGG')
project = get(paste0('project_',prefix,'_',note))
project = subset_samples(project, sample_names(project) %in% tran_sample)
venn.dt.geno = venn_pack(project,rank,label,group)

venn.dt = full_join(venn.dt.geno,venn.dt.tran,by='taxa')
rownames(venn.dt)=venn.dt$taxa
######## overlap ###########
overlap = overlap_pack(venn.dt)
assign(paste0('overlap_',rank,'_',note),overlap) 
save(list = c(paste0('overlap_',rank,'_',note),'rank','note'),file = paste0('overlap_',rank,'_',note))
############################
venn.plot = vennplot_pack(venn.dt)
ggvenn(venn.plot,
       show_percentage = F,text_size = 3,set_name_size = 3,
       stroke_color = "black",stroke_size = 1,
       fill_color = c("#ffb2b2","#b2e7cb"),
       set_name_color = c("#ff0000","#4a9b83"))
ggsave(filename = paste0(folder,'/venn_',group,'_',rank,'_',note,'.jpg'),dpi=900,width=8,height=7,units='cm')
ggsave(filename = paste0(folder,'/venn_',group,'_',rank,'_',note,'.pdf'),dpi=900,width=8,height=7,units='cm')

###############################
rank = 'GO_hit'
label = 'RNA:'
group = 'Patient_Group'
load('project_metatran_GO')
project = get(paste0('project_',prefix,'_',note))
venn.dt.tran = venn_pack(project,rank,label,group)
tran_sample = sample_names(project)

####### match #################
label = 'DNA:'
load('project_metageno_GO')
project = get(paste0('project_',prefix,'_',note))
project = subset_samples(project, sample_names(project) %in% tran_sample)
venn.dt.geno = venn_pack(project,rank,label,group)

venn.dt = full_join(venn.dt.geno,venn.dt.tran,by='taxa')
rownames(venn.dt)=venn.dt$taxa
######## overlap ###########
overlap = overlap_pack(venn.dt)
assign(paste0('overlap_',rank,'_',note),overlap) 
save(list = c(paste0('overlap_',rank,'_',note),'rank','note'),file = paste0('overlap_',rank,'_',note))
############################
venn.plot = vennplot_pack(venn.dt)
ggvenn(venn.plot,
       show_percentage = F,text_size = 3,set_name_size = 3,
       stroke_color = "black",stroke_size = 1,
       fill_color = c("#ffb2b2","#b2e7cb"),
       set_name_color = c("#ff0000","#4a9b83"))
ggsave(filename = paste0(folder,'/venn_',group,'_',rank,'_',note,'.jpg'),dpi=900,width=8,height=7,units='cm')
ggsave(filename = paste0(folder,'/venn_',group,'_',rank,'_',note,'.pdf'),dpi=900,width=8,height=7,units='cm')

###############################
rank = 'eggNOG_hit'
label = 'RNA:'
group = 'Patient_Group'
load('project_metatran_eggNOG')
project = get(paste0('project_',prefix,'_',note))
venn.dt.tran = venn_pack(project,rank,label,group)
tran_sample = sample_names(project)

####### match #################
label = 'DNA:'
load('project_metageno_eggNOG')
project = get(paste0('project_',prefix,'_',note))
project = subset_samples(project, sample_names(project) %in% tran_sample)
venn.dt.geno = venn_pack(project,rank,label,group)

venn.dt = full_join(venn.dt.geno,venn.dt.tran,by='taxa')
rownames(venn.dt)=venn.dt$taxa
######## overlap ###########
overlap = overlap_pack(venn.dt)
assign(paste0('overlap_',rank,'_',note),overlap) 
save(list = c(paste0('overlap_',rank,'_',note),'rank','note'),file = paste0('overlap_',rank,'_',note))
############################
venn.plot = vennplot_pack(venn.dt)
ggvenn(venn.plot,
       show_percentage = F,text_size = 3,set_name_size = 3,
       stroke_color = "black",stroke_size = 1,
       fill_color = c("#ffb2b2","#b2e7cb"),
       set_name_color = c("#ff0000","#4a9b83"))
ggsave(filename = paste0(folder,'/venn_',group,'_',rank,'_',note,'.jpg'),dpi=900,width=8,height=7,units='cm')
ggsave(filename = paste0(folder,'/venn_',group,'_',rank,'_',note,'.pdf'),dpi=900,width=8,height=7,units='cm')
