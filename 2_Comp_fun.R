source('0_function.R')

#################### save dir #######################
folder = paste0('2_Comp_fun')
dir.create(folder)
########################################
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
  dt
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
##########################################
########################################
PTs = theme(#axis.text.x=element_blank(),
  #axis.title=element_blank(),
  legend.position="right",
  legend.title = element_blank(),
  text = element_text(size=10),
  panel.spacing = unit(0.2, "lines"))
scm = scale_color_manual(values = c('orangered','#43a2ca','#018571','#3C5488','#f4a582','slategray','#80cdc1','#ca0020','#b2abd2','#5e3c99'))
sfm = scale_fill_manual(values = c('orangered','#43a2ca','#018571','#3C5488','#f4a582','slategray','#80cdc1','#ca0020','#b2abd2','#5e3c99'))
##########################################
tran_sample =c("A1","A21","A3","A47","A64","A68","A70","A72","A73")
########################################
rank = 'eggNOG_hit'
note = 'eggNOG'
label = 'RNA:'
taxa_list = "COGFunctionalCategoryDescription"
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_eggNOG')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label,tran_sample,taxa_list)

load('project_metageno_eggNOG')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label,tran_sample,taxa_list)

OTU_geno=OTU_geno[rownames(OTU_geno) %>% order(.),]
OTU_tran=OTU_tran[rownames(OTU_tran) %>% order(.),]
### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)

plot_data = table_pack(OTU,tran_sample)
plot_pack(plot_data)

ggsave(filename = paste0(folder,'/Comp_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')
########################################
########################################
rank = 'KEGG_hit'
note = 'KEGG'
label = 'RNA:'
taxa_list = "KEGGLevel2"
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_KEGG')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label,tran_sample,taxa_list)

load('project_metageno_KEGG')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label,tran_sample,taxa_list)

OTU_geno=OTU_geno[rownames(OTU_geno) %>% order(.),]
OTU_tran=OTU_tran[rownames(OTU_tran) %>% order(.),]
### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)

plot_data = table_pack(OTU,tran_sample)
plot_pack(plot_data)

ggsave(filename = paste0(folder,'/Comp_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')
########################################
########################################
rank = 'GO_hit'
note = 'GO'
label = 'RNA:'
taxa_list = "GO_Term"
########################################
load(paste0('overlap_',rank,'_',note))
overlap = get(paste0('overlap_',rank,'_',note))
########################################
load('project_metatran_GO')
project = get(paste0('project_',prefix,'_',note))
OTU_tran = otu_pack(project,rank,label,tran_sample,taxa_list)

load('project_metageno_GO')
project = get(paste0('project_',prefix,'_',note))
OTU_geno = otu_pack(project,rank,label,tran_sample,taxa_list)

OTU_geno=OTU_geno[rownames(OTU_geno) %>% order(.),]
OTU_tran=OTU_tran[rownames(OTU_tran) %>% order(.),]
### must check
all(rownames(OTU_tran) == rownames(OTU_geno))
OTU = cbind(OTU_geno,OTU_tran)

plot_data = table_pack(OTU,tran_sample)
plot_pack(plot_data)

ggsave(filename = paste0(folder,'/Comp_',note,'_',taxa_list,'.jpg'),dpi=900,width=18,height=19,units='cm')
########################################