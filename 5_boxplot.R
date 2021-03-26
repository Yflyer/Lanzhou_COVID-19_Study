source('0_function.R')

#################### save dir #######################
folder = paste0('5_boxplot')
#####################################################
##########################################
tran_sample =c("A1","A21","A3","A47","A64","A68","A70","A72","A73")
########################################
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

### mst plot
st_mr_mst <- ggplot(aes(x=group,y=MST.ij.bray),data = group_tnst_mr)+
  geom_boxplot(aes(color=Study),lwd=0.7)+
  geom_jitter(aes(color=Study),shape=16,alpha=0.3,position=position_jitter(0.2),size=1)
plot = st_mr_mst+PTs+ylim(0,1)+facet_grid(~Study,scales = 'free')
theme_set(theme_bw()+theme(legend.position = "none"))
plot+scale_color_economist()