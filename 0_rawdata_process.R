source('0_function.R')
### Meta information  
meta <- read.csv('0_data/metadata.txt',header = TRUE,row.names = 1,sep='\t')

project_pack = function(dt,taxa,meta){
  ###############################################
  ##### phylo object
  {
    OTU <- otu_table(dt,taxa_are_rows = TRUE) %>% .[,order(colnames(.))]
    TAX <- taxa %>% as.matrix(.) %>% tax_table(.) 
    ENV <- meta  %>% sample_data(.) %>% .[order(rownames(.)),]
  }
  
  ##################################
  project = phyloseq(OTU,TAX,ENV) 
  
  OTU <- otu_table(project)
  ENV = sample_data(project)
  nrow(OTU)
  ### remove site-singleton
  OTU = group_filter(OTU,group = ENV$Patient_Group,freq=2) %>% otu_table(.,taxa_are_rows = T)
  nrow(OTU)
 
  project = phyloseq(OTU,TAX,ENV) 
  
  ### other indexes and add alpha
  {
    sample_data(project)$Shannon=diversity(t(otu_table(project)),index='shannon')
    sample_data(project)$Simpson=diversity(t(otu_table(project)),index='simpson')
    sample_data(project)$Richness=specnumber(t(otu_table(project)))
  }
  project
}
################### load data   ################
prefix = 'metageno'

dt <- read.csv('0_data/metageno_taxa_abund.csv',header = TRUE,row.names = 1,sep=',')
taxa_total <- read.csv('0_data/metageno_taxa_info.csv',header = TRUE,sep=',',row.names = 1)
################### taxa filter ################
note = 'Bacteria'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom=='Bacteria'),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

################### taxa filter ################
note = 'Viruses'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom==note),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

################### taxa filter ################
note = 'Archaea'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom==note),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))


################### taxa filter ################
note = 'Fungi'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom=='Eukaryota'),]
# confirmed fungi on phylum levels and then remove unclassified
taxa$phylum %>% droplevels(.) %>% summary(.)
taxa = taxa[(taxa$phylum!='Eukaryota_unclassified'),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

################################################
prefix = 'metageno'
dt <- read.csv('0_data/metageno_gene_abund.txt',header = TRUE,row.names = 1,sep='\t')

note = 'KEGG'
taxa <- read.csv('0_data/metageno_gene_KEGG.txt',header = TRUE,sep='\t',row.names = NULL)
taxa <-taxa[!duplicated(taxa[,1]),]
rownames(taxa)=taxa$Query
project = project_pack(dt,taxa,meta)
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

note = 'GO'
taxa <- read.csv('0_data/metageno_gene_GO.txt',header = TRUE,sep='\t',quote = "",row.names = NULL)
taxa <-taxa[!duplicated(taxa[,1]),]
rownames(taxa)=taxa$Query
project = project_pack(dt,taxa,meta)
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

note = 'eggNOG'
taxa <- read.csv('0_data/metageno_gene_eggNOG.txt',header = TRUE,sep='\t',quote = "",row.names = NULL)
taxa <-taxa[!duplicated(taxa[,1]),]
rownames(taxa)=taxa$Query
project = project_pack(dt,taxa,meta)
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 


################### load data   ################
prefix = 'metatran'

dt <- read.csv('0_data/metatran_taxa_abund.csv',header = TRUE,row.names = 1,sep=',')
taxa_total <- read.csv('0_data/metatran_taxa_info.csv',header = TRUE,sep=',',row.names = 1)
################### taxa filter ################
note = 'Bacteria'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom=='Bacteria'),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

################### taxa filter ################
note = 'Viruses'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom==note),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

################### taxa filter ################
note = 'Archaea'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom==note),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

################### taxa filter ################
note = 'Fungi'
summary(taxa_total$kingdom)
taxa = taxa_total[(taxa_total$kingdom=='Eukaryota'),]
# confirmed fungi on phylum levels and then remove unclassified
taxa$phylum %>% droplevels(.) %>% summary(.)
taxa = taxa[(taxa$phylum!='Eukaryota_unclassified'),]

project = project_pack(dt,taxa,meta)
##########
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 

##############################################################
prefix = 'metatran'
dt <- read.csv('0_data/metatran_gene_abund.txt',header = TRUE,row.names = 1,sep='\t')

note = 'KEGG'
taxa <- read.csv('0_data/metatran_gene_KEGG.txt',header = TRUE,sep='\t',quote="",row.names = NULL)
taxa <-taxa[!duplicated(taxa[,1]),]
rownames(taxa)=taxa$Query
project = project_pack(dt,taxa,meta)
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

note = 'GO'
taxa <- read.csv('0_data/metatran_gene_GO.txt',header = TRUE,sep='\t',quote="",row.names = NULL)
taxa <-taxa[!duplicated(taxa[,1]),]
rownames(taxa)=taxa$Query
project = project_pack(dt,taxa,meta)
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note))

note = 'eggNOG'
taxa <- read.csv('0_data/metatran_gene_eggNOG.txt',header = TRUE,sep='\t',quote="",row.names = NULL)
taxa <-taxa[!duplicated(taxa[,1]),]
rownames(taxa)=taxa$Query
project = project_pack(dt,taxa,meta)
assign(paste0('project_',prefix,'_',note),project) 
save(list = c(paste0('project_',prefix,'_',note),'prefix','note'),file = paste0('project_',prefix,'_',note)) 
