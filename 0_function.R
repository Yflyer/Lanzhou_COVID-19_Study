########## 0_rawdata_process ######################
# required package
library(ggplot2)
library(ggsci)
library(dplyr)
library(vegan)
library(geosphere)
library(picante)
library(dplyr)
library(phyloseq)

### for quick test, generate a random dataframe:
### table=data.frame(replicate(5,sample(0:1,5,rep=TRUE)))

### for quick test, generate a random matrix:
### table=matrix(rbinom(10*1000, 1, .5), ncol=10)

############ function used #############################
group_count<-function(table,group){
  ### this function is used to count the occurrence frequency of each row in a table in a specific grouping condition
  ### such as: group_count(OTU,treatment)
  ### then we will get each OTU occurrence assigned on treatment levels on rows
  group= group %>% as.factor(.) %>% droplevels(.) # make factor
  count_table = matrix(ncol= nlevels(group),nrow = nrow(table),dimnames = list(rownames(table),levels(group)))
  for (i in 1:ncol(count_table)) {
    dt = table[,group == colnames(count_table)[i]]
    count_table[,i] = apply(dt, 1, function(x) sum(x>0))
  }
  count_table
}

counts_check = function(data_row,factor,level_counts=3,factor_counts=1){
  # this will return the boolean list of  qualified level along data_row which meets requirement.
  factor= factor %>% as.factor(.) %>% droplevels(.) # make factor
  level_check_list = sapply(factor,function(x) data_row[,factor == x]>level_counts) 
  sum(level_check_list)>factor_counts # this will print whether data row meets requirement of level_counts and factor_counts.
}


group_filter = function(table,group,freq=1,cut=1,group_cut=F){
  ### this function is used to filter row data when you group count data table to meet requirement of group count frequency
  ### such as: group_filter(OTU,treatment,freq=3,cut=1,group_cut=F) the otu occured at least 3 times in a level of treatment will be filtered into new datatable.
  ### if set group_cut=T, OTU which occurred lower than cut value in all groups will be discarded
  group= group %>% as.factor(.) %>% droplevels(.) # make factor
  table[is.na(table)]=0 ### must rm NA. since NA will make following judgement NA
  count_table = matrix(ncol= nlevels(group),nrow = nrow(table),dimnames = list(rownames(table),levels(group))) # make count table for grouping
  for (i in 1:ncol(count_table)) {
    dt = table[,group == colnames(count_table)[i]] # make sub-dt for each group
    count_table[,i] = apply(dt, 1, function(x) sum(x>0)) # sum up this sub-dt at each row from this assigned group
  }
  ### get filter index ###
  if (group_cut==FALSE) {
    otu_check_list = apply(count_table, 1, function(table_row) sum(table_row) >freq)
    table = table[otu_check_list,]
  }else{
    otu_check_list = apply(count_table, 1, function(table_row) sum(table_row>=freq) >= cut)
    table = table[otu_check_list,]
  }
  ###
  table
}



group_sum<- function(table,group,margin=1){
  ### we can use margin to cluster samples(row) or otu(col) according to a grouping factor
  group[is.na(group)]='None'
  group = as.factor(group)
  level= unique(group)
  if(margin==1){
    dt = sapply(level,function(one_level) rowSums(table[,group==one_level]))
    colnames(dt)=unique(group)
  }else {
    dt = sapply(level,function(one_level) colSums(table[group==one_level,])) %>% t(.)
    rownames(dt)=unique(group)
  }
  as.data.frame(dt)
}

group_mean<- function(table,group){
  group= group %>% as.factor(.) %>% droplevels(.)
  dt = matrix(ncol= nlevels(group),nrow = nrow(table),
              dimnames = list(rownames(table),levels(group))) %>% as.data.frame(.)
  for (i in levels(group)) {
    dt[[i]] = apply(table,1,function(x) sum(x[group==i]))/summary(group)[i]
  }
  dt
}

###### get three column from dist martix
col3m<-function(m,dist_name = 'dist',diag=F){
  ###
  m = as.matrix(m)
  pair = data.frame(row=rownames(m)[row(m)[upper.tri(m,diag = diag)]], 
                    col=colnames(m)[col(m)[upper.tri(m,diag = diag)]], 
                    dist = m[upper.tri(m)])
  colnames(pair)[3]=dist_name
  pair
}

###### change three column back diag matrix or symmetric matrix 
reshape_col3m<-function(m,name1='name1',name2='name2',value,diag=T,symmetric=F){
  ### this function is very useful to reshpae the col3m data back to a diag matrix or symmetric matrix, which is able to generate dist or run mantel in a proper format
  # m must be dataframe containing pairwise information
  # name1: str, pairwise name1
  # name1: str, pairwise name2
  # value: str, the target index to reshape back
  # format the col3m data
  m = m[,c(name1,name2,value)]
  colnames(m) = c('name1','name2','value')
  # reshape
  reshape.matrix = tidyr::spread(m,key='name2',value='value') %>% as.data.frame(.)
  m.dist =  reshape.matrix[,-1]
  rownames(m.dist) = reshape.matrix[['name1']]
  
  if(diag){
  # add diagonal
  m.dist[,rownames(m.dist)[!rownames(m.dist)%in%colnames(m.dist)]]=NA
  m.dist[colnames(m.dist)[!colnames(m.dist)%in%rownames(m.dist)],]=NA
  m.dist=m.dist[order(rownames(m.dist)),order(colnames(m.dist))]  
  }

  if(symmetric){
  # make symmetric
  m.dist[upper.tri(m.dist)] <- t(m.dist)[upper.tri(m.dist)] 
  }
  as.matrix(m.dist)
}


###### calculate the geospatial distance by lat and lon; need package: dplyr, geoshpere
geo_dist = function(dis,Lon='Lon',Lat='Lat',unit = 'km'){
  dis = dis[,c(Lon,Lat)] %>% as.matrix(.)
  dis = apply(dis, 1, function(y) apply(dis, 1, function(x) distHaversine(y,x)))
  if (unit == 'km'){dis/1000}
}

###### label the pair by factor
label_pair = function(dt,meta,row='row',col='col',factor,sample){
  # for diagonal-pair data to add pair tag
  factor_row = sapply(dt[[row]], function(x) meta[[factor]][which(sample==x)])
  factor_col = sapply(dt[[col]], function(x) meta[[factor]][which(sample==x)])
  dt$pair=paste(factor_row,'-',factor_col)
}

######### stat tools
get_aov_cof = function(formula,dt=NULL,cof='Pr(>F)',var=1){
  # var: string or number, to locate the cof of var in the dataframe
  # y : numeric vector
  # x : treat
  # waiting for extend to be complex formula
  result = aov(as.formula(formula),data = dt) %>% summary(.)
  result[[1]][var,cof]
}

get_mantel = function(var_list,var_name,Env_dt,tar,dist.type='euclid'){
  # var_list = list or vector of the var you want to test or partially test
  # var_name = name of the vars or category to name the dt
  # Env_dt = all the env numeric data for mantel
  # tar = the dist of OTU, etc, bray
  # dist.type = env dist type
  mantel.R = mantel.P = p.mantel.R =p.mantel.P = c()
  for (i in 1:length(var_list)) {
    var = unlist(var_list[i])
    tar.dist = tar
    var.dist = Env_dt[,var] %>% scale(.) %>% vegdist(.,dist.type)
    ctr.dist = Env_dt[,!(all %in% var)] %>% scale(.) %>% vegdist(.,dist.type)
    all.dist = Env_dt[,all] %>% scale(.) %>% vegdist(.,dist.type)
    mantel=mantel(tar.dist,var.dist)
    mantel.R[i] = mantel$statistic
    mantel.P[i] = mantel$signif
    p.mantel =  mantel.partial(tar.dist,all.dist,ctr.dist)
    p.mantel.R[i] = p.mantel$statistic
    p.mantel.P[i] = p.mantel$signif
  }
  mantel = data.frame(var_name,mantel.R,mantel.P,p.mantel.R,p.mantel.P)
}

#################
get_var_name <- function(var) {
  deparse(substitute(var))
}

labelp <- function(x){
  if (x <= 0.05){x = "p < 0.05"} 
  else if (x <= 0.1){x = "0.05 < p < 0.1"} 
  else {x = "p > 0.1"}
  x = factor(x,levels=c('p > 0.1','0.05 < p < 0.1','p < 0.05'))
  }
labelpstar <- function(x){
  if (x <= 0.001){x = "***"}
  else if (x <= 0.01){x = "**"}
  else if (x <= 0.05){x = "*"} 
  else if (x <= 0.1){x = "."}
  else {x = ""}
}

labelpH <- function(x){
  if (x < 5){x = "4 - 5"} 
  else if (x < 6){x = "5 - 6"} 
  else {x = "6 - 7"}
  x = factor(x,levels=c('4 - 5','5 - 6','6 - 7'))
}

labelalt <- function(x){
  if (x < 1000){x = "0-1000m"} 
  else {x = "1000-1300m"}
  x = factor(x,levels=c('0-1000m','1000-1300m'))
}
label_ln_dis <- function(x){
  if (x < 6){x = "Small-scale (<150m)"} 
  else {x = "Mid-scale (100-500km)"}
  x = factor(x,levels=c('Small-scale (<150m)','Mid-scale (100-500km)'))
}

normalization<-function(y){
  x<-y[!is.na(y)]
  x<-(x - min(x)) / (max(x) - min(x))
  y[!is.na(y)]<-x
  return(y)}

############
writeOutput.F <- function(output.title, modelname) {
  sink(file="modeloutputs.txt", append=TRUE) #creates an empty file; if TRUE, it will append to the file
  print("##########################", row.names=FALSE)
  print(output.title, row.names=FALSE)
  print("##########################",  row.names=FALSE)
  summary(modelname) #write summary of the model to text file 
}



