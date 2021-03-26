source('0_function.R')
### 
load('project_metageno')
project = get(paste0('project_',prefix)) 

OTU = otu_table(project)
TAX <- tax_table(project)
ENV = sample_data(project)
##############
Group = ENV$Group
summary(as.factor(Group))

Trend = Taxa = Effect.size = P.value = diff.healthy = diff.pneumonia = c() 
alpha = 0.05
for (i in 1:nrow(OTU)) {
  abund =OTU[i,]%>% as.numeric(.)
  anova.result = aov(abund~Group)
  anova.summary=summary(anova.result)[[1]]
  if (anova.summary[1,5]<alpha) {
    tuk.result = TukeyHSD(anova.result)$Group
    #### both significant on postive control & negative control
    if (tuk.result[1,4]<alpha & tuk.result[2,4]<alpha) {
    Taxa=c(Taxa,paste(TAX[i,2],TAX[i,5],TAX[i,7]))
    P.value = c(P.value,anova.summary[1,5])
    Effect.size=c(Effect.size,anova.summary[1,4])
    diff.healthy = c(diff.healthy,tuk.result[1,1])
    diff.pneumonia = c(diff.pneumonia,tuk.result[2,1])
    # trend
    if (tuk.result[1,1]<0){
      Trend=c(Trend,'Positive')
    } else {
      Trend=c(Trend,'Negative')
    }
  }}}
significant_taxa = data.frame(Taxa,Trend,Effect.size,P.value,diff.healthy,diff.pneumonia)
write.csv(file = 'significant_taxa_group.csv',significant_taxa)

##### check the effect of swab
Group = paste0(ENV$Group,'+',ENV$Swabs)
summary(as.factor(Group))

Trend = Taxa = Effect.size = P.value = diff.healthy = diff.pneumonia = c() 
alpha = 0.05
for (i in 1:nrow(OTU)) {
  abund =OTU[i,]%>% as.numeric(.)
  anova.result = aov(abund~Group)
  anova.summary=summary(anova.result)[[1]]
  if (anova.summary[1,5]<alpha) {
    tuk.result = TukeyHSD(anova.result)$Group
    #### both significant on postive control & negative control
    if (tuk.result[1,4]<alpha & tuk.result[5,4]<alpha) {
      Taxa=c(Taxa,paste(TAX[i,2],TAX[i,5],TAX[i,7]))
      P.value = c(P.value,anova.summary[1,5])
      Effect.size=c(Effect.size,anova.summary[1,4])
      diff.healthy = c(diff.healthy,tuk.result[1,1])
      diff.pneumonia = c(diff.pneumonia,tuk.result[5,1])
      # trend
      if (tuk.result[1,1]<0){
        Trend=c(Trend,'Positive')
      } else {
        Trend=c(Trend,'Negative')
      }
    }}}
significant_taxa_II = data.frame(Taxa,Trend,Effect.size,P.value,diff.healthy,diff.pneumonia)
write.csv(file = 'significant_taxa_group_swab.csv',significant_taxa)

library(ggplot2)