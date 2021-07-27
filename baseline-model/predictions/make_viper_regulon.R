setwd("~/Documents/Projects/DREAMs/d2020_ctd2/predictions/")

make_regulon=function(x){
  result=list()
  result$tfmode=x$DIRECTION
  names(result$tfmode)=x$DRUG
  result$likelihood=rep(1,length(x$DRUG))
  result
}

data=read.csv('regulons/moa_regulon_ctrp.csv',sep = ',',header = T,row.names = 1)
data=split(data[c('DRUG','DIRECTION')],data$GENE)
viper_geneset=lapply(data,make_regulon)
save(viper_geneset,file='regulons//moa_regulon_ctrp.Rdata')
load('regulons/moa_regulon_rep.Rdata')
library(viper)

for (fname in c('sim_ctrp.csv','sim_ctrp_0.csv','sim_ctrp_1.csv','sim_ctrp_2.csv')){
  gex_cv_cor=read.csv(paste0('similarity/',fname),sep=',',header=T,row.names=1)
  tf_activity=viper(eset = gex_cv_cor,regulon = viper_geneset,nes = 1,minsize = 1)
  write.csv(x =tf_activity,file = paste0('enrichments/',fname))
}



