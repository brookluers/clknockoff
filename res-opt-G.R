library(stringr)
library(tidyverse)
fnames <- list.files(pattern='sim-opt-G-.+-nsim1000\\.RData')
N <- as.integer(unique(str_match(fnames, '(-N)([0-9]+)')[,3]))
p <- as.integer(unique(str_match(fnames, '(-p)([0-9]+)')[,3]))
nbeta <- NULL
measures <- c('tpr','fdp','nsel')
resfmt <- vector('list', length(fnames))
sres <- vector('list', length(fnames))
for (i in seq_along(fnames)){
  load(fnames[i])
  rho <- SigmaGen[1,2]
  nbeta <- c(nbeta, simparm$k)
  rtab <- do.call('rbind',lapply(res,function(rr) return(do.call('rbind',rr[measures]))))
  sres[[i]] <- do.call('rbind',lapply(res, function(rr) return(rr$svec)))
  resfmt[[i]] <- 
    as_tibble(rtab) %>% 
    mutate(meas = rownames(rtab),
           sim_ix = rep(1:length(res), each=length(measures)))%>%
    gather(-meas,-sim_ix,equi,sdp,Gdet,key='type',value='val') %>%
    spread(key=meas,value=val)
  resfmt[[i]]$rho <- rho
}
resfmt <- bind_rows(resfmt)

resfmt %>% group_by(type, rho) %>%
  summarise(mfdr = mean(fdp)) %>%
  ggplot(aes(x=rho,y=mfdr,color=type))+
  geom_point()
resfmt %>% group_by(type,rho) %>%
  summarise(mtpr=mean(tpr))%>%
  ggplot(aes(x=rho,y=mtpr,color=type))+
  geom_line()
