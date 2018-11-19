library(stringr)
library(tidyverse)
fnames <- list.files(pattern='sim-opt-G-.+-nsim1000\\.RData')
plist <- as.integer(str_match(fnames,"(-p)([0-9]+)(-)")[,3])
Nlist <- as.integer(str_match(fnames,"(-N)([0-9]+)(-)")[,3])
SigmaList <- str_match(fnames, "(-G)(-)(exch|ar1)(-)")[,4]
rholist <- as.numeric(str_match(fnames, "(-rho)([0-9]\\.[0-9]+)(-)")[,3])
nbeta <- NULL
measures <- c('tpr','fdp','nsel')
resfmt <- vector('list', length(fnames))
sres <- vector('list', length(fnames))
Gres <- vector('list', length(fnames))
for (i in seq_along(fnames)){
  load(fnames[i])
  rho <- rholist[i]
  nbeta <- c(nbeta, simparm$k)
  rtab <- do.call('rbind', lapply(res,function(rr) return(do.call('rbind',rr[measures]))))
  Gres[[i]] <- do.call('rbind', lapply(res, function(rr) return(cbind(j=1:nrow(rr$Geig), rr$Geig))))
  Gres[[i]] <- as_tibble(Gres[[i]]) %>%
    mutate(sim_ix = rep(1:length(res), each = 2*plist[i]),
           p = plist[i], N =Nlist[i], SigmaType = SigmaList[i],
           rho =rholist[i])
  sres[[i]] <- do.call('rbind',lapply(res, function(rr) return(cbind(j=1:plist[i], rr$s))))
  sres[[i]] <- as_tibble(sres[[i]]) %>%
    mutate(sim_ix = rep(1:length(res), each=plist[i]),
           p = plist[i],
           N= Nlist[i],
           SigmaType = SigmaList[i],
           rho=rholist[i])
  
  resfmt[[i]] <- 
    as_tibble(rtab) %>% 
    mutate(meas = rownames(rtab),
           sim_ix = rep(1:length(res), each=length(measures)))%>%
    gather(-meas,-sim_ix,equi,sdp,Gdet,key='type',value='val') %>%
    spread(key=meas,value=val)
  resfmt[[i]]$rho <- rho
  resfmt[[i]]$N <- Nlist[i]
  resfmt[[i]]$p <- plist[i]
  resfmt[[i]]$SigmaType <- SigmaList[i]
}
resfmt <- bind_rows(resfmt)
sres <- bind_rows(sres)
Gres <- bind_rows(Gres)

resfmt %>% group_by(type, rho,SigmaType,p) %>%
  summarise(mfdr = mean(fdp)) %>%
  ggplot(aes(x=rho,y=mfdr,color=type))+
  geom_line() + geom_point()+
  facet_wrap(~SigmaType+p,scales='free')


resfmt %>% group_by(type, rho,SigmaType,p) %>%
  summarise(mtpr = mean(tpr)) %>%
  ggplot(aes(x=rho,y=mtpr,color=type))+
  geom_line() + geom_point()+
  facet_wrap(~SigmaType+p,scales='free')

resfmt %>% 
  ggplot(aes(x=nsel,y=..density..,fill=type))+
  geom_histogram(bins=10,position='dodge')+
  facet_wrap(~SigmaType+p)

Gres%>%
  gather(equi,sdp,Gdet,key='method',value='ev')%>%
  group_by(sim_ix,method,SigmaType,p,rho)%>%
  summarise(maxev=max(ev),
            minev=min(ev)) %>% ungroup %>%
  mutate(minev = ifelse(minev < 1e-8, 0, minev)) ->Gresfmt

Gresfmt %>%
  filter(p==10)%>%
  ggplot(aes(x=minev,y=..density..,fill=method))+
  geom_histogram(position='dodge',bins=20)+
  facet_grid(SigmaType~rho,scales='free')

Gresfmt %>%
  filter(p==10)%>%
  ggplot(aes(x=rho,y=log(maxev/minev),color=method))+
  geom_point(shape=1)+
  facet_wrap(~SigmaType,scales='free')
