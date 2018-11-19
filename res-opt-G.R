library(stringr)
library(tidyverse)
condition_number <- function(Sigma){
  SEig <- eigen(Sigma,symmetric=T)
  return(max(SEig$val) / min(SEig$val))
}
fnames <- list.files(pattern='sim-opt-G-.+-nsim600\\.RData')
plist <- as.integer(str_match(fnames,"(-p)([0-9]+)(-)")[,3])
Nlist <- as.integer(str_match(fnames,"(-N)([0-9]+)(-)")[,3])
SigmaTypeList <- str_match(fnames, "(-G)(-)(exch|ar1)(-)")[,4]
SigmaGenList <- vector('list',length(fnames))
rholist <- as.numeric(str_match(fnames, "(-rho)([0-9]\\.[0-9]+)(-)")[,3])
nbeta <- NULL
measures <- c('tpr','fdp','nsel')
resfmt <- vector('list', length(fnames))
sres <- vector('list', length(fnames))
Gres <- vector('list', length(fnames))
for (i in seq_along(fnames)){
  load(fnames[i])
  rho <- rholist[i]
  SigmaGenList[[i]] <- SigmaGen
  statnames <- sapply(res,function(rr)return(rr$statname))
  nbeta <- c(nbeta, simparm$k)
  rtab <- do.call('rbind', lapply(res,function(rr) return(do.call('rbind',rr[measures]))))
  Gres[[i]] <- do.call('rbind', lapply(res, function(rr) return(cbind(j=1:nrow(rr$Geig), rr$Geig))))
  Gres[[i]] <- as_tibble(Gres[[i]]) %>%
    mutate(sim_ix = rep(1:length(res), each = 2*plist[i]),
           statname = rep(statnames, each=2*plist[i]),
           p = plist[i], N =Nlist[i], SigmaType = SigmaTypeList[i],
           rho =rholist[i])
  sres[[i]] <- do.call('rbind',lapply(res, function(rr) return(cbind(j=1:plist[i], rr$s))))
  sres[[i]] <- as_tibble(sres[[i]]) %>%
    mutate(sim_ix = rep(1:length(res), each=plist[i]),
           statname = rep(statnames, each=plist[i]),
           p = plist[i],
           N= Nlist[i],
           SigmaType = SigmaTypeList[i],
           rho=rholist[i])
  
  resfmt[[i]] <- 
    as_tibble(rtab) %>% 
    mutate(meas = rownames(rtab),
           statname = rep(statnames,each=length(measures)),
           sim_ix = rep(1:length(res), each=length(measures)))%>%
    gather(-meas,-sim_ix, -statname, equi,sdp,Gdet,key='Stype',value='val') %>%
    spread(key=meas,value=val)
  resfmt[[i]]$rho <- rho
  resfmt[[i]]$N <- Nlist[i]
  resfmt[[i]]$p <- plist[i]
  resfmt[[i]]$SigmaType <- SigmaTypeList[i]
}
resfmt <- bind_rows(resfmt)
sres <- bind_rows(sres)
Gres <- bind_rows(Gres)
Gres %>%
  gather(equi, sdp, Gdet, key = 'method', value = 'ev') %>%
  group_by(sim_ix, statname, method, SigmaType, p, rho) %>%
  summarise(maxev = max(ev),
            minev = min(ev),
            Gdim = n()) %>% ungroup %>%
  mutate(minev = ifelse(minev < 1e-8, 0, minev)) -> Gresfmt
rm(Gres)

SigmaGenDat <- 
  tibble(condnum = sapply(SigmaGenList,condition_number),
         SigmaType = SigmaTypeList,
         rho = rholist, p = plist)
resfmt %>% group_by(statname, Stype, p, rho, SigmaType) %>%  do(tibble(
  freq = table(factor(
    .$nsel,
    levels =
      0:(.$p[1]),
    labels =
      0:(.$p[1])
  )),
  nsel =
    as.integer(names(table(
      factor(.$nsel,
             levels =
               0:(.$p[1]),
             labels =
               0:(.$p[1]))
    ))),
  nsim =
    length(.$nsel)
)) %>%
  ungroup %>% mutate(pfreq = freq / nsim) -> nselfmt

ptheme <-
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background = element_rect(fill=NA,color=NA),
        panel.background = element_rect(fill=NA,color='grey'),
        strip.background = element_rect(fill=NA,color='grey'))
scale_Stype <- scale_color_brewer(palette='Set1',
                                  name='S Method',
                                  breaks=c('equi','Gdet','sdp'),
                                  labels=c("Equivariant","Max. det(G)","SDP"))

resfmt %>% group_by(Stype, rho, SigmaType,statname,p) %>%
  summarise(mfdr = mean(fdp)) %>%
  ggplot(aes(x=rho,y=mfdr,color=Stype,linetype=statname))+
  geom_line() + geom_point()+
  scale_Stype + ptheme + 
  geom_hline(aes(yintercept=simparm$FDR))+
  facet_wrap(~SigmaType+p,scales='free')

resfmt %>% group_by(Stype, rho, SigmaType, statname, p) %>%
  summarise(mtpr = mean(tpr)) %>%
  ggplot(aes(x=rho,y=mtpr,color=Stype, linetype=statname))+
  geom_line() + geom_point()+
  scale_Stype + ptheme + 
  facet_wrap(~SigmaType+p,scales='free')

nselfmt %>%
  filter(freq>0) %>%
  ggplot(aes(x=nsel, y=pfreq,color=Stype,shape=statname,linetype=statname))+
  geom_line()+
  geom_point()+
  scale_Stype+ptheme+
  facet_grid(SigmaType~rho+p)


Gresfmt %>%
  filter(p==10)%>%
  left_join(rename(SigmaGenDat, popcond= condnum), by=c('SigmaType','rho','p')) %>%
  ggplot(aes(x=log(popcond), y=log(maxev/minev),color=method))+
  stat_summary(fun.y=mean,geom='point')+
  scale_Stype + ptheme + 
  geom_abline(slope=1,intercept=0)+
  facet_wrap(~SigmaType,scales='free')
