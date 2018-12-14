library(stringr)
library(tidyverse)
condition_number <- function(Sigma){
  SEig <- eigen(Sigma,symmetric=T)
  return(max(SEig$val) / min(SEig$val))
}
dirname <- 'sim-data'
fnames <- list.files(path=dirname, pattern='sim-(exch|ar1|1band)Sigma.+-nsim500\\.RData')
betatypelist <- str_match(fnames, "-(flat|exponential|linear)(BETA)(-)")[,2]
plist <- as.integer(str_match(fnames,"(-p)([0-9]+)(-)")[,3])
Nlist <- as.integer(str_match(fnames,"(-N)([0-9]+)(-)")[,3])
SigmaTypeList <- str_match(fnames, "(-)(exch|ar1|1band)(Sigma)(-)")[,3]
SigmaGenList <- vector('list',length(fnames))
rholist <- as.numeric(str_match(fnames, "(-rho)([0-9]\\.[0-9]+)(-)")[,3])
nbeta <- NULL
measures <- c('tpr','fdp','nsel', 'ppv')
sres <- vector('list', length(fnames))
Gres <- vector('list', length(fnames))
Gdetres <- vector('list',length(fnames))
resfmt <- vector('list', length(fnames))
for (i in seq_along(fnames)){
  load(paste(dirname,fnames[i],sep='/'))
  resfmt_by_stat <- vector('list', length(res[[1]]))
  Gdet_by_stat <- vector('list', length(res[[1]]))
  sres_by_stat <- vector('list', length(res[[1]]))
  names(resfmt_by_stat) <- names(Gdet_by_stat) <-
    names(sres_by_stat) <- 
    names(res[[1]])
  nbeta <- c(nbeta, simparm$k)
  betatype <- betatypelist[i]
  rho <- rholist[i]
  SigmaGenList[[i]] <- SigmaGen
  for (sname in names(resfmt_by_stat)){
    cres <- lapply(res, function(rr) return(rr[[sname]]))
    cGdet <- do.call('rbind', lapply(cres,function(rr)return(exp(rr$lGdet))))
    cGdet <- as_tibble(cGdet) %>% mutate(
      sim_ix = 1:length(res),
      statname = sname,
      betatype = betatype,
      N = Nlist[i],
      p = plist[i],
      rho = rholist[i],
      SigmaType = SigmaTypeList[i]
    )
    Gdet_by_stat[[sname]] <- cGdet
    rtab <- do.call('rbind', lapply(cres, function(rr) return(do.call('rbind',rr[measures]))))
    
    ### Specturm of G, currently not saved in simulation
    ##Gres[[i]] <- do.call('rbind', lapply(cres, function(rr) return(cbind(j=1:nrow(rr$Geig), rr$Geig))))
    #Gres[[i]] <- as_tibble(Gres[[i]]) %>%
    #  mutate(sim_ix = rep(1:length(res), each = nrow(res[[1]]$Geig)),
    #         statname = rep(statnames, each=nrow(res[[1]]$Geig)),
    #         p = plist[i], N =Nlist[i], SigmaType = SigmaTypeList[i],
    #         rho =rholist[i])
    
    c_sres <- do.call('rbind', lapply(cres, function(rr) return(cbind(j=1:plist[i], rr$s))))
    c_sres <- as_tibble(c_sres) %>%
      mutate(sim_ix = rep(1:length(cres), each=plist[i]),
             statname = sname,
             p = plist[i],
             N= Nlist[i],
             betatype=betatype,
             SigmaType = SigmaTypeList[i],
             rho=rholist[i])
    sres_by_stat[[sname]] <- c_sres
    
    c_rfmt <- 
      as_tibble(rtab) %>% 
      mutate(meas = rownames(rtab),
             statname = sname,
             sim_ix = rep(1:length(res), each=length(measures)))%>%
      gather(-meas,-sim_ix, -statname, equi,sdp,Gdet,key='Stype',value='val') %>%
      spread(key=meas,value=val)
    c_rfmt$rho <- rho
    c_rfmt$N <- Nlist[i]
    c_rfmt$p <- plist[i]
    c_rfmt$SigmaType <- SigmaTypeList[i]
    c_rfmt$betatype <- betatype
    resfmt_by_stat[[sname]] <- c_rfmt
  }
  resfmt[[i]] <- bind_rows(resfmt_by_stat)
  sres[[i]] <- bind_rows(sres_by_stat)
  Gdetres[[i]] <- bind_rows(Gdet_by_stat)
}

resfmt <- bind_rows(resfmt)
sres <- bind_rows(sres)
# Gres <- bind_rows(Gres)
Gdetres <- bind_rows(Gdetres) 

Gres %>% 
  gather(equi, sdp, Gdet, key = 'method', value = 'ev') %>%
  group_by(sim_ix, statname, method, SigmaType, p, rho) %>%
  summarise(maxev = max(ev),
            minev = min(ev)) %>% ungroup %>%
  mutate(minev = ifelse(minev < 1e-8, 0, minev)) %>%
  left_join(Gdetres %>% gather(equi,sdp,Gdet,key='method',value='Gdet'),
            by=c('sim_ix','statname','p','rho','SigmaType','method')) -> Gresfmt
rm(Gres)

SigmaGenDat <- 
  tibble(condnum = sapply(SigmaGenList,condition_number),
         SigmaType = SigmaTypeList,
         popdet = sapply(SigmaGenList,function(x) return(exp(sum(log(eigen(x,only.values=T)$val))))),
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
  facet_wrap(~SigmaType+p)->plotfdr
ggsave("sim-fdr.pdf",
       plotfdr,  width=6.5,height=6,units='in')

resfmt %>% 
  group_by(Stype, rho, SigmaType, statname, p) %>%
  summarise(mtpr = mean(tpr)) %>%
  ggplot(aes(x=rho,y=mtpr,color=Stype, linetype=statname))+
  geom_line() + geom_point()+
  scale_Stype + ptheme + 
  facet_wrap(~SigmaType+p)->plottpr
ggsave("sim-tpr.pdf", plottpr, width=6.5,height=6,units='in')

nselfmt %>%
  filter(freq>0) %>%
  ggplot(aes(x=nsel, y=pfreq,color=Stype,shape=statname,linetype=statname))+
  geom_line()+
  geom_point()+
  scale_Stype+ptheme+
  facet_grid(SigmaType~rho+p)

Gresfmt %>%
  filter(p==10)%>%
  left_join(rename(SigmaGenDat, Sigmadet = popdet,popcond= condnum), by=c('SigmaType','rho','p')) %>%
  ggplot(aes(x=Sigmadet, y=Gdet, color=method))+
  stat_summary(fun.y=mean,geom='point')+
  scale_Stype + ptheme + 
  geom_abline(slope=1,intercept=0)+
  facet_wrap(~SigmaType,scales='free')
