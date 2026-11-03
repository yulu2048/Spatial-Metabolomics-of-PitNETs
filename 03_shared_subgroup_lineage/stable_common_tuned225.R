# Descriptive stability in independently screened LC-MS/MS and MSI; not equivalence testing.
base <- '/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
in_dir <- file.path(base,'subtype_common_specific_2B_tuned225/20260914_173710')
out <- file.path(base,'stable_common_2D/20260914')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rd <- function(n) read.csv(file.path(in_dir,n),check.names=FALSE)
wr <- function(x,n) write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='')
ms <- rd('MSI_sample_linear.csv'); det <- rd('MSI_sample_detection.csv'); ct <- rd('MSI_sample_counts.csv')
lc <- read.csv('/home/gaga/huzhan/MSI_ALL_aligned_11samples_consensus0.1Da/reanalysis2B_bundle/previous_lc_endogenous_entity_representatives.csv',check.names=FALSE)
groups <- c('Cushing','GH','NF-GH','NF-ACTH','NF-GnH')
ms_names <- c('Cushing-1','Cushing-2','GH-1','GH-2','GH-3','NF-GH','NF-ACTH','NF-GnH-1','NF-GnH-2','NF-GnH-3')
lc_names <- c('Cushing-ACTH-1051','Cushing-ACTH-1285','GH-1674','GH-1715','GH-1835','NF-GH-1988','NF-ACTH-1962','NF-ACTH-1974','NF-ACTH-2001','NF-CXX-1982','NF-CXX-1984','NF-CXX-1986')
mg <- rep(groups,c(2,3,1,1,3)); lg <- rep(groups,c(2,3,1,3,3))
presets <- list(strict=c(1.25,.20,1.5),primary=c(1.5,.30,2),relaxed=c(2,.40,3))
stats_signal <- function(v,g,detected) {
  gm <- sapply(groups,function(z) rowMeans(v[,g==z,drop=FALSE]))
  d <- data.frame(max_subtype_mean_fold=apply(gm,1,max)/apply(gm,1,min),sample_CV=apply(v,1,sd)/rowMeans(v),max_sample_fold=apply(v,1,max)/apply(v,1,min),all_samples_detected=detected)
  for (n in names(presets)) {p<-presets[[n]];d[[paste0('stable_',n)]]<-detected & d$max_subtype_mean_fold<=p[1]+1e-12 & d$sample_CV<=p[2]+1e-12 & d$max_sample_fold<=p[3]+1e-12}
  d
}
mv <- as.matrix(ms[,ms_names]); lv <- as.matrix(lc[,lc_names]); floor_value <- min(lv)
md <- apply(det[,ms_names]>=.1 & ct[,ms_names]>=20,1,all)
ld <- apply(lv>floor_value,1,all)
mstat <- cbind(ms[,c('feature_index','mz')],stats_signal(mv,mg,md))
lstat <- cbind(lc[,c('entity_key','ID','Metabolites')],stats_signal(2^lv,lg,ld))
wr(mstat,'MSI_stability_all_274');wr(lstat,'LC_stability_all_411')
wr(mstat[mstat$stable_primary,],'MSI_stable_common_primary');wr(lstat[lstat$stable_primary,],'LC_stable_common_primary')
a <- rd('all_retained_LC_mass_matches.csv')
mi <- match(a$feature_index,mstat$feature_index);li<-match(a$entity_key,lstat$entity_key)
stopifnot(!anyNA(mi),!anyNA(li),nrow(a)==2046)
a<-cbind(a,setNames(mstat[mi,-c(1,2)],paste0('MSI_',names(mstat)[-c(1,2)])),setNames(lstat[li,-c(1,2,3)],paste0('LC_',names(lstat)[-c(1,2,3)])))
core <- a$monoisotopic_core & a$fragmentation_score>0
direct <- a$match_path_category=='direct_same_observed_ion' & a$fragmentation_score>0
res<-list();k<-0L
for(n in names(presets)) for(route in c('LC_first_mass_hits','independent_stable_intersection')) for(e in c('all','core_MS2','direct_MS2')) {
  keep<-a[[paste0('LC_stable_',n)]]
  if(route=='independent_stable_intersection') keep<-keep & a[[paste0('MSI_stable_',n)]]
  if(e=='core_MS2')keep<-keep & core
  if(e=='direct_MS2')keep<-keep & direct
  z<-a[keep,];wr(z,paste(n,route,e,sep='__'))
  k<-k+1L;res[[k]]<-data.frame(preset=n,route=route,evidence=e,LC_stable_entities=sum(lstat[[paste0('stable_',n)]]),MSI_stable_features=sum(mstat[[paste0('stable_',n)]]),pairs=nrow(z),LC_entities=length(unique(z$entity_key)),MSI_features=length(unique(z$feature_index)))
}
res<-do.call(rbind,res);wr(res,'threshold_sensitivity_and_route_counts');wr(a,'all_correspondences_with_stability_audit')
stopifnot(sum(mstat$stable_primary)==38,sum(lstat$stable_primary)==23,sum(a$LC_stable_primary & a$MSI_stable_primary)==17,sum(a$LC_stable_primary & a$MSI_stable_primary & core)==0)
png(file.path(out,'01_stability_screen.png'),width=1600,height=700,res=140)
par(mfrow=c(1,2),mar=c(5,5,4,1))
for(n in c('LC','MSI')) {d<-if(n=='LC') lstat else mstat;plot(d$max_subtype_mean_fold,d$sample_CV*100,log='xy',pch=16,col=ifelse(d$stable_primary,'#009E73','#BBBBBB'),xlab='Maximum / minimum subtype mean',ylab='Across-sample CV (%)',main=paste(n,':',sum(d$stable_primary),'stable candidates'));abline(v=1.5,h=30,lty=2,col='firebrick')}
dev.off()
writeLines(c('Primary: subtype mean ratio <=1.5; sample CV <=0.30; sample max/min <=2.', 'MSI uses mean LINEAR signal across ALL tumor pixels in each of 10 physical samples.','LC uses 2^reported log2 signal in 12 samples. All samples must pass detection convention.', 'Strict sensitivity: 1.25/0.20/1.5; relaxed sensitivity: 2/0.40/3.','Primary LC=23, MSI=38; both-stable mass pairs=17; core+LC MS2 pairs=0.','Relaxed both-stable core+LC MS2 pairs=8 (7 MSI features,7 LC entities).','These are descriptive finite-cohort results, not population equivalence or confirmed MSI identities.'),file.path(out,'README.txt'))
file.copy('/home/gaga/huzhan/MSI_NEW/script/stable_common_tuned225.R',file.path(out,'stable_common_tuned225.R'),overwrite=TRUE)
print(res);cat('\nSTABLE_COMMON_R_VERIFIED\n',out,'\n')
rstudioapi::filesPaneNavigate(out)
