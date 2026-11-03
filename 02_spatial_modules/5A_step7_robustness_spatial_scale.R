# Rebuilt 2026-09-14 against tuned225 matrix; source: /home/gaga/huzhan/MSI_ALL_aligned_11samples_consensus0.1Da/spatial_metabolite_modules_10samples/5A_step7_robustness_spatial_scale.R
# 5A Step 7: MSI-only leave-one-sample-out robustness and empirical spatial scales
local({
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(WGCNA);library(dynamicTreeCut);library(matrixStats)})
BASE<-'/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
OUT<-file.path(BASE,'spatial_metabolite_modules_10samples_tuned225')
S7<-file.path(OUT,'step7_robustness_spatial_scale');dir.create(S7,recursive=TRUE,showWarnings=FALSE)
b<-readRDS(file.path(OUT,'step3_modules','MSI_consensus_modules.rds'));reps<-as.data.table(b$representatives);fi<-reps$feature_index;primary<-reps$module;mods<-setdiff(sort(unique(primary)),'grey')
files<-list.files(file.path(OUT,'step2_spatial_consensus','sample_correlations'),pattern='rds$',full.names=TRUE);ck<-lapply(files,readRDS);sn<-vapply(ck,function(z)z$sample,character(1));idx<-match(fi,ck[[1]]$feature_index);R<-lapply(ck,function(z)z$spearman[idx,idx,drop=FALSE])
loso_rows<-list();partitions<-list()
for(o in seq_along(R)){
  rr<-R[-o];flat<-do.call(cbind,lapply(rr,as.vector));nev<-rowSums(is.finite(flat));med<-matrix(matrixStats::rowMedians(flat,na.rm=TRUE),length(fi));pos<-matrix(rowMeans(flat>0,na.rm=TRUE),length(fi));adj<-pmax(med,0)^4;adj[nev<7|pos<.70]<-0;diag(adj)<-1
  tom<-WGCNA::TOMsimilarity(adj,TOMType='unsigned',verbose=0);hc<-hclust(as.dist(1-tom),method='average');cl<-dynamicTreeCut::cutreeDynamic(hc,distM=1-tom,deepSplit=2,pamRespectsDendro=FALSE,minClusterSize=5,verbose=0);partitions[[sn[o]]]<-cl
  for(m in mods){A<-which(primary==m);lev<-setdiff(unique(cl),0);js<-if(length(lev))vapply(lev,function(g){B<-which(cl==g);length(intersect(A,B))/length(union(A,B))},numeric(1)) else NA_real_;k<-if(all(!is.finite(js)))NA_integer_ else which.max(js);loso_rows[[length(loso_rows)+1]]<-data.table(module=m,omitted_sample=sn[o],best_jaccard=if(is.na(k))NA_real_ else js[k],matched_cluster=if(is.na(k))NA_integer_ else lev[k])
  }
  rm(flat,med,pos,adj,tom,hc);gc()
}
loso<-rbindlist(loso_rows);sumlo<-loso[,.(median_LOSO_Jaccard=median(best_jaccard,na.rm=TRUE),min_LOSO_Jaccard=min(best_jaccard,na.rm=TRUE),n_LOSO_ge_0.5=sum(best_jaccard>=.5,na.rm=TRUE)),by=module][order(-median_LOSO_Jaccard)]
fwrite(loso,file.path(S7,'leave_one_sample_out_module_Jaccard.csv'));fwrite(sumlo,file.path(S7,'leave_one_sample_out_module_summary.csv'))
# Empirical semivariograms from frozen unsmoothed pixel module scores, aggregated only for scale estimation.
px<-readRDS(file.path(OUT,'step4_module_activity_spatial_organization','pixel_module_activity_scores.rds'));md<-as.data.table(px$pixel_metadata);sc<-px$scores;vgrows<-list();ranges<-list()
for(si in unique(md$sample_id)){
 ii<-which(md$sample_id==si);d<-copy(md[ii]);d[,bx:=floor((x-min(x))/10)];d[,by:=floor((y-min(y))/10)];key<-paste(d$bx,d$by,sep='_')
 for(m in mods){v<-sc[ii,m];ag<-data.table(block_id=key,bx=d$bx,by=d$by,v=v)[is.finite(v),.(x=mean(bx)*10,y=mean(by)*10,z=mean(v)),by=block_id];if(nrow(ag)<15)next;D<-as.matrix(dist(ag[,.(x,y)]));G<-outer(ag$z,ag$z,function(a,b)(a-b)^2/2);up<-upper.tri(D);dv<-D[up];gv<-G[up];cutmax<-quantile(dv,.8,na.rm=TRUE);br<-seq(0,cutmax,length.out=13);bin<-cut(dv,br,include.lowest=TRUE,labels=FALSE);vv<-data.table(distance=dv,gamma=gv,bin=bin)[!is.na(bin),.(distance=mean(distance),semivariance=mean(gamma),n_pairs=.N),by=bin];vv[,c('sample_id','module'):=list(si,m)];vgrows[[length(vgrows)+1]]<-vv;sill<-mean(tail(vv$semivariance,3));hit<-which(vv$semivariance>=.95*sill);er<-if(length(hit))vv$distance[min(hit)] else max(vv$distance);ranges[[length(ranges)+1]]<-data.table(sample_id=si,module=m,effective_range_coordinate_units=er,empirical_sill=sill,n_blocks=nrow(ag))
 }
}
vg<-rbindlist(vgrows);rg<-rbindlist(ranges);rgs<-rg[,.(median_effective_range=median(effective_range_coordinate_units,na.rm=TRUE),min_effective_range=min(effective_range_coordinate_units,na.rm=TRUE),max_effective_range=max(effective_range_coordinate_units,na.rm=TRUE),n_samples=.N),by=module][order(-median_effective_range)];fwrite(vg,file.path(S7,'empirical_semivariograms.csv'));fwrite(rg,file.path(S7,'per_sample_module_effective_range.csv'));fwrite(rgs,file.path(S7,'module_spatial_scale_summary.csv'))
p1<-ggplot(loso,aes(x=module,y=best_jaccard))+geom_boxplot(fill='#66C2A5')+geom_hline(yintercept=.5,lty=2)+coord_cartesian(ylim=c(0,1))+labs(x=NULL,y='Best Jaccard after omitting one sample',title='MSI module leave-one-sample-out robustness')+theme_bw(base_size=11)
p2<-ggplot(rg,aes(x=reorder(module,effective_range_coordinate_units,median),y=effective_range_coordinate_units))+geom_boxplot(fill='#E78AC3')+coord_flip()+labs(x=NULL,y='Effective range (MSI coordinate units)',title='Empirical spatial scale of frozen MSI module activity')+theme_bw(base_size=11)
ggsave(file.path(S7,'step7_robustness_and_spatial_scale.pdf'),p1/p2,width=11,height=10)
saveRDS(list(loso=loso,loso_summary=sumlo,partitions=partitions,semivariograms=vg,per_sample_ranges=rg,range_summary=rgs,parameters=list(block_width=10,variogram_max_distance_quantile=.8,range_threshold=.95)),file.path(S7,'MSI_module_robustness_and_spatial_scale.rds'))
writeLines(c('All analyses are MSI-only.','LOSO repeats consensus network and dynamic tree cutting after omitting each physical sample; ion-family representatives remain frozen.','Effective range is the first empirical distance bin reaching 95% of the mean semivariance in the last three bins; units are MSI coordinate units and should not be interpreted as micrometers without acquisition calibration.'),file.path(S7,'README.txt'));file.create(file.path(S7,'STEP7_COMPLETE.ok'));cat('STEP7_COMPLETE');print(sumlo);print(rgs)
})
