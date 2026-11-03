# MSI-only subtype enrichment. Shared62 is a frozen user-defined list, historically cross-platform selected.
suppressPackageStartupMessages(library(data.table)); requireNamespace('Matrix')
BASE <- '/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
OUT <- file.path(BASE,'MSI_only_specific_shared62/20260914');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
savecsv <- function(x,n)fwrite(x,file.path(OUT,paste0(n,'.csv')))
groups <- c('Cushing','GH','NF-GH','NF-ACTH','NF-GnH')
samples <- c('Cushing-1','Cushing-2','GH-1','GH-2','GH-3','NF-GH','NF-ACTH','NF-GnH-1','NF-GnH-2','NF-GnH-3')
labs <- sub('-[0-9]+$','',samples)
  groupmean<-function(x,labels) do.call(rbind,lapply(groups,function(g) colMeans(x[labels==g,,drop=FALSE],na.rm=TRUE)))
  auc<-function(x,idx) vapply(seq_len(ncol(x)),function(j){a<-x[idx,j];b<-x[-idx,j];a<-a[is.finite(a)];b<-b[is.finite(b)];if(!length(a)||!length(b))return(NA_real_);mean(outer(a,b,'>'))+.5*mean(outer(a,b,'=='))},numeric(1))
  # Exact one-sided sample-label combination test for target mean vs all other samples.
  # NA positive-signal summaries use available samples only; missing signal is never set to zero.
  exactp<-function(x,idx) {
    vapply(seq_len(ncol(x)),function(j) {
      good<-which(is.finite(x[,j])); a<-intersect(idx,good);b<-setdiff(good,idx)
      if(!length(a)||!length(b))return(NA_real_)
      v<-x[good,j]; k<-length(a);obs<-mean(x[a,j])-mean(x[b,j]); cb<-combn(length(v),k)
      dif<-apply(cb,2,function(z)mean(v[z])-mean(v[-z]))
      mean(dif>=obs-1e-12)
    },numeric(1))
  }
  mat<-readRDS(file.path(BASE,'pixel_by_feature_matrix_aligned_TIC.rds'))
  mz<-as.numeric(readRDS(file.path(BASE,'feature_mz_aligned_filtered.rds')))
  meta<-as.data.table(readRDS(file.path(BASE,'pixel_metadata_sample_source_xy.rds')))
  cl<-as.data.table(readRDS(file.path(BASE,'QC/tumor_pixel_selection_0B_k4/pixel_clusters_12sources.rds')))
  stopifnot(nrow(mat)==nrow(meta),ncol(mat)==length(mz),identical(as.integer(meta$pixel_row),seq_len(nrow(mat))),!anyDuplicated(cl$pixel_row),identical(colnames(mat),sprintf('mz_%.6f',mz)))
  ci<-match(meta$pixel_row,cl$pixel_row);stopifnot(!anyNA(ci),all(meta$sample_id==cl$sample_id[ci]),all(meta$source_id==cl$source_id[ci]),all(meta$x==cl$x[ci]),all(meta$y==cl$y[ci]))
  tumor_file <- file.path(BASE,'intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225/tumor_pixels_within_sample_consensus_domains.rds')
  tumor <- as.data.table(readRDS(tumor_file))
  stopifnot(!anyDuplicated(tumor$original_pixel_row),all(tumor$is_tumor),setequal(as.character(tumor$sample_id),samples))
  ti <- match(tumor$original_pixel_row,meta$pixel_row)
  stopifnot(!anyNA(ti),all(meta$sample_id[ti]==tumor$sample_id),all(meta$source_id[ti]==tumor$source_id),all(meta$x[ti]==tumor$x),all(meta$y[ti]==tumor$y))
  keep <- seq_len(nrow(meta)) %in% ti
  savecsv(tumor[,.(original_pixel_row,sample_id,source_id,x,y,is_tumor)],'tumor_pixel_selection_used')
  pix<-meta[keep]; x<-mat[which(keep),,drop=FALSE];rm(mat);gc()
  stopifnot(all(x@x>=0),all(is.finite(x@x)),setequal(as.character(unique(pix$sample_id)),samples))
  f<-ncol(x);nf<-seq_len(f); sums<-counts<-oldmean<-matrix(0,10,f);np<-integer(10)
  block_sums<-block_counts<-block_n<-vector('list',10)
  for(i in seq_along(samples)) {
    rows<-which(pix$sample_id==samples[i]);xx<-x[rows,,drop=FALSE];np[i]<-length(rows)
    sums[i,]<-Matrix::colSums(xx)*1e6;counts[i,]<-Matrix::colSums(xx>0)
    zz<-xx;zz@x<-log1p(zz@x*1e6);oldmean[i,]<-Matrix::colMeans(zz)
    pp<-pix[rows]; bid<-factor(paste(pp$source_id,floor(pp$x/10),floor(pp$y/10),sep=':'))
    a<-Matrix::sparseMatrix(i=as.integer(bid),j=seq_along(rows),x=1,dims=c(nlevels(bid),length(rows)))
    block_sums[[i]]<-as.matrix(a%*%xx)*1e6;block_counts[[i]]<-as.matrix(a%*%(xx>0));block_n[[i]]<-as.numeric(table(bid))
  }
  linear<-sums/np; detection<-counts/np; primary<-log2(1+linear)
  positive_linear<-sums/counts;positive_linear[counts<20]<-NA_real_; positive<-log2(1+positive_linear)
  stopifnot(max(abs(linear - detection * ifelse(counts > 0, sums / pmax(counts,1), 0))) < 1e-8)
  for(a in c('linear','detection','primary','positive_linear','positive','oldmean','counts')) {v<-get(a);colnames(v)<-as.character(nf);rownames(v)<-samples;assign(a,v)}
  savecsv(data.table(sample_id=samples,class=labs,n_tumor_pixels=np,n_blocks=lengths(block_n)),'MSI_samples')
  for(a in c('linear','detection','primary','positive_linear','positive','oldmean','counts'))savecsv(data.table(feature_index=nf,mz=mz,as.data.table(t(get(a)))) ,paste0('MSI_sample_',a))
  saveRDS(list(primary=primary,detection=detection,positive=positive,linear=linear,positive_linear=positive_linear,oldmean=oldmean,counts=counts,samples=samples,groups=labs,mz=mz),file.path(OUT,'sample_summaries.rds'))
  eligible<-detection>=.1 & counts>=20
  commonpass<-sapply(groups,function(g)colSums(eligible[labs==g,,drop=FALSE])>=ceiling(sum(labs==g)*2/3))
  scommon<-data.table(feature_index=nf,mz=mz,stable_all5=apply(commonpass,1,all),strict_all10=colSums(eligible)==10)
  savecsv(scommon,'MSI_common_all');savecsv(scommon[stable_all5==TRUE],'MSI_common')
  statsfun<-function(v,metric,elig=eligible) {
    gm<-groupmean(v,labs); centered<-sweep(v,1,rowMeans(v,na.rm=TRUE));medcenter<-sweep(v,1,apply(v,1,median,na.rm=TRUE))
    cg<-groupmean(centered,labs);mg<-groupmean(medcenter,labs)
    rbindlist(lapply(seq_along(groups),function(gi){
      g<-groups[gi];idx<-which(labs==g);rest<-setdiff(1:10,idx)
      diff<-colMeans(v[idx,,drop=FALSE],na.rm=TRUE)-colMeans(v[rest,,drop=FALSE],na.rm=TRUE)
      ge<-gm[gi,]-colMeans(gm[-gi,,drop=FALSE]); margin<-gm[gi,]-apply(gm[-gi,,drop=FALSE],2,max)
      sdall<-apply(v,2,sd,na.rm=TRUE); av<-auc(v,idx); p<-exactp(v,idx)
      dr<-colMeans(sweep(v[idx,,drop=FALSE],2,apply(v[rest,,drop=FALSE],2,median,na.rm=TRUE),'>'),na.rm=TRUE)
      ep<-colMeans(elig[idx,,drop=FALSE]);cf<-cg[gi,]-apply(cg[-gi,,drop=FALSE],2,max)>0 & auc(centered,idx)>=.8
      mf<-mg[gi,]-apply(mg[-gi,,drop=FALSE],2,max)>0 & auc(medcenter,idx)>=.8
      lo<-if(length(idx)==1)rep(NA_real_,f) else Reduce('+',lapply(idx,function(h)colMeans(v[setdiff(idx,h),,drop=FALSE],na.rm=TRUE)>apply(gm[-gi,,drop=FALSE],2,max)))/length(idx)
      evalok<-apply(sapply(groups,function(gg)colSums(is.finite(v[labs==gg,,drop=FALSE]))>=ceiling(sum(labs==gg)*2/3)),1,all)
      th<-if(metric=='detection') .1 else .15
      base<-evalok & ep>=(if(length(idx)>1)2/3 else 1)
      base<-base & diff>=th & diff/pmax(sdall,1e-8)>=.8 & margin>0 & av>=.8 & dr>=(if(length(idx)>1)2/3 else 1)
      base<-base & (is.na(lo)|lo>=2/3);base[is.na(base)]<-FALSE
      data.table(metric=metric,class=g,feature_index=nf,mz=mz,n_target=length(idx),target_mean=gm[gi,],effect=diff,equal_group_effect=ge,minimum_pairwise_difference=margin,standardized_effect=diff/pmax(sdall,1e-8),auc=av,direction_fraction=dr,target_eligible_fraction=ep,leave_one_out=lo,mean_center_robust=cf,median_center_robust=mf,evaluable=evalok,p=p,FDR=p.adjust(p,'BH'),base_candidate=base)
    }))
  }
  message('Computing independent MSI statistics: primary, detection, positive intensity, previous metric...')
  st<-rbindlist(list(statsfun(primary,'primary'),statsfun(detection,'detection'),statsfun(positive,'positive',counts>=20)))
  st[,bootstrap_support:=NA_real_]
  # Source-specific 10x10 blocks; resampling never splits a physical sample into independent subjects.
  set.seed(20260914);B<-200L
  candcols<-sort(unique(st[base_candidate==TRUE,feature_index])); nb<-length(candcols)
  wins<-matrix(0,nrow(st),1)
  if(nb) for(b in seq_len(B)) {
    bs<-bc<-matrix(0,10,nb);bn<-numeric(10)
    for(i in 1:10){k<-length(block_n[[i]]);take<-sample.int(k,k,replace=TRUE);bs[i,]<-colSums(block_sums[[i]][take,candcols,drop=FALSE]);bc[i,]<-colSums(block_counts[[i]][take,candcols,drop=FALSE]);bn[i]<-sum(block_n[[i]][take])}
    vals<-list(primary=log2(1+bs/bn),detection=bc/bn,positive=log2(1+bs/bc))
    vals$positive[bc<20]<-NA_real_
    for(metric0 in names(vals)) {
      gm<-groupmean(vals[[metric0]],labs)
      for(gi in seq_along(groups)) {ids<-which(st$metric==metric0 & st$class==groups[gi] & st$base_candidate);if(!length(ids))next;j<-match(st$feature_index[ids],candcols);ok<-gm[gi,j]>apply(gm[-gi,j,drop=FALSE],2,max);ok[is.na(ok)]<-FALSE;wins[ids]<-wins[ids]+ok}
    }
    if(b%%50==0)message('Spatial bootstrap ',b,'/',B)
  }
  st[metric!='previous_meanlog' & base_candidate==TRUE,bootstrap_support:=as.numeric(wins[which(st$metric!='previous_meanlog' & st$base_candidate)])/B]
  st[,candidate:=base_candidate & bootstrap_support>=.8 & (mean_center_robust | median_center_robust)]
  st[is.na(candidate),candidate:=FALSE]
  st[,tier:=ifelse(candidate,ifelse(n_target==1,'singleton_exploratory','candidate'),'not_candidate')]
  savecsv(st,'MSI_statistics_all_metrics');savecsv(st[candidate==TRUE],'MSI_candidates_all_metrics')
  for(mm in unique(st$metric))savecsv(st[metric==mm & candidate==TRUE],paste0('MSI_candidates_',mm))
  # Primary candidate explanations use both supplements without requiring that all three pass.
  expl<-st[metric=='primary' & candidate==TRUE]
  for(mm in c('detection','positive')) {t<-st[metric==mm,.(class,feature_index,equal_group_effect,candidate,evaluable)];setnames(t,c('equal_group_effect','candidate','evaluable'),paste0(mm,c('_effect','_candidate','_evaluable')));expl<-merge(expl,t,by=c('class','feature_index'),all.x=TRUE)}
  expl[,support_pattern:=fifelse(!positive_evaluable,'positive_intensity_not_evaluable',fifelse(detection_candidate & positive_candidate,'both_supplements_support',fifelse(detection_candidate,'detection_supported',fifelse(positive_candidate,'positive_intensity_supported','primary_only_or_subthreshold'))))]
  savecsv(expl,'MSI_primary_candidate_explanations')

shared_ids <- c(93,134,63,121,32,128,141,83,143,65,37,48,54,88,97,31,75,78,89,100,108,136,19,21,30,66,74,76,90,94,96,99,103,104,110,135,12,23,24,28,38,64,105,109,122,123,129,130,133,142,173,174,175,176,205,206,207,208,227,246,261,262)
shared<-data.table(feature_index=shared_ids,mz=mz[shared_ids],set='frozen_shared62',provenance='historical_cross_platform_selection');savecsv(shared,'fixed_shared62')
expl[,in_fixed_shared62:=feature_index %in% shared_ids]
expl[,MSI_FDR_lt_005:=FDR<.05]
savecsv(expl,'MSI_specific_candidates_with_shared_overlap')
summary<-rbindlist(lapply(groups,function(g){z<-expl[class==g];data.table(class=g,n_samples=sum(labs==g),specific_candidates=nrow(z),overlap_shared62=sum(z$in_fixed_shared62),outside_shared62=sum(!z$in_fixed_shared62),FDR_lt_005=sum(z$MSI_FDR_lt_005),detection_supported=sum(z$detection_candidate),positive_supported=sum(z$positive_candidate))}))
savecsv(summary,'MSI_specific_summary')
for(g in groups)savecsv(expl[class==g],paste0('specific_',g))
prev<-fread(file.path(BASE,'subtype_common_specific_2B_tuned225/20260914_173710/MSI_primary_candidate_explanations.csv'))
stopifnot(setequal(paste(expl$class,expl$feature_index),paste(prev$class,prev$feature_index)))
writeLines(c('MSI-only recomputation from full tumor-pixel matrix; no LC input used in subtype selection.','The fixed shared62 list has historical LC-conditioned provenance; it is not a de novo MSI-independent common set.','Specific means subtype-enriched candidate, not exclusive detection. Shared and enriched sets may overlap.','Singleton NF-GH/NF-ACTH are exploratory, not population-validated markers.','Sample-label exact permutation p, BH within each subtype/metric; FDR is reported but not an inclusion gate.','Spatial block bootstrap measures spatial sensitivity, not independent biological replication.','Primary uses log2(1+all-tumor-pixel linear mean *1e6), supplements detection and positive-pixel means.','Primary effect>=0.15, standardized effect>=0.8, AUC>=0.8, highest subtype, target eligibility/direction>=2/3, leave-one-out>=2/3 when replicated, bootstrap>=0.8, mean OR median centering robustness.','Independent MSI criteria match previous 2B; previous LC-dependent subsets are not used.'),file.path(OUT,'METHODS.txt'))
file.copy('/home/gaga/huzhan/MSI_NEW/script/MSI_only_specific_shared62.R',file.path(OUT,'script_executed.R'),overwrite=TRUE)
print(summary);cat('MSI_ONLY_SPECIFIC_RECOMPUTE_VERIFIED\n');rstudioapi::filesPaneNavigate(OUT)
