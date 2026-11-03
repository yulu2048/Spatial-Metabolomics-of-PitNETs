# Independent R recomputation of primary MSI results and verification of delivered tables.
# This script uses only the UPDATED MSI matrix; no LC-MS/MS input or old feature indices.
run_lineage_MSI_tuned225 <- function(BASE='/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225',
                                      OUT=file.path(BASE,'lineage_MSI_tuned225/20260914')) {
  suppressPackageStartupMessages(library(Matrix))
  samples<-c('Cushing-1','Cushing-2','GH-1','GH-2','GH-3','NF-GH','NF-ACTH')
  contrasts<-list(ACTH=list(a=1:2,b=7L),GH=list(a=3:5,b=6L))
  expected<-c('pixel_by_feature_matrix_aligned_TIC.rds'='d708ccb8c9c01b35ee37322528728926',
              'feature_mz_aligned_filtered.rds'='84393ef4812df59aeba46e0590555690',
              'pixel_metadata_sample_source_xy.rds'='eca300f50d856ce11123142e0b04e7f6',
              'intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225/tumor_pixels_within_sample_consensus_domains.rds'='01486d5ab237994bb1c7f4c0de896ea1')
  stopifnot(all(unname(tools::md5sum(file.path(BASE,names(expected))))==unname(expected)))
  m<-readRDS(file.path(BASE,'pixel_by_feature_matrix_aligned_TIC.rds'))
  mz<-as.numeric(readRDS(file.path(BASE,'feature_mz_aligned_filtered.rds')))
  meta<-as.data.frame(readRDS(file.path(BASE,'pixel_metadata_sample_source_xy.rds')))
  t<-as.data.frame(readRDS(file.path(BASE,names(expected)[4])))
  t<-t[t$sample_id%in%samples,];i<-t$original_pixel_row
  stopifnot(ncol(m)==274,nrow(t)==29656,!anyDuplicated(i),all(t$is_tumor),
            all(meta$sample_id[i]==t$sample_id),all(meta$source_id[i]==t$source_id),
            all(meta$x[i]==t$x),all(meta$y[i]==t$y),identical(colnames(m),sprintf('mz_%.6f',mz)))
  m<-m[i,,drop=FALSE];lin<-cnt<-matrix(NA_real_,7,274);n<-integer(7)
  dir.create(file.path(OUT,'input'),recursive=TRUE,showWarnings=FALSE)
  if(!file.exists(file.path(OUT,'input/tumor_matrix.bin.gz'))) {
    con<-gzfile(file.path(OUT,'input/tumor_matrix.bin.gz'),'wb')
    writeBin(as.double(as.matrix(m)),con,size=8,endian='little');close(con)
  }
  for(k in seq_along(samples)) {
    ii<-which(t$sample_id==samples[k]);n[k]<-length(ii)
    lin[k,]<-Matrix::colMeans(m[ii,,drop=FALSE])*1e6
    cnt[k,]<-Matrix::colSums(m[ii,,drop=FALSE]>0)
  }
  pri<-log2(1+lin);det<-cnt/n;eligible<-det>=.1 & cnt>=20
  eff<-function(v,a,b)colMeans(v[a,,drop=FALSE])-colMeans(v[b,,drop=FALSE])
  res<-lapply(names(contrasts),function(con){
    a<-contrasts[[con]]$a;b<-contrasts[[con]]$b;e<-eff(pri,a,b);sg<-sign(e)
    pair<-do.call(rbind,lapply(a,function(i)do.call(rbind,lapply(b,function(j)(pri[i,]-pri[j,])*sg))))
    concord<-colMeans(pair>1e-10)+.5*colMeans(abs(pair)<=1e-10)
    margin<-apply(pair,2,min);loo<-apply(sapply(a,function(i)eff(pri,setdiff(a,i),b)*sg>0),1,all)
    high<-ifelse(e>=0,colMeans(eligible[a,,drop=FALSE]),colMeans(eligible[b,,drop=FALSE]))
    base<-concord>=.8-1e-9 & loo & high>=2/3-1e-9
    v<-pri[c(a,b),,drop=FALSE];cmb<-combn(nrow(v),length(a))
    perm<-apply(cmb,2,function(z)eff(v,z,setdiff(seq_len(nrow(v)),z)))
    p<-rowMeans(abs(perm)>=abs(e)-1e-10)
    center<-eff(sweep(pri,1,apply(pri,1,median)),a,b)
    data.frame(contrast=con,feature_index=seq_along(mz),mz=mz,effect=e,minimum_pairwise_margin=margin,
               direction_concordance=concord,leave_one_functional_direction=loo,
               high_group_eligible_fraction=high,candidate=base & abs(e)>=.5,
               median_centered_effect=center,p_exact_two_sided=p,FDR=p.adjust(p,'BH'))
  })
  res<-do.call(rbind,res);dir.create(file.path(OUT,'tables'),recursive=TRUE,showWarnings=FALSE)
  write.csv(res,file.path(OUT,'tables/R_independent_primary_results.csv'),row.names=FALSE)
  delivered<-file.path(OUT,'tables/MSI_differential_all.csv')
  if(file.exists(delivered)) {
    z<-read.csv(delivered);j<-match(paste(res$contrast,res$feature_index),paste(z$contrast,z$feature_index))
    stopifnot(!anyNA(j),identical(res$candidate,as.logical(z$candidate[j])))
    for(v in c('mz','effect','minimum_pairwise_margin','direction_concordance','high_group_eligible_fraction','median_centered_effect','p_exact_two_sided','FDR'))
      stopifnot(max(abs(res[[v]]-z[[v]][j]))<1e-8)
    writeLines(c('PASS: fresh R calculation from original updated matrix matches all 548 delivered rows.',
                 'Matched: candidate flags, mz, effect, minimum margin, concordance, eligibility, centering, exact P, BH FDR.',
                 'Cushing 2 vs NF-ACTH 1; GH 3 vs NF-GH 1. No LC data.'),file.path(OUT,'R_verification.txt'))
  }
  print(with(res,table(contrast,candidate)));cat('R_PRIMARY_REANALYSIS_VERIFIED\n');invisible(res)
}
run_lineage_MSI_tuned225()
