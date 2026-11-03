## 2B: Independent LC and MSI discovery, followed by frozen 0.1 Da mass matching.
## Primary MSI: log2(1 + mean(TIC-normalized signal over ALL tumor pixels)*1e6).
## Supplements: detection fraction; log2(1 + mean positive TIC signal*1e6).
## Original mean(pixel log1p CPM) is a sensitivity analysis, not the primary endpoint.
## The biological unit is a physical sample; spatial bootstrap measures technical/spatial stability only.
.script2B <- '/home/gaga/huzhan/MSI_NEW/script/2B_subtype_common_specific_tuned225.R'
run_subtype_tuned225 <- function() {
  suppressPackageStartupMessages(library(data.table)); requireNamespace('Matrix')
  OLD <- '/home/gaga/huzhan/MSI_ALL_aligned_11samples_consensus0.1Da'
  BASE <- '/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
  INPUT <- file.path(OLD,'reanalysis2B_bundle')
  OUT <- file.path(BASE,'subtype_common_specific_2B_tuned225',format(Sys.time(),'%Y%m%d_%H%M%S'))
  dir.create(OUT,recursive=TRUE)
  file.copy(.script2B,file.path(OUT,'script_executed.R'))
  savecsv <- function(x,n) fwrite(x,file.path(OUT,paste0(n,'.csv')))
  groups <- c('Cushing','GH','NF-GH','NF-ACTH','NF-GnH')
  samples <- c('Cushing-1','Cushing-2','GH-1','GH-2','GH-3','NF-GH','NF-ACTH','NF-GnH-1','NF-GnH-2','NF-GnH-3')
  labs <- sub('-[0-9]+$','',samples)
  lc_cols <- c('Cushing-ACTH-1051','Cushing-ACTH-1285','GH-1674','GH-1715','GH-1835','NF-ACTH-1962','NF-ACTH-1974','NF-ACTH-2001','NF-CXX-1982','NF-CXX-1984','NF-CXX-1986','NF-GH-1988')
  lc_labs <- c(rep('Cushing',2),rep('GH',3),rep('NF-ACTH',3),rep('NF-GnH',3),'NF-GH')
  reps <- fread(file.path(INPUT,'previous_lc_endogenous_entity_representatives.csv'))
  lc <- fread(file.path(OLD,'lc_ms_ms.csv'),select=c('ID','m/z','Metabolites',lc_cols))
  lc <- lc[!is.na(ID)&nzchar(ID)]
  stopifnot(nrow(reps)==411,!anyDuplicated(lc$ID),!anyDuplicated(reps$entity_key))
  ix<-match(reps$ID,lc$ID);stopifnot(!anyNA(ix),all(lc$Metabolites[ix]==reps$Metabolites),all(abs(lc[['m/z']][ix]-reps[['m/z']])<1e-6))
  lx<-t(as.matrix(lc[ix,..lc_cols]));storage.mode(lx)<-'double'
  stopifnot(all(is.finite(lx)),max(abs(lx-t(as.matrix(reps[,..lc_cols]))))<1e-6)
  floor_lc<-min(lx); lp<-lx>floor_lc
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
  lg<-groupmean(lx,lc_labs); lpres<-groupmean(lp,lc_labs)
  lcstats<-rbindlist(lapply(seq_along(groups),function(gi){
    g<-groups[gi];idx<-which(lc_labs==g);rest<-colMeans(lg[-gi,,drop=FALSE]);eff<-lg[gi,]-rest; margin<-lg[gi,]-apply(lg[-gi,,drop=FALSE],2,max)
    av<-auc(lx,idx);pv<-exactp(lx,idx);dr<-colMeans(sweep(lx[idx,,drop=FALSE],2,rest,'>'));det<-lpres[gi,]
    if(length(idx)==1) tier<-ifelse(eff>=.25 & margin>0 & av>=.8 & det>=1,'singleton_candidate','not_candidate') else tier<-ifelse(eff>=.25 & margin>0 & av>=.8 & det>=2/3 & dr>=2/3,'tier_A_candidate',ifelse(eff>=.15 & margin>0 & av>=.7 & det>=2/3 & dr>=2/3,'tier_B_candidate','not_candidate'))
    data.table(class=g,entity_key=reps$entity_key,lc_ID=reps$ID,lc_row=reps$lc_row,metabolite=reps$Metabolites,fragmentation_score=reps[['Fragmentation Score']],n_target=length(idx),effect=eff,minimum_pairwise_difference=margin,auc=av,direction_fraction=dr,nonimputed_fraction=det,p=pv,FDR=p.adjust(pv,'BH'),tier=tier)
  }))
  lcc<-lcstats[tier!='not_candidate'];lcommon<-reps[apply(lpres>=2/3-1e-9,2,all)]
  priorlc<-fread(file.path(INPUT,'previous_lc_class_specific_candidates.csv'))
  stopifnot(setequal(paste(lcc$class,lcc$entity_key),paste(priorlc$class,priorlc$entity_key)))
  savecsv(lcstats,'LC_statistics_all');savecsv(lcc,'LC_candidates');savecsv(lcommon,'LC_common')
  message('LC independently recomputed and agrees with previous candidate set: ',nrow(lcc),' candidates. Loading MSI TIC matrix...')
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
  # Recompute the complete previous ion-path search against NEW m/z values.
  proton<-1.007276466621;electron<-0.000548579909065;water<-18.01056468403
  fa<-46.00547930360;ac<-60.02112936806;chloride<-34.968852682+electron
  na<-22.9897692820-electron;potassium<-38.9637064864-electron
  ladd<-data.table(adduct=c('M-H','M+FA-H','M+Cl','M-H2O-H','2M-H','M-'),mult=c(1,1,1,1,2,1),offset=c(-proton,fa-proton,chloride,-water-proton,-proton,electron),canonical=c('[M-H]-','[M+HCOO]-','[M+Cl]-','[M-H2O-H]-','[2M-H]-','[M]-'))
  ions<-data.table(ion=c('[M-H]-','[M+Na-2H]-','[M+K-2H]-','[M+Cl]-','[M+HCOO]-','[M+CH3COO]-','[M-H2O-H]-','[2M-H]-','[M]-'),mult=c(rep(1,7),2,1),offset=c(-proton,na-2*proton,potassium-2*proton,chloride,fa-proton,ac-proton,-water-proton,-proton,electron),family=c(rep('core_adduct',6),'neutral_loss','dimer','radical_ion'))
  paths<-list(); unsupported<-reps[!Adducts%in%ladd$adduct];savecsv(unsupported,'LC_unsupported_adducts')
  for(i in seq_len(nrow(reps))) {
    r<-reps[i];la<-ladd[adduct==r$Adducts];if(!nrow(la))next
    neutral<-(r[['m/z']]-la$offset)/la$mult
    for(j in seq_len(nrow(ions)))for(iso in 0:2) {
      io<-ions[j];expected<-io$mult*neutral+io$offset+iso*1.00335483507
      ids<-which(abs(mz-expected)<=.1+1e-12);if(!length(ids)||expected<=0)next
      direct<-iso==0 && io$ion==la$canonical
      cat0<-if(iso>0)'isotope_exploratory' else if(direct)'direct_same_observed_ion' else if(io$family=='core_adduct')'neutral_mass_alternative_adduct' else paste0('neutral_mass_',io$family,'_exploratory')
      pr<-if(iso>0)3L else if(direct)0L else if(io$family=='core_adduct')1L else 2L
      er<-abs(mz[ids]-expected)
      paths[[length(paths)+1L]]<-data.table(feature_index=ids,msi_mz=mz[ids],lc_row=r$lc_row,entity_key=r$entity_key,metabolite=r$Metabolites,formula=r$Formula,lc_ID=r$ID,lc_mz=r[['m/z']],retention_time_min=r[['Retention time (min)']],declared_lc_adduct=r$Adducts,declared_lc_adduct_canonical=la$canonical,neutral_mass_from_lc_observed=neutral,predicted_msi_ion=io$ion,predicted_ion_family=io$family,predicted_isotope_number=iso,predicted_msi_mz=expected,msi_minus_predicted_da=mz[ids]-expected,absolute_mass_difference_da=er,mass_difference_ppm=er/expected*1e6,absolute_da_tier=ifelse(er<=.01,'A_<=0.01Da',ifelse(er<=.05,'B_0.01-0.05Da','C_0.05-0.10Da')),match_path_category=cat0,path_rank=pr,monoisotopic_core=iso==0 & io$family=='core_adduct',fragmentation_score=r[['Fragmentation Score']],MS2_fragment_supported=r[['Fragmentation Score']]>0,annotation_mass_error_ppm=r[['Mass Error (ppm)']],HMDB=r$HMDB,PubChem=r$PubChem,KEGG=r$KEGG,absolute_da_tier_rank=ifelse(er<=.01,0L,ifelse(er<=.05,1L,2L)))
    }
  }
  allpaths<-rbindlist(paths);setorder(allpaths,feature_index,lc_row,path_rank,absolute_da_tier_rank,absolute_mass_difference_da,predicted_isotope_number)
  allpaths[,best_path_for_feature_lc_entity:=!duplicated(paste(feature_index,entity_key))]
  savecsv(allpaths,'all_mass_match_paths_new_matrix');matches<-allpaths[best_path_for_feature_lc_entity==TRUE]
  stopifnot(!anyDuplicated(matches[,.(entity_key,feature_index)]),all(matches$lc_ID%in%reps$ID),all(abs(mz[matches$feature_index]-matches$msi_mz)<1e-6),all(matches$absolute_mass_difference_da<=.1+1e-9))
  stopifnot(all(abs(matches$lc_mz-reps[['m/z']][match(matches$lc_ID,reps$ID)])<1e-6))
  savecsv(matches,'all_retained_LC_mass_matches')
  cross<-merge(matches,lcc[,.(class,entity_key,LC_effect=effect,LC_auc=auc,LC_FDR=FDR)],by='entity_key',allow.cartesian=TRUE)
  countsout<-list();byclass<-list();addcount<-function(name,d) {
    countsout[[length(countsout)+1]]<<-data.table(result_set=name,pairs=nrow(d),LC_entities=uniqueN(d$entity_key),MSI_features=uniqueN(d$feature_index),LC_with_MS2=uniqueN(d[fragmentation_score>0,entity_key]))
    if('class'%in%names(d))byclass[[length(byclass)+1]]<<-rbindlist(lapply(groups,function(g)d[class==g,.(result_set=name,class=g,pairs=.N,LC_entities=uniqueN(entity_key),MSI_features=uniqueN(feature_index))]))
  }
  layers<-function(prefix,d) {
    addcount(prefix,d)
    addcount(paste0(prefix,'__core'),d[monoisotopic_core==TRUE])
    addcount(paste0(prefix,'__core_MS2'),d[monoisotopic_core==TRUE & fragmentation_score>0])
    addcount(paste0(prefix,'__direct_MS2'),d[match_path_category=='direct_same_observed_ion' & fragmentation_score>0])
  }
  layers('all_retained_mass_only',matches);layers('LC_candidates_mass_only',cross)
  for(mm in c('primary','detection','positive')) {
    v<-switch(mm,primary=primary,detection=detection,positive=positive,previous_meanlog=oldmean)
    gm<-groupmean(v,labs)
    ss<-st[metric==mm,.(class,feature_index,MSI_effect=equal_group_effect,MSI_margin=minimum_pairwise_difference,MSI_auc=auc,MSI_FDR=FDR,MSI_candidate=candidate,MSI_evaluable=evaluable)]
    z<-merge(cross,ss,by=c('class','feature_index'))
    z[,profile_spearman:=vapply(seq_len(.N),function(i){a<-lg[,match(entity_key[i],reps$entity_key)];bb<-gm[,feature_index[i]];if(any(!is.finite(bb))||sd(bb)==0)return(NA_real_);suppressWarnings(cor(a,bb,method='spearman'))},numeric(1))]
    z[,direction_consistent:=MSI_evaluable & MSI_margin>0 & is.finite(profile_spearman) & profile_spearman>=-1e-12]
    savecsv(z,paste0('LC_to_MSI_',mm,'_all_mass_matches'))
    dr<-z[direction_consistent==TRUE];inter<-z[MSI_candidate==TRUE]
    layers(paste0(mm,'__LCfirst_direction'),dr);layers(paste0(mm,'__independent_intersection'),inter)
    addcount(paste0(mm,'__LCfirst_A'),dr[monoisotopic_core==TRUE & fragmentation_score>0 & absolute_mass_difference_da<=.05 & profile_spearman>=.5-1e-12])
    savecsv(dr,paste0('LCfirst_direction_',mm));savecsv(inter,paste0('independent_intersection_',mm))
    savecsv(inter[monoisotopic_core==TRUE & fragmentation_score>0],paste0('independent_core_MS2_',mm))
    savecsv(inter[match_path_category=='direct_same_observed_ion' & fragmentation_score>0],paste0('independent_direct_MS2_',mm))
  }
  priority<-merge(fread(file.path(OUT,'independent_core_MS2_primary.csv')),expl[,.(class,feature_index,detection_effect,detection_candidate,positive_effect,positive_candidate,positive_evaluable,support_pattern)],by=c('class','feature_index'),all.x=TRUE)
  savecsv(priority,'priority_primary_core_MS2_with_supplements')
  savecsv(priority[support_pattern=='both_supplements_support'],'priority_primary_core_MS2_both_supplements')
  common_lcfirst<-matches[entity_key%in%lcommon$entity_key]
  common_lcfirst[,MSI_common_independent:=feature_index%in%scommon[stable_all5==TRUE,feature_index]]
  layers('common_LCfirst_mass_only',common_lcfirst);savecsv(common_lcfirst,'common_LCfirst_mass_only')
  commoncross<-matches[entity_key%in%lcommon$entity_key & feature_index%in%scommon[stable_all5==TRUE,feature_index]]
  layers('common_independent_intersection',commoncross);savecsv(commoncross,'common_independent_intersection')
  savecsv(rbindlist(countsout),'matching_counts');savecsv(rbindlist(byclass),'matching_counts_by_class')
  cs<-rbindlist(lapply(groups,function(g) {
    rbindlist(list(data.table(platform='LC',metric='log2_reported',class=g,n_features=uniqueN(lcc[class==g,entity_key]),n_MS2=uniqueN(lcc[class==g & fragmentation_score>0,entity_key]),n_FDR05=lcstats[class==g,sum(FDR<.05,na.rm=TRUE)]),st[class==g,.(platform='MSI',n_features=sum(candidate),n_MS2=NA_integer_,n_FDR05=sum(FDR<.05,na.rm=TRUE)),by=.(metric,class)]),fill=TRUE)
  }))
  savecsv(cs,'independent_candidate_counts')
  savecsv(expl[,.N,by=.(class,support_pattern)],'primary_support_pattern_counts')
  writeLines(c(
    'New tuned225 aligned matrix; actual feature count is read from data. Original 2B statistical thresholds retained.',
    'Tumor mask: ALL is_tumor rows from completed tuned225 1A output; includes low-confidence domain pixels; no domain/core subset.',
    'Primary: log2(1 + mean saved TIC-normalized signal across ALL tumor pixels *1e6). No renormalization. Relative signal, not concentration.',
    'Supplements: positive-pixel fraction and log2(1+positive-pixel mean*1e6). Positive intensity NA when <20 positive pixels.',
    'Biological unit is physical sample; the upstream MSI input contains 10 samples and does not include Cushing-3; pool all three NF-GnH-3 exports. MSI n=2/3/1/1/3; LC n=2/3/1/3/3 in Cushing/GH/NF-GH/NF-ACTH/NF-GnH order.',
    'Common MSI: >=10% detected and >=20 positive pixels in >=ceil(2n/3) samples in EACH of five classes. strict_all10 additionally requires all 10.',
    'Common LC: above repeated global floor in >=2/3 samples in every class. Floor is inherited operational imputation convention.',
    'LC candidates: original tier A/B and singleton effect, AUC, direction, detection, highest-class gates; recomputed and validated against previous list.',
    'MSI candidates: eligible target samples >=2/3 (singleton all); highest class; effect vs all other samples >=0.15 (detection >=0.10); standardized effect>=0.8; AUC>=0.8; direction>=2/3; target leave-one-out superiority>=2/3; spatial bootstrap>=0.8; mean OR median sample-centering sensitivity.',
    'Positive analyses require finite values for >=ceil(2n/3) samples in EVERY class. Supplements explain primary candidates, not additional primary gates.',
    'Exact one-sided sample-label permutation p; BH correction within class and metric. P-values use sample-weighted one-vs-rest; reported equal_group_effect weights other class means equally.',
    'Spatial bootstrap: 200 source-specific 10x10-block resamples within physical sample; spatial sensitivity, not biological replication.',
    'LC-first: LC candidate to any MSI feature mass matches. Direction-supported subset: target MSI class highest and five-class Spearman>=0.',
    'Independent intersection: LC and MSI candidates pass independently in SAME target class; no extra Spearman selection.',
    'Common LC-first mass matches and common independent intersection are separately exported. Common does not mean equal abundance.',
    'Fresh mass matching against new m/z: previous six LC adducts, nine MSI ion forms, isotope numbers 0-2; absolute error<=0.1Da, no ppm cap. One best path per feature-LC entity, path rank before mass tier/error.',
    'Core: monoisotopic core adduct. direct_MS2: same observed ion plus LC fragmentation score>0; direct can include non-core dehydrated ions.',
    'LC MS2 supports LC annotation, not MSI chemical identity. Candidate correspondences and singleton subtype enrichments remain exploratory.',
    'Old MSI feature indices and candidate lists are not reused. previous_meanlog is exported only as sample values, not as stale candidate membership.'
  ),file.path(OUT,'METHODS.txt'))
  saveRDS(list(seed=20260914,bootstrap_B=B,block_width=10,lc_floor=floor_lc,source_script=.script2B),file.path(OUT,'parameters.rds'))
  writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
  savecsv(data.table(file=c(.script2B,file.path(INPUT,list.files(INPUT))),md5=unname(tools::md5sum(c(.script2B,file.path(INPUT,list.files(INPUT)))))), 'provenance_md5')
  writeLines(OUT,file.path(dirname(OUT),'LATEST_RUN.txt'))
  exportset<-function(d,name) {
    savecsv(d,name)
    for(k in c('primary','detection','positive')) {
      fi<-sort(unique(d$feature_index));vv<-get(k)
      savecsv(data.table(feature_index=fi,mz=mz[fi],as.data.table(t(vv[,fi,drop=FALSE]))),paste0(name,'__MSI_',k,'_matrix'))
    }
  }
  exportset(scommon[stable_all5==TRUE],'MSI_shared_features')
  exportset(expl,'MSI_subtype_specific_features')
  for(g in groups)exportset(expl[class==g],paste0('MSI_specific_',g))
  exportset(common_lcfirst,'route1_common_LCfirst')
  exportset(fread(file.path(OUT,'LCfirst_direction_primary.csv')),'route1_specific_direction_supported')
  exportset(commoncross,'route2_common_independent')
  exportset(fread(file.path(OUT,'independent_intersection_primary.csv')),'route2_specific_independent')
  oldcounts<-fread(file.path(OLD,'reanalysis_all_tumor_pixels_2B/20260908_094509/matching_counts.csv'))
  nc<-rbindlist(countsout);compare<-merge(oldcounts,nc,by='result_set',all=TRUE,suffixes=c('_old1430','_new274'));savecsv(compare,'old_new_matching_count_comparison')
  input_files<-c(file.path(BASE,c('pixel_by_feature_matrix_aligned_TIC.rds','feature_mz_aligned_filtered.rds','pixel_metadata_sample_source_xy.rds')),tumor_file,file.path(OLD,'lc_ms_ms.csv'))
  savecsv(data.table(file=input_files,size=file.info(input_files)$size,md5=unname(tools::md5sum(input_files))),'actual_data_input_md5')
  message('REANALYSIS_TUNED225_COMPLETE: ',OUT)
  print(cs);print(rbindlist(countsout))
  invisible(OUT)
}
result_subtype_tuned225 <- run_subtype_tuned225()
