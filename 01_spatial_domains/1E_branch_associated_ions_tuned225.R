# Exploratory ion interpretation of the fixed four-branch view from 1D.
# No refitting, no pixel-level p-values, no claim of independent tree validation.
local({
  suppressPackageStartupMessages({library(data.table); library(Matrix); library(ggplot2)})
  root<-"/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225"
  ddir<-file.path(root,"domain_relationships_21domains_1D_tuned225")
  adir<-file.path(root,"intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225")
  out<-file.path(ddir,"branch_associated_ions_1E")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  set.seed(20260915); nboot<-2000L; min_samples<-3L; min_effect<-.35; min_direction<-.8
  r<-readRDS(file.path(ddir,"domain_relationships.rds")); dm<-copy(r$domain_metadata)
  cuts<-fread(file.path(ddir,"exploratory_dendrogram_cuts_not_common_domain_labels.csv"))
  dm[,branch:=cuts$exploratory_cut_k4[match(domain_id,cuts$domain_id)]]
  core<-fread(file.path(ddir,"spatial_block_bootstrap_qc.csv"))
  dm[,n_core:=core$n_core[match(domain_id,core$domain_id)]]
  roles<-fread(file.path(adir,"within_sample_discovery_validation_feature_qc.csv"))
  nf<-ncol(r$effect); nd<-nrow(dm)
  stopifnot(nd==21L,!anyNA(dm$branch),!anyDuplicated(dm[,.(sample_id,branch)]))
  passed<-validation<-matrix(FALSE,nd,nf)
  for(i in seq_len(nd)) {
    rr<-roles[sample_id==dm$sample_id[i]]
    passed[i,rr$matrix_feature_index]<-rr$passed_strict_filter
    validation[i,rr$matrix_feature_index]<-rr$feature_role=="validation"
  }
  med<-function(v) if(length(v)) median(v) else NA_real_
  frac<-function(v) if(length(v)) mean(v) else NA_real_
  ci<-function(v) {
    if(length(v)<min_samples) return(c(NA_real_,NA_real_))
    bm<-matrix(v[sample.int(length(v),length(v)*nboot,replace=TRUE)],nrow=nboot)
    as.numeric(quantile(apply(bm,1,median),c(.025,.975)))
  }
  spatial_pos<-spatial_neg<-matrix(NA_real_,nd,nf)
  for(i in seq_len(nd)) {spatial_pos[i,]<-colMeans(r$bootstrap_effect[,i,]>0,na.rm=TRUE); spatial_neg[i,]<-colMeans(r$bootstrap_effect[,i,]<0,na.rm=TRUE)}
  message("[1E A] Branch-versus-own-tumor contrasts; physical-sample bootstrap.")
  rows<-list(); ix<-0L
  for(g in 1:4) for(f in seq_len(nf)) {
    ids<-which(dm$branch==g & passed[,f]); vids<-ids[validation[ids,f]]
    v<-r$effect[ids,f]; raw<-r$unresidualized_effect[ids,f]
    m<-med(v); sign_m<-if(is.finite(m)&&m<0) -1 else 1
    bound<-ci(v); vv<-r$effect[vids,f]; vb<-ci(vv)
    sp<-if(sign_m>0) spatial_pos[ids,f] else spatial_neg[ids,f]
    loo<-if(length(v)>=min_samples) vapply(seq_along(v),function(j) median(v[-j]),numeric(1)) else numeric()
    ix<-ix+1L
    rows[[ix]]<-data.table(branch=g,matrix_feature_index=f,mz=r$mz[f],n_branch_samples=sum(dm$branch==g),n_qc_samples=length(ids),median_effect=m,direction=if(sign_m>0)"enriched"else"depleted",direction_fraction=frac(sign_m*v>0),sample_bootstrap_ci_low=bound[1],sample_bootstrap_ci_high=bound[2],leave_one_sample_out_direction_fraction=frac(sign_m*loo>0),spatial_supported_sample_fraction=frac(sp>=.95),median_raw_logCPM_difference=med(raw),raw_direction_fraction=frac(sign_m*raw>0),n_nonmodel_samples=length(vids),nonmodel_median_effect=med(vv),nonmodel_direction_fraction=frac(sign_m*vv>0),nonmodel_ci_low=vb[1],nonmodel_ci_high=vb[2],used_in_tree_feature_set=f %in% r$common_features)
  }
  stats<-rbindlist(rows)
  stats[,low_mz_review:=mz<200]
  stats[,candidate:=n_qc_samples>=pmax(min_samples,ceiling(.8*n_branch_samples)) & abs(median_effect)>=min_effect & direction_fraction>=min_direction & ((direction=="enriched" & sample_bootstrap_ci_low>0)|(direction=="depleted" & sample_bootstrap_ci_high<0)) & leave_one_sample_out_direction_fraction==1]
  stats[,nonmodel_supported:=candidate & n_nonmodel_samples>=min_samples & abs(nonmodel_median_effect)>=min_effect & nonmodel_direction_fraction>=min_direction & ((direction=="enriched" & nonmodel_ci_low>0)|(direction=="depleted" & nonmodel_ci_high<0))]
  stats[,robust_association:=nonmodel_supported & raw_direction_fraction>=min_direction & spatial_supported_sample_fraction>=min_direction]
  for(nm in c("candidate","nonmodel_supported","robust_association")) set(stats,which(is.na(stats[[nm]])),nm,FALSE)
  stats[,specificity_status:="not_established_insufficient_shared_samples_for_all_pairs"]
  fwrite(stats,file.path(out,"branch_ion_statistics_all.csv"))
  fwrite(stats[candidate==TRUE],file.path(out,"branch_associated_candidates.csv"))
  fwrite(stats[robust_association==TRUE],file.path(out,"branch_ions_with_nonmodel_and_sensitivity_support.csv"))
  message("[1E B] Paired branch comparisons; no unpaired substitution.")
  # Recover each domain's mean relative to the sample's pooled stable-core mean:
  # (mu_d - mu_rest) * (N - n_d)/N = mu_d - mu_all.
  dm[,n_sample_core:=sum(n_core),by=sample_id]
  centered<-r$effect*((dm$n_sample_core-dm$n_core)/dm$n_sample_core)
  design<-list(); comparisons<-list(); evidence<-list(); jx<-0L; ex<-0L
  for(a in 1:3) for(b in (a+1):4) {
    ss<-intersect(dm$sample_id[dm$branch==a],dm$sample_id[dm$branch==b])
    design[[length(design)+1L]]<-data.table(branch_a=a,branch_b=b,n_shared_samples=length(ss),shared_samples=paste(ss,collapse="; "),comparison_status=if(length(ss)>=min_samples)"estimable_exploratory"else"insufficient_shared_samples")
    ia<-match(paste(ss,a),paste(dm$sample_id,dm$branch)); ib<-match(paste(ss,b),paste(dm$sample_id,dm$branch))
    for(f in seq_len(nf)) {
      ok<-passed[ia,f]&passed[ib,f]; ai<-ia[ok]; bi<-ib[ok]; ids<-ss[ok]
      v<-centered[ai,f]-centered[bi,f]; rawv<-r$absolute_mean[ai,f]-r$absolute_mean[bi,f]
      vi<-validation[ai,f]&validation[bi,f]; m<-med(v); sgn<-if(is.finite(m)&&m<0)-1 else 1; bd<-ci(v)
      jx<-jx+1L
      comparisons[[jx]]<-data.table(branch_a=a,branch_b=b,matrix_feature_index=f,mz=r$mz[f],n_qc_shared_samples=length(v),n_nonmodel_samples=sum(vi),median_paired_difference=m,direction_fraction=frac(sgn*v>0),ci_low=bd[1],ci_high=bd[2],raw_direction_fraction=frac(sgn*rawv>0),nonmodel_median_difference=med(v[vi]),nonmodel_direction_fraction=frac(sgn*v[vi]>0))
      if(length(ids)) {ex<-ex+1L; evidence[[ex]]<-data.table(branch_a=a,branch_b=b,sample_id=ids,matrix_feature_index=f,mz=r$mz[f],paired_difference=v,raw_logCPM_difference=rawv,not_used_in_pixel_domain_model=vi)}
    }
  }
  pw<-rbindlist(comparisons); ds<-rbindlist(design)
  pw[,paired_candidate:=n_qc_shared_samples>=min_samples & abs(median_paired_difference)>=min_effect & direction_fraction>=min_direction & ((median_paired_difference>0 & ci_low>0)|(median_paired_difference<0 & ci_high<0))]
  pw[,paired_nonmodel_supported:=paired_candidate & n_nonmodel_samples>=min_samples & nonmodel_direction_fraction>=min_direction & abs(nonmodel_median_difference)>=min_effect & raw_direction_fraction>=min_direction]
  for(nm in c("paired_candidate","paired_nonmodel_supported")) set(pw,which(is.na(pw[[nm]])),nm,FALSE)
  fwrite(ds,file.path(out,"branch_pair_comparison_design.csv")); fwrite(pw,file.path(out,"paired_branch_ion_statistics_all.csv"))
  fwrite(pw[paired_nonmodel_supported==TRUE],file.path(out,"paired_branch_ions_with_nonmodel_support.csv"))
  fwrite(rbindlist(evidence),file.path(out,"physical_sample_paired_branch_ion_effects.csv"))
  counts<-stats[,.(n_candidates=sum(candidate),n_nonmodel_supported=sum(nonmodel_supported),n_robust=sum(robust_association)),by=.(branch,direction)]
  fwrite(counts,file.path(out,"branch_ion_counts.csv")); fwrite(dm,file.path(out,"branch_domain_metadata.csv"))
  # Per-sample evidence, never treating domains or export partitions as independent patients.
  long<-data.table(domain_id=rep(dm$domain_id,each=nf),sample_id=rep(dm$sample_id,each=nf),branch=rep(dm$branch,each=nf),matrix_feature_index=rep(seq_len(nf),nd),mz=rep(r$mz,nd),effect=as.vector(t(r$effect)),passed_qc=as.vector(t(passed)),not_used_in_pixel_domain_model=as.vector(t(validation)))
  fwrite(long,file.path(out,"domain_ion_effects_with_feature_roles.csv"))
  message("[1E C] Effect heatmaps, sample-level plots and representative ion maps.")
  chosen<-stats[candidate==TRUE][order(-robust_association,-nonmodel_supported,-abs(median_effect)),head(.SD,6L),by=.(branch,direction)]
  fwrite(chosen,file.path(out,"heatmap_selected_ions.csv")); features<-unique(chosen$matrix_feature_index)
  annotation<-data.frame(Branch=paste0("B",dm$branch),Diagnosis=as.character(dm$diagnosis),row.names=dm$domain_id)
  if(length(features)>=2L) {
    hm<-t(r$effect[,features,drop=FALSE]); rownames(hm)<-sprintf("%.4f",r$mz[features]); colnames(hm)<-dm$domain_id
    dots<-t(ifelse(validation[,features,drop=FALSE],".","")); dimnames(dots)<-dimnames(hm)
    hm[hm>2]<-2; hm[hm< -2]<- -2
    for(ext in c("pdf","png")) pheatmap::pheatmap(hm,cluster_cols=r$hclust,cluster_rows=TRUE,annotation_col=annotation,display_numbers=dots,fontsize_number=8,fontsize=9,color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),breaks=seq(-2,2,length.out=101),border_color=NA,main="Branch-associated ion effects\nDots: not used for pixel-domain modelling",filename=file.path(out,paste0("01_ion_by_21domain_effect_heatmap.",ext)),width=13,height=max(9,length(features)*.18+3))
  }
  saveplot<-function(g,name,w=10,h=8) {for(ext in c("pdf","png")) ggsave(file.path(out,paste0(name,".",ext)),g,width=w,height=h,dpi=200)}
  display<-stats[candidate==TRUE][order(-robust_association,-nonmodel_supported,-abs(median_effect)),head(.SD,5),by=.(branch,direction)]
  display[,ion_label:=sprintf("%.4f",mz)]
  saveplot(ggplot(display,aes(median_effect,ion_label,color=robust_association))+geom_segment(aes(x=sample_bootstrap_ci_low,xend=sample_bootstrap_ci_high,yend=ion_label))+geom_point()+geom_vline(xintercept=0,lty=2)+facet_wrap(~branch,scales="free_y",labeller=label_both)+labs(x="Median within-tumor effect (sample-bootstrap 95% interval)",y="m/z",color="Nonmodel +\nsensitivity support")+theme_bw(),"02_branch_ion_effect_intervals",12,10)
  reps<-stats[robust_association==TRUE][order(direction!="enriched",-abs(median_effect)),head(.SD,2),by=branch]
  fwrite(reps,file.path(out,"representative_ions_for_spatial_maps.csv"))
  if(nrow(reps)) {
    points<-merge(long,reps[,.(branch,matrix_feature_index)],by=c("branch","matrix_feature_index"))
    points[,panel:=paste0("B",branch," | m/z ",sprintf("%.4f",mz))]
    saveplot(ggplot(points,aes(sample_id,effect,color=not_used_in_pixel_domain_model))+geom_hline(yintercept=0,lty=2)+geom_point(size=2)+facet_wrap(~panel,scales="free_x",ncol=2)+labs(x=NULL,y="Within-tumor effect",color="Not used in\npixel model")+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1)),"03_representative_ion_sample_consistency",12,11)
    pm<-as.data.table(readRDS(file.path(adir,"tumor_pixels_within_sample_consensus_domains.rds")))
    mat<-readRDS(file.path(root,"pixel_by_feature_matrix_aligned_logCPM.rds"))
    dir.create(file.path(out,"ion_maps"),showWarnings=FALSE)
    for(f in unique(reps$matrix_feature_index)) {
      mp<-pm[,.(sample_id,source_id,x,y)]; mp[,value:=as.numeric(mat[pm$original_pixel_row,f])]
      mp[,panel:=paste(sample_id,source_id,sep=" / ")]; limits<-as.numeric(quantile(mp$value,c(.02,.98)))
      if(diff(limits)==0) limits<-range(mp$value)
      g<-ggplot(mp,aes(x,y,fill=value))+geom_tile(width=1,height=1)+coord_fixed()+scale_y_reverse()+facet_wrap(~panel,ncol=3)+scale_fill_viridis_c(limits=limits,oob=scales::squish)+labs(title=paste("m/z",sprintf("%.4f",r$mz[f])),subtitle="Tumor pixels; original logCPM; shared color scale across sources",x="X",y="Y",fill="logCPM")+theme_bw(base_size=9)+theme(panel.grid=element_blank())
      saveplot(g,paste0("ion_maps/mz_",sprintf("%.4f",r$mz[f])),12,13)
    }
  }
  paired_top<-pw[paired_nonmodel_supported==TRUE][order(-abs(median_paired_difference)),head(.SD,10L),by=.(branch_a,branch_b)]
  fwrite(paired_top,file.path(out,"top_paired_branch_ions.csv"))
  pe<-merge(rbindlist(evidence),paired_top[,.(branch_a,branch_b,matrix_feature_index)],by=c("branch_a","branch_b","matrix_feature_index"))
  if(nrow(pe)) {
    pe[,contrast:=paste0("B",branch_a," - B",branch_b)]; pe[,ion:=sprintf("%.4f",mz)]
    saveplot(ggplot(pe,aes(sample_id,ion,fill=paired_difference))+geom_tile()+geom_point(data=pe[not_used_in_pixel_domain_model==TRUE],size=.7)+facet_wrap(~contrast,scales="free",ncol=2)+scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-2,2),oob=scales::squish)+labs(x=NULL,y="m/z",fill="Paired effect",title="Within-sample branch contrasts",subtitle="Dots: feature not used in that sample's pixel-domain model")+theme_bw()+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank()),"04_paired_branch_ion_effects",11,7)
  }
  writeLines(c("Fixed exploratory k=4 tree branches; no refitting.","Each domain effect is stable core versus remaining stable core within the same physical sample; sample contributions are equally weighted.","Candidate: QC in >=80% of branch samples and >=3; |median effect|>=0.35; >=80% same direction; sample-bootstrap 95% interval excludes zero; all leave-one-sample-out medians retain direction.","Nonmodel support: >=3 samples in which feature was not used in pixel-domain modelling; same magnitude/direction thresholds; sample-bootstrap interval excludes zero.","Robust association additionally requires >=80% samples to retain direction in raw logCPM contrasts and >=80% samples to show >=95% same-direction conditional spatial-bootstrap support.","Bootstrap of physical samples: 2000 replicates, median statistic; only 5-6 patients per branch; intervals are exploratory and not multiplicity-adjusted tests.","Nonmodel features are NOT independent tree validation: some contributed to the tree and other ions may be correlated.","Pairwise comparisons require >=3 physical samples containing both branches; missing branches are never assigned zero.","No ion can currently be established as specific against all three other branches: B1/B3 and B1/B4 have zero shared samples; B2/B3 and B2/B4 have one.","Pairwise residual-domain differences reconstruct sample-centered domain means from one-versus-rest effects and core counts. Their uncertainty uses sample bootstrap, not the one-versus-rest spatial-bootstrap arrays.","Ions not chemically annotated or deduplicated here. m/z<200 is flagged by mass alone, not removed or identified as matrix signal.","Heatmap: selected exploratory candidates, standardized residual effects clipped to +/-2 for display only; dots indicate not used in corresponding sample pixel model.","Ion maps: tumor pixels only, original logCPM, shared 2nd-98th percentile color limits across sources within each ion; not an absolute concentration measurement."),file.path(out,"analysis_rules_and_limitations.txt"))
  capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  print(counts); print(ds); print(pw[,.(n_paired_candidates=sum(paired_candidate),n_nonmodel_supported=sum(paired_nonmodel_supported)),by=.(branch_a,branch_b)])
  print(reps[,.(branch,mz,median_effect,n_nonmodel_samples,nonmodel_direction_fraction)])
  message("[1E DONE] ",out)
})
