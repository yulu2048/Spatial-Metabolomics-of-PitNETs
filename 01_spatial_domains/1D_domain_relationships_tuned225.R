# Compare existing domains; does not refit labels or change 1A/1B/1C outputs.
local({
  suppressPackageStartupMessages({library(data.table); library(Matrix); library(ggplot2)})
  root <- "/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225"
  indir <- file.path(root, "intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225")
  out <- file.path(root, "domain_relationships_21domains_1D_tuned225")
  dir.create(out, recursive=TRUE, showWarnings=FALSE)
  set.seed(20260914); nboot <- 200L; block_width <- 10
  meta <- as.data.table(readRDS(file.path(indir,"tumor_pixels_within_sample_consensus_domains.rds")))
  meta[, sample_id := as.character(sample_id)]
  roles <- fread(file.path(indir,"within_sample_discovery_validation_feature_qc.csv"))
  mz <- readRDS(file.path(root,"feature_mz_aligned_filtered.rds"))
  x <- readRDS(file.path(root,"pixel_by_feature_matrix_aligned_logCPM.rds"))
  stopifnot(!anyDuplicated(meta$original_pixel_row),all(meta$sample_id != "Cushing-3"),ncol(x)==length(mz))
  samples <- unique(meta$sample_id)
  common <- roles[,.(n=sum(passed_strict_filter)),by=matrix_feature_index][n==length(samples),matrix_feature_index]
  stopifnot(length(common)>=30L)
  domains <- unique(meta[,.(sample_id,diagnosis,consensus_domain,spatial_organization_label)])
  setorder(domains,sample_id,consensus_domain)
  domains[, domain_id := paste(sample_id,consensus_domain,sep=" / ")]
  nd <- nrow(domains); nf <- ncol(x)
  ef <- raw <- absolute <- matrix(NA_real_,nd,nf,dimnames=list(domains$domain_id,NULL))
  boot <- array(NA_real_,c(nboot,nd,nf)); block_qc <- list()
  safe_z <- function(v) {s<-sd(v); if(!is.finite(s)||s<1e-12) return(rep(0,length(v))); (v-mean(v))/s}
  rho <- function(a,b) {ok<-is.finite(a)&is.finite(b); if(sum(ok)<10L||sd(a[ok])==0||sd(b[ok])==0) return(NA_real_); as.numeric(stats::cor(a[ok],b[ok],method="spearman"))}
  for(si in seq_along(samples)) {
    sid<-samples[si]; message("[domain profiles / spatial bootstrap] ",sid)
    m<-copy(meta[sample_id==sid]); z0<-as.matrix(x[m$original_pixel_row,,drop=FALSE])
    zp<-safe_z(log1p(m$detected_peak_count)); za<-safe_z(m$alignment_nmatch)
    design<-cbind(1,zp,safe_z(lm.fit(cbind(1,zp),za)$residuals),safe_z(log1p(m$raw_tic)))
    qq<-qr(design); q<-qr.Q(qq)[,seq_len(qq$rank),drop=FALSE]
    z<-z0-q%*%crossprod(q,z0); z<-sweep(z,2,colMeans(z),"-")
    s<-sqrt(colSums(z^2)/(nrow(z)-1)); s[!is.finite(s)|s<1e-3]<-1
    z<-sweep(z,2,s,"/"); z[!is.finite(z)]<-0; z[z>5]<-5; z[z< -5]<- -5
    # Source-specific coordinates; statistics pooled across all exports of the physical sample.
    block<-paste(m$source_id,floor(m$x/block_width),floor(m$y/block_width),sep=":")
    keys<-unique(block); bid<-match(block,keys); nb<-length(keys)
    source_by_block<-as.character(m$source_id[match(keys,block)])
    strata<-split(seq_len(nb),source_by_block)
    w<-matrix(0,nboot,nb)
    for(b in seq_len(nboot)) for(v in strata) w[b,]<-w[b,]+tabulate(v[sample.int(length(v),length(v),replace=TRUE)],nbins=nb)
    core<-m$stable_core_pixel %in% TRUE
    for(j in which(domains$sample_id==sid)) {
      inside<-core & m$consensus_domain==domains$consensus_domain[j]; outside<-core & !inside
      stopifnot(sum(inside)>=20,sum(outside)>=20)
      ef[j,]<-colMeans(z[inside,,drop=FALSE])-colMeans(z[outside,,drop=FALSE])
      raw[j,]<-colMeans(z0[inside,,drop=FALSE])-colMeans(z0[outside,,drop=FALSE])
      absolute[j,]<-colMeans(z0[inside,,drop=FALSE])
      aggregate_blocks<-function(sel) {ans<-matrix(0,nb,nf); rs<-rowsum(z[sel,,drop=FALSE],bid[sel],reorder=FALSE); ans[as.integer(rownames(rs)),]<-rs; ans}
      ni<-as.vector(w%*%tabulate(bid[inside],nbins=nb)); no<-as.vector(w%*%tabulate(bid[outside],nbins=nb))
      bi<-(w%*%aggregate_blocks(inside))/ifelse(ni>0,ni,NA_real_)
      bo<-(w%*%aggregate_blocks(outside))/ifelse(no>0,no,NA_real_)
      boot[,j,]<-bi-bo
      block_qc[[length(block_qc)+1L]]<-data.table(domain_id=domains$domain_id[j],n_core=sum(inside),n_other_core=sum(outside),n_blocks=nb,n_domain_blocks=uniqueN(bid[inside]),n_valid_bootstrap=sum(ni>0 & no>0))
    }
  }
  rm(x); gc()
  sim<-stats::cor(t(ef[,common,drop=FALSE]),method="spearman")
  abs_sim<-stats::cor(t(absolute[,common,drop=FALSE]),method="spearman")
  hc<-hclust(as.dist(1-sim),method="average")
  val_sim<-matrix(NA_real_,nd,nd,dimnames=dimnames(sim)); diag(val_sim)<-1
  prs<-list(); contributors<-list(); ix<-0L
  for(i in seq_len(nd-1L)) for(j in (i+1L):nd) {
    ix<-ix+1L; sa<-domains$sample_id[i]; sb<-domains$sample_id[j]
    val<-intersect(roles[sample_id==sa & feature_role=="validation",matrix_feature_index],roles[sample_id==sb & feature_role=="validation",matrix_feature_index])
    r<-rho(ef[i,val],ef[j,val]); val_sim[i,j]<-val_sim[j,i]<-r
    br<-vapply(seq_len(nboot),function(b) rho(boot[b,i,val],boot[b,j,val]),numeric(1))
    ci<-if(sum(is.finite(br))>=.9*nboot) quantile(br,c(.025,.975),na.rm=TRUE) else c(NA_real_,NA_real_)
    high<-common[mz[common]>=200]
    prs[[ix]]<-data.table(domain_a=domains$domain_id[i],domain_b=domains$domain_id[j],sample_a=sa,sample_b=sb,same_sample=sa==sb,n_common_features=length(common),rho_common=sim[i,j],rho_absolute=abs_sim[i,j],rho_unresidualized_contrast=rho(raw[i,common],raw[j,common]),rho_mz_ge200=rho(ef[i,high],ef[j,high]),n_shared_validation=length(val),rho_validation=r,validation_ci_low=ci[1],validation_ci_high=ci[2],bootstrap_positive_fraction=mean(br>0,na.rm=TRUE),n_valid_bootstrap=sum(is.finite(br)))
    if(sa!=sb) contributors[[length(contributors)+1L]]<-data.table(domain_a=domains$domain_id[i],domain_b=domains$domain_id[j],matrix_feature_index=val,mz=mz[val],effect_a=ef[i,val],effect_b=ef[j,val],concordant=ef[i,val]*ef[j,val]>0,shared_effect_strength=pmin(abs(ef[i,val]),abs(ef[j,val])))
  }
  pairs<-rbindlist(prs)
  # Exploratory effect/support thresholds, not biological replication or inferential p-values.
  pairs[, supported_positive_pair := !same_sample & n_shared_validation>=30 & rho_common>=.5 & rho_validation>=.5 & validation_ci_low>0 & n_valid_bootstrap>=.9*nboot]
  pairs[, sensitivity_consistent := supported_positive_pair & rho_unresidualized_contrast>0 & rho_mz_ge200>0]
  fwrite(pairs,file.path(out,"domain_pair_similarity_and_validation.csv"))
  fwrite(pairs[supported_positive_pair==TRUE][order(-rho_validation)],file.path(out,"supported_cross_sample_pairs.csv"))
  fwrite(rbindlist(block_qc),file.path(out,"spatial_block_bootstrap_qc.csv")); fwrite(domains,file.path(out,"domain_metadata.csv"))
  ct<-rbindlist(contributors); fwrite(ct,file.path(out,"shared_validation_feature_effects_by_pair.csv"))
  top<-merge(ct,pairs[supported_positive_pair==TRUE,.(domain_a,domain_b)],by=c("domain_a","domain_b"))
  fwrite(top[concordant==TRUE][order(-shared_effect_strength),head(.SD,15),by=.(domain_a,domain_b)],file.path(out,"top15_concordant_validation_features_supported_pairs.csv"))
  saveRDS(list(domain_metadata=domains,common_features=common,mz=mz,effect=ef,unresidualized_effect=raw,absolute_mean=absolute,similarity=sim,validation_similarity=val_sim,hclust=hc,bootstrap_effect=boot),file.path(out,"domain_relationships.rds"))
  save_plot<-function(p,name,w=10,h=9) {ggsave(file.path(out,paste0(name,".pdf")),p,width=w,height=h); ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=200)}
  heat<-function(mat,name,title,ord=hc$order) {
    dd<-as.data.table(as.table(mat)); setnames(dd,c("Domain_A","Domain_B","Correlation"))
    dd[,Domain_A:=factor(Domain_A,levels=rownames(mat)[ord])]; dd[,Domain_B:=factor(Domain_B,levels=rev(rownames(mat)[ord]))]
    p<-ggplot(dd,aes(Domain_A,Domain_B,fill=Correlation))+geom_tile()+scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-1,1),na.value="grey85")+coord_equal()+labs(x=NULL,y=NULL,title=title,fill="Spearman r")+theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank())
    save_plot(p,name); p
  }
  mainplot<-heat(sim,"01_domain_relative_effect_similarity","Within-tumor domain effect similarity (194 shared QC features)")
  heat(val_sim,"02_pairwise_heldout_validation_similarity","Held-out feature similarity (pair-specific feature sets)")
  heat(abs_sim,"03_absolute_logCPM_similarity","Absolute domain mean logCPM: sensitivity comparison",hclust(as.dist(1-abs_sim),method="average")$order)
  pdf(file.path(out,"04_domain_relationship_dendrogram.pdf"),width=12,height=7)
  plot(hc,main="Relationships among 21 within-tumor domains",xlab="",sub="Average linkage; 1 - Spearman correlation of relative effects",ylab="Dissimilarity",cex=.8); dev.off()
  png(file.path(out,"04_domain_relationship_dendrogram.png"),width=2400,height=1400,res=200)
  plot(hc,main="Relationships among 21 within-tumor domains",xlab="",sub="Exploratory branches, not predefined common states",ylab="Dissimilarity",cex=.8); dev.off()
  save_plot(ggplot(pairs[same_sample==FALSE],aes(rho_common,rho_validation,color=supported_positive_pair))+geom_hline(yintercept=.5,lty=2,color="grey70")+geom_vline(xintercept=.5,lty=2,color="grey70")+geom_point(size=2)+scale_color_manual(values=c("grey65","#B2182B"))+labs(x="Similarity: common QC features",y="Similarity: shared held-out features",color="Supported pair")+theme_classic(),"05_discovery_validation_concordance",7,6)
  # Do not infer a preferred number of meta-domains: export descriptive cuts only.
  cuts<-copy(domains); for(k in 2:6) cuts[,(paste0("exploratory_cut_k",k)):=as.integer(cutree(hc,k))]
  fwrite(cuts,file.path(out,"exploratory_dendrogram_cuts_not_common_domain_labels.csv"))
  long<-data.table(domain_id=rep(domains$domain_id,each=nf),matrix_feature_index=rep(seq_len(nf),nd),mz=rep(mz,nd),relative_effect=as.vector(t(ef)),raw_logCPM_difference=as.vector(t(raw)),absolute_logCPM_mean=as.vector(t(absolute)))
  fwrite(long,file.path(out,"domain_feature_profiles.csv"))
  settings<-c("Primary profiles reproduce 1A: within-physical-sample residualization of logCPM against detection count, alignment matches and log raw TIC; SD scaling; cap +/-5.","Domain contrasts use stable core versus other stable cores of the SAME sample. NF-GnH-3 exports pooled.",paste("Primary fixed feature set passes strict 1A QC in all samples; n =",length(common)),"Held-out validation uses pairwise intersection of validation roles (not necessarily primary feature set).",paste("Bootstrap:",nboot,"resamples of",block_width,"coordinate-unit square blocks, stratified by export; fixed labels, fixed core selection, fixed residualization."),"Bootstrap CIs quantify conditional spatial sampling sensitivity, NOT patient replication or a fresh clustering validation. They depend on block size.","Supported cross-sample pair: >=30 jointly held-out features; common and held-out rho >=0.5; held-out bootstrap lower 95% bound >0; >=90% valid resamples. Exploratory thresholds.","Two domains within one sample have exactly opposite contrasts; two dendrogram branches cannot establish two universal biological states.","Held-out ions are not necessarily chemically independent; correlations do not identify metabolites or prove shared mechanisms.")
  writeLines(settings,file.path(out,"analysis_settings_and_limitations.txt")); capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
  print(pairs[supported_positive_pair==TRUE][order(-rho_validation),.(domain_a,domain_b,rho_common,n_shared_validation,rho_validation,validation_ci_low,validation_ci_high,sensitivity_consistent)])
  print(pairs[,.(n_pairs=.N,n_supported=sum(supported_positive_pair),n_sensitivity_consistent=sum(sensitivity_consistent)),by=same_sample])
  print(mainplot); message("[DONE] ",out)
})

# Additional tree stability and an audit against the original 1A effect tables.
local({
  library(data.table); library(ggplot2)
  root<-"/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225"
  out<-file.path(root,"domain_relationships_21domains_1D_tuned225")
  a<-file.path(root,"intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225")
  r<-readRDS(file.path(out,"domain_relationships.rds")); dm<-r$domain_metadata
  p<-fread(file.path(out,"domain_pair_similarity_and_validation.csv"))
  old<-rbindlist(list(fread(file.path(a,"independent_validation_feature_effects_by_domain.csv")),fread(file.path(a,"descriptive_discovery_feature_effects_by_domain.csv"))))
  old[,domain_id:=paste(sample_id,domain,sep=" / ")]
  old[,recomputed_effect:=r$effect[cbind(match(domain_id,dm$domain_id),matrix_feature_index)]]
  discrepancy<-max(abs(old$recomputed_effect-old$mean_residualized_z_difference))
  fwrite(data.table(n_compared=nrow(old),maximum_absolute_difference=discrepancy),file.path(out,"effect_reproduction_audit.csv"))
  stopifnot(discrepancy<1e-8)
  leaves<-function(h) {v<-vector("list",nrow(h$merge)); for(i in seq_len(nrow(h$merge))) v[[i]]<-sort(unlist(lapply(h$merge[i,],function(j) if(j<0) -j else v[[j]]))); v}
  key<-function(v) paste(sort(v),collapse=",")
  ref<-leaves(r$hclust); support<-numeric(length(ref)); co<-lapply(2:6,function(k) matrix(0,nrow(dm),nrow(dm),dimnames=list(dm$domain_id,dm$domain_id)))
  valid<-0L
  for(b in seq_len(dim(r$bootstrap_effect)[1])) {
    eb<-r$bootstrap_effect[b,,r$common_features]
    if(any(!is.finite(eb))) next
    h<-hclust(as.dist(1-stats::cor(t(eb),method="spearman")),method="average")
    bk<-vapply(leaves(h),key,character(1)); support<-support+(vapply(ref,key,character(1)) %in% bk)
    valid<-valid+1L
    for(k in 2:6) {cl<-cutree(h,k); co[[k-1L]]<-co[[k-1L]]+outer(cl,cl,"==")}
  }
  node<-data.table(node=seq_along(ref),n_domains=lengths(ref),height=r$hclust$height,conditional_bootstrap_support=support/valid,members=vapply(ref,function(v) paste(dm$domain_id[v],collapse="; "),character(1)))
  fwrite(node,file.path(out,"dendrogram_branch_bootstrap_support.csv"))
  stability<-rbindlist(lapply(2:6,function(k) {cl<-cutree(r$hclust,k); rbindlist(lapply(unique(cl),function(g) {ii<-which(cl==g); cc<-co[[k-1L]][ii,ii,drop=FALSE]/valid; data.table(k=k,branch=g,n_domains=length(ii),n_samples=uniqueN(dm$sample_id[ii]),mean_pair_coassignment=if(length(ii)>1)mean(cc[upper.tri(cc)])else NA_real_,minimum_pair_coassignment=if(length(ii)>1)min(cc[upper.tri(cc)])else NA_real_,members=paste(dm$domain_id[ii],collapse="; "))}))}))
  fwrite(stability,file.path(out,"exploratory_cut_bootstrap_stability.csv"))
  cp<-rbindlist(lapply(2:6,function(k) {v<-as.data.table(as.table(co[[k-1L]]/valid)); setnames(v,c("i","j","coassignment")); v[,c("domain_a","domain_b","k"):=list(as.character(i),as.character(j),k)]; v[,.(domain_a,domain_b,k,coassignment)]}))
  fwrite(cp,file.path(out,"bootstrap_coassignment_all_cuts.csv"))
  dd<-cp[k==4]; ord<-dm$domain_id[r$hclust$order]
  dd[,domain_a:=factor(domain_a,levels=ord)]; dd[,domain_b:=factor(domain_b,levels=rev(ord))]
  g<-ggplot(dd,aes(domain_a,domain_b,fill=coassignment))+geom_tile()+scale_fill_gradient(low="white",high="#54278F",limits=c(0,1))+coord_equal()+labs(x=NULL,y=NULL,title="Exploratory four-branch coassignment",subtitle="Conditional spatial bootstrap; not refitted pixel domains",fill="Proportion")+theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=60,hjust=1),panel.grid=element_blank())
  ggsave(file.path(out,"06_exploratory_k4_coassignment.pdf"),g,width=10,height=9); ggsave(file.path(out,"06_exploratory_k4_coassignment.png"),g,width=10,height=9,dpi=200)
  # Each physical sample pair counted once here; complementary domain pairs are not independent evidence.
  sp<-p[same_sample==FALSE,.(n_domain_pairs=.N,n_supported=sum(supported_positive_pair),n_sensitivity_consistent=sum(sensitivity_consistent),maximum_validation_rho=max(rho_validation)),by=.(sample_a,sample_b)]
  fwrite(sp,file.path(out,"physical_sample_pair_summary.csv"))
  near<-rbindlist(list(p[same_sample==FALSE,.(domain_id=domain_a,partner=domain_b,rho_validation,rho_common,supported_positive_pair)],p[same_sample==FALSE,.(domain_id=domain_b,partner=domain_a,rho_validation,rho_common,supported_positive_pair)]))
  fwrite(near[order(-rho_validation),head(.SD,1),by=domain_id],file.path(out,"nearest_cross_sample_domain.csv"))
  if(requireNamespace("pheatmap",quietly=TRUE)) {
    annotation<-data.frame(Diagnosis=as.character(dm$diagnosis),Organization=as.character(dm$spatial_organization_label),row.names=dm$domain_id)
    for(ext in c("pdf","png")) pheatmap::pheatmap(r$similarity,cluster_rows=r$hclust,cluster_cols=r$hclust,annotation_row=annotation,annotation_col=annotation,color=colorRampPalette(c("#2166AC","white","#B2182B"))(100),breaks=seq(-1,1,length.out=101),border_color=NA,fontsize=9,main="Relative domain-effect similarity",filename=file.path(out,paste0("07_annotated_domain_similarity.",ext)),width=13,height=11)
  }
  print(stability[k==4]); print(node[n_domains>=3 & n_domains<21]); print(data.table(n_valid_tree_bootstrap=valid,maximum_1A_effect_discrepancy=discrepancy))
  message("[TREE STABILITY DONE] ",out)
})
