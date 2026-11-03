# Rebuilt 2026-09-14 against tuned225 matrix; source: /home/gaga/huzhan/MSI_ALL_aligned_11samples_consensus0.1Da/spatial_metabolite_modules_10samples/5A_step1_feature_universe.R
# 5A Step 1: freeze MSI-only feature universe and mass-relation candidates
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
BASE <- '/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
OUT <- file.path(BASE,'spatial_metabolite_modules_10samples_tuned225')
S1 <- file.path(OUT,'step1_feature_universe')
dir.create(S1,recursive=TRUE,showWarnings=FALSE)
qfile <- file.path(BASE,'cross_sample_common_domains_10samples_1B_tuned225','cross_sample_feature_qc_physical_sample_level.csv')
q <- fread(qfile)
stopifnot(all(q$matrix_feature_index==seq_len(nrow(q))))
q[, module_discovery_eligible := passed_strict_technical_filter]
q[, exclusion_reason := fifelse(module_discovery_eligible,'eligible',fifelse(n_physical_samples_passing_prevalence < min_physical_samples_required,'prevalence_lt_8_samples',fifelse(median_abs_detection_qc_correlation > 0.35,'technical_correlation','other_existing_QC')))]
features <- q[module_discovery_eligible==TRUE, .(feature_index=matrix_feature_index,mz,n_samples_prevalent=n_physical_samples_passing_prevalence,median_logCPM_variance=median_physical_sample_logCPM_variance,fraction_samples_variable=fraction_physical_samples_variable,median_abs_detection_qc_correlation)]
stopifnot(nrow(features)>=50L)
# Candidate ion-family relations only; no feature is collapsed until spatial concordance is evaluated.
deltas <- data.table(relation=c('isotope_M+1','isotope_M+2','Na_minus_H','K_minus_H','K_minus_Na','NH4_minus_H_or_NH3','neutral_loss_H2O'),target_delta=c(1.003355,2.006710,21.981943,37.955882,15.973939,17.026549,18.010565),tolerance_Da=c(.05,.05,.05,.05,.05,.05,.05))
m <- features$mz; idx <- features$feature_index
pairs <- rbindlist(lapply(seq_len(nrow(deltas)),function(k){d<-deltas[k]; a<-outer(m,m,function(a,b) abs(abs(a-b)-d$target_delta)<=d$tolerance_Da); z<-which(upper.tri(a)&a,arr.ind=TRUE); if(!nrow(z)) return(NULL); data.table(feature_i=idx[z[,1]],mz_i=m[z[,1]],feature_j=idx[z[,2]],mz_j=m[z[,2]],observed_delta=abs(m[z[,1]]-m[z[,2]]),relation=d$relation,target_delta=d$target_delta,delta_error_Da=abs(abs(m[z[,1]]-m[z[,2]])-d$target_delta))}))
setorder(pairs,feature_i,feature_j,delta_error_Da)
fwrite(q,file.path(S1,'all_1430_features_with_frozen_module_eligibility.csv'))
fwrite(features,file.path(S1,'eligible_876_MSI_features.csv'))
fwrite(pairs,file.path(S1,'mass_relation_candidate_pairs_before_spatial_check.csv'))
summary <- rbindlist(list(data.table(metric=c('all_aligned_MSI_features','eligible_after_existing_prevalence_variance_technical_QC','excluded','eligible_fraction'),value=c(nrow(q),nrow(features),sum(!q$module_discovery_eligible),nrow(features)/nrow(q))),pairs[,.(value=.N),by=.(metric=paste0('candidate_pairs_',relation))]))
fwrite(summary,file.path(S1,'step1_summary.csv'))
p1<-ggplot(q,aes(x=mz,fill=module_discovery_eligible))+geom_histogram(binwidth=25,boundary=0,color='white',linewidth=.15)+scale_fill_manual(values=c('FALSE'='#B8B8B8','TRUE'='#2A6FBB'),labels=c('Excluded','Eligible'))+labs(x='m/z',y='Feature count',fill=NULL,title='Step 1: frozen MSI feature universe')+theme_bw(base_size=11)
p2<-ggplot(q,aes(x=n_physical_samples_passing_prevalence,fill=module_discovery_eligible))+geom_histogram(binwidth=1,boundary=-.5,color='white')+scale_fill_manual(values=c('FALSE'='#B8B8B8','TRUE'='#2A6FBB'),guide='none')+scale_x_continuous(breaks=0:10)+labs(x='Physical samples passing prevalence',y='Feature count')+theme_bw(base_size=11)
ggsave(file.path(S1,'step1_feature_universe.pdf'),p1/p2,width=8,height=7)
params <- list(analysis_principle='MSI-only discovery; domain labels and LC-MS/MS are prohibited from changing eligibility, correlations, ion-family collapse, or modules',eligible_rule='Reuse the pre-domain 1B strict technical filter: prevalence in >=8/10 physical samples plus variance and detection-QC safeguards',n_all=nrow(q),n_eligible=nrow(features),candidate_mass_relations=deltas,input_md5=tools::md5sum(c(qfile,file.path(BASE,'feature_mz_aligned_filtered.rds'),file.path(BASE,'pixel_by_feature_matrix_aligned_logCPM.rds'))),session_info=utils::sessionInfo())
saveRDS(params,file.path(S1,'step1_parameters_and_sessionInfo.rds'))
writeLines(c('Step 1 freezes the MSI-only feature universe before any module analysis.','The 1B strict feature QC is reused because it was computed before domain clustering and contains only prevalence, variance and technical-correlation criteria.','Mass relations are only candidates. They will be collapsed only after very high, recurrent spatial correlation is demonstrated.','LC-MS/MS is not read in this step.'),file.path(S1,'README.txt'))
cat('STEP1_COMPLETE\n'); print(summary); cat('OUT=',S1,'\n')
