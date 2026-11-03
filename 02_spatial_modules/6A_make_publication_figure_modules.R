# Rebuilt 2026-09-14 against tuned225 matrix; source: /home/gaga/huzhan/MSI_ALL_aligned_11samples_consensus0.1Da/spatial_metabolite_modules_10samples/6A_make_publication_figure_modules.R
# Publication figure: MSI metabolic modules, spatial organization, domain validation and LC-MS/MS support
local({
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(ggrepel);library(viridis)})
BASE <- '/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
OUT <- file.path(BASE,'spatial_metabolite_modules_10samples_tuned225')
FD <- file.path(OUT,'publication_figure_metabolic_modules')
dir.create(FD,recursive=TRUE,showWarnings=FALSE)
mods <- sort(unique(fread(file.path(OUT,'MASTER_module_evidence_table.csv'))$module))
sample_order <- c('GH-1','GH-2','GH-3','NF-GH','NF-ACTH','Cushing-1','Cushing-2','NF-GnH-1','NF-GnH-2','NF-GnH-3')
cols_flag <- c(high_stability='#2166AC',moderate_stability='#FDAE61',exploratory='#BDBDBD')
theme_set(theme_bw(base_size=9)+theme(panel.grid=element_blank(),plot.title=element_text(face='bold',size=10),plot.subtitle=element_text(size=8,color='grey30'),legend.title=element_text(size=8),legend.text=element_text(size=7)))
save_both <- function(p,name,w,h){ggsave(file.path(FD,paste0(name,'.pdf')),p,width=w,height=h,units='in',device=cairo_pdf);ggsave(file.path(FD,paste0(name,'.png')),p,width=w,height=h,units='in',dpi=300,bg='white')}
master <- fread(file.path(OUT,'MASTER_module_evidence_table.csv'))
master[,module:=factor(module,levels=mods)]
# B: consensus module structure and robustness
mem <- fread(file.path(OUT,'step3_modules','representative_feature_module_membership.csv'))[module!='grey',.(feature_index,module)]
pairs <- fread(file.path(OUT,'step2_spatial_consensus','all_consensus_feature_pairs.csv'))
pp <- merge(pairs,mem,by.x='feature_i',by.y='feature_index');setnames(pp,'module','mi')
pp <- merge(pp,mem,by.x='feature_j',by.y='feature_index');setnames(pp,'module','mj')
bm <- pp[,.(rho=median(median_spearman,na.rm=TRUE),n_pairs=.N),by=.(mi,mj)]
bm <- rbind(bm,bm[,.(mi=mj,mj=mi,rho,n_pairs)])[,.(rho=mean(rho),n_pairs=sum(n_pairs)),by=.(mi,mj)]
diagv <- master[,.(mi=as.character(module),mj=as.character(module),rho=median_within_consensus_rho,n_pairs=NA_integer_)]
bm <- rbind(bm[mi!=mj],diagv,fill=TRUE)
bm[,mi:=factor(mi,levels=mods)];bm[,mj:=factor(mj,levels=rev(mods))]
pB1 <- ggplot(bm,aes(mi,mj,fill=rho))+geom_tile()+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-.2,.6),oob=scales::squish,na.value='grey92')+coord_equal()+labs(x=NULL,y=NULL,fill='Median rho',title='Consensus co-variation structure')+theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1,size=6),axis.text.y=element_text(size=6))
rob <- copy(master);rob[,module_f:=factor(module,levels=rev(mods))]
pB2 <- ggplot(rob,aes(y=module_f))+geom_vline(xintercept=.5,lty=3,color='grey45')+geom_segment(aes(x=min_LOSO_Jaccard,xend=median_LOSO_Jaccard,yend=module_f,color=evidence_flag),linewidth=1)+geom_point(aes(x=median_LOSO_Jaccard,color=evidence_flag,size=n_representative_features))+geom_point(aes(x=min_alt_cut_jaccard),shape=17,size=1.7,color='black')+scale_color_manual(values=cols_flag)+scale_size(range=c(2,5))+scale_x_continuous(limits=c(0,1))+labs(x='Robustness / Jaccard',y=NULL,color=NULL,size='Features',title='Leave-one-sample-out robustness',subtitle='Circle: median LOSO; line: minimum-to-median; triangle: alternative cut minimum')+theme(legend.position='bottom',axis.text.y=element_text(size=6))
pB <- pB1+pB2+plot_layout(widths=c(1,1.15))+plot_annotation(title='B | MSI-only module construction and stability')
# C: representative spatial metabolic programs
pix <- fread(file.path(OUT,'step4_module_activity_spatial_organization','pixel_module_activity_scores.csv'))
spm <- fread(file.path(OUT,'step4_module_activity_spatial_organization','module_spatial_metrics_by_sample.csv'))
m14s <- spm[module=='M14'][which.max(moran_I),sample]
examples <- data.table(sample=c('GH-1','NF-GnH-3','GH-1',m14s),module=c('M02','M03','M04','M14'),type=c('Recurrent domain program','Locally reconfigured','Common-domain program','Trans-domain / continuous'))
pl <- rbindlist(lapply(seq_len(nrow(examples)),function(i){ss<-examples$sample[i];mm<-examples$module[i];d<-pix[sample_id==ss,.(x,y,score=get(mm))];sc<-sd(d$score,na.rm=TRUE);if(!is.finite(sc)||sc<1e-8)sc<-1;d[,score_z:=pmax(-2.5,pmin(2.5,(score-mean(score,na.rm=TRUE))/sc))];d[,panel:=paste0(ss,' | ',mm,'\n',examples$type[i])];d}))
pl[,panel:=factor(panel,levels=unique(panel))]
pC <- ggplot(pl,aes(x,y,color=score_z))+geom_point(shape=15,size=.55)+facet_wrap(~panel,nrow=1)+coord_equal()+scale_y_reverse()+scale_color_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-2.5,2.5),name='Module z')+labs(x=NULL,y=NULL,title='C | Representative MSI spatial metabolic programs')+theme(axis.text=element_blank(),axis.ticks=element_blank(),strip.text=element_text(size=7),legend.position='right')
# D: cross-sample spatial organization
spm[,sample:=factor(sample,levels=rev(sample_order))];spm[,module:=factor(module,levels=mods)]
pD <- ggplot(spm,aes(module,sample,fill=moran_I))+geom_tile(color='white',linewidth=.15)+scale_fill_viridis_c(option='C',limits=c(0,max(spm$moran_I,na.rm=TRUE)),oob=scales::squish)+labs(x=NULL,y=NULL,fill="Moran's I",title='D | Cross-sample spatial organization',subtitle='MSI module scores only')+theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1,size=6),axis.text.y=element_text(size=7))
# E: validation against sample-specific and common domains
rel <- fread(file.path(OUT,'step5B_sample_specific_domain_relation','sample_module_specific_domain_relations.csv'))
rel[,sample:=factor(sample,levels=rev(sample_order))];rel[,module:=factor(module,levels=mods)]
pE1 <- ggplot(rel,aes(module,sample,fill=eta2_validation_stable))+geom_tile(color='white',linewidth=.15)+geom_point(data=rel[is.finite(eta2_validation_stable)&eta2_validation_stable>=.25],shape=21,size=1,fill='black',color='black')+scale_fill_gradient(low='white',high='#B2182B',limits=c(0,.62),oob=scales::squish,na.value='grey85')+labs(x=NULL,y=NULL,fill=expression(eta^2),title='Sample-specific domains',subtitle='Validation-only stable cores; dot: eta² >= 0.25')+theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1,size=6),axis.text.y=element_text(size=7))
es <- master[is.finite(median_common_eta2_validation_only)&is.finite(median_sample_specific_eta2_validation)]
labelmods <- c('M01','M02','M03','M04','M05','M07','M08','M12','M14')
pE2 <- ggplot(es,aes(median_common_eta2_validation_only,median_sample_specific_eta2_validation))+geom_abline(slope=1,intercept=0,lty=2,color='grey50')+geom_point(aes(size=n_samples_evaluable,color=evidence_flag),alpha=.9)+ggrepel::geom_text_repel(data=es[as.character(module)%in%labelmods],aes(label=module),size=2.6,max.overlaps=Inf,box.padding=.25)+scale_color_manual(values=cols_flag)+scale_size(range=c(2.5,6))+coord_equal(xlim=c(0,.48),ylim=c(0,.48))+labs(x='Common-domain median eta²',y='Sample-specific median eta²',color=NULL,size='Evaluable\nsamples',title='Common versus local domains',subtitle='Above diagonal: locally reconfigured')+theme(legend.position='bottom')
pE <- pE1+pE2+plot_layout(widths=c(1.3,1))+plot_annotation(title='E | Independent domain validation of frozen MSI modules')
# F: LC-MS/MS post-hoc orthogonal support
lc <- copy(master);lc[,fragment_fraction:=fifelse(n_representative_families>0,n_families_with_fragment_supported/n_representative_families,NA_real_)]
lf <- melt(lc,id.vars='module',measure.vars=c('fraction_families_annotated','fragment_fraction','class_profile_spearman','n_unique_candidate_metabolite_names'),variable.name='metric',value.name='value')
metric_labels <- c(fraction_families_annotated='Candidate coverage',fragment_fraction='Fragment-supported fraction',class_profile_spearman='LC-MSI class concordance',n_unique_candidate_metabolite_names='Unique candidate metabolites')
lf[,metric_label:=factor(metric_labels[as.character(metric)],levels=unname(metric_labels))]
lf[,score01:={v<-value;if(all(!is.finite(v)))rep(NA_real_,.N) else {rg<-range(v,na.rm=TRUE);if(diff(rg)<1e-10)rep(.5,.N) else (v-rg[1])/diff(rg)}},by=metric]
lf[,module:=factor(module,levels=rev(mods))]
pF <- ggplot(lf[is.finite(score01)],aes(metric_label,module,size=score01,fill=score01))+geom_point(shape=21,color='grey20')+scale_size(range=c(1.5,7),limits=c(0,1))+scale_fill_viridis_c(option='B',limits=c(0,1))+labs(x=NULL,y=NULL,size='Within-metric\nrelative support',fill='Relative\nsupport',title='F | LC-MS/MS orthogonal support',subtitle='Post-hoc annotation and class-profile concordance only')+theme(axis.text.x=element_text(angle=30,hjust=1,size=7),axis.text.y=element_text(size=6),legend.position='bottom')
# G: integrated evidence landscape
ge <- master[,.(module,Module_cohesion=median_within_consensus_rho,LOSO_stability=median_LOSO_Jaccard,Spatial_autocorrelation=median_moran_I,Sample_specific_domain=median_sample_specific_eta2_validation,Common_domain=median_common_eta2_validation_only,LC_candidate_coverage=fraction_families_annotated,Fragment_support=fifelse(n_representative_families>0,n_families_with_fragment_supported/n_representative_families,NA_real_),LC_MSI_concordance=class_profile_spearman,evidence_flag)]
gl <- melt(ge,id.vars=c('module','evidence_flag'),variable.name='evidence',value.name='raw')
gl[,scaled:={v<-raw;if(all(!is.finite(v)))rep(NA_real_,.N) else {rg<-range(v,na.rm=TRUE);if(diff(rg)<1e-10)rep(.5,.N) else (v-rg[1])/diff(rg)}},by=evidence]
ord <- gl[,.(combined=mean(scaled,na.rm=TRUE)),by=module][order(combined),module]
gl[,module:=factor(module,levels=ord)]
gl[,evidence:=factor(evidence,levels=c('Module_cohesion','LOSO_stability','Spatial_autocorrelation','Sample_specific_domain','Common_domain','LC_candidate_coverage','Fragment_support','LC_MSI_concordance'))]
pG <- ggplot(gl,aes(evidence,module,fill=scaled))+geom_tile(color='white',linewidth=.2)+scale_fill_viridis_c(option='C',limits=c(0,1),na.value='grey85')+labs(x=NULL,y=NULL,fill='Within-column\nrelative evidence',title='G | Integrated module evidence landscape',subtitle='Column-wise scaling; domain and LC-MS/MS fields are post-hoc')+theme(axis.text.x=element_text(angle=40,hjust=1,size=7),axis.text.y=element_text(size=6))
# Save individual panels
save_both(pB,'Panel_B_module_structure_stability',12,5.2)
save_both(pC,'Panel_C_representative_spatial_programs',13,4.2)
save_both(pD,'Panel_D_cross_sample_spatial_organization',8.2,5.4)
save_both(pE,'Panel_E_domain_validation',13,5.4)
save_both(pF,'Panel_F_LCMS_support',8.5,7)
save_both(pG,'Panel_G_integrated_evidence',9.5,7)
main <- (pB|pC)/(pD|pE)/(pF|pG)+plot_layout(heights=c(1,1,1.15),widths=c(1,1.55))+plot_annotation(title='Spatial metabolic programs in pituitary tumors',subtitle='MSI-derived frozen modules with spatial, domain and LC-MS/MS post-hoc support',theme=theme(plot.title=element_text(face='bold',size=16),plot.subtitle=element_text(size=10,color='grey30')))
ggsave(file.path(FD,'Main_Figure_modules_spatial_domain_LCMS.pdf'),main,width=22,height=18,units='in',device=cairo_pdf)
ggsave(file.path(FD,'Main_Figure_modules_spatial_domain_LCMS.png'),main,width=22,height=18,units='in',dpi=300,bg='white')
ggsave(file.path(FD,'Main_Figure_modules_spatial_domain_LCMS.tiff'),main,width=22,height=18,units='in',dpi=300,compression='lzw',bg='white')
writeLines(c('Main figure excludes the study-design panel.','B: MSI-only consensus modules and robustness.','C: representative module activity maps.','D: cross-sample spatial autocorrelation.','E: independent sample-specific and common-domain validation.','F: LC-MS/MS post-hoc candidate, fragment and concordance evidence.','G: integrated evidence landscape.','Domain and LC-MS/MS information never changes module discovery, membership or MSI spatial metrics.'),file.path(FD,'README.txt'))
file.create(file.path(FD,'FIGURE_COMPLETE.ok'))
cat('FIGURE_COMPLETE\n',FD,'\n',paste(list.files(FD),collapse='\n'))
})
