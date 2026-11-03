# Figure ABC: MSI only. Current shared62 is a frozen historical cross-platform list.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
BASE <- '/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
OUT <- file.path(BASE,'figure_ABC_MSI/20260914_v2');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
groups <- c('GH','NF-GH','NF-ACTH','Cushing','NF-GnH')
map_samples <- c('GH-1','NF-GH','NF-ACTH','Cushing-1','NF-GnH-1')
all_samples <- c('GH-1','GH-2','GH-3','NF-GH','NF-ACTH','Cushing-1','Cushing-2','NF-GnH-1','NF-GnH-2','NF-GnH-3')
labels <- sub('-[0-9]+$','',all_samples)
A_ids <- c(12,38,63,100,121,143,173,205,227,262)
BC_ids <- c(158,210,117,29,1,95,4,62,106,20,71,98,260,261,251,17,208,262,209,207,10,18,3,16,12,42,2,197,5,196,204,150,229,242,11,231,235,47,234,241,37,81,92,91,132,139,35,90)
BC_class <- c("GH","GH","GH","GH","GH","GH","GH","GH","GH","GH","NF-GH","NF-GH","NF-GH","NF-GH","NF-GH","NF-GH","NF-GH","NF-GH","NF-GH","NF-GH","NF-ACTH","NF-ACTH","NF-ACTH","NF-ACTH","NF-ACTH","NF-ACTH","NF-ACTH","NF-ACTH","NF-ACTH","NF-ACTH","Cushing","Cushing","Cushing","Cushing","Cushing","Cushing","Cushing","Cushing","Cushing","Cushing","NF-GnH","NF-GnH","NF-GnH","NF-GnH","NF-GnH","NF-GnH","NF-GnH","NF-GnH")
shared_ids <- c(12,19,21,23,24,28,30,31,32,37,38,48,54,63,64,65,66,74,75,76,78,83,88,89,90,93,94,96,97,99,100,103,104,105,108,109,110,121,122,123,128,129,130,133,134,135,136,141,142,143,173,174,175,176,205,206,207,208,227,246,261,262)

mat <- readRDS(file.path(BASE,'pixel_by_feature_matrix_aligned_TIC.rds'))
mz <- as.numeric(readRDS(file.path(BASE,'feature_mz_aligned_filtered.rds')))
meta <- as.data.table(readRDS(file.path(BASE,'pixel_metadata_sample_source_xy.rds')))
tumor <- as.data.table(readRDS(file.path(BASE,'intratumor_spatial_heterogeneity_10samples_part1_consensus_tuned225/tumor_pixels_within_sample_consensus_domains.rds')))
ix <- match(tumor$original_pixel_row,meta$pixel_row)
stopifnot(!anyNA(ix),all(tumor$is_tumor),!anyDuplicated(ix),ncol(mat)==274,all(meta$sample_id[ix]==tumor$sample_id),all(meta$x[ix]==tumor$x),all(meta$y[ix]==tumor$y))
need <- sort(unique(c(A_ids,BC_ids,shared_ids)))
x <- as.matrix(mat[ix,need,drop=FALSE]);colnames(x)<-as.character(need);rm(mat);gc()
pix <- meta[ix];pix[,class:=sub('-[0-9]+$','',sample_id)]
stopifnot(nrow(x)==46574,all(is.finite(x)),all(x>=0),setequal(unique(pix$sample_id),all_samples))

# Select from the frozen shared62 by robust signal in all FIVE displayed samples.
map_rows<-lapply(map_samples,function(s)which(pix$sample_id==s))
aq<-rbindlist(lapply(shared_ids,function(fid){
 v<-x[,as.character(fid)];q95<-as.numeric(quantile(v[unlist(map_rows)],.95))
 ss<-sapply(map_rows,function(j)c(mean=mean(v[j]),median=median(v[j]),q25=as.numeric(quantile(v[j],.25)),detection=mean(v[j]>0)))
 data.table(feature_index=fid,mz=mz[fid],min_mean=min(ss['mean',]),min_median=min(ss['median',]),min_q25_rel=min(ss['q25',])/q95,min_median_rel=min(ss['median',])/q95,min_detection=min(ss['detection',]),fold_means=max(ss['mean',])/min(ss['mean',]))
}))
aq[,eligible:=min_detection>=.95 & min_median_rel>=.50]
aq[,strength_score:=min_median*min_median_rel]
aq[,mz_bin_50Da:=floor(mz/50)]
pool<-aq[eligible==TRUE][order(-strength_score,mz)];selected<-integer(0)
for(b in sort(unique(pool$mz_bin_50Da))){
 opts<-pool[mz_bin_50Da==b]
 for(fid in opts$feature_index)if(!length(selected)||all(abs(mz[fid]-mz[selected])>=2)){selected<-c(selected,fid);break}
}
for(fid in pool$feature_index){
 if(all(abs(mz[fid]-mz[selected])>=2))selected<-c(selected,fid)
 if(length(selected)==20)break
}
stopifnot(length(selected)==20,!anyDuplicated(selected))
A_ids<-selected[order(mz[selected])];aq[,selected_A20:=feature_index %in% A_ids]
fwrite(aq[order(mz)],file.path(OUT,'A_shared62_selection_screen.csv'))
fwrite(aq[selected_A20==TRUE][order(mz)],file.path(OUT,'A_selected_shared20.csv'))
fwrite(data.table(class=BC_class,feature_index=BC_ids,mz=mz[BC_ids]),file.path(OUT,'C_specific48_unchanged.csv'))
stat_ids<-sort(unique(c(A_ids,BC_ids)))
per_sample<-rbindlist(lapply(stat_ids,function(fid)rbindlist(lapply(seq_along(all_samples),function(i){
 v<-x[pix$sample_id==all_samples[i],as.character(fid)]
 data.table(feature_index=fid,mz=mz[fid],sample_id=all_samples[i],class=labels[i],n_pixels=length(v),mean_TIC=mean(v),median_TIC=median(v),Q25_TIC=as.numeric(quantile(v,.25)),Q75_TIC=as.numeric(quantile(v,.75)),detection_fraction=mean(v>0))
}))))
fwrite(per_sample,file.path(OUT,'AC_per_sample_signal_statistics.csv'))
cat('A20_SELECTED\n');print(aq[selected_A20==TRUE,.(feature_index,mz,min_mean,min_detection,min_median_rel)][order(mz)])
box_stats<-list();k<-0L
boxpage <- function(fid) {
 v<-log1p(x[,as.character(fid)]);yr<-range(v);pad<-max(diff(yr)*.07,1e-8)
 par(mar=c(0,0,0,0),fig=c(0,1,0,1),new=FALSE,bg='white');plot.new()
 text(.08,.955,sprintf('A | Shared ion f%d | m/z %.4f',fid,mz[fid]),adj=0,cex=1.25)
 text(.08,.915,'All tumor pixels pooled within subtype; no trimming or subsampling',adj=0,cex=.78)
 par(fig=c(.08,.92,.20,.86),new=TRUE,mar=c(0,0,0,0),xaxs='i',yaxs='i')
 plot(NA,xlim=c(.5,5.5),ylim=c(max(0,yr[1]-pad),yr[2]+pad),axes=FALSE,xlab='',ylab='')
 abline(h=axTicks(2),col='#E6E6E6');abline(v=1:5,col='#EFEFEF')
 set.seed(20260914+fid)
 for(i in 1:5){u<-v[pix$class==groups[i]];q<-quantile(u,c(.25,.5,.75));mu<-mean(u);ss<-sd(u);iqr<-q[3]-q[1];lo<-min(u[u>=q[1]-1.5*iqr]);hi<-max(u[u<=q[3]+1.5*iqr]);xb<-i-.18;xp<-i+.16
  points(xp+runif(length(u),-.09,.09),u,pch=16,cex=.15,col=ifelse(u>mu+ss,adjustcolor('#D7301F',.5),adjustcolor('#2C7FB8',.4)))
  rect(xb-.10,q[1],xb+.10,q[3],col='white',border='#333333',lwd=.8)
  segments(xb-.10,mu,xb+.10,mu,lwd=1.2)
  segments(xb,lo,xb,q[1]);segments(xb,q[3],xb,hi);segments(xb-.05,lo,xb+.05,lo);segments(xb-.05,hi,xb+.05,hi)
  segments(xp,max(0,mu-ss),xp,mu+ss,lwd=.7);points(xp,mu,pch=23,bg='black',col='white',cex=.75)
  k<<-k+1L;box_stats[[k]]<<-data.table(feature_index=fid,class=groups[i],n_pixels=length(u),mean=mu,median=q[2],Q25=q[1],Q75=q[3],sd=ss,mean_outside_IQR=mu<q[1]||mu>q[3])
 }
 axis(1,at=1:5,labels=groups,tick=FALSE,line=.3,cex.axis=.95);axis(2,las=1,cex.axis=.8);mtext('log1p(TIC-normalized pixel intensity)',2,line=3.4,cex=.85)
 par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0));plot.window(c(0,1),c(0,1))
 text(.08,.06,'Box = Q25-Q75; horizontal line = mean; diamond/bar = mean +/- SD; red pixels > mean + SD.',adj=0,cex=.70)
}
# Coordinates centered inside equal-size canvases. Common scale and pixel size across samples.
map_rows<-lapply(map_samples,function(s)which(pix$sample_id==s))
stopifnot(all(lengths(map_rows)>0))
coords<-lapply(map_rows,function(j){p<-pix[j];stopifnot(uniqueN(p$source_id)==1,!anyDuplicated(paste(p$x,p$y)));data.frame(x=as.integer(p$x-min(p$x)),y=as.integer(p$y-min(p$y)))})
W<-max(vapply(coords,function(p)max(p$x)+1,numeric(1)));H<-max(vapply(coords,function(p)max(p$y)+1,numeric(1)))
for(i in 1:5){coords[[i]]$x<-coords[[i]]$x+floor((W-max(coords[[i]]$x)-1)/2);coords[[i]]$y<-coords[[i]]$y+floor((H-max(coords[[i]]$y)-1)/2)}
palette_spatial<-colorRampPalette(c('blue','cyan','green','yellow','magenta','red'))(256)

spatial_stats<-list();sk<-0L
spatialpage<-function(fid,panel='C',target=''){
 rawv<-x[,as.character(fid)]
 vals<-if(panel=='A')log1p(rawv)else rawv*1e6
 lower<-if(panel=='A')0 else as.numeric(quantile(vals[unlist(map_rows)],.50))
 upper<-as.numeric(quantile(vals[unlist(map_rows)],.95));stopifnot(is.finite(upper),upper>lower)
 mean5<-vapply(map_rows,function(j)mean(rawv[j])*1e6,numeric(1))
 target_i<-match(target,groups)
 fc<-if(panel=='C')mean5[target_i]/max(mean5[-target_i])else NA_real_
 par(fig=c(0,1,0,1),new=FALSE,mar=c(0,0,0,0),bg='white',xaxs='i',yaxs='i')
 plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab='',ylab='')
 text(.08,.95,sprintf('%s | %s f%d | m/z %.4f',panel,if(panel=='A')'Shared'else paste0(target,'-enriched'),fid,mz[fid]),adj=0,cex=1.15)
 if(panel=='C')text(.08,.89,sprintf('Selected target sample mean / highest other sample mean = %.2fx%s',fc,if(fc<=1)' (target sample is not the highest)'else ''),adj=0,cex=.78)
 iw<-min(.160,.50*3.8/12*W/H);ih<-iw*12/3.8*H/W
 for(i in 1:5){cx<-.164+(i-1)*.168
  rect(cx-.080,.25,cx+.080,.75,col='black',border=NA)
  pp<-coords[[i]];v<-pmax(pmin((vals[map_rows[[i]]]-lower)/(upper-lower),1),0)
  img<-matrix('#000000',H,W);img[cbind(pp$y+1,pp$x+1)]<-palette_spatial[1+floor(v*255)]
  rasterImage(as.raster(img),cx-iw/2,.50-ih/2,cx+iw/2,.50+ih/2,interpolate=FALSE)
  text(cx,.815,map_samples[i],cex=.92,font=if(panel=='C'&&i==target_i)2 else 1)
  if(panel=='C')text(cx,.215,sprintf('Mean: %.0f',mean5[i]),cex=.78)
  sk<<-sk+1L;spatial_stats[[sk]]<<-data.table(panel=panel,target=target,feature_index=fid,mz=mz[fid],sample_id=map_samples[i],lower=lower,upper=upper,scale=if(panel=='A')'0-Q95 log1p(TIC)'else 'Q50-Q95 linear TICx1e6',fraction_at_or_below_lower=mean(vals[map_rows[[i]]]<=lower),fraction_at_or_above_upper=mean(vals[map_rows[[i]]]>=upper),mean_TICx1e6=mean5[i],target_vs_max_other_fold=fc)
 }
 rasterImage(as.raster(matrix(palette_spatial,nrow=1)),.36,.12,.64,.15,interpolate=FALSE)
 text(.36,.09,if(panel=='A')'0'else sprintf('<= Q50: %.0f',lower),cex=.73)
 text(.64,.09,if(panel=='A')sprintf('>= Q95: %.4g',upper)else sprintf('>= Q95: %.0f',upper),cex=.73)
 text(.5,.035,if(panel=='A')'log1p(TIC signal); shared 0-Q95 scale across five samples'else 'TIC signal x 1e6; shared Q50-Q95 contrast window; values outside limits are clipped',cex=.73)
}
pdf(file.path(OUT,'A_common20_pixel_boxplots.pdf'),width=12,height=5.4,useDingbats=FALSE)
for(fid in A_ids)boxpage(fid)
dev.off();fwrite(rbindlist(box_stats),file.path(OUT,'A_pixel_box_statistics.csv'))
pdf(file.path(OUT,'A_common20_spatial_maps.pdf'),width=12,height=3.8,useDingbats=FALSE)
for(fid in A_ids)spatialpage(fid,'A','Shared')
dev.off()
for(g in groups){
 pdf(file.path(OUT,paste0('C_specific_',g,'_spatial_Q50_Q95.pdf')),width=12,height=3.8,useDingbats=FALSE)
 for(fid in BC_ids[BC_class==g])spatialpage(fid,'C',g)
 dev.off()
}
fwrite(rbindlist(spatial_stats),file.path(OUT,'AC_spatial_scales_and_clipping.csv'))
contrast<-rbindlist(spatial_stats)[panel=='C',.(target_vs_max_other_fold=first(target_vs_max_other_fold)),by=.(target,feature_index,mz)]
contrast[,display_priority:=fifelse(target_vs_max_other_fold>=1.5,'stronger_selected_sample_contrast',fifelse(target_vs_max_other_fold>=1.25,'moderate_selected_sample_contrast',fifelse(target_vs_max_other_fold>1,'small_selected_sample_difference','selected_target_not_highest')))]
fwrite(contrast,file.path(OUT,'C_actual_five_sample_contrast.csv'))
writeLines(c(
'A20 selection from fixed shared62: all five displayed samples detection >=95%, and each sample median >=50% of five-sample pooled Q95 (linear TIC signal).',
'Strength score = minimum five-sample median * minimum five-sample median / pooled Q95. Prioritize the weakest sample signal, not the strongest sample.',
'Ensure one eligible feature per occupied 50-Da bin, then fill by decreasing score with at least 2-Da display spacing to 20. This spacing is not chemical de-isotoping or a metabolite identity claim.',
'A pixel boxes still pool all 10 physical samples by subtype, every tumor pixel retained; Q25-Q75 box, mean line, mean +/- SD.',
'A maps: unchanged 0-Q95 log1p(TIC) display scale, shared across five samples.',
'C: same 48 features and order as B; B itself is unchanged. Use pooled five-sample Q50-Q95 linear TIC x1e6 display window for every feature, identical mapping across samples; below/above endpoints clipped to endpoint colors.',
'The user reference 3C_specific_feature_plot.R uses min to min+0.8*(max-min). Here robust pooled quantiles replace the extreme-based endpoints and raise the floor to Q50 to emphasize high signals. This is a contrast adjustment, not a change in measured abundance.',
'C page reports actual all-pixel sample means and target mean / maximum other mean, unaffected by clipping. Weak or reversed contrasts are retained and identified, not made stronger through sample-specific scaling.',
'Five sample order: GH-1, NF-GH, NF-ACTH, Cushing-1, NF-GnH-1. Equal spacing and common coordinate extents; no rotations.',
'Input is already TIC normalized. No second TIC normalization. Source base is MSI_ALL_aligned_10samples_consensus_tuned225.',
'Cushing and NF-GnH keep relaxed exploratory feature lists; NF-GH and NF-ACTH have a single physical sample. Spatial appearance is descriptive.'
),file.path(OUT,'METHODS_v2.txt'))
file.copy('/home/gaga/huzhan/MSI_NEW/script/figure_AC_MSI_v2.R',file.path(OUT,'figure_AC_MSI_v2.R'),overwrite=TRUE)
oldwd<-getwd();setwd(OUT);utils::zip('Figure_AC_MSI_v2.zip',files=list.files()[!grepl('zip$',list.files())],flags='-q');setwd(oldwd)
cat('AC_V2_COMPLETE\n',OUT,'\n');print(contrast[,.(n=.N),by=.(target,display_priority)]);rstudioapi::filesPaneNavigate(OUT)
