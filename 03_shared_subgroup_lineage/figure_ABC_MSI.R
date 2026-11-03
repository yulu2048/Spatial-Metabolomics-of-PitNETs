# Figure ABC: MSI only. Current shared62 is a frozen historical cross-platform list.
suppressPackageStartupMessages({library(data.table);library(Matrix)})
BASE <- '/home/gaga/huzhan/MSI_ALL_aligned_10samples_consensus_tuned225'
OUT <- file.path(BASE,'figure_ABC_MSI/20260914');dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
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
fwrite(data.table(feature_index=A_ids,mz=mz[A_ids]),file.path(OUT,'A_selected_shared10.csv'))
fwrite(data.table(class=BC_class,feature_index=BC_ids,mz=mz[BC_ids]),file.path(OUT,'BC_selected_specific48.csv'))
fwrite(pix[,.(n_pixels=.N),by=.(class,sample_id)],file.path(OUT,'pixel_counts.csv'))
sample_means<-sapply(all_samples,function(s)colMeans(x[pix$sample_id==s,,drop=FALSE]))
rownames(sample_means)<-as.character(need)
group_means<-sapply(groups,function(g)rowMeans(sample_means[,labels==g,drop=FALSE]))
rownames(group_means)<-as.character(need)
fwrite(data.table(feature_index=need,mz=mz[need],as.data.table(sample_means)),file.path(OUT,'all_selected_sample_linear_TIC_means.csv'))
# Box and scatter: 5 equal-width class columns, pooling ALL tumor pixels within class.
# Q25-Q75 box with MEAN line (requested), Tukey whiskers computed from quartiles.
# Separate diamond and bar: mean +/- SD in log1p(TIC intensity) display space.
box_stats<-list();k<-0L
boxpage <- function(fid) {
 v<-log1p(x[,as.character(fid)]);yr<-range(v);pad<-max(diff(yr)*.07,1e-8)
 par(mar=c(0,0,0,0),fig=c(0,1,0,1));plot.new()
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
spatial_scale<-list();sk<-0L
spatialpage<-function(fid,panel='C',target=''){
 vals<-log1p(x[,as.character(fid)]);q95<-as.numeric(quantile(vals[unlist(map_rows)],.95));stopifnot(is.finite(q95),q95>0)
 par(fig=c(0,1,0,1),new=FALSE,mar=c(0,0,0,0),bg='white',xaxs='i',yaxs='i')
 plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab='',ylab='')
 text(.08,.93,sprintf('%s | %s f%d | m/z %.4f',panel,target,fid,mz[fid]),adj=0,cex=1.2)
 iw<-min(.160,.50*3.8/12*W/H);ih<-iw*12/3.8*H/W
 for(i in 1:5){cx<-.164+(i-1)*.168
  rect(cx-.080,.25,cx+.080,.75,col='black',border=NA)
  pp<-coords[[i]];v<-pmin(vals[map_rows[[i]]]/q95,1);img<-matrix('#000000',H,W);img[cbind(pp$y+1,pp$x+1)]<-palette_spatial[1+floor(v*255)]
  rasterImage(as.raster(img),cx-iw/2,.50-ih/2,cx+iw/2,.50+ih/2,interpolate=FALSE)
  text(cx,.815,map_samples[i],cex=.92)
 }
 rasterImage(as.raster(matrix(palette_spatial,nrow=1)),.36,.13,.64,.16,interpolate=FALSE)
 text(.36,.10,'0',cex=.75);text(.64,.10,sprintf('Q95 = %.4g',q95),cex=.75)
 text(.5,.045,'log1p(TIC intensity); one shared 0-Q95 scale per ion across the five samples',cex=.75)
 sk<<-sk+1L;spatial_scale[[sk]]<<-data.table(panel=panel,target=target,feature_index=fid,Q95=q95)
}
pdf(file.path(OUT,'A_common10_pixel_boxplots.pdf'),width=12,height=5.4,useDingbats=FALSE)
for(fid in A_ids)boxpage(fid)
dev.off();fwrite(rbindlist(box_stats),file.path(OUT,'A_pixel_box_statistics.csv'))
pdf(file.path(OUT,'A_common10_spatial_maps.pdf'),width=12,height=3.8,useDingbats=FALSE)
for(fid in A_ids)spatialpage(fid,'A','Shared')
dev.off()
# B has 5 subtype rows and 48 feature columns, subtype means weight physical samples equally.
raw<-t(group_means[as.character(BC_ids),,drop=FALSE])*1e6
zscore<-scale(raw,center=TRUE,scale=TRUE)
fwrite(data.table(class=groups,as.data.table(raw)),file.path(OUT,'B_raw_TICx1e6_matrix.csv'))
fwrite(data.table(class=groups,as.data.table(zscore)),file.path(OUT,'B_zscore_matrix.csv'))
heatpage<-function(v,standardized=FALSE){
 lim<-if(standardized)c(-2,2)else c(0,max(v));pal<-if(standardized)colorRampPalette(c('#2166AC','white','#B2182B'))(256) else hcl.colors(256,'YlOrRd',rev=TRUE)
 par(fig=c(0,1,0,1),new=FALSE,mar=c(0,0,0,0),bg='white',xaxs='i',yaxs='i')
 plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab='',ylab='')
 text(.05,.94,if(standardized)'B | Subtype-enriched ions: Z-score'else 'B | Subtype-enriched ions: linear TIC-normalized intensity',adj=0,cex=1.25)
 xx<-seq(.09,.94,length.out=length(BC_ids)+1);yy<-seq(.80,.30,length.out=6)
 for(j in seq_along(BC_ids))for(i in 1:5){t<-min(max((v[i,j]-lim[1])/diff(lim),0),1);rect(xx[j],yy[i+1],xx[j+1],yy[i],col=pal[1+floor(t*255)],border='#BBBBBB',lwd=.4)}
 for(i in 1:5)text(.078,mean(yy[c(i,i+1)]),groups[i],adj=1,cex=.8)
 for(j in seq_along(BC_ids))text(mean(xx[c(j,j+1)]),.27,sprintf('%.4f',mz[BC_ids[j]]),srt=90,adj=1,cex=.59)
 starts<-cumsum(c(0,head(as.numeric(table(factor(BC_class,levels=groups))),-1)));ends<-cumsum(as.numeric(table(factor(BC_class,levels=groups))))
 for(i in 1:5){segments(xx[starts[i]+1],.825,xx[ends[i]+1],.825,lwd=2);text(mean(c(xx[starts[i]+1],xx[ends[i]+1])),.86,groups[i],cex=.8)}
 rasterImage(as.raster(matrix(pal,nrow=1)),.33,.10,.67,.135,interpolate=FALSE)
 text(.33,.075,format(lim[1],digits=3),cex=.8);text(.67,.075,format(lim[2],digits=3),cex=.8)
 text(.5,.035,if(standardized)'Z-score across five subtype means for each feature; sample means weighted equally'else 'Mean TIC-normalized signal x 1e6; no Z-score or log transform; not absolute concentration',cex=.8)
}
for(mode in c('raw','zscore')){pdf(file.path(OUT,paste0('B_specific48_heatmap_',mode,'.pdf')),width=15,height=5.5,useDingbats=FALSE);heatpage(if(mode=='raw')raw else zscore,mode=='zscore');dev.off()}
for(g in groups){pdf(file.path(OUT,paste0('C_specific_',g,'_spatial_maps.pdf')),width=12,height=3.8,useDingbats=FALSE);for(fid in BC_ids[BC_class==g])spatialpage(fid,'C',g);dev.off()}
fwrite(rbindlist(spatial_scale),file.path(OUT,'spatial_Q95_scales.csv'))
# Supplement: fixed62, MSI only, every point is a physical sample. 4 features/page.
pdf(file.path(OUT,'S_shared62_MSI_sample_distributions.pdf'),width=12,height=9,useDingbats=FALSE)
for(start in seq(1,length(shared_ids),by=4)){par(mfrow=c(2,2),mar=c(5,5,3,1));ids<-shared_ids[start:min(start+3,length(shared_ids))];for(fid in ids){v<-sample_means[as.character(fid),];ref<-mean(group_means[as.character(fid),]);v<-v/ref;plot(NA,xlim=c(.5,5.5),ylim=c(0,max(2.5,max(v)*1.1)),xaxt='n',xlab='',ylab='Relative linear sample mean',main=sprintf('f%d | m/z %.4f',fid,mz[fid]));axis(1,1:5,groups,las=2,cex.axis=.8);abline(h=1,lty=3,col='grey');for(i in 1:5){u<-v[labels==groups[i]];xx<-if(length(u)>1)seq(i-.13,i+.13,length.out=length(u))else i;points(xx,u,pch=16,col=c('#0072B2','#009E73','#CC79A7','#D55E00','#B68B00')[i]);segments(i-.2,mean(u),i+.2,mean(u))}};if(length(ids)<4)for(i in seq_len(4-length(ids)))plot.new()}
dev.off()
file.copy('/home/gaga/huzhan/MSI_NEW/script/figure_ABC_MSI.R',file.path(OUT,'figure_ABC_MSI.R'),overwrite=TRUE)
writeLines(c('Five-class order: GH, NF-GH, NF-ACTH, Cushing, NF-GnH.','Spatial sample order: GH-1, NF-GH, NF-ACTH, Cushing-1, NF-GnH-1.','A: ten shared ions chosen to span m/z, from frozen shared62. ALL tumor pixels from all 10 samples pooled by subtype; no trimming/downsampling. Box Q25-Q75, mean line, Tukey whiskers; separate mean +/- SD, all in log1p(TIC) space.','Boxplot and spatial column centers use normalized page x coordinates 0.164,0.332,0.500,0.668,0.836 on 12-inch pages. Spatial tiles have identical widths and common coordinate extents.','Input matrix already TIC-normalized; do not normalize it again. Spatial maps include selected tumor pixels only, backgrounds black; zero signal blue.','Spatial colors 0 to per-feature Q95 of log1p(TIC intensity) over all selected five-sample tumor pixels. Same scale for five samples, clipped above Q95.','B: sample-wise all-tumor-pixel linear mean, then arithmetic mean across physical samples per subtype. Raw=linear TIC signal x1e6, not absolute concentration. Z score is across five subtype means for each ion, display clipped at +/-2.','Specific selection uses original criteria for GH/NF-GH/NF-ACTH; exploratory broad criteria for Cushing/NF-GnH. Take up to10 per subtype by AUC, standardized effect, minimum pairwise difference; NF-GnH all8. These groups do not share identical discovery thresholds.','All names are m/z only, no LC annotations. Frozen shared62 retains historical cross-platform provenance. Some subtype-enriched features overlap shared62.','Singleton NF-GH/NF-ACTH are exploratory. Pixel plots are descriptive and do not provide biological replication.','S:62 shared features, MSI sample points only, no LC panels.'),file.path(OUT,'METHODS.txt'))
files<-list.files(OUT,full.names=TRUE);files<-files[!grepl('zip$',files)];oldwd<-getwd();setwd(OUT);utils::zip('Figure_ABC_MSI_complete.zip',files=basename(files),flags='-q');setwd(oldwd)
cat('FIGURE_ABC_COMPLETE\n',OUT,'\n');print(list.files(OUT,pattern='pdf$'));rstudioapi::filesPaneNavigate(OUT)
