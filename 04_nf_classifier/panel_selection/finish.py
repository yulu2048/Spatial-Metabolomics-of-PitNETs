from pathlib import Path
import sys,json,gzip
import numpy as np,pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;R=H.parent;P=R/'result3_endogenous_constrained_20260915';sys.path.insert(0,str(P));import analysis as old
rm=pd.read_csv(P/'geometry/ROI_manifest.csv');train=rm[rm.width_pixels.eq(10)];x=np.load(P/'cache/ROI_arrays.npz')['response'];inc=sparse.load_npz(P/'cache/ROI_incidence.npz');sel=json.loads((H/'final_selection.json').read_text());assoc=pd.read_csv(H/'final_association.csv');allfid=sorted(assoc.feature_id,key=lambda f:int(f[1:]));fidmap={f:i for i,f in enumerate(allfid)};X=x[:,[int(f[1:])-1 for f in allfid]];checks=[]
for f in range(4):
 te=rm[rm.fold.eq(f)];tr=train[~train.fold.eq(f)];tr=tr[old.purge(tr,te,2)]
 shared=np.intersect1d(inc[tr.index].indices,inc[te.index].indices).size;assert shared==0
 checks.append(dict(check='outer_shared_pixels',fold=f,value=shared,passed=True))
 for ff in sorted(tr.fold.unique()):
  it=tr[tr.fold.eq(ff)];itr=tr[~tr.fold.eq(ff)];itr=itr[old.purge(itr,it,2)]
  shared=np.intersect1d(inc[itr.index].indices,inc[it.index].indices).size;assert shared==0;checks.append(dict(check='inner_shared_pixels',fold=f,inner_fold=ff,value=shared,passed=True))
for k in [5]:
 m=np.load(H/f'fixed_{k}_MSI_model.npz');ids=[fidmap[f] for f in m['features']];co,it=old.fit(old.stats(X[train.index],train),{'indices':ids,'lambda':.75});np.testing.assert_allclose(X[:,ids]@m['coefficient']+m['intercept'],X@co+it,rtol=1e-10,atol=1e-9);checks.append(dict(check='exported_model_score_equivalence',k=k,passed=True))
nd=pd.read_csv(H/'nested_predictions.csv');assert len(nd)==len(rm) and not nd.duplicated(['k','roi_id']).any();pd.DataFrame(checks).to_csv(H/'verification.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
SA=old.SA;colors=['#B04478','#367AB0','#268E79','#268E79','#268E79'];selected=sorted(set(sel['5']['features']),key=lambda f:int(f[1:]))
fig,axes=plt.subplots(2,3,figsize=(15,9));fig.subplots_adjust(left=.065,right=.98,bottom=.14,top=.85,wspace=.30,hspace=.50)
for ax,f in zip(axes.flat,selected):
 vals=[x[train[train.sample_id.eq(s)].index,int(f[1:])-1] for s in SA];vi=ax.violinplot(vals,positions=range(5),showextrema=False,widths=.75,points=250)
 for i,(b,v) in enumerate(zip(vi['bodies'],vals)):
  b.set_facecolor(colors[i]);b.set_edgecolor(colors[i]);b.set_alpha(.4);q1,med,q3=np.quantile(v,[.25,.5,.75]);ax.plot([i,i],[q1,q3],color=colors[i],lw=4);ax.scatter(i,med,s=20,facecolors='white',edgecolors=colors[i],zorder=3);ax.scatter(i,np.mean(v),marker='D',s=15,c='#17354A',zorder=4)
 r=assoc.set_index('feature_id').loc[f];which='/'.join(str(k) for k in [5] if f in sel[str(k)]['features']);ax.set_title(f'{f} | {r.target_class} {r.direction} | {which}-feature panel',fontsize=11,loc='left');ax.set_xticks(range(5),['GH','ACTH','GnH-1','GnH-2','GnH-3']);ax.set_ylabel('log2(1 + ROI mean TIC x 10^6)');ax.grid(axis='y',alpha=.15);ax.margins(y=.10)
fig.suptitle('Five-feature MSI panel',fontsize=18,weight='bold',x=.065,ha='left',y=.97)
fig.text(.065,.91,'One violin per specimen | 10 x 10-pixel ROIs | GnH-3 pieces pooled | Diamond: mean; circle: median; line: IQR',color='#555')
fig.text(.065,.045,'LC abundance and LC group trends are not used. Candidate chemistry remains putative.\nROIs within a specimen overlap; these distributions do not constitute independent patient replicates.',fontsize=10,color='#555',linespacing=1.6)
axes.flat[-1].axis('off')
fig.savefig(H/'selected_feature_ROI_distributions.png',dpi=180);fig.savefig(H/'selected_feature_ROI_distributions.svg');plt.close(fig)
a=assoc[assoc.shortlisted].copy();mat=a[['ROI_mean_'+s for s in SA]].to_numpy();z=(mat-mat.mean(1,keepdims=True))/mat.std(1,keepdims=True)
fig,ax=plt.subplots(figsize=(10,10));fig.subplots_adjust(left=.35,right=.86,top=.90,bottom=.12);im=ax.imshow(z,cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto');ax.set_xticks(range(5),SA,rotation=20,ha='right');ax.set_yticks(range(len(a)),[f'{r.feature_id} | {r.target_class} {r.direction}' for _,r in a.iterrows()]);ax.tick_params(length=0)
for i,j in np.ndindex(mat.shape):ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white' if abs(z[i,j])>1.4 else '#17354A')
fig.colorbar(im,ax=ax,fraction=.04,pad=.05,label='Row-standardized specimen mean');fig.suptitle(f'{len(a)} eligible features',fontsize=15,y=.96,weight='bold');fig.text(.06,.025,'Numbers: specimen mean log2 ROI signal. Color: row z score. Ranking uses the worst specimen-pair AUC.\nScreen: worst AUC >= 0.75 and minimum signed specimen mean gap >= 0.25; no LC expression criterion.',fontsize=9,linespacing=1.6)
fig.savefig(H/'eligible18_specimen_heatmap.png',dpi=170);fig.savefig(H/'eligible18_specimen_heatmap.svg');plt.close(fig)
# Additional pixel-level descriptive means: distinct from ROI-mean log signals.
pix=pd.read_csv(P/'geometry/tumor_pixel_metadata.csv');raw=np.frombuffer(gzip.decompress((P/'input/tumor_matrix.bin.gz').read_bytes()),dtype='<f8').reshape((len(pix),274),order='F');rows=[]
for f in selected:
 for s in SA:
  v=np.log2(1+raw[pix.sample_id.eq(s),int(f[1:])-1]*1e6);rows.append(dict(feature_id=f,sample_id=s,n_pixels=len(v),pixel_log_mean=v.mean(),pixel_log_median=np.median(v),zero_fraction=(v==0).mean()))
pd.DataFrame(rows).to_csv(H/'selected_pixel_summary.csv',index=False)
# Full candidate table includes all competing endogenous annotations; never choose by LC performance.
pool=pd.read_csv(R/'result3_joint55_20260915/candidate_pool.csv');assoc.merge(pool,on='feature_id').to_csv(H/'all55_association_and_candidate_annotations.csv',index=False)
print('Validation checks passed; figures and complete audit tables created.')
