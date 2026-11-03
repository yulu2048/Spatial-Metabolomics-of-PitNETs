from pathlib import Path
import gzip,json,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from analyze_msi import OUT,TABLE,SAMPLES,CONTRASTS

FIG=OUT/'output/pdf';FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
    'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
COL={'functional_higher':'#bc3b37','nonfunctional_higher':'#2668a2'}
PAT={'detection_dominant':'#3a9b87','intensity_dominant':'#a65ba6','mixed':'#d5a338',
    'opposing_components':'#52657a','insufficient_positive_pixels':'#cccccc'}
D=pd.read_csv(TABLE/'MSI_differential_all.csv');C=D[D.candidate].copy()
A=np.load(OUT/'input/sample_arrays.npz'); mz=A['mz']
pix=pd.read_csv(TABLE/'tumor_pixel_metadata.csv')
x=np.frombuffer(gzip.decompress((OUT/'input/tumor_matrix.bin.gz').read_bytes()),dtype='<f8').reshape((len(pix),274),order='F')*1e6
ixs={s:np.flatnonzero(pix.sample_id==s) for s in SAMPLES}
TITLE={'ACTH':'Cushing vs NF-ACTH','GH':'GH vs NF-GH'}
def footer(fig,page,text='Exploratory MSI associations. Biological unit: patient; nonfunctional group n = 1 in each contrast.'):
    fig.text(.04,.025,text,fontsize=8,color='#555555');fig.text(.97,.025,str(page),ha='right',fontsize=8)
def tidy(ax):ax.spines[['top','right']].set_visible(False)
def highlabel(row):return ('Cushing' if row.contrast=='ACTH' else 'GH') if row.direction=='functional_higher' else ('NF-ACTH' if row.contrast=='ACTH' else 'NF-GH')
def choose(con,n=8):
    d=C[C.contrast==con];return pd.concat([d[d.direction==di].head(n) for di in COL])
def atlas():
    index=[];limits=[];pdfpath=FIG/'MSI_all_candidates_spatial_atlas_269pages.pdf'
    with PdfPages(pdfpath) as pdf:
        for page,(_,row) in enumerate(C.iterrows(),1):
            con=row.contrast;j=int(row.feature_index)-1;si=sum(CONTRASTS[con],[]);samples=[SAMPLES[k] for k in si]
            vals=[];weights=[]
            for s in samples:
                z=np.log2(1+x[ixs[s],j]);z=z[z>0]
                if len(z):vals.append(z);weights.append(np.full(len(z),1/len(z)))
            v=np.concatenate(vals);w=np.concatenate(weights);order=np.argsort(v);cs=np.cumsum(w[order]);cs/=cs[-1]
            lo,hi=np.interp([.02,.98],cs,v[order]);hi=max(hi,lo+1e-6);norm=Normalize(lo,hi,clip=True)
            fig=plt.figure(figsize=(13,8.5))
            fig.suptitle(f'{TITLE[con]} | F{j+1:03d} | m/z {mz[j]:.6f} | {highlabel(row)} higher',x=.04,y=.96,ha='left',fontsize=16,weight='bold')
            fig.text(.04,.906,f'Effect {row.effect:+.3f} | minimum pairwise margin {row.minimum_pairwise_margin:.3f} | robust tier: {"yes" if row.robust_candidate else "no"} | exact P {row.p_exact_two_sided:.3f} | FDR {row.FDR:.3f}',fontsize=10)
            gs=fig.add_gridspec(2,len(si),left=.05,right=.94,bottom=.47,top=.86,height_ratios=[1,.06],hspace=.25,wspace=.08)
            for k,s in enumerate(samples):
                p=pix.iloc[ixs[s]];assert p.source_id.nunique()==1
                xx=(p.x-p.x.min()).to_numpy(int);yy=(p.y-p.y.min()).to_numpy(int)
                arr=np.full((yy.max()+1,xx.max()+1),np.nan);z=np.log2(1+x[ixs[s],j]);arr[yy,xx]=z
                ax=fig.add_subplot(gs[0,k]);gray=np.where(np.isfinite(arr),.84,np.nan)
                ax.imshow(gray,cmap='gray',vmin=0,vmax=1,origin='upper',interpolation='nearest')
                ax.imshow(np.ma.masked_where(~np.isfinite(arr)|(arr==0),arr),cmap='RdBu_r',norm=norm,origin='upper',interpolation='nearest')
                ax.set_title(f'{s}\n{len(p):,} tumor pixels',fontsize=10);ax.set_axis_off()
                positive=z[z>0];low=float((positive<lo).mean()) if len(positive) else 0;high=float((positive>hi).mean()) if len(positive) else 0
                limits.append(dict(page=page,contrast=con,feature_index=j+1,sample=s,lower_limit=lo,upper_limit=hi,positive_clipped_low=low,positive_clipped_high=high))
            cax=fig.add_subplot(gs[1,:]);fig.colorbar(ScalarMappable(norm=norm,cmap='RdBu_r'),cax=cax,orientation='horizontal')
            cax.set_xlabel('Pixel log2(1 + TIC signal x 10^6): equal-specimen weighted positive P2-P98; same range within page',fontsize=9)
            gs2=fig.add_gridspec(1,3,left=.07,right=.94,bottom=.18,top=.36,wspace=.4)
            for k,(key,label,trans) in enumerate([('primary','All-tumor mean\nlog2(1 + mean signal x 10^6)',False),('detection','Detected tumor pixels (%)',False),('positive','Positive-pixel mean\nlog2(1 + mean signal x 10^6)',True)]):
                ax=fig.add_subplot(gs2[0,k]);y=A[key][si,j].copy()
                if trans:y=np.log2(1+y)
                if key=='detection':y*=100
                ax.scatter(range(len(si)),y,c=[COL['functional_higher']]*len(CONTRASTS[con][0])+[COL['nonfunctional_higher']],s=42,zorder=3)
                for h,val in enumerate(y):
                    if np.isfinite(val):ax.annotate(f'{val:.2f}',(h,val),xytext=(0,6),textcoords='offset points',ha='center',fontsize=8)
                    else:ax.text(h,.1,'NA',transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
                ax.set_xticks(range(len(si)),samples,rotation=20,ha='right',fontsize=8);ax.set_title(label,fontsize=10)
                ax.set_xlim(-.4,len(si)-.6);ax.margins(y=.35);ax.grid(axis='y',alpha=.2)
                if key=='detection':ax.set_ylim(0,115)
            component_text=(f'detection contribution {row.detection_contribution:+.1f}; positive-intensity contribution {row.positive_intensity_contribution:+.1f} (linear signal x 10^6)' if row.decomposition_evaluable else 'components not estimated: fewer than 20 positive pixels in at least one specimen')
            fig.text(.05,.105,f'Difference source: {row.decomposition_pattern.replace("_"," ")} | {component_text}',fontsize=9)
            fig.text(.05,.079,f'Conditional direction support: 10-pixel blocks {row.conditional_direction_support:.1%}; minimum over 5/10/20-pixel blocks {row.minimum_multiscale_support:.1%}. Grey = zero; white = outside tumor.',fontsize=9)
            footer(fig,page);pdf.savefig(fig,dpi=150);plt.close(fig)
            index.append(dict(page=page,contrast=con,feature_index=j+1,mz=mz[j],higher_group=highlabel(row),robust_candidate=row.robust_candidate))
            if page%40==0:print('Atlas pages:',page,flush=True)
    pd.DataFrame(index).to_csv(TABLE/'atlas_page_index.csv',index=False)
    pd.DataFrame(limits).to_csv(TABLE/'atlas_scale_limits.csv',index=False)
    print('Atlas complete:',len(index),flush=True)

def overview():
    with PdfPages(FIG/'MSI_lineage_analysis_summary.pdf') as pdf:
        fig=plt.figure(figsize=(13,8.5));fig.suptitle('Within-lineage hormone secretion phenotypes',x=.05,y=.95,ha='left',fontsize=21,weight='bold')
        fig.text(.05,.89,'MSI only | consensus_tuned225 | 274 aligned peaks | 7 patients | 29,656 tumor pixels',fontsize=12)
        ax=fig.add_axes([.05,.53,.9,.3]);ax.axis('off')
        rows=[['Cushing vs NF-ACTH','2 vs 1','38 / 89','127','21 / 36'],['GH vs NF-GH','3 vs 1','126 / 16','142','81 / 6']]
        t=ax.table(cellText=rows,colLabels=['Contrast','Patients','F high / NF high','Candidates','Robust: F / NF'],loc='center',cellLoc='center');t.auto_set_font_size(False);t.set_fontsize(12);t.scale(1,2.8)
        for (r,c),cell in t.get_celld().items():cell.set_edgecolor('white');cell.set_facecolor('#dbe7ef' if r==0 else '#f0f3f6')
        fig.text(.06,.48,'Candidate definition',fontsize=14,weight='bold')
        texts=['Absolute difference >= 0.5 in log2(1 + all-tumor mean signal x 10^6).',
          'Pairwise directional concordance >= 0.8; direction retained after leaving out each functional patient.',
          'At least 2/3 of high-group patients: detection >= 10% and >= 20 positive tumor pixels.',
          'Robust subset additionally: minimum pairwise margin >= 0.5, median-centering direction preserved,',
          'and conditional spatial direction support >= 80% at all three block widths (5 / 10 / 20).']
        for k,tx in enumerate(texts):fig.text(.06,.435-k*.042,tx,fontsize=11)
        fig.text(.06,.15,'No peak reaches patient-level FDR < 0.05.',fontsize=14,weight='bold',color='#9c3030')
        fig.text(.06,.11,'Exact label assignments: 3 (ACTH) and 4 (GH). These are exploratory candidates, not validated biomarkers.',fontsize=11)
        footer(fig,1);pdf.savefig(fig);plt.close(fig)
        for page,con in enumerate(CONTRASTS,2):
            d=choose(con,10);js=d.feature_index.to_numpy(int)-1;si=sum(CONTRASTS[con],[])
            mat=A['primary'][si][:,js].T;center=mat-mat.mean(1,keepdims=True);v=max(abs(center).max(),.5)
            fig=plt.figure(figsize=(13,8.5));fig.suptitle(f'{TITLE[con]} | representative candidates in both directions',x=.05,y=.95,ha='left',fontsize=18,weight='bold')
            ax=fig.add_axes([.15,.17,.3,.65]);im=ax.imshow(center,cmap='RdBu_r',vmin=-v,vmax=v,aspect='auto')
            ax.set_xticks(range(len(si)),[SAMPLES[k] for k in si],rotation=25,ha='right');ax.set_yticks(range(len(js)),[f'F{j+1:03d} | {mz[j]:.4f}' for j in js],fontsize=9)
            ax.set_title('Sample summary, centered within peak');cb=fig.colorbar(im,ax=ax,shrink=.6);cb.set_label('Deviation from peak mean')
            ax2=fig.add_axes([.57,.17,.37,.65]);yy=np.arange(len(js));ax2.barh(yy,d.effect,color=[COL[x] for x in d.direction]);ax2.set_ylim(len(js)-.5,-.5);ax2.set_yticks(yy,['*' if x else '' for x in d.robust_candidate]);ax2.axvline(0,color='black',lw=.8);ax2.set_xlabel('Functional minus nonfunctional\nlog2(1 + mean signal x 10^6)');ax2.set_title('Effect; * robust subset');ax2.grid(axis='x',alpha=.2)
            fig.text(.05,.083,'Top 10 per direction, ordered by robust subset and minimum pairwise margin. Full candidates are in the spatial atlas.',fontsize=9)
            footer(fig,page);pdf.savefig(fig);plt.close(fig)
        fig,axs=plt.subplots(1,2,figsize=(13,8.5));fig.subplots_adjust(left=.09,right=.97,top=.82,bottom=.23,wspace=.35)
        fig.suptitle('Where does the difference arise?',x=.05,y=.95,ha='left',fontsize=20,weight='bold')
        fig.text(.05,.88,'Exact symmetric decomposition: total mean difference = detection contribution + positive-intensity contribution',fontsize=11)
        for ax,con in zip(axs,CONTRASTS):
            d=C[C.contrast==con];counts=d.groupby(['direction','decomposition_pattern']).size()
            bottom=np.zeros(2)
            for pat,col in PAT.items():
                vv=np.array([counts.get((di,pat),0) for di in COL]);ax.bar([0,1],vv,bottom=bottom,color=col,label=pat.replace('_',' '))
                for k,n in enumerate(vv):
                    if n>=5:ax.text(k,bottom[k]+n/2,str(n),ha='center',va='center',color='white',weight='bold')
                bottom+=vv
            ax.set_xticks([0,1],['Functional higher','Nonfunctional higher']);ax.set_ylabel('Candidate peaks');ax.set_title(TITLE[con]);ax.set_ylim(0,max(bottom)*1.08)
        handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.075),ncol=3,frameon=False)
        footer(fig,4,'Dominance uses absolute component shares >= 2/3 with concordant signs. Opposing components are reported separately.');pdf.savefig(fig);plt.close(fig)
        fig=plt.figure(figsize=(13,8.5));fig.suptitle('Across lineages and sensitivity checks',x=.05,y=.95,ha='left',fontsize=20,weight='bold')
        cross=pd.read_csv(TABLE/'cross_lineage_effects.csv');ax=fig.add_axes([.08,.25,.39,.56])
        ax.scatter(cross.effect_ACTH,cross.effect_GH,s=14,c='#c5ccd3',label='All other peaks')
        for pat,col in [('both_functional_higher','#bc3b37'),('both_nonfunctional_higher','#2668a2'),('opposite_or_zero','#765798')]:
            z=cross[cross.both_candidates&(cross.pattern==pat)];ax.scatter(z.effect_ACTH,z.effect_GH,s=26,c=col,label=f'{pat.replace("_"," ")} ({len(z)})')
        ax.axhline(0,color='gray',lw=.8);ax.axvline(0,color='gray',lw=.8);ax.set_xlabel('ACTH: functional minus nonfunctional');ax.set_ylabel('GH: functional minus nonfunctional');ax.set_title('Same aligned peak; each contrast analyzed separately');ax.legend(fontsize=8,loc='upper left')
        ax=fig.add_axes([.61,.25,.34,.56]);sen=pd.read_csv(TABLE/'threshold_sensitivity.csv');sen=sen[sen.minimum_margin_threshold==0]
        for con,ls in [('ACTH','-'),('GH','--')]:
            for di,col in COL.items():
                z=sen[sen.contrast==con];ax.plot(z.abs_effect_threshold,z[di],marker='o',ls=ls,color=col,label=f'{con}: {di.replace("_"," ")}')
        ax.set_xlabel('Absolute primary effect threshold');ax.set_ylabel('Candidate peaks');ax.set_xticks([.25,.5,1]);ax.set_title('Effect-threshold sensitivity');ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.text(.07,.13,'56 peaks are candidates in both contrasts: 10 higher in both functional groups, 2 higher in both nonfunctional groups,',fontsize=10)
        fig.text(.07,.095,'and 44 in opposite directions. Shared m/z does not establish molecular identity or a causal secretion mechanism.',fontsize=10)
        footer(fig,5);pdf.savefig(fig);plt.close(fig)
    print('Summary complete: 5 pages',flush=True)

if __name__=='__main__':overview();atlas()
