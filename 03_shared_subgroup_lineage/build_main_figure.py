"""Single manuscript figure: within-lineage differences, decomposition, spatial examples.
No cross-lineage comparison panel. Frozen 20260914 MSI-only results are not reselected.
"""
from pathlib import Path
import gzip,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch

BASE=Path(__file__).resolve().parent/'20260914'
OUT=BASE/'main_figure';OUT.mkdir(exist_ok=True)
PDF=OUT/'output/pdf';PDF.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,
 'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'axes.linewidth':.5,
 'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none',
 'savefig.facecolor':'white'})
SAMPLES=['Cushing-1','Cushing-2','GH-1','GH-2','GH-3','NF-GH','NF-ACTH']
CONS={'ACTH':([0,1],[6]),'GH':([2,3,4],[5])}
NAMES={'ACTH':('Cushing','NF-ACTH'),'GH':('GH','NF-GH')}
COL={'functional_higher':'#B43B3B','nonfunctional_higher':'#286DA8'}
PATS={'intensity_dominant':('#9966A9','Intensity-dominant'),
      'detection_dominant':('#359584','Coverage-dominant'),
      'mixed':('#D1A641','Mixed'),
      'opposing_components':('#667684','Opposing'),
      'insufficient_positive_pixels':('#CDD2D6','Not estimable')}
d=pd.read_csv(BASE/'tables/MSI_differential_all.csv')
c=d[d.candidate].copy();a=np.load(BASE/'input/sample_arrays.npz');mz=a['mz']
pix=pd.read_csv(BASE/'tables/tumor_pixel_metadata.csv')
x=np.frombuffer(gzip.decompress((BASE/'input/tumor_matrix.bin.gz').read_bytes()),dtype='<f8').reshape((len(pix),274),order='F')*1e6
ix={s:np.flatnonzero(pix.sample_id==s) for s in SAMPLES}
selected=[]
for con in CONS:
    for pat,di in [('intensity_dominant','functional_higher'),('intensity_dominant','nonfunctional_higher'),('detection_dominant','functional_higher')]:
        eligible=c[(c.contrast==con)&(c.decomposition_pattern==pat)&(c.direction==di)&c.robust_candidate]
        row=eligible.sort_values(['minimum_pairwise_margin','abs_effect','feature_index'],ascending=[False,False,True]).iloc[0].copy()
        row['selection_rule']='Highest minimum pairwise margin within predeclared contrast/source/direction stratum among robust candidates'
        selected.append(row)
selected=pd.DataFrame(selected)
assert selected.feature_index.tolist()==[245,16,139,70,71,252]
selected.to_csv(OUT/'selected_spatial_features.csv',index=False)

fig=plt.figure(figsize=(190/25.4,250/25.4))
def txt(x,y,s,**kw):return fig.text(x,y,s,**kw)
def panel(letter,x,y,title):
    txt(x,y,letter,fontsize=12,weight='bold',va='top')
    txt(x+.027,y-.003,title,fontsize=9,weight='bold',va='top')
heatmaps=[];scales=[];spatialdata=[]
panel('A',.035,.978,'Cushing vs NF-ACTH')
panel('B',.535,.978,'GH vs NF-GH')

for con,left in [('ACTH',.035),('GH',.535)]:
    sub=pd.concat([c[(c.contrast==con)&(c.direction==di)].sort_values('abs_effect',ascending=False) for di in COL])
    js=sub.feature_index.to_numpy(int)-1;si=sum(CONS[con],[]);v=a['primary'][si][:,js].T
    zz=(v-v.mean(1,keepdims=True))/v.std(1,keepdims=True,ddof=0)
    assert np.isfinite(zz).all() and abs(zz).max()<2
    nf=(sub.direction=='functional_higher').sum();n=len(sub)
    txt(left+.028,.946,f'{len(CONS[con][0])} vs {len(CONS[con][1])} patients | {n} candidate peaks',fontsize=7,color='#444444')
    # No clustering: deterministic order by direction and absolute effect.
    hx=left+.105;hy=.713;hh=.213;hw=.185
    ax=fig.add_axes([hx,hy,hw,hh]);ax.imshow(zz,aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2,interpolation='nearest')
    ax.set_xticks(range(len(si)),[SAMPLES[s] for s in si],rotation=35,ha='right',fontsize=6.3)
    ax.tick_params(axis='x',length=0,pad=3);ax.set_yticks([]);ax.spines[:].set_visible(False)
    ax.axhline(nf-.5,color='white',lw=1.2)
    ax.text(-.14,1-nf/n/2,f'{NAMES[con][0]} high\nn = {nf}',transform=ax.transAxes,ha='right',va='center',fontsize=6.5,color=COL['functional_higher'])
    ax.text(-.14,(n-nf)/n/2,f'{NAMES[con][1]} high\nn = {n-nf}',transform=ax.transAxes,ha='right',va='center',fontsize=6.5,color=COL['nonfunctional_higher'])
    ax2=fig.add_axes([hx+hw+.020,hy,.127,hh]);yy=np.arange(n)
    ax2.barh(yy,sub.effect,height=1,color=[COL[z] for z in sub.direction],linewidth=0)
    ax2.set_ylim(n-.5,-.5);ax2.set_xlim(-8,8);ax2.axvline(0,color='#404040',lw=.5)
    ax2.set_xticks([-8,0,8]);ax2.set_yticks([]);ax2.set_xlabel('Effect',labelpad=1,fontsize=7)
    ax2.spines[['left','top','right']].set_visible(False)
    for order,(_,r) in enumerate(sub.iterrows(),1):
        for s in si:heatmaps.append({'panel':'A' if con=='ACTH' else 'B','contrast':con,'row':order,'feature_index':int(r.feature_index),'mz':r.mz,'sample':SAMPLES[s],'primary':a['primary'][s,int(r.feature_index)-1],'row_zscore':zz[order-1,si.index(s)],'effect':r.effect})

cax=fig.add_axes([.27,.657,.15,.007]);fig.colorbar(ScalarMappable(norm=Normalize(-2,2),cmap='RdBu_r'),cax=cax,orientation='horizontal',ticks=[-2,0,2]);cax.tick_params(length=2,pad=1,labelsize=6)
txt(.035,.658,'Heatmaps: within-peak z-score',fontsize=6.8)
txt(.535,.658,'Effect: functional minus nonfunctional',fontsize=6.8)

panel('C',.035,.626,'Sources of the whole-tumor signal difference')
ax=fig.add_axes([.235,.502,.475,.096])
summary=[]
for k,(con,di) in enumerate([(con,di) for con in CONS for di in COL]):
    sub=c[(c.contrast==con)&(c.direction==di)];n=len(sub);start=0
    for pat,(color,label) in PATS.items():
        count=int((sub.decomposition_pattern==pat).sum());width=count/n*100
        ax.barh(k,width,left=start,color=color,height=.66,edgecolor='white',linewidth=.35)
        if width>10:ax.text(start+width/2,k,str(count),ha='center',va='center',color='white' if pat!='mixed' else '#333333',fontsize=7,weight='bold')
        summary.append({'contrast':con,'direction':di,'pattern':pat,'n':count,'total':n,'percent':width});start+=width
    label=f'{NAMES[con][0 if di=="functional_higher" else 1]} high (n = {n})'
    ax.text(-.025,k,label,transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=7,color=COL[di])
ax.set_xlim(0,100);ax.set_ylim(3.6,-.6);ax.set_xticks([0,50,100],['0','50','100%']);ax.set_yticks([])
ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='x',length=2,pad=2)
ax.legend(handles=[Patch(facecolor=color,label=label) for color,label in PATS.values()],loc='upper left',bbox_to_anchor=(1.03,1.09),frameon=False,fontsize=6.7,handlelength=1,handletextpad=.5,labelspacing=.6)
txt(.235,.477,'Numbers within segments denote peaks; widths denote proportions.',fontsize=6.3,color='#555555')

panel('D',.035,.45,'ACTH: spatial examples')
panel('E',.535,.45,'GH: spatial examples')

for con,left in [('ACTH',.035),('GH',.535)]:
    si=sum(CONS[con],[]);ss=[SAMPLES[i] for i in si];w=.432;gap=.008;tilew=(w-gap*(len(si)-1))/len(si)
    # Within each contrast, all maps use a shared acquisition-coordinate frame.
    xmax=max(int(pix.iloc[ix[s]].x.max()-pix.iloc[ix[s]].x.min()+1) for s in ss)
    ymax=max(int(pix.iloc[ix[s]].y.max()-pix.iloc[ix[s]].y.min()+1) for s in ss)
    rows=selected[selected.contrast==con]
    for rownum,(_,r) in enumerate(rows.iterrows()):
        j=int(r.feature_index)-1;top=.423-rownum*.128
        typ='Intensity' if r.decomposition_pattern=='intensity_dominant' else 'Coverage'
        higher=NAMES[con][0 if r.direction=='functional_higher' else 1]
        txt(left+.025,top,f'F{j+1:03d}  m/z {mz[j]:.4f} | {higher} high | {typ}',fontsize=7.1,weight='bold',va='bottom')
        values=[];weights=[]
        for s in ss:
            v=np.log2(1+x[ix[s],j]);v=v[v>0];values.append(v);weights.append(np.full(len(v),1/len(v)))
        v=np.concatenate(values);wt=np.concatenate(weights);order=np.argsort(v);cdf=np.cumsum(wt[order]);cdf/=cdf[-1]
        lo,hi=np.interp([.02,.98],cdf,v[order]);norm=Normalize(lo,hi,clip=True)
        for k,s in enumerate(ss):
            p=pix.iloc[ix[s]];xx=(p.x-p.x.min()).to_numpy(int);yy=(p.y-p.y.min()).to_numpy(int)
            arr=np.full((ymax,xmax),np.nan);xx=xx+(xmax-xx.max()-1)//2;yy=yy+(ymax-yy.max()-1)//2
            z=np.log2(1+x[ix[s],j]);arr[yy,xx]=z
            ax=fig.add_axes([left+.025+k*(tilew+gap),top-.075,tilew,.065])
            ax.imshow(np.where(np.isfinite(arr),.85,np.nan),cmap='gray',vmin=0,vmax=1,interpolation='nearest')
            ax.imshow(np.ma.masked_where(~np.isfinite(arr)|(arr==0),arr),cmap='RdBu_r',norm=norm,interpolation='nearest');ax.axis('off')
            ax.set_title(s,fontsize=6.2,pad=1)
            sval=f'{a["primary"][si[k],j]:.2f} / {100*a["detection"][si[k],j]:.1f}%'
            fig.text(left+.025+k*(tilew+gap)+tilew/2,top-.083,sval,ha='center',va='top',fontsize=6.2)
            scales.append({'panel':'D' if con=='ACTH' else 'E','contrast':con,'feature_index':j+1,'sample':s,'lower':lo,'upper':hi,'positive_clipped_low':float((z[z>0]<lo).mean()),'positive_clipped_high':float((z[z>0]>hi).mean())})
            spatialdata.append({'contrast':con,'feature_index':j+1,'mz':mz[j],'sample':s,'primary':a['primary'][si[k],j],'detected_percent':100*a['detection'][si[k],j],'positive_mean_log2p1':np.log2(1+a['positive'][si[k],j]),'tumor_pixels':len(p)})
        cb=fig.add_axes([left+.05,top-.098,.245,.005]);fig.colorbar(ScalarMappable(norm=norm,cmap='RdBu_r'),cax=cb,orientation='horizontal',ticks=[lo,hi]);cb.set_xticklabels([f'{lo:.1f}',f'{hi:.1f}']);cb.tick_params(length=1.5,pad=1,labelsize=5.5)
        txt(left+.314,top-.102,f'Effect {r.effect:+.2f}',fontsize=6.2)

txt(.06,.043,'Below maps: whole-tumor log-mean / detected pixels (%). Grey: zero; white: outside tumor.',fontsize=6.5)
txt(.06,.028,'Map colors: pixel log2(1 + TIC signal x 10^6), shared positive P2-P98 limits within each peak and comparison.',fontsize=6.3)
txt(.06,.013,'Exploratory candidates; patient-level FDR >= 0.05 throughout. NF-ACTH and NF-GH each have one patient.',fontsize=6.3,color='#555555')

pd.DataFrame(heatmaps).to_csv(OUT/'panels_AB_source_data.csv',index=False)
pd.DataFrame(summary).to_csv(OUT/'panel_C_source_data.csv',index=False)
pd.DataFrame(spatialdata).to_csv(OUT/'panels_DE_source_data.csv',index=False)
pd.DataFrame(scales).to_csv(OUT/'panels_DE_scale_limits.csv',index=False)
fig.savefig(PDF/'Figure_lineage_secretion_MSI.pdf',dpi=600)
fig.savefig(OUT/'Figure_lineage_secretion_MSI.svg',dpi=600)
fig.savefig(OUT/'Figure_lineage_secretion_MSI.png',dpi=600)
fig.savefig(OUT/'Figure_lineage_secretion_MSI.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
plt.close(fig)
print('Main figure created. Selected peaks:',selected[['contrast','feature_index','mz','decomposition_pattern','direction']].to_dict('records'))
