"""MSI-only within-lineage reanalysis. No LC inputs, annotations or legacy features.
Run with Python and installed numpy/pandas. Plotting additionally needs matplotlib.
"""
from pathlib import Path
import itertools, json, gzip, hashlib, warnings
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=Path(__file__).resolve().parent/'20260914'
SOURCE=OUT/'input/reference'
TABLE=OUT/'tables'; TABLE.mkdir(parents=True,exist_ok=True)
SAMPLES=['Cushing-1','Cushing-2','GH-1','GH-2','GH-3','NF-GH','NF-ACTH']
CONTRASTS={'ACTH':([0,1],[6]),'GH':([2,3,4],[5])}
SEED=20260914
def save(d,n): d.to_csv(TABLE/(n+'.csv'),index=False)
def effect(v,a,b):return v[a].mean(axis=0)-v[b].mean(axis=0)
def bh(p):
    order=np.argsort(p); v=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.zeros_like(p);out[order]=np.minimum(v,1);return out
def prepare():
    pix=pd.read_csv(SOURCE/'tumor_pixel_selection_used.csv')
    pix=pix[pix.sample_id.isin(SAMPLES)].reset_index(drop=True)
    raw=(OUT/'input/tumor_matrix.bin.gz').read_bytes()
    assert hashlib.md5(raw).hexdigest()=='114ae80bd86c47038f7ead3b43789f1a'
    x=np.frombuffer(gzip.decompress(raw),dtype='<f8').reshape((len(pix),274),order='F')*1e6
    assert len(pix)==29656 and not pix.original_pixel_row.duplicated().any()
    assert np.isfinite(x).all() and (x>=0).all()
    mz=pd.read_csv(SOURCE/'MSI_sample_primary.csv').mz.to_numpy()
    n=np.array([(pix.sample_id==s).sum() for s in SAMPLES])
    sums=np.stack([x[pix.sample_id==s].sum(0) for s in SAMPLES])
    counts=np.stack([(x[pix.sample_id==s]>0).sum(0) for s in SAMPLES])
    linear=sums/n[:,None];det=counts/n[:,None]
    positive=np.divide(sums,counts,out=np.full_like(sums,np.nan),where=counts>0)
    positive[counts<20]=np.nan
    primary=np.log2(1+linear)
    oldmean=np.stack([np.log1p(x[pix.sample_id==s]).mean(0) for s in SAMPLES])
    errors={}
    for tag,v in [('primary',primary),('linear',linear),('detection',det),('positive_linear',positive),('counts',counts),('oldmean',oldmean)]:
        ref=pd.read_csv(SOURCE/f'MSI_sample_{tag}.csv')[SAMPLES].to_numpy().T
        errors[tag]=float(np.nanmax(abs(v-ref)))
        assert np.allclose(v,ref,rtol=1e-10,atol=1e-8,equal_nan=True),tag
        save(pd.DataFrame(v.T,columns=SAMPLES).assign(feature_index=np.arange(1,275),mz=mz),f'MSI_sample_{tag}')
    save(pix,'tumor_pixel_metadata')
    save(pd.DataFrame({'sample':SAMPLES,'n_tumor_pixels':n}), 'sample_manifest')
    (OUT/'input/pixels.csv').write_text(pix.to_csv(index=False))
    (OUT/'input/summary_validation.json').write_text(json.dumps(errors,indent=2))
    print('Fresh pixel summaries agree with all 6 current source tables:',errors,flush=True)
    return pix,x,mz,linear,det,positive,primary,counts,oldmean

def differential(mz,primary,det,counts,oldmean):
    results=[];sensitivity=[]
    for con,(aa,bb) in CONTRASTS.items():
        a=np.array(aa);b=np.array(bb);e=effect(primary,a,b);sg=np.sign(e)
        dif=primary[a,None,:]-primary[None,b,:]
        concord=(dif*sg>1e-10).mean((0,1))+.5*(abs(dif)<=1e-10).mean((0,1))
        strict=(dif*sg).min((0,1))
        v=primary[np.r_[a,b]]
        perm=np.stack([effect(v,list(c),[i for i in range(len(v)) if i not in c]) for c in itertools.combinations(range(len(v)),len(a))])
        p=(abs(perm)>=abs(e)-1e-10).mean(0)
        loo=np.stack([effect(primary,a[a!=i],b)*sg>0 for i in a]).all(0)
        eligible=(det>=.1)&(counts>=20)
        high=np.where(e>=0,eligible[a].mean(0),eligible[b].mean(0))
        centered=effect(primary-np.median(primary,axis=1,keepdims=True),a,b)
        base=(concord>=.8-1e-9)&loo&(high>=2/3-1e-9)
        d=pd.DataFrame(dict(contrast=con,feature_index=np.arange(1,275),mz=mz,effect=e,
          direction=np.where(e>0,'functional_higher',np.where(e<0,'nonfunctional_higher','equal')),
          functional_mean=primary[a].mean(0),nonfunctional_mean=primary[b].mean(0),
          direction_concordance=concord,minimum_pairwise_margin=strict,strict_separation=strict>1e-10,
          leave_one_functional_direction=loo,high_group_eligible_fraction=high,base_eligible=base,
          candidate=base&(abs(e)>=.5),median_centered_effect=centered,centering_same_direction=centered*e>0,
          mean_pixel_log1p_effect=effect(oldmean,a,b),mean_pixel_log1p_same_direction=effect(oldmean,a,b)*e>0,
          detection_difference=effect(det,a,b),p_exact_two_sided=p,FDR=bh(p),
          n_label_assignments=len(perm),n_functional=len(a),n_nonfunctional=len(b)))
        for s in np.r_[a,b]:d[SAMPLES[s]]=primary[s]
        results.append(d)
        for th in [.25,.5,1.]:
            for margin in [0,.25,.5]:
                m=base&(abs(e)>=th)&(strict>margin)
                sensitivity.append(dict(contrast=con,abs_effect_threshold=th,minimum_margin_threshold=margin,
                    functional_higher=int((m&(e>0)).sum()),nonfunctional_higher=int((m&(e<0)).sum())))
    save(pd.DataFrame(sensitivity),'threshold_sensitivity')
    return pd.concat(results,ignore_index=True)

def components(lin,p,q,a,b):
    d=effect(lin,a,b)
    cross=p[a].mean(0)*q[b].mean(0)-p[b].mean(0)*q[a].mean(0)
    return d,(d+cross)/2,(d-cross)/2

def spatial(pix,x,lin,det,pos,counts,stats):
    dec=[];boots=[];manifest=[]
    for con,(a,b) in CONTRASTS.items():
        d,ec,ic=components(lin,det,pos,a,b);ok=(counts[a+b]>=20).all(0)
        ec[~ok]=np.nan;ic[~ok]=np.nan
        share=abs(ec)/(abs(ec)+abs(ic)+1e-30)
        pattern=np.where(~ok,'insufficient_positive_pixels',np.where(ec*ic<0,'opposing_components',
            np.where(share>=2/3,'detection_dominant',np.where(share<=1/3,'intensity_dominant','mixed'))))
        dec.append(pd.DataFrame(dict(contrast=con,feature_index=np.arange(1,275),linear_difference=d,
            detection_contribution=ec,positive_intensity_contribution=ic,absolute_detection_share=share,
            decomposition_evaluable=ok,decomposition_pattern=pattern,reconstruction_error=d-ec-ic)))
    dec=pd.concat(dec,ignore_index=True)
    for width in [5,10,20]:
        B=300 if width==10 else 100;rng=np.random.default_rng(SEED+width)
        sampleboot=[]
        for s in SAMPLES:
            ix=np.flatnonzero(pix.sample_id==s); pp=pix.iloc[ix]
            keys=pp.source_id.astype(str)+':'+(pp.x//width).astype(str)+':'+(pp.y//width).astype(str)
            codes,uniq=pd.factorize(keys)
            bs=np.zeros((len(uniq),274));bc=np.zeros_like(bs)
            np.add.at(bs,codes,x[ix]);np.add.at(bc,codes,x[ix]>0);bn=np.bincount(codes)
            # Resample blocks within each source, preserving source identity and block counts.
            weights=np.zeros((B,len(uniq)))
            sources=np.array([key.rsplit(':',2)[0] for key in uniq])
            for src in np.unique(sources):
                si=np.flatnonzero(sources==src)
                weights[:,si]=rng.multinomial(len(si),np.full(len(si),1/len(si)),size=B)
            su=weights@bs;ct=weights@bc;nn=weights@bn
            ml=su/nn[:,None];mp=ct/nn[:,None];mq=np.divide(su,ct,out=np.full_like(su,np.nan),where=ct>0)
            sampleboot.append((ml,mp,mq,ct))
            manifest.append(dict(sample=s,width_pixels=width,blocks=len(uniq),sources=len(np.unique(sources)),pixels=len(ix),B=B))
        ml,mp,mq,mc=[np.stack([s[i] for s in sampleboot]) for i in range(4)]
        for con,(a,b) in CONTRASTS.items():
            e=effect(np.log2(1+ml),a,b);d,ec,ic=components(ml,mp,mq,a,b)
            valid=(mc[a+b]>=20).all(0);ec[~valid]=np.nan;ic[~valid]=np.nan
            obs=stats[stats.contrast==con].effect.to_numpy()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                qs=np.nanquantile(np.stack([e,ec,ic]),[.025,.975],axis=1)
            boots.append(pd.DataFrame(dict(contrast=con,feature_index=np.arange(1,275),block_width=width,B=B,
                conditional_direction_support=(e*obs>0).mean(0),primary_low=qs[0,0],primary_high=qs[1,0],
                detection_contribution_low=qs[0,1],detection_contribution_high=qs[1,1],
                intensity_contribution_low=qs[0,2],intensity_contribution_high=qs[1,2],
                decomposition_valid_fraction=valid.mean(0))))
        print('Spatial bootstrap complete:',width,'pixels',flush=True)
    boot=pd.concat(boots,ignore_index=True)
    save(dec,'spatial_decomposition_all');save(boot,'spatial_bootstrap_sensitivity');save(pd.DataFrame(manifest),'block_manifest')
    return dec,boot

def main():
    pix,x,mz,lin,det,pos,pri,counts,oldmean=prepare()
    stats=differential(mz,pri,det,counts,oldmean)
    dec,boot=spatial(pix,x,lin,det,pos,counts,stats)
    stats=stats.merge(dec,on=['contrast','feature_index'],validate='one_to_one').merge(
        boot[boot.block_width==10],on=['contrast','feature_index'],validate='one_to_one')
    minboot=boot.groupby(['contrast','feature_index']).conditional_direction_support.min().rename('minimum_multiscale_support').reset_index()
    stats=stats.merge(minboot,on=['contrast','feature_index'],validate='one_to_one')
    stats['robust_candidate']=stats.candidate&(stats.minimum_pairwise_margin>=.5)&stats.centering_same_direction&(stats.minimum_multiscale_support>=.8)
    stats['abs_effect']=abs(stats.effect)
    stats=stats.sort_values(['contrast','direction','robust_candidate','minimum_pairwise_margin','abs_effect'],ascending=[True,True,False,False,False])
    save(stats,'MSI_differential_all');save(stats[stats.candidate],'MSI_candidates');save(stats[stats.robust_candidate],'MSI_robust_candidates')
    summary=stats.groupby(['contrast','direction']).agg(n_tested=('feature_index','size'),candidates=('candidate','sum'),
        robust_candidates=('robust_candidate','sum'),FDR_lt_005=('FDR',lambda s:(s<.05).sum())).reset_index()
    save(summary,'candidate_counts')
    save(stats[stats.candidate].groupby(['contrast','direction','decomposition_pattern']).size().rename('n').reset_index(),'candidate_decomposition_counts')
    a=stats[stats.contrast=='ACTH'];b=stats[stats.contrast=='GH']
    cols=['feature_index','mz','effect','candidate','robust_candidate','direction']
    cross=a[cols].merge(b[cols],on=['feature_index','mz'],suffixes=('_ACTH','_GH'),validate='one_to_one')
    cross['both_candidates']=cross.candidate_ACTH&cross.candidate_GH
    cross['pattern']=np.where(cross.effect_ACTH*cross.effect_GH>0,np.where(cross.effect_ACTH>0,'both_functional_higher','both_nonfunctional_higher'),'opposite_or_zero')
    save(cross,'cross_lineage_effects')
    params={'matrix_version':'consensus_tuned225','features':274,'patients':7,'pixels':29656,'LC_used':False,
       'seed':SEED,'primary':'log2(1 + all-tumor-pixel mean(TIC signal)*1e6); equal sample weights',
       'candidate':{'abs_effect':.5,'oriented_pairwise_concordance':.8,'leave_one_functional_direction':True,'high_group_eligible_fraction':2/3,'eligible_detection':.1,'eligible_positive_pixels':20},
       'robust_candidate_additional':{'minimum_pairwise_margin':.5,'median_centered_same_direction':True,'minimum_multiscale_spatial_support':.8},
       'bootstrap':{'widths':[5,10,20],'replicates':[100,300,100],'resampling':'blocks within source within sample; conditional spatial uncertainty, not biological replication'},
       'FDR':'exact two-sided sample label permutation, BH per contrast across 274 peaks',
       'singleton_limitation':'cannot estimate population variance or assess leave-one-out for NF samples',
       'input_md5':pd.read_csv(SOURCE/'actual_data_input_md5.csv').to_dict('records')}
    (OUT/'parameters.json').write_text(json.dumps(params,indent=2,ensure_ascii=False))
    np.savez_compressed(OUT/'input/sample_arrays.npz',primary=pri,linear=lin,detection=det,positive=pos,counts=counts,mz=mz)
    print(summary.to_string(index=False),flush=True)
    print('Cross-lineage candidates:',cross[cross.both_candidates].pattern.value_counts().to_dict(),flush=True)

if __name__=='__main__':main()
