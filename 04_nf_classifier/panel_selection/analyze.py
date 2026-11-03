from pathlib import Path
import sys,ast,itertools,json,time,hashlib
import numpy as np,pandas as pd
from scipy.stats import rankdata
H=Path(__file__).resolve().parent;R=H.parent;P=R/'result3_endogenous_constrained_20260915';sys.path.insert(0,str(P));import analysis as old
CL=old.CL;SA=old.SA
pool=pd.read_csv(R/'result3_joint55_20260915/candidate_pool.csv');fid=np.array(sorted(pool.feature_id.unique(),key=lambda s:int(s[1:])));assert len(fid)==55
rm=pd.read_csv(P/'geometry/ROI_manifest.csv');assert np.array_equal(rm.roi_id,rm.index)
train=rm[rm.width_pixels.eq(10)];X=np.load(P/'cache/ROI_arrays.npz')['response'][:,[int(s[1:])-1 for s in fid]]
assert np.isfinite(X).all()
# Reuse only pure numerical definitions. Do not read LC abundance or prediction tables.
tree=ast.parse((R/'result3_joint55_20260915/run.py').read_text());nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name in ['stats','scores','metrics','Spatial']];exec(compile(ast.Module(body=nodes,type_ignores=[]),'shared_numerics','exec'))
def say(s):print(time.strftime('%H:%M:%S'),s,flush=True)
def association(meta):
 vals={s:X[g.index] for s,g in meta.groupby('sample_id')};assert set(vals)==set(SA)
 means={s:v.mean(0) for s,v in vals.items()};pair={}
 for a,b in itertools.combinations(SA,2):
  na,nb=len(vals[a]),len(vals[b]);ranks=rankdata(np.concatenate([vals[a],vals[b]]),axis=0);auc=(ranks[:na].sum(0)-na*(na+1)/2)/(na*nb);pair[a,b]=auc;pair[b,a]=1-auc
 w=old.weights(meta);mu=w@X[meta.index];total=w@((X[meta.index]-mu)**2);between=np.zeros(len(fid))
 for c in range(3):
  ix=meta.class_index.to_numpy()==c;cm=w[ix]@X[meta.index[ix]]/w[ix].sum();between+=w[ix].sum()*(cm-mu)**2
 rows=[]
 for j,f in enumerate(fid):
  opts=[]
  for c in range(3):
   aa=[s for s in SA if old.LABEL[s]==c];bb=[s for s in SA if old.LABEL[s]!=c]
   for direction in [1,-1]:
    auc=min((pair[a,b][j] if direction==1 else 1-pair[a,b][j]) for a in aa for b in bb)
    gap=min(direction*(means[a][j]-means[b][j]) for a in aa for b in bb)
    opts.append((auc,gap,c,direction))
  auc,gap,c,direction=max(opts,key=lambda z:(z[0],z[1],-z[2],z[3]))
  rows.append(dict(feature_id=f,index=j,target_class=CL[c],direction='high' if direction==1 else 'low',worst_specimen_pair_AUC=auc,min_signed_mean_gap_log2=gap,weighted_eta_squared=between[j]/max(total[j],1e-12),eligible=bool(auc>=.75 and gap>=.25),**{'ROI_mean_'+s:means[s][j] for s in SA}))
 return pd.DataFrame(rows).sort_values(['worst_specimen_pair_AUC','min_signed_mean_gap_log2','weighted_eta_squared','index'],ascending=[False,False,False,True])
def search(meta,tag):
 a=association(meta);eligible=a[a.eligible].copy();a['shortlisted']=a.index.isin(eligible.index);a.to_csv(H/f'{tag}_association.csv',index=False)
 ids=sorted(eligible['index'].tolist());assert len(ids)>=5, f'{tag}: only {len(ids)} eligible features'
 space=Spatial(meta);out={};rows=[]
 for k in [5]:
  pp=np.array(list(itertools.combinations(ids,k)),dtype=int);assert len(pp)>0;acc=[];loss=[]
  for st in range(0,len(pp),128):
   aa,ll=space.get(pp[st:st+128]);acc.extend(aa);loss.extend(ll)
  acc=np.array(acc);loss=np.array(loss);j=np.lexsort((np.arange(len(pp)),loss,-np.round(acc,12)))[0]
  out[k]=dict(indices=pp[j].tolist(),features=fid[pp[j]].tolist(),inner_BA=float(acc[j]),inner_logloss=float(loss[j]),n_combinations=len(pp),n_eligible=int(a.eligible.sum()),n_shortlist=len(ids))
  rows.extend(dict(k=k,features=';'.join(fid[p]),inner_BA=aa,inner_logloss=ll,selected=i==j) for i,(p,aa,ll) in enumerate(zip(pp,acc,loss)))
  say(tag+' '+str(k)+' '+str(out[k]))
 pd.DataFrame(rows).to_csv(H/f'{tag}_all_combinations.csv',index=False);return out
(H/'DESIGN.json').write_text(json.dumps(dict(pool='55 unique MSI features with prior reasonable-host LC candidate annotations; LC abundance never loaded',association='On training 10x10 ROIs: all specimen-pair rank AUCs; choose target class and high/low direction with greatest WORST pair AUC; require worst AUC>=0.75 and smallest signed specimen mean log2 gap>=0.25',shortlist='All eligible features; ranking for reporting only, no top-N truncation',search='Exhaustive search of every five-feature subset among all features passing the training-set association screen; no feature-identity exclusion, mutual-exclusion rule, or m/z-composition rule; maximize inner spatial class/specimen weighted recall then minimize logloss',classifier='Weighted standardized shrinkage LDA: lambda0.75 diagonal ridge1e-4',validation='4 outer spatial folds; re-screen all55 within each outer train; repeat association filtering and five-feature selection inside each outer training set; overlap+2pixel buffer purged. The full-data selected panel then receives a clearly labeled conditional spatial assessment with fold-specific refitting.',weights='Equal classes, equal physical specimens within class, equal ROIs within specimen',scope='Spatial performance in existing specimens. GH and ACTH each n1, GnH n3. No new-patient generalization or class significance claim.',input_sha256=hashlib.sha256((P/'cache/ROI_arrays.npz').read_bytes()).hexdigest()),indent=2))
final=search(train,'final');(H/'final_selection.json').write_text(json.dumps(final,indent=2));nested=[];folds=[]
for f in range(4):
 te=rm[rm.fold.eq(f)];tr=train[~train.fold.eq(f)];tr=tr[old.purge(tr,te,2)];out=search(tr,f'outer_{f}');st=stats(X[tr.index],tr.class_index.to_numpy(),old.weights(tr));folds.append(dict(fold=f,n_train=len(tr),n_test=len(te),panels=out))
 for k in [5]:
  sc=scores(st,X[te.index],np.array([out[k]['indices']]))[0];g=te.copy();g['k']=k;g['prediction']=sc.argmax(1);g['correct']=g.prediction==g.class_index;g['features']=';'.join(out[k]['features']);nested.append(g)
nd=pd.concat(nested);nd.to_csv(H/'nested_predictions.csv',index=False);(H/'outer_selection.json').write_text(json.dumps(folds,indent=2))
# Fixed-panel spatial CV is explicitly conditional on full-development selection.
fixed=[]
for f in range(4):
 te=rm[rm.fold.eq(f)];tr=train[~train.fold.eq(f)];tr=tr[old.purge(tr,te,2)];st=stats(X[tr.index],tr.class_index.to_numpy(),old.weights(tr))
 for k in [5]:
  sc=scores(st,X[te.index],np.array([final[k]['indices']]))[0];g=te.copy();g['k']=k;g['prediction']=sc.argmax(1);g['correct']=g.prediction==g.class_index;fixed.append(g)
fd=pd.concat(fixed);fd.to_csv(H/'fixed_panel_conditional_predictions.csv',index=False)
summ=[]
for method,d in [('nested_selection',nd),('fixed_panel_conditional',fd)]:
 for (k,w),g in d.groupby(['k','width_pixels']):summ.append(dict(method=method,k=k,width_pixels=w,n=len(g),correct=int(g.correct.sum()),raw_accuracy=g.correct.mean(),weighted_recall=old.weights(g)@g.correct.to_numpy()))
pd.DataFrame(summ).to_csv(H/'summary.csv',index=False)
nd.groupby(['k','width_pixels','sample_id','source_id']).agg(n=('correct','size'),correct=('correct','sum'),accuracy=('correct','mean')).reset_index().to_csv(H/'nested_piece_performance.csv',index=False)
a=association(train).set_index('feature_id');st=old.stats(X[train.index],train)
for k in [5]:
 ids=final[k]['indices'];co,it=old.fit(st,dict(indices=ids,lambda_=.75,**{'lambda':.75}));np.savez(H/f'fixed_{k}_MSI_model.npz',coefficient=co[ids],intercept=it,features=fid[ids],classes=np.array(CL),transform='log2(1+mean(TIC_fraction)*1e6)')
 pd.DataFrame({'feature_id':fid[ids],**{c:co[ids,j] for j,c in enumerate(CL)}}).to_csv(H/f'fixed_{k}_coefficients.csv',index=False)
 a.loc[fid[ids]].to_csv(H/f'panel_{k}_association.csv');pool[pool.feature_id.isin(fid[ids])].to_csv(H/f'panel_{k}_candidate_annotations.csv',index=False)
(H/'model_intercepts.json').write_text(json.dumps({str(k):np.load(H/f'fixed_{k}_MSI_model.npz')['intercept'].tolist() for k in [5]},indent=2))
print(pd.DataFrame(summ).to_string(index=False),flush=True);say('DONE')
