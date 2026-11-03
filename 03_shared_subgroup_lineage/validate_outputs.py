from pathlib import Path
import json
import numpy as np
import pandas as pd
import pypdfium2 as pdfium
from pypdf import PdfReader
from PIL import Image,ImageDraw
from analyze_msi import OUT,TABLE,CONTRASTS
qa=OUT/'qa';qa.mkdir(exist_ok=True)
d=pd.read_csv(TABLE/'MSI_differential_all.csv');c=d[d.candidate];ix=pd.read_csv(TABLE/'atlas_page_index.csv')
assert len(d)==548 and len(c)==269 and len(ix)==269
assert not d[['contrast','feature_index']].duplicated().any()
assert c.strict_separation.all() and (c.minimum_pairwise_margin>0).all()
assert (c.abs_effect>=.5).all() and (c.direction_concordance==1).all()
assert c.centering_same_direction.all() and (c.minimum_multiscale_support==1).all()
assert (d.FDR>=d.p_exact_two_sided-1e-12).all() and (d.FDR>=.05).all()
a=np.load(OUT/'input/sample_arrays.npz');maxerr=0
for con,(aa,bb) in CONTRASTS.items():
    z=d[d.contrast==con].sort_values('feature_index');dif=a['primary'][aa,None,:]-a['primary'][None,bb,:]
    e=dif.mean((0,1));assert np.allclose(e,z.effect)
    m=(dif*np.sign(e)).min((0,1));assert np.allclose(m,z.minimum_pairwise_margin)
    cov=[];inte=[]
    for i in aa:
        for j in bb:
            p=a['detection'];q=a['positive']
            cov.append((p[i]-p[j])*(q[i]+q[j])/2)
            inte.append((q[i]-q[j])*(p[i]+p[j])/2)
    assert np.allclose(np.mean(cov,axis=0),z.detection_contribution,equal_nan=True)
    assert np.allclose(np.mean(inte,axis=0),z.positive_intensity_contribution,equal_nan=True)
    maxerr=max(maxerr,np.nanmax(abs(z.reconstruction_error)))
assert maxerr<1e-8
boots=pd.read_csv(TABLE/'spatial_bootstrap_sensitivity.csv')
assert len(boots)==548*3 and set(boots.block_width)=={5,10,20}
manifest=[]
for name,expected in [('MSI_lineage_analysis_summary.pdf',5),('MSI_all_candidates_spatial_atlas_269pages.pdf',269)]:
    path=OUT/'output/pdf'/name;reader=PdfReader(path);assert len(reader.pages)==expected
    doc=pdfium.PdfDocument(str(path))
    thumbs=[]
    for i in range(expected):
        tx=reader.pages[i].extract_text()
        assert 'Exploratory' in tx or 'Dominance uses' in tx
        if expected==269:
            row=ix.iloc[i];assert f'F{int(row.feature_index):03d}' in tx and f'{row.mz:.6f}' in tx
        im=doc[i].render(scale=.3).to_pil().convert('RGB')
        assert np.asarray(im).std()>8
        if expected==5 or i in [0,37,38,126,127,252,253,268]:
            doc[i].render(scale=1.35).to_pil().save(qa/f'{"summary" if expected==5 else "atlas"}_{i+1:03d}.png')
        im.thumbnail((312,210));tile=Image.new('RGB',(322,235),'#dddddd');tile.paste(im,((322-im.width)//2,0));ImageDraw.Draw(tile).text((5,215),str(i+1),fill='black');thumbs.append(tile)
    for start in range(0,expected,36):
        part=thumbs[start:start+36];sheet=Image.new('RGB',(322*6,235*((len(part)+5)//6)),'white')
        for k,im in enumerate(part):sheet.paste(im,((k%6)*322,(k//6)*235))
        sheet.save(qa/f'{"summary" if expected==5 else "atlas"}_contact_{start+1:03d}.png')
    manifest.append(dict(file=name,pages=expected,size_bytes=path.stat().st_size))
pd.DataFrame(manifest).to_csv(OUT/'PDF_manifest.csv',index=False)
(qa/'validation.json').write_text(json.dumps(dict(status='PASS',differential_rows=len(d),candidate_pages=len(c),
    decomposition_max_reconstruction_error=float(maxerr),all_pdf_pages_rendered=True,all_atlas_page_identities_verified=True),indent=2))
print(json.dumps(manifest,indent=2));print('VALIDATION PASS')
