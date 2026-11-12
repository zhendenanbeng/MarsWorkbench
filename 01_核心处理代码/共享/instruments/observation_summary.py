"""Observation-level descriptive summaries; no new mineral or uncertainty inference."""
import numpy as np

class MeanSpectrum:
    def __init__(self):
        self.x=None;self.total=None;self.count=None;self.signature=None;self.points=0;self.aligned=False;self.incompatible=False
    def add(self,data,backend):
        x,y=np.asarray(data['x'],float),np.asarray(data['y'],float)
        signature=(data.get('xlabel') or backend.xlabel,data.get('ylabel') or backend.ylabel,(data.get('meta') or {}).get('product'),(data.get('meta') or {}).get('spectral_column'))
        if x.ndim!=1 or y.shape!=x.shape or len(x)<3 or not np.isfinite(x).all() or not np.all(np.diff(x)>=0) or not np.any(np.diff(x)>0):self.incompatible=True;return
        if self.x is None:
            self.x=x.copy();self.total=np.zeros_like(x);self.count=np.zeros_like(x,dtype=int);self.signature=signature
        if signature!=self.signature:self.incompatible=True;return
        if np.array_equal(x,self.x):values=y
        else:
            if np.any(np.diff(x)==0) or np.any(np.diff(self.x)==0):self.incompatible=True;return
            self.aligned=True;values=np.full(self.x.shape,np.nan);indices=np.flatnonzero(np.isfinite(y));step=float(np.median(np.diff(x)))
            groups=np.split(indices,np.flatnonzero((np.diff(indices)>1)|(np.diff(x[indices])>5*step))+1)
            for group in groups:
                if len(group)<2:continue
                inside=(self.x>=x[group[0]])&(self.x<=x[group[-1]])
                values[inside]=np.interp(self.x[inside],x[group],y[group])
        finite=np.isfinite(values);self.total[finite]+=values[finite];self.count[finite]+=1;self.points+=1
    def result(self):
        if self.incompatible:return None,'谱轴、单位或谱列不一致，未混合平均'
        if self.x is None or not self.count.any():return None,'没有可平均的有效谱'
        y=np.divide(self.total,self.count,out=np.full(self.x.shape,np.nan),where=self.count>0)
        message='有效点等权平均：%d点；各通道贡献%d–%d点'%(self.points,int(self.count.min()),int(self.count.max()))
        if self.aligned:message+='；不同轴分段插值到首点网格，不跨空档或外推'
        return (self.x,y),message

def summarize(rows,backend_id,passive=False):
    valid=[r for r in rows if not r.get('error')];elements=sorted({e for r in valid for e in r['elements']})
    rates=[]
    for element in elements:
        if backend_id in ('chemcam','supercam'):
            denominator=sum(int(r.get('n_shots') or 0) for r in valid)
            numerator=sum(r.get('all_shot_support',{}).get(element,0) for r in valid);unit='发'
        else:denominator=len(valid);numerator=sum(element in r['elements'] for r in valid);unit='点'
        if denominator:rates.append({'element':element,'numerator':numerator,'denominator':denominator,'fraction':numerator/denominator,'unit':unit})
    rates.sort(key=lambda r:(-r['fraction'],r['element']))
    official={}
    source='moc' if backend_id=='supercam' else 'apxs'
    for oxide in (('SiO2','FeOT') if backend_id=='supercam' else ('FeO','SiO2') if backend_id=='apxs' else ()):
        values=[v['wt'] for r in rows for v in (r.get(source) or {}).get('rows',[]) if v['oxide']==oxide and v.get('wt') is not None and np.isfinite(v['wt'])]
        official[oxide]={'mean':float(np.mean(values)) if values else None,'available_points':len(values),'total_points':len(rows)}
    return {'element_count':None if passive or not valid else len(elements),'top':rates[:5],'applicable':not passive,
            'official':official,'valid_points':len(valid),'failed_points':len(rows)-len(valid)}

