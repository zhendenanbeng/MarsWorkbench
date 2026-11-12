"""Lazy metadata/element search; failures remain unresolved, not negative detections."""
from collections import OrderedDict
from copy import deepcopy
import hashlib
from pathlib import Path
import re
import threading
from . import lines,evidence
_CACHE=OrderedDict();_LOCK=threading.RLock()
ELEMENTS={r[0].casefold():r[0] for r in lines.LIBS_LINES+lines.APXS_LINES if r[0] not in ('C2','CN')}

def point_support(backend,obs,point,steps):
    digest=hashlib.sha256(repr((backend.id,obs.product,sorted(steps.items()),point.meta,lines.LIBS_LINES,lines.APXS_LINES)).encode())
    for path in sorted(point.files.values()):digest.update(Path(path).read_bytes())
    key=digest.hexdigest()
    with _LOCK:
        if key in _CACHE:_CACHE.move_to_end(key);return deepcopy(_CACHE[key])
    data=backend.load_point(point,None,steps);marks,elements=backend.annotate_data(data)
    result={'elements':elements,'rates':evidence.element_statistics(marks,elements)}
    with _LOCK:
        _CACHE[key]=result
        while len(_CACHE)>256:_CACHE.popitem(last=False)
    return deepcopy(result)

def select(backend,observations,query,steps,cancelled=lambda:False):
    query=re.sub(r'(?i)\bsol\s+(\d+)\b',r'sol:\1',query.strip())
    tokens=query.replace(',',' ').casefold().split();wanted={ELEMENTS[t] for t in tokens if t in ELEMENTS}
    kept=[];unresolved=[];errors=[]
    for i,obs in enumerate(observations):
        if cancelled():return {'indices':[],'cancelled':True,'unresolved':[],'errors':[]}
        text=' '.join([str(obs.sol),obs.seq,obs.product,getattr(obs,'target','') or '']).casefold()
        valid=True
        for token in tokens:
            if token in ELEMENTS:continue
            sol=re.fullmatch(r'(?:sol[:_]?)?(\d+)',token)
            if (sol and int(sol.group(1))!=int(obs.sol)) or (not sol and token not in text):valid=False;break
        if not valid:continue
        if not wanted:kept.append(i);continue
        if getattr(obs,'passive',False) or (backend.id=='supercam' and obs.product in ('CP2','CP3')):continue
        seen=set();failed=False
        for point in obs.points:
            if cancelled():return {'indices':[],'cancelled':True,'unresolved':[],'errors':[]}
            try:seen.update(point_support(backend,obs,point,steps)['elements'])
            except Exception as exc:failed=True;errors.append({'index':i,'error':str(exc)[:180]});continue
            if wanted.issubset(seen):break
        if wanted.issubset(seen):kept.append(i)
        elif failed:unresolved.append(i)
    return {'indices':kept,'cancelled':False,'unresolved':unresolved,'errors':errors[:20]}

