"""PIXL PDS acquisition: inventory -> trusted registry refs -> verified local originals."""
import copy
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import threading
import urllib.parse
import urllib.request
import uuid
import xml.etree.ElementTree as ET
import appenv
from instruments import pixl_browser as browser
from instruments.download_queue import check_length
from sherloc import pixl

HOSTS={'pds.nasa.gov','pds.mcp.nasa.gov','pds-geosciences.wustl.edu','planetarydata.jpl.nasa.gov','pds-imaging.jpl.nasa.gov'}
MANIFEST='_pixl_download_manifest.json'

class Cancelled(Exception):pass

def utc():return datetime.now(timezone.utc).isoformat()

def destination():
    return appenv.empty_data('pixl')

def safe_url(url):
    parsed=urllib.parse.urlparse(url)
    if parsed.scheme!='https' or parsed.hostname not in HOSTS or parsed.username or parsed.password or parsed.port not in (None,443):raise ValueError('文件直链不在已核验的PDS HTTPS主机上')
    return url

def read_url(url):
    request=urllib.request.Request(safe_url(url),headers={'User-Agent':'MarsWorkbench/'+appenv.VERSION})
    with urllib.request.urlopen(request,timeout=30) as response:
        safe_url(response.geturl());raw=response.read();check_length(response,len(raw));return raw

def parse_inventory(raw,collection,kind):
    products={}
    for row in csv.reader(io.StringIO(raw.decode('utf-8-sig'))):
        if len(row)<2 or row[0].strip()!='P':continue
        lid=row[1].strip().split('::')[0]
        if not lid.startswith('urn:nasa:pds:mars2020_pixl:'+collection+':'):continue
        meta=pixl.parse_lid(lid)
        if not meta or meta['type']!=kind:continue
        if not re.fullmatch(r'[a-z0-9_]+',meta['stem']):continue
        key=meta['stem'][:-2];version=int(meta['stem'][-2:])
        if key not in products or version>products[key][0]:products[key]=(version,meta)
    return [value[1] for value in products.values()]

def inventory(collection,kind,refresh=False):
    cache=appenv.cache_write('pixl_download/inv_'+collection+'.csv')
    legacy=Path(pixl.CACHE)/('inv_'+collection+'.csv')
    fallback=cache if cache.exists() else legacy if legacy.exists() else None
    warning=None;source='本地缓存'
    if refresh or fallback is None:
        url=pixl.BASE+'/'+collection+'/collection_'+collection+'_inventory.csv'
        try:
            raw=read_url(url)
            if not parse_inventory(raw,collection,kind):raise ValueError('官方清单未包含当前支持的'+kind.upper()+'产品')
            temp=cache.with_name(cache.name+'.'+uuid.uuid4().hex+'.part');temp.write_bytes(raw);os.replace(temp,cache)
            fallback=cache;source='在线PDS清单'
        except Exception as exc:
            if fallback is None:raise ValueError('无法取得PIXL清单；本地文件仍可浏览。'+str(exc)) from exc
            warning='联网刷新未完成，使用已有缓存：'+str(exc)
    raw=fallback.read_bytes()
    return parse_inventory(raw,collection,kind),{'collection':collection,'source':source,'cache_time_utc':datetime.fromtimestamp(fallback.stat().st_mtime,timezone.utc).isoformat(),'sha256':hashlib.sha256(raw).hexdigest(),'warning':warning}

def catalog(refresh=False):
    chem,cm=inventory('data_oxides_pmc','rqc',refresh);positions,pm=inventory('data_processed','rxl',refresh)
    scans=[]
    for item in chem:
        matches=[p for p in positions if p['pos']==item['pos'] and p['seq'][:9]==item['seq'][:9]]
        scans.append({'id':hashlib.sha256(item['lid'].encode()).hexdigest()[:20],'sol':int(item['sol']),'placement':item['pos'],
            'rtt':item['seq'][:9],'file':item['stem']+'.csv','lid':item['lid'],'rxl_lids':[p['lid'] for p in matches],
            'coordinate_sols':sorted({int(p['sol']) for p in matches}),'coordinate_candidates':len(matches)})
    return {'scans':sorted(scans,key=lambda v:(v['sol'],v['file'])),'inventories':[cm,pm],
        'note':'按RQC成分处理sol列出；仅列已验证RQC格式，不代表PIXL全部观测。RXL按同定位+RTT配套，图像按官方引用下载。'}

def refs_from_record(record,lid):
    if not (record.get('summary') or {}).get('hits'):raise ValueError('PDS注册表没有该产品：'+lid)
    candidates=[item['properties'] for item in record.get('data',[]) if lid in item.get('properties',{}).get('lid',[])]
    if len(candidates)!=1:raise ValueError('注册表产品标识不唯一或不匹配')
    props=candidates[0];result={}
    def values(key):
        value=props.get(key,[]);return value if isinstance(value,list) else [value]
    for area in ('Data','Label'):
        prefix='ops:'+area+'_File_Info.ops:'
        names=values(prefix+'file_name');urls=values(prefix+'file_ref');md5s=values(prefix+'md5_checksum')
        if len(names)!=len(urls):raise ValueError('PDS文件名与文件引用未对齐')
        for index,(name,url) in enumerate(zip(names,urls)):
            if not isinstance(name,str) or Path(name).name!=name or not re.fullmatch(r'[A-Za-z0-9_.-]+',name):raise ValueError('PDS文件名无效')
            if urllib.parse.unquote(Path(urllib.parse.urlparse(url).path).name).lower()!=name.lower():raise ValueError('PDS直链文件名与产品引用不同')
            digest=md5s[index].lower() if index<len(md5s) and isinstance(md5s[index],str) and re.fullmatch(r'[a-fA-F0-9]{32}',md5s[index]) else None
            if name.lower() in result:raise ValueError('PDS文件引用重名')
            result[name.lower()]={'name':name,'url':safe_url(url),'md5':digest,'lid':lid}
    return result

def resolve(lid):
    query=urllib.parse.quote('lid eq "'+lid+'"',safe='')
    return refs_from_record(json.loads(read_url(pixl.IMGOPS_API+'?limit=1&q='+query)),lid)

def validate(path,entry,kind):
    raw=path.read_bytes()
    if not raw or b'<html' in raw[:256].lower() or b'<!doctype html' in raw[:256].lower():raise ValueError('下载返回空文件或HTML错误页')
    if entry.get('md5') and hashlib.md5(raw).hexdigest()!=entry['md5']:raise ValueError('官方MD5校验失败：'+entry['name'])
    if entry['name'].lower().endswith('.xml'):
        root=ET.fromstring(raw);lids=[n.text for n in root.iter() if n.tag.rsplit('}',1)[-1]=='logical_identifier']
        if lids!=[entry['lid']]:raise ValueError('XML产品标识与选择的LID不一致')
    elif kind=='rqc':
        values,_=browser.chemistry(path)
        if not values:raise ValueError('成分CSV没有有效PMC')
    elif kind=='rxl':
        if not browser.read_coordinates(path)['pts']:raise ValueError('坐标CSV没有有效PMC')
    elif kind=='image':
        image=pixl.mcc_gray(str(path))
        if image.ndim!=2 or not image.size:raise ValueError('MCC不是支持的二维图像')
    return {'file':entry['name'],'url':entry['url'],'lid':entry['lid'],'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'published_md5':entry.get('md5'),'verified_utc':utc()}

def get_file(folder,entry,kind,manifest,cancel,progress):
    folder=Path(folder).resolve();path=(folder/entry['name']).resolve()
    if not path.is_relative_to(folder) or Path(entry['name']).name!=entry['name']:raise ValueError('下载文件不能越出目录')
    def reuse():
        row=validate(path,entry,kind);old=manifest.get(path.name)
        if not entry.get('md5') and (not old or old.get('sha256')!=row['sha256'] or old.get('url')!=entry['url']):raise ValueError('同名本地文件未能核验，保留原件；请选择空目录：'+path.name)
        return dict(row,reused=True)
    if cancel():raise Cancelled()
    if path.exists():return reuse()
    temp=folder/(path.name+'.'+uuid.uuid4().hex+'.part')
    try:
        request=urllib.request.Request(safe_url(entry['url']),headers={'User-Agent':'MarsWorkbench/'+appenv.VERSION})
        count=0
        with urllib.request.urlopen(request,timeout=30) as response,temp.open('xb') as stream:
            safe_url(response.geturl())
            total=int(response.headers.get('Content-Length') or 0)
            while True:
                if cancel():raise Cancelled()
                chunk=response.read(65536)
                if not chunk:break
                stream.write(chunk);count+=len(chunk);progress(count,total)
            check_length(response,count)
        if cancel():raise Cancelled()
        row=validate(temp,entry,kind)
        try:
            if os.name=='nt':os.rename(temp,path)  # Windows rename refuses an existing destination.
            else:os.link(temp,path)
        except FileExistsError:return reuse()
        return row
    finally:
        temp.unlink(missing_ok=True)

class Controller:
    def __init__(self):
        self.lock=threading.RLock();self.catalog=None;self.token=None;self.job={'running':False,'status':'idle'}
    def query(self,refresh=False):
        data=catalog(refresh)
        with self.lock:self.catalog=data;self.token=uuid.uuid4().hex
        return dict(data,catalog_id=self.token,destination=destination())
    def snapshot(self):
        with self.lock:return copy.deepcopy(self.job)
    def update(self,**values):
        with self.lock:self.job.update(values)
    def log(self,text):
        with self.lock:self.job.setdefault('log',[]).append(text);self.job['log']=self.job['log'][-300:]
    def cancelled(self):
        with self.lock:return bool(self.job.get('cancel'))
    def cancel(self):
        with self.lock:
            if self.job.get('running'):self.job['cancel']=True
        return self.snapshot()
    def start(self,token,ids,destination,images=True):
        root=Path(destination).expanduser().resolve() if isinstance(destination,str) and destination.strip() else None
        if root is None:raise ValueError('请选择下载目录')
        if root.is_relative_to(appenv.app_dir().resolve()) and not root.is_relative_to(appenv.data_home().resolve()):raise PermissionError('请选择软件目录外的可写数据目录')
        with self.lock:
            if self.job.get('running'):raise ValueError('已有下载任务运行中，请等待或取消')
            if token!=self.token or not self.catalog:raise ValueError('清单已变化，请重新选择扫描')
            if not isinstance(ids,list) or not 1<=len(ids)<=50 or any(not isinstance(i,str) for i in ids) or len(set(ids))!=len(ids):raise ValueError('请选择1至50个不同扫描')
            by_id={row['id']:row for row in self.catalog['scans']}
            if any(i not in by_id for i in ids):raise ValueError('扫描不在当前官方清单内')
            scans=[copy.deepcopy(by_id[i]) for i in ids]
            root.mkdir(parents=True,exist_ok=True)
            self.job={'running':True,'status':'running','cancel':False,'folder':str(root),'scan_total':len(scans),'scan_done':0,
                'files_done':0,'reused':0,'bytes_done':0,'bytes_total':0,'log':[],'results':[],'ready':False,'started_utc':utc()}
            self.thread=threading.Thread(target=self.run,args=(scans,root,images),daemon=True);self.thread.start()
        return self.snapshot()
    def run(self,scans,root,images):
        try:
            for index,scan in enumerate(scans):
                if self.cancelled():break
                folder=(root/('sol_%05d'%scan['sol'])/('rtt_'+scan['rtt']+'_'+scan['id'][:8])).resolve()
                if not folder.is_relative_to(root):raise ValueError('下载目录越界')
                folder.mkdir(parents=True,exist_ok=True);manifest_path=folder/MANIFEST;manifest={}
                if manifest_path.exists():
                    prior=json.loads(manifest_path.read_text(encoding='utf-8'))
                    if prior.get('instrument')!='pixl':raise ValueError('目标目录已有不兼容的下载记录，未覆盖')
                    manifest=prior.get('files',{})
                self.update(scan_index=index+1,current=scan['file']);self.log('扫描 %d/%d：%s'%(index+1,len(scans),scan['file']))
                failures=[];plan=[]
                try:
                    for lid in [scan['lid'],*scan['rxl_lids']]:
                        if self.cancelled():raise Cancelled()
                        try:
                            meta=pixl.parse_lid(lid);refs=resolve(lid)
                            for ext in ('.csv','.xml'):
                                entry=refs.get((meta['stem']+ext).lower())
                                if not entry:raise ValueError('官方API缺文件直链：'+meta['stem']+ext)
                                plan.append((entry,meta['type']))
                            if images and meta['type']=='rxl':
                                plan.extend((entry,'image') for name,entry in refs.items() if name.startswith(('pcw_','pcd_','pc__')) and name.endswith(('.tif','.img','.tiff')))
                        except Exception as exc:failures.append(str(exc));self.log('文件引用失败：'+str(exc))
                    seen=set()
                    for entry,kind in plan:
                        if entry['name'].lower() in seen:continue
                        seen.add(entry['name'].lower())
                        if self.cancelled():raise Cancelled()
                        self.update(current=entry['name'],bytes_done=0,bytes_total=0)
                        try:
                            row=get_file(folder,entry,kind,manifest,self.cancelled,lambda n,t:self.update(bytes_done=n,bytes_total=t))
                            manifest[entry['name']]=row
                            appenv.atomic_json(manifest_path,{'instrument':'pixl','version':1,'scan':scan,'files':manifest,'updated_utc':utc()})
                            with self.lock:
                                self.job['files_done']+=1;self.job['reused']+=bool(row.get('reused'))
                                if kind=='rqc' and entry['name'].lower().endswith('.csv'):self.job['ready']=True
                            self.log(('已核验复用：' if row.get('reused') else '已下载：')+entry['name'])
                        except Cancelled:raise
                        except Exception as exc:failures.append(entry['name']+'：'+str(exc));self.log('失败：'+failures[-1])
                    _,obs,_=browser.scan(folder)
                    if images:
                        for ob in obs:
                            for rxl in ob['rxl']:
                                raw=browser.read_coordinates(rxl)
                                missing=[f['name'] for f in browser.frames(folder,rxl,raw) if f['missing']]
                                for name in missing:
                                    failure='缺引用MCC图像：'+name
                                    if failure not in failures:failures.append(failure);self.log(failure)
                    # File validation is independent; only a valid RQC unlocks opening this dataset.
                    ready=any(browser.chemistry(ob['path'])[0] for ob in obs)
                    with self.lock:self.job['ready']=self.job['ready'] or ready
                    if not scan['rxl_lids']:self.log('该成分产品未匹配RXL；已保留成分表，不猜坐标')
                except Cancelled:break
                except Exception as exc:failures.append(str(exc));self.log('扫描失败：'+str(exc))
                with self.lock:
                    self.job['results'].append({'scan':scan['id'],'folder':str(folder),'status':'failed' if failures else 'completed','errors':failures})
                    self.job['scan_done']+=1
            status='cancelled' if self.cancelled() else 'completed_with_errors' if any(r['errors'] for r in self.job['results']) else 'completed'
            self.update(status=status,finished_utc=utc())
        except Exception as exc:self.update(status='failed',error=str(exc));self.log('下载任务失败：'+str(exc))
        finally:self.update(running=False)
