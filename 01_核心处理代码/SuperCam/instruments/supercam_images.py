"""SuperCam CI and WATSON context images. Display only; no spatial association."""
import io
import json
import os
from pathlib import Path
import re
import threading
import urllib.request
from datetime import datetime, timezone

import numpy as np
from PIL import Image
from sherloc import watson
from sherloc.vicar import _parse_label
from .fitsio import _read_header
from . import pds_listing
from .io_integrity import check_length
import appenv

ROOT = Path(__file__).resolve().parents[1]
CACHE = appenv.cache_write('scam_images')
CI_ROOT = 'https://pds-geosciences.wustl.edu/m2020/urn-nasa-pds-mars2020_supercam/data_radcal_rmi'
CI_RX = re.compile(r'^scam_(\d{4})_(\d{10})_(\d{3})_ci__(scam\d{5})_(.+)_(\d{2})p(\d{2})\.fits$', re.I)
SCOPE = '仅浏览：不叠加、不投影、不做光斑判定；同Sol不表示共靶。WATSON来自另一台仪器。'
LOCK = threading.RLock()


def sol_number(value):
    s = str(value)
    if not re.fullmatch(r'(?:sol_?)?\d{1,5}', s, re.I):
        raise ValueError('火星日格式无效')
    n = int(re.sub(r'^sol_?', '', s, flags=re.I))
    if n > 9999:
        raise ValueError('图像火星日超出产品命名范围')
    return n


def parse_ci(name):
    m = CI_RX.fullmatch(name)
    if not m:
        return None
    return dict(name=name, sol=int(m[1]), sclk=m[2]+'.'+m[3], seq=m[4].lower(),
                target=m[5].rstrip('_'), point=int(m[6]), version=int(m[7]), kind='CI')


def _atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    import tempfile
    fd, temp = tempfile.mkstemp(dir=path.parent, suffix='.part')
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(data)
        os.replace(temp, path)
    finally:
        if os.path.exists(temp):
            os.unlink(temp)


def ci_frames(sol, refresh=False):
    n = sol_number(sol)
    cache = CACHE / ('sol%05d' % n) / 'listing.json'
    with LOCK:
        try:
            legacy = appenv.SOURCE/'.cache/scam_images'/('sol%05d'%n)/'listing.json'
            if CACHE == appenv.cache_write('scam_images') and not cache.exists() and legacy.exists() and not refresh:
                cache=legacy
            if cache.exists() and not refresh:
                saved = json.loads(cache.read_text(encoding='utf-8'))
                rows = saved['frames']
                if not all(parse_ci(r['name']) and parse_ci(r['name'])['sol'] == n for r in rows):
                    raise ValueError('缓存清单内容无效，请刷新')
                return rows
            entries = pds_listing.list_files(CI_ROOT, 'sol_%05d' % n)
            rows = []
            for e in entries:
                r = parse_ci(e['name'])
                if r and r['sol'] == n:
                    rows.append(dict(r, size=e['size']))
            rows.sort(key=lambda r: (r['seq'], r['sclk'], r['name']))
            _atomic(cache, json.dumps({'fetched_utc': datetime.now(timezone.utc).isoformat(), 'frames': rows}, ensure_ascii=False).encode())
            return rows
        except Exception as exc:
            raise ValueError('CI目录列表失败（可稍后重试；已有图像缓存仍可使用）：%s' % exc) from exc


def watson_groups(sol):
    """Local inventory only. Never trigger load_inventory's network fallback."""
    n = sol_number(sol)
    if not Path(watson._inventory_path()).exists():
        raise ValueError('WATSON本地清单未缓存；请先准备watson_inventory.csv，当前列表不自动联网抓取全任务清单')
    groups = {}
    for stem in sorted(set(watson.load_inventory())):
        r = watson.parse_name(stem)
        if not r or int(r['sol']) != n or r['prefix'] != 'SIF':
            continue
        m = re.fullmatch(r'([A-Z]\d{7})(SRLC\d+)', r['att'])
        if not m:
            continue
        key = m[2].lower() + ' / ' + m[1].lower()
        group = groups.setdefault(key, {'id': key, 'seq': m[2].lower(), 'placement': m[1].lower(), 'frames': []})
        group['frames'].append({'stem': stem, 'sclk': r['sclk'], 'kind': r['kind'], 'encoding': r['enc']})
    return [dict(g, count=len(g['frames']), representative=g['frames'][0]['stem']) for _, g in sorted(groups.items())]


def listing(sol, refresh=False):
    n = sol_number(sol)
    result = {'sol': n, 'scope': SCOPE, 'ci': [], 'watson': []}
    for key, fetch in [('ci', lambda: ci_frames(n, refresh)), ('watson', lambda: watson_groups(n))]:
        try:
            result[key] = fetch()
        except Exception as exc:
            result[key+'_error'] = str(exc)
    return result


def fits_image(path):
    """First 2D image HDU, band zero if 3D; native scaling and BLANK respected."""
    raw = Path(path).read_bytes()
    off = 0
    types = {8: '>u1', 16: '>i2', 32: '>i4', 64: '>i8', -32: '>f4', -64: '>f8'}
    while off < len(raw):
        h, start = _read_header(raw, off)
        ndim = int(h.get('NAXIS', 0))
        axes = [int(h.get('NAXIS%d' % i, 0)) for i in range(1, ndim+1)]
        bitpix = int(h.get('BITPIX', 8))
        size = (int(np.prod(axes)) if axes else 0) * abs(bitpix)//8 + int(h.get('PCOUNT', 0))
        if ndim in (2, 3) and all(a > 0 for a in axes) and 'TABLE' not in h.get('XTENSION', ''):
            if bitpix not in types or start + size > len(raw):
                raise ValueError('FITS类型不支持或图像栅格截断')
            original = np.frombuffer(raw, types[bitpix], count=axes[0]*axes[1], offset=start).reshape(axes[1], axes[0])
            number = lambda v: float(str(v).replace('D', 'E'))
            a = original.astype(float)*number(h.get('BSCALE', 1))+number(h.get('BZERO', 0))
            if bitpix > 0 and 'BLANK' in h:
                a[original == int(h['BLANK'])] = np.nan
            return a
        off = start + ((size + 2879)//2880)*2880
    raise ValueError('FITS中没有可浏览的二维图像')


def vicar_image(path):
    """Read IMAGE object's band zero, attached record/byte pointer, BSQ only."""
    raw = Path(path).read_bytes()
    end = re.search(rb'(?m)^END\s*\r?$' , raw[:1 << 20])
    if not end:
        raise ValueError('未找到完整ODL标签')
    text = raw[:end.end()].decode('ascii', 'replace')
    block = re.search(r'\bOBJECT\s*=\s*IMAGE\s+(.*?)END_OBJECT\s*=\s*IMAGE\b', text, re.S)
    if not block:
        raise ValueError('ODL没有IMAGE对象')
    label = _parse_label(block[1])
    head = _parse_label(text[:block.start()])
    lines, samples = int(label['LINES']), int(label['LINE_SAMPLES'])
    bands, bits = int(label.get('BANDS', 1)), int(label['SAMPLE_BITS'])
    storage = str(label.get('BAND_STORAGE_TYPE', 'BAND_SEQUENTIAL'))
    if storage not in ('BAND_SEQUENTIAL', 'BSQ'):
        raise ValueError('仅支持BSQ波段排列，未猜测其他排列')
    st = str(label['SAMPLE_TYPE']).upper()
    if st in ('IEEE_REAL', 'MSB_REAL'):
        dtype = '>f%d' % (bits//8)
    elif st in ('PC_REAL', 'LSB_REAL'):
        dtype = '<f%d' % (bits//8)
    elif st in ('MSB_INTEGER', 'MSB_UNSIGNED_INTEGER', 'LSB_INTEGER', 'LSB_UNSIGNED_INTEGER', 'UNSIGNED_INTEGER'):
        dtype = ('<' if 'LSB' in st else '>') + ('u' if 'UNSIGNED' in st else 'i') + str(bits//8)
    else:
        raise ValueError('不支持的SAMPLE_TYPE：'+st)
    if bits not in (8, 16, 32, 64) or lines <= 0 or samples <= 0 or bands <= 0:
        raise ValueError('ODL图像维度/位深无效')
    pointer = re.search(r'(?m)^\^IMAGE\s*=\s*(\d+)\s*(<BYTES>)?', text)
    if not pointer:
        raise ValueError('没有可核验的附加IMAGE指针，未用文件尾猜偏移')
    offset = (int(pointer[1])-1)*(1 if pointer[2] else int(head['RECORD_BYTES']))
    need = lines*samples*bands*(bits//8)
    if offset < end.end() or offset+need > len(raw):
        raise ValueError('ODL图像指针越界或栅格截断')
    a = np.frombuffer(raw, dtype=dtype, count=lines*samples, offset=offset).reshape(lines, samples).astype(float)
    for key in ('INVALID_CONSTANT', 'MISSING_CONSTANT'):
        if key in label:
            a[a == float(label[key])] = np.nan
    return a


def stretch_png(array):
    a = np.asarray(array, dtype=float)
    valid = np.isfinite(a)
    if a.ndim != 2 or not valid.any():
        raise ValueError('没有可显示的有限灰度像素')
    lo, hi = np.percentile(a[valid], [1, 99])
    out = np.zeros(a.shape, dtype=np.uint8)
    out[valid] = np.clip((a[valid]-lo)/(hi-lo)*255, 0, 255).astype(np.uint8) if hi > lo else 128
    b = io.BytesIO()
    Image.fromarray(out).save(b, format='PNG')
    return b.getvalue()


def _download(url, path):
    req = urllib.request.Request(url, headers={'User-Agent': 'sherloc-workbench/1.0'})
    # File stays private until HTTP length and format conversion both succeed.
    with urllib.request.urlopen(req, timeout=120) as response:
        data = response.read()
        check_length(response, len(data))
    _atomic(path, data)


def ci_png(name, sol):
    n = sol_number(sol)
    row = parse_ci(name)
    if not row or row['sol'] != n:
        raise ValueError('CI产品名与火星日不匹配')
    with LOCK:
        folder = CACHE / ('sol%05d' % n)
        png = folder / (name+'.png')
        old_png=appenv.SOURCE/'.cache/scam_images'/('sol%05d'%n)/(name+'.png')
        if not png.exists() and old_png.exists():return old_png.read_bytes()
        if png.exists():
            return png.read_bytes()
        path = folder / name
        old_raw=appenv.SOURCE/'.cache/scam_images'/('sol%05d'%n)/name
        if not path.exists() and old_raw.exists():path=old_raw
        old = ROOT / '.cache/phase8_probe' / ('sol_%05d_ci' % n) / name
        if old.exists() and not path.exists():
            path = old
        if not path.exists():
            try:
                _download(CI_ROOT+'/sol_%05d/' % n+name, path)
            except Exception as exc:
                raise ValueError('CI下载失败：%s' % exc) from exc
        try:
            data = stretch_png(fits_image(path))
        except Exception as exc:
            raise ValueError('CI转换失败（缓存原件保留）：%s' % exc) from exc
        _atomic(png, data)
        return data


def watson_png(stem):
    if not re.fullmatch(r'sif_[a-z0-9_\-]+', stem) or not watson.parse_name(stem):
        raise ValueError('WATSON产品名无效')
    with LOCK:
        path = Path(watson.CACHE) / 'img' / (stem+'.IMG')
        png = CACHE / 'watson_png' / (stem+'.png')
        if png.exists():
            return png.read_bytes()
        if not path.exists():
            url = watson.lookup_url(stem)
            if not url:
                raise ValueError('WATSON API查询失败或无直链，请稍后重试（可能限速或网络不可用）')
            try:
                _download(url, path)
            except Exception as exc:
                raise ValueError('WATSON下载失败：%s' % exc) from exc
        try:
            data = stretch_png(vicar_image(path))
        except Exception as exc:
            raise ValueError('WATSON转换失败（原件保留）：%s' % exc) from exc
        _atomic(png, data)
        return data



