# -*- coding: utf-8 -*-
"""极简 FITS 读取器（只服务 SuperCam 的 .fits 光谱产品）。

FITS 结构：若干 HDU，每个 = 2880 字节整数倍的 ASCII 头 + 数据块。
头是 80 字节一行的卡片，``END`` 卡收尾（头本身可跨多个 2880 块——
SuperCam 的主头就有 7 块，第一版解析器没读全块导致整个文件错位）。

只支持二进制表扩展（XTENSION='BINTABLE'）里的定长列：
B/I/J/K/E/D/A。SuperCam 光谱产品用到 E/D/A/K/B，够了；
遇到变长数组（P/Q）抛明确错误，提示装 astropy。

如果机器上装了 astropy，read_bintable 会优先用它兜底。
"""

from __future__ import annotations

import numpy as np

_BLOCK = 2880


class FitsError(Exception):
    pass


def _read_header(raw, off):
    """从 off 读一个 FITS 头。返回 (关键字 dict, 数据起始偏移)。"""
    cards = {}
    while True:
        block = raw[off:off + _BLOCK]
        if len(block) < _BLOCK:
            raise FitsError("FITS 文件在头中间就结束了（偏移 %d）" % off)
        done = False
        for i in range(36):
            card = block[i * 80:(i + 1) * 80].decode("ascii", "replace")
            key = card[:8].strip()
            if key == "END":
                done = True
                break
            if not key or key.startswith(("COMMENT", "HISTORY")):
                continue
            if card[8:10].rstrip() == "=":
                val = card[10:].split("/")[0].strip()
                if val.startswith("'") and val.count("'") >= 2:
                    val = val.split("'")[1].strip()
                cards[key] = val
        off += _BLOCK
        if done:
            return cards, off


_TFORM_DTYPES = {"B": ">u1", "I": ">i2", "J": ">i4", "K": ">i8",
                 "E": ">f4", "D": ">f8"}


def _column_dtype(tform):
    m = (tform or "").strip()
    count = ""
    while m[:1].isdigit():
        count += m[0]
        m = m[1:]
    code = m[:1].upper()
    repeat = int(count) if count else 1
    if code == "A":
        return ("S%d" % repeat).replace("S1", "S1"), 1 if repeat == 1 else repeat
    if code in ("P", "Q"):
        raise FitsError("变长数组列（%s）不支持，请 pip install astropy"
                        " 后重试（read_bintable 会自动改用 astropy）" % tform)
    if code not in _TFORM_DTYPES:
        raise FitsError("不认识的列类型 %r" % tform)
    if repeat == 1:
        return np.dtype(_TFORM_DTYPES[code]), 1
    # 重复定长：按结构化子列展开
    return np.dtype([("", _TFORM_DTYPES[code])] * repeat), repeat


def read_bintable(path, extname):
    """读指定 EXTNAME 的二进制表扩展。

    返回 dict：{列名(TTYPE): 一维 numpy 数组}。
    多重子列（罕见）会摊平成 名.子序号。
    """
    try:
        import astropy.io.fits as afits        # noqa: F401
    except ImportError:
        afits = None

    if afits is not None:
        with afits.open(path) as hdul:
            for hdu in hdul:
                if hdu.header.get("EXTNAME", "").strip().upper() == extname.upper():
                    data = hdu.data
                    if data is None:
                        return {}
                    return {str(name): np.asarray(data[name])
                            for name in hdu.columns.names}
        raise FitsError("找不到扩展 %r" % extname)

    raw = open(path, "rb").read()
    off = 0
    while off < len(raw):
        hdr, data_off = _read_header(raw, off)
        name = str(hdr.get("EXTNAME", "")).strip()
        if name.upper() == extname.upper() and "NAXIS1" in hdr:
            naxis1 = int(hdr["NAXIS1"])
            naxis2 = int(hdr["NAXIS2"])
            tfields = int(hdr.get("TFIELDS", "0") or "0")
            fields = []
            for i in range(1, tfields + 1):
                tname = str(hdr.get("TTYPE%d" % i, "col%d" % i)).strip()
                dt, _rep = _column_dtype(str(hdr.get("TFORM%d" % i, "E")))
                fields.append((tname, dt))
            dtype = np.dtype(fields)
            if dtype.itemsize != naxis1:
                raise FitsError("列宽 %d 与 NAXIS1=%d 对不上（%s）"
                                % (dtype.itemsize, naxis1, path))
            table = np.frombuffer(raw, dtype=dtype, count=naxis2,
                                  offset=data_off)
            out = {}
            for field_i, fname in enumerate(dtype.names, 1):
                col = table[fname]
                def physical(values):
                    # FITS 物理值 = 存储值*TSCAL + TZERO；无 Astropy 时也必须一致。
                    # 常见的无符号整数用有符号存储 + 2^(位数-1) 零点。
                    if values.dtype.kind not in "iuf":
                        return values
                    scale = float(str(hdr.get("TSCAL%d" % field_i, "1")).replace("D", "E"))
                    zero = float(str(hdr.get("TZERO%d" % field_i, "0")).replace("D", "E"))
                    if scale == 1 and zero == 0:
                        return values
                    bits = values.dtype.itemsize * 8
                    if values.dtype.kind == "i" and scale == 1 and zero == 2 ** (bits - 1):
                        unsigned = values.dtype.str.replace("i", "u")
                        return np.bitwise_xor(values.view(unsigned), np.array(2 ** (bits - 1), dtype=unsigned))
                    return values.astype(np.float64) * scale + zero
                if col.dtype.names:               # 多重子列摊平
                    for j, sub in enumerate(col.dtype.names):
                        out["%s.%d" % (fname, j + 1)] = physical(col[sub])
                else:
                    out[fname] = physical(col)
            return out
        # 跳过这个 HDU 的数据块
        naxis1 = int(hdr.get("NAXIS1", "0") or "0")
        naxis2 = int(hdr.get("NAXIS2", "0") or "0")
        pcount = int(hdr.get("PCOUNT", "0") or "0")
        gcount = int(hdr.get("GCOUNT", "1") or "1")
        ext = str(hdr.get("XTENSION", "")).upper()
        if "TABLE" in ext:
            data_size = (naxis1 * naxis2 + pcount) * gcount
        elif naxis1 and naxis2:
            data_size = naxis1 * naxis2 * max(abs(int(hdr.get("BITPIX", "8"))) // 8, 1)
        else:
            data_size = 0
        off = data_off + ((data_size + _BLOCK - 1) // _BLOCK) * _BLOCK
    raise FitsError("找不到扩展 %r（%s）" % (extname, path))
