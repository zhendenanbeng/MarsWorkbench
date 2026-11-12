# -*- coding: utf-8 -*-
"""覆盖真实文件名、发号对应和启动接口；不下载数据、不打开窗口。"""
from __future__ import annotations

import sys as _path_sys
from pathlib import Path as _BootPath
_path_sys.path.insert(0,str(_BootPath(__file__).resolve().parents[3]/'02_图形界面代码'))
import source_paths as _layout

import json
import os
import sys
import threading
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from unittest.mock import patch

import numpy as np

ROOT = str(_layout.project_layout.APP)
sys.path.insert(0, ROOT)

import console
import launch
import view
from instruments import REGISTRY
from sherloc import qc
from sherloc.organize import scan_folder

FIX = str(_layout.project_layout.VERIFY/'fixtures')
TOTAL = FAIL = 0


def check(condition, message):
    global TOTAL, FAIL
    TOTAL += 1
    if not condition:
        FAIL += 1
        print("  FAIL " + message)


def check_download_types():
    backend = REGISTRY["supercam"]
    names = [
        "scam_0998_0755522849_085_cl1_scam02998_aegis_0998a__________06p01.fits",
        "scam_0998_0755522849_085_cl1_scam02998_aegis_0998a__________06p01.xml",
        "SCAM_1859_0830000000_838_CP2_SCAM01859_TARGET_WITH_UNDERSCORES__01P01.FITS",
        "not_a_supercam_product.xml",
    ]
    result = backend.parse_listing([{"name": name, "size": 10} for name in names])
    groups = {g["key"]: g for g in result["groups"]}
    check(set(groups) == {"CL1", "CP2", "OTHER"}, "产品类型必须排除时间小数段")
    check(groups.get("CL1", {}).get("count") == 2, "CL1 必须包含成对 FITS/XML")
    check(result["preset"] == ["CL1"], "看图预选必须选中 CL1")
    check(result["total_files"] == 4, "未知产品也必须保留在下载清单")


def check_shots():
    for bid in ("chemcam", "supercam"):
        backend = REGISTRY[bid]
        observations = backend.scan_folder(os.path.join(FIX, bid))
        observations = [o for o in observations if not getattr(o,'passive',False)]
        if bid == "supercam":
            # fixtures 里也有无发数的被动 CP2/CP3 观测，发号检查只对 LIBS 生效
            observations = [o for o in observations if o.product == "CL1"]
        point = observations[0].points[0]
        dust = backend.load_point(point, 0, {"dust": True})
        raw6 = backend.load_point(point, 5, {})
        raw1 = backend.load_point(point, 0, {})
        check(np.array_equal(dust["y"], raw6["y"]), bid + " 除尘后首发谱值仍是原始第6发")
        check(not np.array_equal(dust["y"], raw1["y"]), bid + " 原始第1发与第6发不能混淆")
        check(dust["shot_numbers"] == list(range(6, 31)), bid + " 除尘后原始发号应为6–30")
        check(raw1["shot_numbers"] == list(range(1, 31)), bid + " 未除尘应保留1–30")
        plotted = console.plot_for({"observations": observations}, backend, 0, 0, 0, {"dust": True})
        check(plotted["series"][0]["label"] == "原始第 6 发", bid + " 曲线标签必须指向原始第6发")
        check(plotted["shots"][1] == {"i": 0, "label": "原始第 6 发"}, bid + " 下拉与曲线发号一致")
        check(plotted["log"][-1]["text"] == "原始第 6 发", bid + " 处理日志也保留原始发号")
        # 少于九发时，现有除尘策略不会删除；不能只按开关硬加5。
        short = backend._process(np.arange(4), np.arange(32).reshape(4, 8),
                                 0, {"dust": True}, [], point)
        check(short["shot_numbers"] == list(range(1, 9)), bid + " 除尘未执行时不得偏移发号")


def check_apxs():
    backend = REGISTRY["apxs"]
    observations = backend.scan_folder(os.path.join(FIX, "apxs"))
    d = console.plot_for({"observations": observations}, backend, 0, 0, -1, {})
    check(d["series"][0]["label"] == "全谱", "APXS 曲线不得写激光平均发数")
    check(d["shots"] == [{"i": -1, "label": "全谱"}], "APXS 下拉只表示全谱")
    check(len(d["oxides"]["rows"]) == 16, "APXS 成分表保持16项")


def check_mineral_summary():
    folder = os.path.join(os.path.dirname(ROOT), "03_下载数据/SHERLOC/练习数据",
                          "processed", "sol_00293")
    _, scans = scan_folder(folder)
    scan = next(sc for sc in scans if sc.spectra is not None)
    expected = qc.analyse_spectra(*scan.read_regions())
    actual = view.scan_minerals({"mineral_cache": {}}, scan)
    check(actual is not None, "界面必须产生矿物摘要")
    counts = [(g["classes"][0]["key"], g["n_points"])
              for g in actual["groups"] if g["classes"]]
    want = [(g["classes"][0]["key"], g["n_points"])
            for g in expected["mineral_groups"] if g["classes"]]
    check(counts == want, "界面矿物类别与点数必须与质控/地图一致")
    check(any(k == "carbonate" and 26 <= n <= 62 for k, n in counts),
          "Sol293 界面摘要必须保留已验证的碳酸盐信号")
    check(actual["n_points_with_band"] == expected["n_points_with_band"],
          "界面合格波段点数必须与质控一致")
    want_points = [sorted({c["class"] for c in p["classes"]})
                   for p in expected["point_detail"]]
    check(actual["point_classes"] == want_points, "全部逐点标签必须与质控一致")


def check_launcher():
    class ReadyHandler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def do_GET(self):
            self.send_response(200)
            self.end_headers()
            self.wfile.write(json.dumps({"build": view.BUILD}).encode())

    service = ThreadingHTTPServer(("127.0.0.1", 0), ReadyHandler)
    threading.Thread(target=service.serve_forever, daemon=True).start()
    port = service.server_port
    expected = "http://127.0.0.1:%d/" % port
    server = ThreadingHTTPServer(("127.0.0.1", 0), launch.make_handler({
        "launcher_port": 8790, "console_port": 8766, "sherloc_port": port}))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = "http://127.0.0.1:%d/" % server.server_port
    try:
        with patch.object(launch.subprocess, "Popen") as spawn:
            with urllib.request.urlopen(base + "api/open/sherloc") as response:
                result = json.load(response)
            check(result.get("url") == expected, "启动接口必须返回实际 SHERLOC 端口")
            check(not spawn.called, "当前版本服务已就绪时不得重复启动")
        with urllib.request.urlopen(base) as response:
            page = response.read().decode()
        check('href="' + expected + '"' in page, "启动卡片地址必须与接口端口一致")
        check('href="http://127.0.0.1:8790/"' not in page, "SHERLOC 卡片不能跳回启动器")
        with patch.object(launch, "ensure_sherloc", side_effect=RuntimeError("测试启动失败")):
            try:
                urllib.request.urlopen(base + "api/open/sherloc")
                check(False, "启动失败必须返回非成功状态")
            except urllib.error.HTTPError as exc:
                check(exc.code == 503, "启动失败返回503")
                check("测试启动失败" in json.load(exc)["error"], "启动错误必须传到前端")
        with urllib.request.urlopen(base + "api/open/sherloc") as response:
            check(json.load(response).get("url") == expected, "启动失败后仍可重试")
    finally:
        server.shutdown()
        server.server_close()
        service.shutdown()
        service.server_close()


def main():
    check_download_types()
    check_shots()
    check_apxs()
    check_mineral_summary()
    check_launcher()
    print("共 %d 项断言，%d 项失败" % (TOTAL, FAIL))
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
