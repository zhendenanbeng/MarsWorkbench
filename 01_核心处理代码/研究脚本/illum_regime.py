# -*- coding: utf-8 -*-
"""对比两台相机的照明/曝光配置（判断谁在被拍时的光源是什么）。"""
import os
import re

X = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "xml")


def info(name, tag):
    p = os.path.join(X, name)
    if not os.path.exists(p):
        print("%-52s （标签不在缓存里）" % tag)
        return
    t = open(p, encoding="utf-8", errors="replace").read()
    exp = re.search(r'<img:exposure_duration unit="ms">([\d.]+)', t)
    cnt = re.search(r"<img:exposure_count>\s*(\d+)", t)
    led = dict(re.findall(r"<img:device_name>([^<]+)</img:device_name>\s*"
                          r"<img:device_state>([^<]+)</img:device_state>", t))
    leds = {k: v for k, v in led.items() if "LED" in k}
    ms = float(exp.group(1)) if exp else float("nan")
    n = int(cnt.group(1)) if cnt else 1
    print("%-46s 曝光 %7.1f ms × %-5d = %8.1f s   LED %s"
          % (tag, ms, n, ms * n / 1000.0, leds))


info("SIF_0185_0683368184_652FDR_N0070000SRLC00720_0000LMJ04.xml",
     "WATSON 论文那张 (sol185)")
info("SIF_0186_0683466024_335FDR_N0070000SRLC00621_0000LMJ04.xml",
     "WATSON srlc00621 全画幅 (sol186)")
info("SIF_0186_0683466037_832FDR_N0070000SRLC00621_0000LUJ05.xml",
     "WATSON srlc00621 LUJ (sol186)")
info("SC3_0186_0683479674_054FDR_N0070000SRLC11420_0000LMJ05.xml", "ACI (sol186)")
info("SC3_0186_0683478583_214FDR_N0070000SRLC11420_0000LMJ05.xml",
     "ACI 打激光前那张 (sol186)")
