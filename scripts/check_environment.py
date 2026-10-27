#!/usr/bin/env python3
import os,sys
print("Python:",sys.version.split()[0])
print("CPU count:",os.cpu_count())
try:
    import lhapdf
    print("LHAPDF: OK",getattr(lhapdf,"version",lambda:"unknown")())
    try:
        print("LHAPDF paths:",lhapdf.paths())
    except Exception: pass
except Exception as e:
    raise SystemExit("LHAPDF import FAILED. Run ./setup_macos.sh first.\n"+str(e))
for p in ["numpy","scipy","yaml","pandas","matplotlib"]:
    __import__(p)
    print(p,": OK")
