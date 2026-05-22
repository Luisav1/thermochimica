import re, sys

path = "finalCaseBdata.txt"  # filename
msfl_re = re.compile(r'^\s*([0-9.]+(?:E[+-]\d+)?)\s+mol\s+MSFL')
p3_re   = re.compile(r'^\s*\+\s*([0-9.]+(?:E[+-]\d+)?)\s+mol\s+P3c1')
T_re    = re.compile(r'^\s*Temperature\s*=\s*([0-9.]+)\s*\[K\]')

rows = []
msfl = None
p3 = None

with open(path, "r", encoding="utf-8", errors="ignore") as f:
    for line in f:
        m = msfl_re.match(line)
        if m:
            msfl = float(m.group(1))
            continue
        m = p3_re.match(line)
        if m:
            p3 = float(m.group(1))
            continue
        m = T_re.match(line)
        if m:
            T = float(m.group(1))
            Tc = T - 273.15
            if p3 is None:
                p3 = 0.0
            rows.append((T, Tc, msfl, p3))
            msfl = None
            p3 = None

print("T_K T_C MSFL_mol P3c1_mol")
for T, Tc, msfl, p3 in rows:
    # T_C is integer for *.15 K temps
    print(f"{T:.2f} {int(round(Tc))} {msfl:.0f} {p3:.4f}")