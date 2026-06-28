"""Mine ORCA's DFT manual -> {functional keyword: [references]} (offline)."""
import json
import re
from bs4 import BeautifulSoup

dft = BeautifulSoup(open("dft.html").read(), "html.parser")
idx = BeautifulSoup(open("manual_index.html").read(), "html.parser")

# ref-id (e.g. 'id924') -> {"number","text","doi"} from the global bibliography
def clean(s):
    return " ".join(s.split())

refs = {}
for a in dft.find_all("a", class_=re.compile("reference")):
    href = a.get("href", "")
    m = re.search(r"#(id\d+)$", href)
    if not m:
        continue
    rid = m.group(1)
    if rid in refs:
        continue
    el = idx.find(id=rid)
    if el is None:
        continue
    text = clean(el.get_text(" "))
    text = re.sub(r"^\[\s*\d+\s*\]\s*", "", text)  # drop leading [N]
    doi = re.search(r"10\.\d{4,}/\S+?(?=[\s.]*$|\s)", text)
    refs[rid] = {"text": text, "doi": (doi.group(0).rstrip(".") if doi else None)}

# functional keyword -> [ref-ids], from the tables
func_refs = {}
for table in dft.find_all("table"):
    header = [clean(th.get_text()) for th in table.find_all("th")]
    if not header:
        continue
    kw_hdr = "Keyword" if "Keyword" in header else ("Keywords" if "Keywords" in header else None)
    if kw_hdr is None or "Functional" not in header:
        continue
    ki = header.index(kw_hdr)
    fi = header.index("Functional")
    for tr in table.find_all("tr"):
        cells = tr.find_all("td")
        if len(cells) <= max(ki, fi):
            continue
        rids = []
        for a in cells[fi].find_all("a", class_=re.compile("reference")):
            m = re.search(r"#(id\d+)$", a.get("href", ""))
            if m and m.group(1) in refs:
                rids.append(m.group(1))
        if not rids:
            continue
        kw_text = clean(cells[ki].get_text(" "))
        for kw in re.split(r"[,\s]+", kw_text):
            kw = kw.strip()
            if kw:
                func_refs.setdefault(kw, [])
                for r in rids:
                    if r not in func_refs[kw]:
                        func_refs[kw].append(r)

json.dump({"refs": refs, "func_refs": func_refs}, open("dft_mined.json", "w"), indent=1)

print(f"functionals with refs: {len(func_refs)} | unique references: {len(refs)}")
print(f"references with a DOI: {sum(1 for r in refs.values() if r['doi'])}/{len(refs)}")
print("\n=== sample functionals ===")
for kw in ["B3LYP", "B3LYP/G", "PBE0", "wB97X-D3", "M062X", "M06-2X", "TPSSh", "PBE"]:
    if kw in func_refs:
        for rid in func_refs[kw]:
            r = refs[rid]
            print(f"{kw:10s} -> doi={r['doi']}  {r['text'][:75]}")
    else:
        print(f"{kw:10s} -> (keyword not found in tables)")
