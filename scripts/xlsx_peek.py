"""Read-only peek at an .xlsx with the standard library: sheet names, header row, row counts, and
optionally the rows whose first column matches a pattern. For reading downloaded supplements.

uv run python -I scripts/xlsx_peek.py <file.xlsx> [regex_on_any_cell] [max_rows]
"""

import re
import sys
import xml.etree.ElementTree as ET
import zipfile

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def main(path, pattern=None, max_rows=40):
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            ss.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}
    for sh in wb.find("m:sheets", NS):
        name, rid = sh.get("name"), sh.get(RNS)
        t = target[rid].lstrip("/")
        t = t if t.startswith("xl/") else "xl/" + t
        rows = ET.fromstring(z.read(t)).find("m:sheetData", NS).findall("m:row", NS)

        def vals(row):
            out = []
            for c in row.findall("m:c", NS):
                v = c.find("m:v", NS)
                x = v.text if v is not None else ""
                out.append(ss[int(x)] if c.get("t") == "s" and x else x)
            return out

        print(f"## sheet {name!r}: {len(rows)} rows; header: {vals(rows[0]) if rows else []}")
        if pattern == "--clusters":
            hdr = vals(rows[0])
            if "cluster" in hdr:
                j = hdr.index("cluster")
                from collections import Counter

                cnt = Counter(vals(r)[j] for r in rows[1:] if len(vals(r)) > j)
                print("   clusters:", dict(cnt))
        elif pattern:
            n = 0
            for r in rows[1:]:
                v = vals(r)
                if any(re.search(pattern, str(x)) for x in v):
                    print("  ", v)
                    n += 1
                    if n >= int(max_rows):
                        break


if __name__ == "__main__":
    main(*sys.argv[1:])
