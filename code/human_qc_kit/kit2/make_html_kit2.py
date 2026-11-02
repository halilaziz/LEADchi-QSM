#!/usr/bin/env python3
"""make_html_kit2.py - (kit 2; copy of kit-1 make_html.py with kit-2 paths; template read from kit 1 unchanged)
embed the rater item list (kit_id, order, forced_review; no AI information) into rater_kit2/rate.html.

Input : analysis/items_for_html.json (written by build_kit.py), analysis/rate_template.html
Output: rater_kit2/rate.html (single self-contained offline file; images are loaded by relative path images/...)
"""
import hashlib
import json
import os
import sys

K1 = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
KIT = K1 + "/kit2_2026-10-28"
items = json.load(open(os.path.join(KIT, "analysis", "items_for_html.json")))
assert [it["order"] for it in items] == list(range(1, len(items) + 1))
assert set(items[0]) == {"order", "kit_id", "forced_review"}, "item list must carry no AI information"
js = json.dumps(items, separators=(",", ":"))
tag = hashlib.sha256(js.encode()).hexdigest()[:12]
tpl = open(os.path.join(K1, "analysis", "rate_template.html"), encoding="utf-8").read()   # kit-1 template, unchanged
assert tpl.count("/*__ITEMS__*/[]") == 1 and tpl.count("/*__KITTAG__*/") == 1
RUBRIC = "/path/to/Penn_LEAD/LEADchi/prereg/qc_ai/RUBRIC_sections1-3.md"


def md2html(md):
    """Minimal converter for the rubric (headings, tables, lists, bold, code, paragraphs). Text is kept verbatim."""
    import html as H
    import re

    def inline(t):
        t = H.escape(t, quote=False)
        t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
        return re.sub(r"`(.+?)`", r"<code>\1</code>", t)
    out, lines, i = [], md.splitlines(), 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("## "):
            out.append(f"<h2>{inline(ln[3:])}</h2>")
        elif ln.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                if not re.fullmatch(r"\|[-|: ]+\|", lines[i].strip()):
                    rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            i -= 1
            out.append("<table><tr>" + "".join(f"<th>{inline(c)}</th>" for c in rows[0]) + "</tr>" +
                       "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in rows[1:]) +
                       "</table>")
        elif re.match(r"^(- |\d+\. )", ln):
            tag = "ol" if re.match(r"^\d+\. ", ln) else "ul"
            items = []
            while i < len(lines) and re.match(r"^(- |\d+\. |  +\S)", lines[i]):
                if re.match(r"^(- |\d+\. )", lines[i]):
                    items.append(re.sub(r"^(- |\d+\. )", "", lines[i]))
                else:
                    items[-1] += " " + lines[i].strip()
                i += 1
            i -= 1
            out.append(f"<{tag}>" + "".join(f"<li>{inline(x)}</li>" for x in items) + f"</{tag}>")
        elif ln.strip():
            para = [ln]
            while i + 1 < len(lines) and lines[i + 1].strip() and not re.match(r"^(## |\||- |\d+\. )", lines[i + 1]):
                i += 1
                para.append(lines[i])
            out.append("<p>" + inline(" ".join(para)) + "</p>")
        i += 1
    return "\n".join(out)


md = open(RUBRIC, encoding="utf-8").read()
sec23 = md[md.index("## 2."):]
rub = ("<h3>Registered rubric (qc_rubric.md), sections 2 and 3, verbatim</h3>"
       "<p class=\"hint\">Section 1 (blinding of the original procedure) is reproduced in README_rater.md.</p>" + md2html(sec23))
assert tpl.count("<!--__RUBRIC__-->") == 1
html = tpl.replace("/*__ITEMS__*/[]", js).replace("/*__KITTAG__*/", tag).replace("<!--__RUBRIC__-->", rub)
out = os.path.join(KIT, "rater_kit2", "rate.html")
if os.path.exists(out) and "--force" not in sys.argv:
    sys.exit(f"{out} exists (use --force to regenerate the HTML only)")
open(out, "w", encoding="utf-8").write(html)
print(f"wrote {out}: {len(items)} items, storage tag {tag}")
