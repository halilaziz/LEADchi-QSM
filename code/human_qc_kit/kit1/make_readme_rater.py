#!/usr/bin/env python3
"""make_readme_rater.py - assemble rater_kit/README_rater.md: instructions (README_rater_head.md) + the registered
rubric copied byte-for-byte (sections 2-3 as the operative rubric; section 1 in an appendix for context)."""
import os
import sys

KIT = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
RUB = "/path/to/Penn_LEAD/LEADchi/prereg/qc_ai/RUBRIC_sections1-3.md"
md = open(RUB, encoding="utf-8").read()
i2 = md.index("## 2.")
sec1, sec23 = md[:i2], md[i2:]
head = open(os.path.join(KIT, "analysis", "README_rater_head.md"), encoding="utf-8").read()
out = (head + "\n" + sec23.rstrip("\n") + "\n\n<!-- END VERBATIM -->\n\n---\n\n"
       "## Appendix: registered rubric section 1 (context only; verbatim)\n\n"
       "This section describes the blinding of the original grading round. A later amendment replaced its single human "
       "rater by the automated rater that this re-rating evaluates. Its procedural rules (who grades, where grades are "
       "entered) do not apply to you; sections 1 to 6 of this page do.\n\n"
       "<!-- BEGIN VERBATIM: qc_rubric.md section 1 -->\n\n" + sec1.rstrip("\n") + "\n\n<!-- END VERBATIM -->\n")
assert sec23 in out and sec1.rstrip("\n") in out
dst = os.path.join(KIT, "rater_kit", "README_rater.md")
if os.path.exists(dst) and "--force" not in sys.argv:          # exists-guard added 2026-10-07 (review item L10)
    sys.exit(f"{dst} exists (use --force to regenerate README_rater.md only)")
open(dst, "w", encoding="utf-8").write(out)
print("wrote README_rater.md", len(out), "chars")
