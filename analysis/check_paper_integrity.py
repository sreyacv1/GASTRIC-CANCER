#!/usr/bin/env python3
"""Guard against silent content loss in PAPER.md.

Written after a caption-conversion pass (commit 1c9b0bc) deleted Results
sections 3.2-3.11 and 17 figure embeds, which then survived three commits
unnoticed. Run before committing PAPER.md.
"""
import re, os, sys
T = open("PAPER.md").read()
fail = []
n_sub = len(re.findall(r'^### 3\.\d+ ', T, re.M))
if n_sub != 11: fail.append(f"Results subsections: {n_sub}, expected 11")
n_main = len(re.findall(r'!\[Figure [1-8]\]', T))
if n_main != 8: fail.append(f"main figure embeds: {n_main}, expected 8")
n_supp = len(re.findall(r'!\[Supplementary Figure S\d+\]', T))
if n_supp != 11: fail.append(f"supplementary embeds: {n_supp}, expected 11")
for lab, p in re.findall(r'!\[((?:Supplementary )?Figure [^\]]+)\]\(([^)]+)\)', T):
    if not os.path.exists(p): fail.append(f"{lab}: missing file {p}")
for n in range(1, 38):
    if not re.search(rf'^{n}\. ', T, re.M): fail.append(f"reference {n} absent")
if "most immune-infiltrated" in T:
    fail.append("pre-correction EBV/MSI ranking wording present")
if "StromaScore" in T:
    fail.append("phantom StromaScore facet referenced")
print("\n".join("FAIL: " + f for f in fail) if fail else
      f"OK  {n_sub} Results subsections, {n_main} main + {n_supp} supplementary embeds, all files present")
sys.exit(1 if fail else 0)
