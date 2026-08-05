#!/usr/bin/env python3
"""Guard against silent content loss in PAPER.md.

Written after a caption-conversion pass (commit 1c9b0bc) deleted Results
sections 3.2-3.11 and 17 figure embeds, which then survived three commits
unnoticed. Run before committing PAPER.md.
"""
import re, os, sys
T = open("PAPER.md").read()
fail = []
note = []
n_sub = len(re.findall(r'^### 3\.\d+ ', T, re.M))
if n_sub != 11: fail.append(f"Results subsections: {n_sub}, expected 11")
n_main = len(re.findall(r'!\[Figure [1-8]\]', T))
if n_main != 8: fail.append(f"main figure embeds: {n_main}, expected 8")
n_supp = len(re.findall(r'!\[Supplementary Figure S\d+\]', T))
if n_supp != 9: fail.append(f"supplementary embeds: {n_supp}, expected 9")
for lab, p in re.findall(r'!\[((?:Supplementary )?Figure [^\]]+)\]\(([^)]+)\)', T):
    if not os.path.exists(p): fail.append(f"{lab}: missing file {p}")
for n in range(1, 38):
    if not re.search(rf'^{n}\. ', T, re.M): fail.append(f"reference {n} absent")
if "most immune-infiltrated" in T:
    fail.append("pre-correction EBV/MSI ranking wording present")
if "StromaScore" in T:
    fail.append("phantom StromaScore facet referenced")

# Split-package consistency: the section files and their .docx are generated from
# PAPER.md, so a PAPER.md edit that is not re-split ships a stale deliverable.
SPLIT = "submission_package/split_manuscript"
if os.path.isdir(SPLIT):
    import subprocess
    src_mtime = os.path.getmtime("PAPER.md")
    for g in sorted(os.listdir(SPLIT)):
        if re.match(r'^\d\d_.*\.(md|docx)$', g):
            if os.path.getmtime(os.path.join(SPLIT, g)) < src_mtime - 1:
                fail.append(f"{g} is older than PAPER.md - re-split the package")
    rm = os.path.join(SPLIT, "README.md")
    if os.path.exists(rm) and "referenced, not embedded" in open(rm).read():
        fail.append("README claims figures are not embedded, but they are")

# Interior legibility SCREEN (advisory, not a gate).
# Edge-ink and aspect checks pass on a figure whose text is too small to read, so
# this reports the median glyph height of each embedded figure converted to points
# at 180 mm print width, and flags anything below a reference value for a human to
# look at.
# CALIBRATION AND ITS LIMIT. The metric orders a known real pair correctly: the
# 7-facet s15_immune montage at git 4a2c90d scored 1.89 pt and was illegible, its
# replacement 2.27 pt and reads cleanly. It also caught two genuine defects - S9
# and S10 scored 1.53-1.56 pt from over-narrow 3-across layouts and were
# relayouted. But it CANNOT separate bad from acceptable at the boundary: the
# relayouted S9/S10 also score 1.89 pt and are legible by eye, the same value as
# the known-illegible montage. So the threshold is set at 2.30 pt, which flags
# everything at or below the legible-but-dense band and accepts nothing on trust.
# Expect standing notes on figures whose small type is correct by convention
# (Fig 5's WGCNA dendrogram leaf labels, Fig 6's heatmap gene rows) and on the
# relayouted S9/S10. A note means LOOK AT IT, never a failure.
LEGIBILITY_REF_PT = 2.30  # see calibration note above
try:
    import numpy as np
    from PIL import Image
    from scipy import ndimage
    PRINT_IN = 7.09  # 180 mm double-column
    low = []
    for lab, p in re.findall(r'!\[((?:Supplementary )?Figure [^\]]+)\]\(([^)]+)\)', T):
        if not os.path.exists(p):
            continue
        im = Image.open(p)
        arr = np.array(Image.alpha_composite(
            Image.new("RGBA", im.size, (255, 255, 255, 255)),
            im.convert("RGBA")).convert("L"))
        lb, _ = ndimage.label(arr < 128)
        hs = [(s0.stop - s0.start) for s0, s1 in ndimage.find_objects(lb)
              if 5 <= s0.stop - s0.start <= 60 and 3 <= s1.stop - s1.start <= 45
              and 0.6 <= (s0.stop - s0.start) / max(s1.stop - s1.start, 1) <= 6]
        if len(hs) < 40:
            continue
        pt = float(np.median(hs)) / (im.size[0] / PRINT_IN) * 72
        if pt <= LEGIBILITY_REF_PT:
            low.append((pt, lab, im.size))
    for pt, lab, size in sorted(low):
        note.append(f"legibility screen: {lab} median text {pt:.2f} pt at 180 mm "
                    f"({size[0]}x{size[1]}) - at/below the {LEGIBILITY_REF_PT} pt "
                    f"reference; view it and confirm the tick labels read")
except ImportError:
    pass

for n in note:
    print("NOTE:", n)
print("\n".join("FAIL: " + f for f in fail) if fail else
      f"OK  {n_sub} Results subsections, {n_main} main + {n_supp} supplementary embeds, all files present")
sys.exit(1 if fail else 0)
