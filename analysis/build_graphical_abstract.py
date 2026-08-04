import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Ellipse, Rectangle

BLUE, ORANGE, GREEN, GREENL = "#2E86C1", "#E8721C", "#4B9B3F", "#D9EDCB"

# Times New Roman is not installed here; Tinos is metric-compatible with it.
from matplotlib import font_manager as fm
import glob
for _f in glob.glob("assets/fonts/Tinos-*.ttf"):
    fm.fontManager.addfont(_f)
TNR = "Tinos"

fig, ax = plt.subplots(figsize=(13.0, 7.6), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

def box(x, y, w, h, txt, fc, fs=12.5, tc="white", bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.6",
                                fc=fc, ec="none", zorder=3))
    ax.text(x+w/2, y+h/2, txt, ha="center", va="center", color=tc, zorder=4,
            fontsize=fs, fontweight="bold" if bold else "normal",
            fontfamily=TNR, linespacing=1.45)

def arrow(x, y0, y1):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>",
                mutation_scale=17, lw=1.9, color="#222", zorder=5,
                shrinkA=0, shrinkB=0))

def band(label, y):
    ax.text(2.0, y, label, ha="center", va="center", rotation=90, fontsize=13,
            fontweight="bold", color=GREEN, fontfamily=TNR, zorder=4)

ax.text(50, 95.5, "A multi-omics analysis of gastric-cancer prognosis",
        ha="center", va="center", fontsize=17.5, fontweight="bold", fontfamily=TNR)

# ---- row 1: datasets ----
band("DATASETS", 82)
box(9,  76, 26.5, 11, "TCGA + GTEx\n(transcriptome)",            BLUE)
box(37.5,76, 26.5, 11, "Three GEO cohorts\n(n = 922)",            ORANGE)
box(66, 76, 27.5, 11, "scRNA-seq · 16S microbiome\n· GWAS",       GREEN)

for x in (22.2, 50.7, 79.7): arrow(x, 75.0, 68.0)

# ---- row 2: analysis (frame hugs the boxes) ----
band("ANALYSIS", 58)
ax.add_patch(Rectangle((7.0, 51.0), 88.5, 16.0, fill=False, ec=ORANGE, lw=2.0, zorder=2))
box(9,  53.5, 26.5, 11, "Integrated differential\nexpression + WGCNA",  BLUE, fs=12)
box(37.5,53.5, 26.5, 11, "LASSO-Cox signature\n+ nested cross-validation", ORANGE, fs=11.5)
box(66, 53.5, 27.5, 11, "Module preservation ·\nMR · batch audit",       GREEN, fs=12)

arrow(50.7, 50.0, 44.0)

# ---- middle synthesis ----
box(24, 36.0, 53, 8.0, "Multi-omics triangulation & external validation", GREEN, fs=13)

arrow(29.5, 35.0, 28.0)
arrow(71.5, 35.0, 28.0)

# ---- row 3: results (separated, no overlap) ----
band("RESULTS", 16)
for cx, txt in [(29.5, "Stromal / CAF prognostic\nprogramme\n(externally validated)"),
                (71.5, "Microbiome signal =\nbatch artefact;\nMR null")]:
    ax.add_patch(Ellipse((cx, 16.0), 38, 21, fc=GREENL, ec=GREEN, lw=1.8, zorder=3))
    ax.text(cx, 16.0, txt, ha="center", va="center", fontsize=12.5, fontweight="bold",
            color="#1B5E20", fontfamily=TNR, linespacing=1.55, zorder=4)

fig.savefig("graphical_abstract.png", dpi=300, bbox_inches="tight",
            facecolor="white", pad_inches=0.22)
print("saved")
