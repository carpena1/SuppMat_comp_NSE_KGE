"""
Figure A3 (proposed): diagnostic pathway for simulations relative to the
+/-5% k_s uncertainty band (Eq. 11) around the simplified NSE-KGE* relation (Eq. 7).

Cases inside UB (Eq. 7 applies) are not classified. Each case outside UB falls in exactly one of:
  1   outside UB, |1-k_s| <= 0.05        bias-driven
  2a  outside UB, 1-k_s < -0.05              imbalance: correlation > variability (alpha < 0.952)
  2b  outside UB, 1-k_s >  0.05, alpha < 1   imbalance: variability underestimated
  2c  outside UB, 1-k_s >  0.05, alpha >= 1  imbalance: variability overestimated
(k_s = rho/alpha; all conditions written in terms of 1 - k_s)

(a) Decision chart with CAMELS counts.
(b-e) Theoretical class regions in the (1 - k_s, |beta_n|) plane at fixed NSE
      (large-n limit; alpha = sqrt((NSE + beta_n^2)/(2k_s - 1)), rho = k_s alpha,
      KGE* from Eq. 4), with CAMELS cases overlaid in the panel nearest their NSE.

Usage:  python figA5_diagnostic_map.py [NSEvsKGE_datasets_camels_48.txt | metrics.csv] [out.png | out.pdf]   (.pdf/.svg/.eps -> fully vector, no raster)
        (CSV columns: NSE, KGEs, rho, alpha, beta_n, N)
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import FancyBboxPatch, Patch

plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8,
                     "mathtext.fontset": "custom", "mathtext.it": "Liberation Sans:italic",
                     "mathtext.rm": "Liberation Sans", "hatch.linewidth": 0.6})
EPS = 0.05

# ================================================================ core logic
def band(nse):
    """Eq. (11): (lower, upper) KGE* bounds for +/-5% k_s."""
    b1 = 1 - np.sqrt(2.114 * nse - 4.111 * np.sqrt(nse) + 2)   # k_s = 0.95
    b2 = 1 - np.sqrt(1.911 * nse - 3.909 * np.sqrt(nse) + 2)   # k_s = 1.05
    return np.minimum(b1, b2), np.maximum(b1, b2)

def theory(nse, ks, bn):
    alpha = np.sqrt((nse + bn**2) / (2 * ks - 1))
    rho = ks * alpha
    return 1 - np.sqrt((rho - 1)**2 + (alpha - 1)**2 + bn**2), rho, alpha

def classify(nse, kge, rho, alpha):
    """1:'0' | 2:"0'" | 3:'1' | 4:'2a' | 5:'2b' | 6:'2c'"""
    ks = rho / alpha
    lo, hi = band(nse)
    t = 1e-9
    inside = (kge >= lo - t) & (kge <= hi + t)
    bal = np.abs(1 - ks) <= EPS + t
    c = np.zeros(np.shape(ks), int)
    c = np.where(inside & bal, 1, c)
    c = np.where(inside & ~bal, 2, c)
    c = np.where(~inside & bal, 3, c)
    c = np.where(~inside & ~bal & (ks > 1 + EPS), 4, c)
    c = np.where(~inside & ~bal & (ks < 1 - EPS) & (alpha < 1), 5, c)
    c = np.where(~inside & ~bal & (ks < 1 - EPS) & (alpha >= 1), 6, c)
    return c

CODE = {1: "0", 2: "0???", 3: "1", 4: "2a", 5: "2b", 6: "2c"}

# ================================================================ data / counts
# Defaults = CAMELS values (overwritten when the data file is given)
COUNTS = {1: 236, 2: 31, 3: 11, 4: 132, 5: 164, 6: 39}
SIDE = dict(above=202, below=144, below4=132, below3=11, below_other=1, above56=202, n56=203)
DK_MED, BACK = 0.003, 1

cases = None
if len(sys.argv) > 1 and sys.argv[1].endswith((".csv", ".txt")):
    import pandas as pd
    if sys.argv[1].endswith(".txt"):          # raw FITEVAL table (NSEvsKGE_datasets_camels_48.txt)
        df = pd.read_csv(sys.argv[1], sep="\t")
        df.columns = [c.strip() for c in df.columns]
        df = df[df.NSE > 0.2].copy()
        df["alpha"] = df.sd_Yp / df.sd_Yo
        df["rho"] = df.r
        df["beta_n"] = (df.mean_Yp - df.mean_Yo) / df.sd_Yo
        df["KGEs"] = df["KGE*"]
    else:
        df = pd.read_csv(sys.argv[1])
    df["KGEs"] = df.KGEs.round(3)                  
    df["ks"] = df.rho / df.alpha
    df["cls"] = classify(df.NSE.values, df.KGEs.values, df.rho.values, df.alpha.values)
    hi = (df.KGEs >= 0.95) & (df.cls >= 3)         # KGE* >= 0.95 is never counted outside
    df.loc[hi, "cls"] = np.where((df.ks[hi] - 1).abs() <= EPS, 1, 2)
    lo_, hi_ = band(df.NSE)
    df["side"] = np.where(df.cls <= 2, "in", np.where(df.KGEs > hi_, "above", "below"))
    # bias-removal check: beta_n = 0, rho and alpha unchanged (Eqs. 4-5)
    fac = df.N / (df.N - 1) if "N" in df else 1.0
    df["NSE0"] = df.NSE + df.beta_n**2 * fac
    df["K0"] = (1 - np.sqrt((df.rho - 1)**2 + (df.alpha - 1)**2)).round(3)
    l0, h0 = band(df.NSE0)
    df["in0"] = (df.K0 >= l0) & (df.K0 <= h0)
    imb = df[df.cls >= 4]
    COUNTS = {k: int((df.cls == k).sum()) for k in range(1, 7)}
    out = df[df.cls >= 3]
    SIDE = dict(above=int((out.side == "above").sum()), below=int((out.side == "below").sum()),
                below4=int(((out.side == "below") & (out.cls == 4)).sum()),
                below3=int(((out.side == "below") & (out.cls == 3)).sum()),
                below_other=int(((out.side == "below") & out.cls.isin([5, 6])).sum()),
                above56=int(((out.side == "above") & out.cls.isin([5, 6])).sum()),
                n56=int(out.cls.isin([5, 6]).sum()))
    DK_MED = float((imb.K0 - imb.KGEs).median())
    BACK = int(imb.in0.sum())
    print(f"inside UB: {COUNTS[1] + COUNTS[2]} ({COUNTS[1]} with |1-k_s| <= 0.05, {COUNTS[2]} with |1-k_s| > 0.05)")
    print(f"outside UB: {sum(COUNTS[k] for k in (3, 4, 5, 6))} -> 1: {COUNTS[3]}, 2a: {COUNTS[4]}, 2b: {COUNTS[5]}, 2c: {COUNTS[6]}")
    print("side:", SIDE,
          f"\nimbalance cases: median dKGE*(bias removed) = {DK_MED:.3f}; back in band = {BACK}")
    cases = df

# ================================================================ styling
FILL = {1: "#ffffff", 2: "#ffffff", 3: "#f7b98f", 4: "#8fd3d0", 5: "#9cbcf0", 6: "#c7b6ee"}
INFEAS, INK, MUTED = "#f2f2f2", "#1a1a1a", "#555555"
cmap = ListedColormap([INFEAS] + [FILL[i] for i in range(1, 7)])
C = COUNTS
N_IN, N_OUT = C[1] + C[2], C[3] + C[4] + C[5] + C[6]

OUT = sys.argv[-1] if sys.argv[-1].endswith((".png", ".pdf", ".svg", ".eps")) else "/mnt/user-data/outputs/FigA5_diagnostic_guide.png"
VECTOR = not OUT.endswith(".png")
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(7.48, 7.45), dpi=300)
gs = fig.add_gridspec(2, 4, height_ratios=[1.24, 1.0], hspace=0.17, wspace=0.10,
                      left=0.07, right=0.985, top=0.99, bottom=0.175)

# ================================================================ (a) decision chart (cases outside UB only)
axf = fig.add_subplot(gs[0, :]); axf.set_xlim(0, 100); axf.set_ylim(0, 65); axf.axis("off")

def box(x, y, w, h, txt="", fc="white", bold=False, fs=7.2, ec=INK, lw=0.8, tc=INK):
    axf.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, lw=lw, fc=fc, ec=ec,
                                 boxstyle="round,pad=0.25,rounding_size=1.2"))
    if txt:
        axf.text(x, y, txt, ha="center", va="center", fontsize=fs, color=tc,
                 fontweight="bold" if bold else "normal", linespacing=1.25)

def arrow(x0, y0, x1, y1, lab=None, lx=0, ly=0, ha="center"):
    axf.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.8,
                 color=INK, shrinkA=0, shrinkB=0, mutation_scale=7))
    if lab:
        axf.text((x0 + x1)/2 + lx, (y0 + y1)/2 + ly, lab, ha=ha, va="center",
                 fontsize=6.8, color=MUTED, style="italic", linespacing=1.15)

W, H, YL = 23.0, 13.2, 6.8
TOP = YL + H / 2
X = [12.0, 37.0, 62.0, 87.0]
def leaf(i, k, title, body):
    x = X[i]
    box(x, YL, W, H, fc=FILL[k])
    axf.text(x, YL + H/2 - 2.0, f"Class {CODE[k]}   {title}", ha="center", va="center",
             fontsize=7.4, fontweight="bold", color=INK)
    axf.text(x, YL + 0.2, body, ha="center", va="center", fontsize=6.5, color=INK, linespacing=1.2)
    axf.text(x, YL - H/2 + 1.6, f"CAMELS: {C[k]}", ha="center", va="center",
             fontsize=6.8, color=INK, fontweight="bold")

D = "1 ??? k$_s$"
axf.text(0, 64.5, "a", fontsize=10, fontweight="bold", va="top")
box(50, 61.0, 30, 4.2, "Compute NSE and KGE*", bold=True)
box(50, 53.0, 44, 4.2, "Is KGE* inside the ??5% band around Eq. 7 (Eq. 11)?")
arrow(50, 58.9, 50, 55.1)
# inside: stop (not analysed in Appendix A)
box(13.5, 53.0, 20, 4.2, f"Eq. 7 applies ({N_IN} cases)", fc="#f4f4f2", ec="#9a9a96", fs=6.8, tc=MUTED)
arrow(28, 53.0, 23.6, 53.0, "yes", 0, 1.3)

# outside
box(40.0, 37.0, 17, 4.0, "|1 ??? k$_s$| ??? 0.05 ?")
box(43.0, 45.0, 34, 4.2, "Compute ??, ??, ??$_n$ and k$_s$ = ??/??", bold=True)
arrow(48.0, 50.9, 45.5, 47.1, f"no: {N_OUT}", 4.0, 0.4)
arrow(41.5, 42.9, 40.5, 39.0)
arrow(36.0, 35.0, X[0] + 1.0, TOP, "yes", -1.8, 0)
box(57.5, 27.5, 15, 4.0, f"{D} < ???0.05 ?")
arrow(44.5, 35.0, 55.0, 29.5, "no", 1.6, 0.8)
arrow(54.5, 25.5, X[1] + 1.0, TOP, "yes", -1.8, 0)
box(74.5, 19.5, 10, 3.8, "?? < 1 ?")
arrow(61.5, 25.5, 72.0, 21.4, "no", 1.5, 0.8)
arrow(72.0, 17.6, X[2] + 1.0, TOP, "yes", -1.8, 0)
arrow(77.0, 17.6, X[3] - 1.0, TOP, "no", 1.8, 0)

# (side notes on band position and role of bias are given in the figure caption)

leaf(0, 3, "Bias-driven", f"?? ??? ??, but |??$_n$| exceeds its\ntolerance, which is narrow near\n{D} ??? ???0.05 (panels b???e)")
leaf(1, 4, "??/?? > 1.05", "correlation > variability; since ?? ??? 1,\n?? < 0.952: variability\nunderestimated; Fig. A4a")
leaf(2, 5, "??/?? < 0.95, ?? < 1", "variability relatively larger than\ncorrelation, but underestimated;\nFig. A4b")
leaf(3, 6, "??/?? < 0.95, ?? ??? 1", "variability relatively larger than\ncorrelation, and overestimated\nwhen ?? > 1; Fig. A4c")

# ================================================================ (b-e) regime maps
ks = np.linspace(0.6, 1.5, 901)
bn = np.linspace(0, 0.30, 451)
KS, BN = np.meshgrid(ks, bn)
nse_panels = [0.50, 0.65, 0.80, 0.90]
# NSE ranges centred on each panel value (regions/lines are computed at the centre value);
# half-widths 0.075, 0.075, 0.05, 0.05 -- the largest symmetric ranges that do not overlap
HALF = [0.075, 0.075, 0.05, 0.05]
edges = [(c - h, c + h) for c, h in zip(nse_panels, HALF)]
if cases is not None:
    _o = cases[cases.cls >= 3]
    _shown = sum(int(((_o.NSE >= a) & (_o.NSE < b)).sum()) for a, b in edges)
    print(f"cases outside UB shown in panels b-e: {_shown} of {len(_o)}")
for j, (nse, (e0, e1)) in enumerate(zip(nse_panels, edges)):
    ax = fig.add_subplot(gs[1, j])
    kge, rho, alpha = theory(nse, KS, BN)
    Cm = np.where(rho > 1, 0, classify(nse, kge, rho, alpha))
    if VECTOR:   # one vector polygon set per class (no raster image in the PDF)
        for k in range(0, 7):
            ind = (Cm == k).astype(float)
            if not ind.any(): continue
            col = INFEAS if k == 0 else FILL[k]
            cs = ax.contourf(1 - KS, BN, ind, levels=[0.5, 1.5], colors=[col], antialiased=True)
            try: cs.set_edgecolor("face"); cs.set_linewidth(0.3)   # hide hairline seams between regions
            except Exception: pass
    else:
        ax.pcolormesh(1 - KS, BN, Cm, cmap=cmap, vmin=-0.5, vmax=6.5, shading="auto", rasterized=True)
    hat = ax.contourf(1 - KS, BN, (rho > 1).astype(float), levels=[0.5, 1.5], colors="none",
                      hatches=["////"], zorder=2)
    try: hat.set_edgecolor("#b5b5b5"); hat.set_linewidth(0)
    except Exception: pass
    inside = ((Cm == 1) | (Cm == 2)).astype(float)
    ax.contour(1 - KS, BN, inside, levels=[0.5], colors=INK, linewidths=0.9)
    # alpha = 1 only matters where it separates 2b from 2c (outside UB, 1 - k_s > 0.05)
    a_imb = np.where((Cm == 5) | (Cm == 6), alpha, np.nan)
    ax.contour(1 - KS, BN, a_imb, levels=[1.0], colors=INK,
               linewidths=0.6, linestyles=[(0, (3, 2))])
    # band edge at 1 - k_s = 0: |beta_n| a balanced simulation can sustain before leaving UB
    bb = np.arange(0, 0.5, 0.0005)
    kk, rr, aa = theory(nse, np.ones_like(bb), bb)
    b_edge = bb[np.argmax(classify(nse, kk, rr, aa) >= 3)]
    if b_edge <= 0.30:
        ax.plot([0], [b_edge], marker="D", ms=4.2, mfc=INK, mec="white", mew=0.7, zorder=10)
        ax.annotate(f"{b_edge:.2f}", (0, b_edge), (-0.035, b_edge + 0.006), ha="right", va="bottom",
                    fontsize=7, fontweight="bold", color=INK, zorder=10,
                    bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.8))
    # 1 - k_s = +/-0.05 class boundaries, drawn only outside UB (not across the white inside area)
    for x0 in (-0.05, 0.05):
        ks0 = 1 - x0
        kk, rr, aa = theory(nse, np.full_like(bn, ks0), bn)
        cc = classify(nse, kk, rr, aa)
        yy = np.where((cc >= 3) & (rr <= 1), bn, np.nan)
        ax.plot(np.full_like(bn, x0), yy, color=INK, lw=0.5, ls=":", zorder=3)
    if cases is not None:
        sub = cases[(cases.NSE >= e0) & (cases.NSE < e1)]
        # markers filled with the case's own class (computed at its own NSE)
        for k in (4, 5, 6):
            g = sub[sub.cls == k]
            ax.scatter(1 - g.ks, g.beta_n.abs(), s=8, marker="o", fc=FILL[k], ec=INK, lw=0.5, zorder=8)
        g = sub[sub.cls == 3]
        ax.scatter(1 - g.ks, g.beta_n.abs(), s=11, marker="^", fc=FILL[3], ec=INK, lw=0.6, zorder=8)
    ax.set_xlim(-0.5, 0.4); ax.set_ylim(0, 0.30)
    ax.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    ax.tick_params(direction="in", length=2.5, width=0.6, labelsize=7.2, top=True, right=True)
    for s in ax.spines.values(): s.set_linewidth(0.6)
    ax.set_xlabel("1 ??? k$_s$   (k$_s$ = ??/??)", fontsize=7.6)
    if j == 0: ax.set_ylabel("|??$_n$|", fontsize=8)
    else: ax.set_yticklabels([])
    ax.text(0.0, 1.015, "bcde"[j], transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
    rng = f"[{e0:.3g} ??? {e1:.3g})"   # regions/lines at nse; symbols = cases with NSE in rng
    ax.text(0.5, 1.02, f"NSE = {nse:.2f}  {rng}", transform=ax.transAxes, fontsize=7.6, va="bottom", ha="center")
    mi = (Cm == 1) | (Cm == 2)          # inside UB: one neutral region, not classified here
    xs, ys = KS[mi], BN[mi]
    if mi.any():
        i = np.argmin((xs - np.median(xs))**2 / 0.01 + (ys - 0.02)**2 / 0.0025)
        pass
    for k in range(3, 7):
        m = Cm == k
        if m.mean() < 0.006: continue
        if k == 4: x, y = 1.30, 0.10
        elif k == 3: x, y = 1.025, 0.255 if nse < 0.6 else 0.22
        else:
            xs, ys = KS[m], BN[m]; x, y = np.median(xs), np.median(ys)
            if k == 1: y = min(y, 0.035)
            if not m[np.argmin(abs(bn - y)), np.argmin(abs(ks - x))]:
                i = np.argmin((xs - x)**2 / 0.01 + (ys - y)**2 / 0.0025); x, y = xs[i], ys[i]
        ax.text(1 - x, y, CODE[k], fontsize=7.4, fontweight="bold", ha="center", va="center",
                color=INK, zorder=9, bbox=dict(fc=FILL[k], ec="none", pad=0.4, alpha=0.85))

lab = {1: "0   Inside UB, |1 ??? k$_s$| ??? 0.05: Eq. 7 applies",
       2: "0???  Inside UB, |1 ??? k$_s$| > 0.05: match by compensation",
       3: "Class 1   Outside UB, |1 ??? k$_s$| ??? 0.05: bias-driven",
       4: "Class 2a  Outside UB, 1 ??? k$_s$ < ???0.05: correlation > variability",
       5: "Class 2b  Outside UB, 1 ??? k$_s$ > 0.05, ?? < 1: variability underestimated",
       6: "Class 2c  Outside UB, 1 ??? k$_s$ > 0.05, ?? ??? 1: variability overestimated"}
handles = [Patch(fc=FILL[k], ec="#888", lw=0.4, label=lab[k]) for k in range(3, 7)]
handles.insert(0, Patch(fc="#ffffff", ec="#888", lw=0.4, label="Inside UB: Eq. 7 applies (not classified here)"))
handles += [Patch(fc=INFEAS, ec="#b5b5b5", hatch="////", lw=0.4, label="Infeasible (?? > 1)"),
            plt.Line2D([], [], color=INK, lw=0.9, label="Edge of ??5% band (Eq. 11)"),
            plt.Line2D([], [], color=INK, lw=0.6, ls=(0, (3, 2)),
                       label="?? = 1 (class 2b | 2c boundary);   dotted: 1 ??? k$_s$ = ??0.05")]
if cases is not None:
    handles += [plt.Line2D([], [], ls="", marker="o", ms=3.5, mfc="white", mec=INK, mew=0.6,
                           label="CAMELS cases in classes 2a???2c, filled with own class (NSE within panel range)"),
                plt.Line2D([], [], ls="", marker="^", ms=3.9, mfc=FILL[3], mec=INK, mew=0.6, label="CAMELS cases in class 1")]
handles.append(plt.Line2D([], [], ls="", marker="D", ms=4.2, mfc=INK, mec="white", mew=0.7,
                          label="Band edge at 1 ??? k$_s$ = 0: |??$_n$| a balanced simulation can sustain"))
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=6.6, frameon=False,
           bbox_to_anchor=(0.52, 0.0), handlelength=2.2, columnspacing=1.6, labelspacing=0.42)

fig.savefig(OUT, dpi=300)
print("saved", OUT)
