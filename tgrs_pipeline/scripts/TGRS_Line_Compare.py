#!/usr/bin/env python3
"""Per-line comparison of two TGRS simulations (A = reference, B = test) and the TGRS data,
using the line list of Weidenspointner et al. 2005 (ApJS 156, 69, Tables 5 and 6).

For every line: net count rate (counts/s) in a window around the line energy, minus a linear
continuum fitted to side bands, for Sim A, Sim B and the data. Errors are propagated from the
Poisson errors stored in the ROOT histograms (TGRS.cxx uses Sumw2).

Usage:
  python3 TGRS_Line_Compare.py A/TGRS_output.root B/TGRS_output.root \
      --lines tgrs_lines.csv --out analysis/lines --ref-label "10.2" --test-label "dc6" \
      [--sigma 3] [--ntop 12] [--fwhm0 0.9 --fwhm1 0.0013] [--component 2]

Outputs in --out: line_comparison.csv, line_comparison.json, worst_lines.png/.pdf, line_summary.png/.pdf
--component N additionally reports the line in Component_N of both files (0 photons, 1 prompt p, 2 activation).
"""
import argparse, csv, json, os, sys
import numpy as np
# load matplotlib/scipy before uproot: on Savio this pulls in the newer libstdc++ that pyarrow needs
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import scipy.stats  # noqa
except Exception:
    pass


def read_hist(path, name):
    import uproot
    with uproot.open(path) as f:
        keys = [k.split(";")[0] for k in f.keys()]
        if name not in keys:
            raise KeyError(f"'{name}' not in {path}. Available: {keys}")
        h = f[name]
        return (np.asarray(h.values(), float), np.asarray(h.errors(), float),
                np.asarray(h.axis().edges(), float))


def integrate(vals, errs, edges, lo, hi):
    """Integral of a per-keV density over [lo,hi] (partial bins by overlap) -> (counts/s, err, width)."""
    ov = np.clip(np.minimum(edges[1:], hi) - np.maximum(edges[:-1], lo), 0, None)
    return float(np.sum(vals * ov)), float(np.sqrt(np.sum((errs * ov) ** 2))), float(ov.sum())


def net_rate(vals, errs, edges, e0, fwhm):
    """Net line rate: window +-1.5 FWHM, continuum = linear from two side bands (2-4 FWHM away).
    The window is widened to >=1.5 bins each side so coarse bins at high energy still work."""
    centers = 0.5 * (edges[1:] + edges[:-1])
    i0 = np.searchsorted(edges, e0) - 1
    if i0 < 0 or i0 >= len(vals):
        return None
    bw = edges[i0 + 1] - edges[i0]
    w = max(1.5 * fwhm, 1.5 * bw)
    sb0, sb1 = w + 0.5 * fwhm, w + 2.5 * fwhm + 2 * bw
    if e0 - sb1 < edges[0] or e0 + sb1 > edges[-1]:
        return None
    S, eS, wS = integrate(vals, errs, edges, e0 - w, e0 + w)
    Ls, eLs, wLs = integrate(vals, errs, edges, e0 - sb1, e0 - sb0)
    Rs, eRs, wRs = integrate(vals, errs, edges, e0 + sb0, e0 + sb1)
    if wLs <= 0 or wRs <= 0 or wS <= 0:
        return None
    cont_density = 0.5 * (Ls / wLs + Rs / wRs)
    e_cd = 0.5 * np.hypot(eLs / wLs, eRs / wRs)
    net = S - cont_density * wS
    err = float(np.hypot(eS, e_cd * wS))
    return net, err, cont_density * wS


def ratio(a, ea, b, eb):
    if b is None or b == 0 or a is None:
        return np.nan, np.nan
    r = a / b
    e = abs(r) * np.hypot(ea / a if a else 0.0, eb / b)
    return r, e


def split_id(idstr):
    """'73mGe(67-13) | 53.4 | a, s, ...' -> (transition, lab energy, comments)."""
    parts = [p.strip() for p in idstr.split(" | ")]
    trans = parts[0] if parts else ""
    lab = parts[1] if len(parts) > 1 else ""
    com = ", ".join(parts[2:]) if len(parts) > 2 else ""
    if lab and not any(ch.isdigit() for ch in lab):
        com = ", ".join(parts[1:]); lab = ""
    return trans, lab, com


def fmt(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    ax = abs(x)
    return f"{x:.2e}" if (ax < 1e-3 and ax > 0) else f"{x:.3g}"


def write_table_pdf(path, rows, ref_label, test_label, sigma, with_comp, title, scope_note=""):
    """Table in the style of Table 5 of Weidenspointner et al. 2005, with simulation columns added."""
    from reportlab.lib.pagesizes import letter, landscape
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
    from xml.sax.saxutils import escape
    ss = getSampleStyleSheet()
    cell = ParagraphStyle("c", parent=ss["BodyText"], fontSize=6.5, leading=7.6)
    hdr = ParagraphStyle("h", parent=cell, fontName="Helvetica-Bold")
    doc = SimpleDocTemplate(path, pagesize=landscape(letter), leftMargin=30, rightMargin=30,
                            topMargin=30, bottomMargin=30, title=title)
    head = ["TGRS line (keV)", "Paper rate (1/s)", "Transition (levels, keV)", "Lab E (keV)",
            "Data (1/s)", f"A {ref_label} (1/s)", f"B {test_label} (1/s)", "A/B", "(A-B) /sigma"]
    if with_comp:
        head += ["A act. only", "B act. only"]
    head += ["Paper comments / flag"]
    data = [[Paragraph(h, hdr) for h in head]]
    style = [("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
             ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey), ("VALIGN", (0, 0), (-1, -1), "TOP")]
    for i, r in enumerate(rows, start=1):
        trans, lab, com = split_id(r["id"])
        sg = r["sig_AB"]
        flag = ""
        if abs(sg) > sigma:
            flag = f"<b>{'B higher' if sg < 0 else 'B lower'}</b> ({abs(sg):.0f} sigma). "
            style.append(("BACKGROUND", (0, i), (-1, i), colors.Color(1, .93, .8)))
        ratio = f"{r['ratio_AB']:.2f}" if np.isfinite(r["ratio_AB"]) and r["A"] > 0 and r["B"] > 0 else "-"
        vals = [f"{r['energy']:.1f}", fmt(r["paper_rate"]), escape(trans), escape(lab),
                fmt(r["data"]), fmt(r["A"]), fmt(r["B"]), ratio, f"{sg:.1f}"]
        if with_comp:
            vals += [fmt(r.get("compA")), fmt(r.get("compB"))]
        vals += [flag + escape(com)]
        data.append([Paragraph(v, cell) for v in vals])
    cw = [36, 34, 112, 32, 42, 48, 48, 28, 38] + ([46, 46] if with_comp else [])
    cw.append(max(120, 750 - sum(cw)))
    t = Table(data, colWidths=cw, repeatRows=1)
    t.setStyle(TableStyle(style))
    body = ParagraphStyle("b", parent=ss["BodyText"], fontSize=8.5, leading=10.5, spaceAfter=3)
    story = [Paragraph(f"<b>{escape(title)}</b>", ss["Heading2"]),
             Paragraph("<b>What this table shows.</b> Each row is one gamma-ray line from the TGRS line list of "
                       "Weidenspointner et al. 2005 (ApJS 156, 69, Tables 5 and 6), in the same order and with the same "
                       "identification as the paper (only the first identification of a blend is given). For every line we "
                       f"measure the same quantity three times: in the TGRS data, in simulation A ({escape(ref_label)}, "
                       f"the reference build) and in simulation B ({escape(test_label)}, the test build). The point is to "
                       "see which lines differ between the two builds, so that the cause can be debugged, and whether "
                       "either build reproduces the measured line.", body),
             Paragraph("<b>Which lines are included.</b> " + escape(scope_note), body),
             Paragraph("<b>How the numbers are made.</b> The line strength is the net count rate in counts per second: "
                       "the rate inside a window of +-1.5 FWHM around the line energy, minus a straight-line continuum "
                       "estimated from two side bands next to the window. The FWHM is an approximate detector resolution. "
                       "Errors are Poisson errors from the ROOT histograms. <i>A/B</i> is the ratio of the two simulated "
                       "line strengths (1 means the builds agree); <i>(A-B)/sigma</i> is their difference in units of its "
                       "statistical error. "
                       + ("The 'act. only' columns give the same net rate for the activation (delayed decay) component "
                          "alone, which shows whether a difference comes from the isotope decay or from the photon and "
                          "prompt components. " if with_comp else ""), body),
             Paragraph("<b>How to read it.</b> Shaded rows differ between A and B by more than "
                       f"{sigma:g} sigma; 'B higher' means the test build produces the stronger line. Compare both "
                       "simulations with the Data column to see which one is closer to the measurement. "
                       "<b>Caution:</b> a negative rate means the line sits next to a stronger line or on a steep continuum, "
                       "so the sideband estimate failed; such rows (and the A/B shown as '-') are unreliable. Lines above "
                       "about 800 keV fall in coarse simulation bins and are only roughly resolved. Weak lines depend on the "
                       "assumed resolution.", body),
             Spacer(1, 6), t]
    doc.build(story)


def append_pdfs(base, extra, out):
    try:
        from pypdf import PdfWriter, PdfReader
    except ImportError:
        try:
            from PyPDF2 import PdfWriter, PdfReader
        except ImportError:
            print("pypdf/PyPDF2 not installed (pip install --user pypdf): wrote the line pages separately only")
            return False
    w = PdfWriter()
    for f in [base] + extra:
        for pg in PdfReader(f).pages:
            w.add_page(pg)
    with open(out, "wb") as fh:
        w.write(fh)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref"); ap.add_argument("test")
    ap.add_argument("--lines", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "tgrs_lines.csv"))
    ap.add_argument("--out", default="line_comparison")
    ap.add_argument("--ref-label", default="A (reference)"); ap.add_argument("--test-label", default="B (test)")
    ap.add_argument("--sigma", type=float, default=3.0, help="flag |A-B| above this many sigma")
    ap.add_argument("--min-rate", type=float, default=0.0, help="ignore lines weaker than this in the paper (counts/s)")
    ap.add_argument("--ntop", type=int, default=12)
    ap.add_argument("--fwhm0", type=float, default=0.9, help="FWHM(E)=fwhm0+fwhm1*E in keV (approximate TGRS resolution)")
    ap.add_argument("--fwhm1", type=float, default=0.0013)
    ap.add_argument("--component", type=int, default=None)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--append-to", default=None,
                    help="existing validation report PDF: the line table and plots are appended as new pages "
                         "(a copy <name>_with_lines.pdf is written; the original is not changed)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    H = {}
    for tag, p in (("A", a.ref), ("B", a.test)):
        H[tag] = read_hist(p, "Sim")
        if a.component is not None:
            H[tag + "c"] = read_hist(p, f"Component_{a.component}")
    H["D"] = read_hist(a.ref, "TGRS")

    rows = []
    n_in_list = n_weak = n_range = 0
    with open(a.lines) as f:
        for r in csv.DictReader(f):
            n_in_list += 1
            e0, rate = float(r["energy_keV"]), float(r["rate_per_s"])
            weak = rate < a.min_rate
            fw = a.fwhm0 + a.fwhm1 * e0
            res = {k: net_rate(*H[k], e0, fw) for k in H}
            if any(v is None for v in res.values()):
                n_range += 1
                continue
            (nA, eA, _), (nB, eB, _), (nD, eD, _) = res["A"], res["B"], res["D"]
            d = nA - nB; sd = np.hypot(eA, eB)
            rAB, erAB = ratio(nA, eA, nB, eB)
            rAD, _ = ratio(nA, eA, nD, eD); rBD, _ = ratio(nB, eB, nD, eD)
            row = dict(table=r["table"], energy=e0, paper_rate=rate, id=r["id"], fwhm=fw,
                       A=nA, eA=eA, B=nB, eB=eB, data=nD, edata=eD,
                       diff_AB=d, sig_AB=d / sd if sd > 0 else 0.0, ratio_AB=rAB, eratio_AB=erAB,
                       ratio_A_data=rAD, ratio_B_data=rBD,
                       sig_A_data=(nA - nD) / np.hypot(eA, eD), sig_B_data=(nB - nD) / np.hypot(eB, eD))
            if a.component is not None:
                cA, cB = res["Ac"], res["Bc"]
                row.update(compA=cA[0], compB=cB[0],
                           comp_sig=(cA[0] - cB[0]) / max(np.hypot(cA[1], cB[1]), 1e-30))
            row["in_report"] = not weak
            rows.append(row)
    if not rows:
        sys.exit("no lines inside the histogram range")
    rows.sort(key=lambda r: r["energy"])
    rows_all = rows                              # every analysable line -> CSV
    rows = [r for r in rows_all if r["in_report"]]   # lines at or above --min-rate -> table, plots, summary
    n_weak = len(rows_all) - len(rows)
    if not rows:
        sys.exit("no lines at or above --min-rate")

    keys = list(rows_all[0].keys())
    with open(os.path.join(a.out, "line_comparison.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(rows_all)

    flagged = [r for r in rows if abs(r["sig_AB"]) > a.sigma]
    flagged.sort(key=lambda r: -abs(r["sig_AB"]))
    sA = np.array([r["sig_A_data"] for r in rows]); sB = np.array([r["sig_B_data"] for r in rows])
    summ = dict(n_lines=len(rows), n_flagged_AB=len(flagged), sigma=a.sigma,
                median_ratio_AB=float(np.nanmedian([r["ratio_AB"] for r in rows])),
                median_abs_sig_AB=float(np.median([abs(r["sig_AB"]) for r in rows])),
                median_ratio_A_data=float(np.nanmedian([r["ratio_A_data"] for r in rows])),
                median_ratio_B_data=float(np.nanmedian([r["ratio_B_data"] for r in rows])),
                n_A_data_off=int(np.sum(np.abs(sA) > a.sigma)), n_B_data_off=int(np.sum(np.abs(sB) > a.sigma)),
                flagged=[dict(energy=r["energy"], id=r["id"], A=r["A"], B=r["B"], data=r["data"],
                              ratio_AB=r["ratio_AB"], sig_AB=r["sig_AB"]) for r in flagged])
    json.dump(summ, open(os.path.join(a.out, "line_comparison.json"), "w"), indent=1, default=float)

    print(f"\n{len(rows)} lines compared; {len(flagged)} differ between A and B by > {a.sigma} sigma")
    print(f"median A/B = {summ['median_ratio_AB']:.3f}; median sim/data: A {summ['median_ratio_A_data']:.2f}, B {summ['median_ratio_B_data']:.2f}")
    print(f"{'E(keV)':>8} {'A':>10} {'B':>10} {'data':>10} {'A/B':>6} {'sigAB':>6}  ID")
    for r in flagged[:40]:
        print(f"{r['energy']:8.1f} {r['A']:10.4g} {r['B']:10.4g} {r['data']:10.4g} {r['ratio_AB']:6.2f} {r['sig_AB']:6.1f}  {r['id'][:50]}")

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # summary: A/B per line and each sim vs data
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    E = np.array([r["energy"] for r in rows])
    ax[0].errorbar(E, [r["ratio_AB"] for r in rows], [r["eratio_AB"] for r in rows], fmt="o", ms=3, color="k", lw=.6)
    ax[0].axhline(1, color="gray"); ax[0].set_ylim(0, 3); 
    ax[0].set_ylabel("A / B  (1 = builds agree)")
    ax[0].set_title("Each point is one line: ratio of the net line rate in simulation A to B (red labels: largest differences, keV)", fontsize=9)
    for r in flagged[:8]:
        ax[0].annotate(f"{r['energy']:.0f}", (r["energy"], min(r["ratio_AB"], 2.9)), fontsize=7, color="r")
    ax[1].plot(E, [r["ratio_A_data"] for r in rows], "o", ms=3, label=a.ref_label)
    ax[1].plot(E, [r["ratio_B_data"] for r in rows], "s", ms=3, label=a.test_label)
    ax[1].axhline(1, color="gray"); ax[1].set_ylim(-1, 4); ax[1].set_xscale("log"); ax[1].set_xlabel("Line energy (keV)")
    ax[1].set_ylabel("simulation / data  (1 = reproduces the measurement)")
    ax[1].set_title("Each point is one line: net line rate in a simulation divided by the net line rate in the TGRS data", fontsize=9); ax[1].legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, f"line_summary.{ext}"), dpi=a.dpi)
    plt.close(fig)

    # zoom plots of the worst lines
    top = flagged[:a.ntop] if flagged else sorted(rows, key=lambda r: -abs(r["sig_AB"]))[:a.ntop]
    if top:
        n = len(top); nc = 3; nr = (n + nc - 1) // nc
        fig, axs = plt.subplots(nr, nc, figsize=(13, 3.2 * nr), squeeze=False)
        for axx, r in zip(axs.ravel(), top):
            e0, fw = r["energy"], r["fwhm"]
            span = max(6 * fw, 12)
            for k, lab, st in (("D", "data", dict(color="k", drawstyle="steps-mid", lw=.8)),
                               ("A", a.ref_label, dict(color="C0", drawstyle="steps-mid", lw=1)),
                               ("B", a.test_label, dict(color="C3", drawstyle="steps-mid", lw=1))):
                v, _, ed = H[k]; c = .5 * (ed[1:] + ed[:-1]); m = (c > e0 - span) & (c < e0 + span)
                axx.plot(c[m], v[m], label=lab, **st)
            axx.axvline(e0, color="gray", ls=":")
            trans, _lab, _c = split_id(r["id"])
            rt = (f"A/B = {r['ratio_AB']:.2f}" if np.isfinite(r["ratio_AB"]) and r["A"] > 0 and r["B"] > 0
                  else "A/B not defined (negative net rate)")
            axx.set_title(f"{e0:.1f} keV: {trans[:34]}\n{rt}, (A-B)/σ = {r['sig_AB']:.1f}", fontsize=8)
            axx.set_xlabel("Energy (keV)"); axx.set_ylabel("Rate (counts/s/keV)")
        for axx in axs.ravel()[n:]:
            axx.axis("off")
        axs[0, 0].legend(fontsize=7)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(a.out, f"worst_lines.{ext}"), dpi=a.dpi)
        plt.close(fig)
    # Table 5 style PDF (+ plots) and optional append to the validation report
    tpdf = os.path.join(a.out, "line_table.pdf")
    scope = (f"The line list of the paper has {n_in_list} entries; this table contains {len(rows)} of them. "
             f"{n_weak} weaker than {a.min_rate:g} counts/s in the paper were left out of this table (they are analysed and listed in line_comparison.csv, column in_report = False)" if a.min_rate > 0 else
             f"The line list of the paper has {n_in_list} entries; this table contains {len(rows)} of them")
    scope += (f", and {n_range} could not be analysed because their window lies outside the simulated energy range. "
              f"{len(flagged)} of the {len(rows)} lines differ between A and B by more than {a.sigma:g} sigma. "
              "This is not a list of all real differences: some flagged rows are artifacts of a neighbouring line "
              "(see the caution below), and weak lines that were left out may differ too. The complete results are in "
              "line_comparison.csv.")
    write_table_pdf(tpdf, rows, a.ref_label, a.test_label, a.sigma, a.component is not None,
                    "Line-by-line comparison of the two simulations and the TGRS data (Table 5 style)", scope)
    extra = [tpdf]
    try:
        from reportlab.platypus import SimpleDocTemplate, Image, Paragraph
        from reportlab.lib.pagesizes import landscape, letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from xml.sax.saxutils import escape
        pp = os.path.join(a.out, "line_plots.pdf")
        ss = getSampleStyleSheet()
        d = SimpleDocTemplate(pp, pagesize=landscape(letter), leftMargin=30, rightMargin=30, topMargin=30, bottomMargin=30)
        cap = ParagraphStyle("cap", parent=ss["BodyText"], fontSize=8.5, leading=10.5)
        st = [Paragraph("<b>Line-by-line comparison: summary of all lines</b>", ss["Heading2"]),
              Paragraph("Top: for every line in the table, the ratio of the net line rate in simulation A "
                        f"({escape(a.ref_label)}) to B ({escape(a.test_label)}); points away from 1 are lines where the "
                        "builds differ, vertical bars are statistical errors, and the red numbers mark the most "
                        "different lines (energy in keV). Bottom: the net line rate of each simulation divided by the "
                        "measured rate; 1 means the simulation reproduces the TGRS line, below 1 means the line is too weak. "
                        "Negative or very noisy rates (weak lines, blends) are not reliable.", cap)]
        for png, h in (("line_summary.png", 330), ("worst_lines.png", None)):
            fp = os.path.join(a.out, png)
            if os.path.exists(fp):
                from PIL import Image as PI
                w0, h0 = PI.open(fp).size
                maxh = 400 if png == "line_summary.png" else 440
                sc = min(720 / w0, maxh / h0)  # page is 752 x 552 pt usable
                st.append(Image(fp, w0 * sc, h0 * sc))
                if png == "line_summary.png":
                    from reportlab.platypus import PageBreak
                    st.append(PageBreak())
                    st.append(Paragraph("<b>Line-by-line comparison: the lines with the largest A-B differences</b>", ss["Heading2"]))
                    st.append(Paragraph("Each panel is a zoom on one line. Black: TGRS data; blue: simulation A "
                                        f"({escape(a.ref_label)}); red: simulation B ({escape(a.test_label)}); dotted line: "
                                        "line energy from the paper. The title gives the energy, the identification and how "
                                        "far apart A and B are. Lines are ordered by the significance of the A-B difference "
                                        "(largest first). Panels of a line next to a stronger line show the neighbour, not the "
                                        "line itself.", cap))
        d.build(st); extra.append(pp)
    except Exception as e:
        print(f"[note] plot pages for the report not built: {e}")
    if a.append_to:
        out = a.append_to[:-4] + "_with_lines.pdf" if a.append_to.endswith(".pdf") else a.append_to + "_with_lines.pdf"
        if append_pdfs(a.append_to, extra, out):
            print(f"Report with line pages: {out}")
    print(f"\nWrote {a.out}/line_comparison.csv, .json, line_table.pdf, line_summary.*, worst_lines.*")


if __name__ == "__main__":
    main()
