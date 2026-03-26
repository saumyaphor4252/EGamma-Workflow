"""
Cut optimization for a single variable (e.g. sigma2vv) in pT bins.

- Loads signal and background ROOT trees.
- For each pT bin: finds the cut that gives a target signal efficiency,
  reports background rejection at that cut, and builds ROC curves.
- Fits cut vs pT center with cut = a + b*pT and prints CMSSW-style parameters.
"""
import argparse
import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cut optimization: target signal efficiency per pT bin, background rejection, ROC curves."
    )
    parser.add_argument("--signal", "-s", default="signal.root", help="Signal ROOT file")
    parser.add_argument("--background", "-b", default="background.root", help="Background ROOT file")
    parser.add_argument("--tree", "-t", default="Events", help="TTree name")
    parser.add_argument("--var", "-v", default="eg_sigma2vv", help="Variable to cut on (e.g. eg_sigma2vv)")
    parser.add_argument("--pt-var", default="eg_et", help="pT branch name (e.g. eg_et)")
    parser.add_argument("--weight-var", default=None, help="Optional event weight branch (default: unweighted)")
    parser.add_argument("--target-eff", type=float, default=0.70, help="Target signal efficiency (e.g. 0.70)")
    parser.add_argument(
        "--var-min",
        type=float,
        default=None,
        help="Ignore objects with variable < var_min (e.g. 1e-9 to exclude zeros). Default: use all.",
    )
    parser.add_argument(
        "--pt-bins",
        default="30,40 40,50 50,70 70,100 100,200",
        help="pT bin edges as space-separated 'min,max' pairs or plain pairs (e.g. '30,40 40,50' or '30 40 40 50').",
    )
    parser.add_argument("--roc-xmin", type=float, default=0.0, help="ROC plot x-axis (signal efficiency) minimum")
    parser.add_argument("--roc-xmax", type=float, default=1.0, help="ROC plot x-axis (signal efficiency) maximum")
    args = parser.parse_args()

    # Robust parsing for pt bins.
    tokens = [t.strip() for t in args.pt_bins.replace(";", " ").split() if t.strip()]
    pt_bins = []
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if "," in tok:
            a, b = tok.split(",", 1)
            pt_bins.append((float(a.strip()), float(b.strip())))
            i += 1
            continue
        if i + 1 >= len(tokens):
            parser.error(f"--pt-bins token '{tok}' is missing its pair. Use 'min,max' or provide an even number of values.")
        pt_bins.append((float(tok), float(tokens[i + 1])))
        i += 2
    if not pt_bins:
        parser.error("--pt-bins must contain at least one bin")
    args.pt_bins = pt_bins

    # Strip accidental whitespace/newlines in branch names
    args.tree = args.tree.strip()
    args.var = args.var.strip()
    args.pt_var = args.pt_var.strip()
    if args.weight_var is not None:
        args.weight_var = args.weight_var.strip()
    return args


_args = parse_args()
file_sig = _args.signal
file_bkg = _args.background
tree_name = _args.tree
var = _args.var
pt_var = _args.pt_var
weight_var = _args.weight_var
target_eff = _args.target_eff
var_min = _args.var_min
pt_bins = _args.pt_bins
roc_xmin = _args.roc_xmin
roc_xmax = _args.roc_xmax

f_sig = ROOT.TFile.Open(file_sig)
f_bkg = ROOT.TFile.Open(file_bkg)
t_sig = f_sig.Get(tree_name)
t_bkg = f_bkg.Get(tree_name)


def find_cut(values, weights, target_eff):
    if len(values) == 0:
        return np.nan
    idx = np.argsort(values)
    values = values[idx]
    weights = weights[idx]
    cum_weights = np.cumsum(weights)
    total = cum_weights[-1]
    if total <= 0:
        return np.nan
    eff = cum_weights / total
    i = np.searchsorted(eff, target_eff)
    return values[min(i, len(values) - 1)]


def summarize_values(label, values, weights, round_digits=6, topk=5):
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    n = v.size
    if n == 0:
        print(f"   {label}: n=0")
        return
    finite = np.isfinite(v) & np.isfinite(w)
    n_fin = int(np.sum(finite))
    n_nonfin = n - n_fin
    v_fin = v[finite]
    n_neg = int(np.sum(v_fin < 0)) if n_fin else 0
    n_zero = int(np.sum(v_fin == 0)) if n_fin else 0
    print(
        f"   {label}: n={n}  finite={n_fin} ({n_fin/n:.3f})  nonfinite={n_nonfin} ({n_nonfin/n:.3f})"
        f"  neg={n_neg} ({(n_neg/max(1,n_fin)):.3f})  zero={n_zero} ({(n_zero/max(1,n_fin)):.3f})"
    )
    if n_fin:
        r = np.round(v_fin, round_digits)
        uniq, cnt = np.unique(r, return_counts=True)
        order = np.argsort(cnt)[::-1][:topk]
        common = ", ".join([f"{uniq[i]:g}({cnt[i]})" for i in order])
        print(f"   {label}: most common (rounded to {round_digits}): {common}")


def filter_finite(values, weights):
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    m = np.isfinite(v) & np.isfinite(w)
    return v[m], w[m]


def filter_var_min(values, weights, vmin):
    if vmin is None:
        return values, weights
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    m = v >= vmin
    return v[m], w[m]


def weighted_eff_below_cut(values, weights, cut):
    if len(values) == 0:
        return 0.0
    vals = np.asarray(values, dtype=float)
    wgts = np.asarray(weights, dtype=float)
    total = np.sum(wgts)
    if total <= 0:
        return 0.0
    return np.sum(wgts[vals <= cut]) / total


def compute_roc(sig_vals, sig_wgts, bkg_vals, bkg_wgts, target_sig_eff, max_points=2000):
    sig_vals, sig_wgts = filter_finite(sig_vals, sig_wgts)
    bkg_vals, bkg_wgts = filter_finite(bkg_vals, bkg_wgts)
    if sig_vals.size == 0 or bkg_vals.size == 0:
        return {"fpr": np.array([]), "tpr": np.array([]), "auc": np.nan, "bkg_eff_at_target": np.nan}
    sig_tot = np.sum(sig_wgts)
    bkg_tot = np.sum(bkg_wgts)
    if sig_tot <= 0 or bkg_tot <= 0:
        return {"fpr": np.array([]), "tpr": np.array([]), "auc": np.nan, "bkg_eff_at_target": np.nan}

    sig_idx = np.argsort(sig_vals)
    sig_v = sig_vals[sig_idx]
    sig_cw = np.cumsum(sig_wgts[sig_idx])
    bkg_idx = np.argsort(bkg_vals)
    bkg_v = bkg_vals[bkg_idx]
    bkg_cw = np.cumsum(bkg_wgts[bkg_idx])
    all_v = np.concatenate([sig_v, bkg_v])
    if all_v.size > max_points:
        cuts = np.quantile(all_v, np.linspace(0.0, 1.0, max_points, dtype=float))
    else:
        cuts = np.sort(np.unique(all_v))

    def eff_from_sorted(v_sorted, cw_sorted, total, cut_values):
        pos = np.searchsorted(v_sorted, cut_values, side="right") - 1
        eff = np.zeros_like(cut_values, dtype=float)
        m = pos >= 0
        eff[m] = cw_sorted[pos[m]] / total
        return eff

    tpr = eff_from_sorted(sig_v, sig_cw, sig_tot, cuts)
    fpr = eff_from_sorted(bkg_v, bkg_cw, bkg_tot, cuts)
    tpr = np.clip(tpr, 0.0, 1.0)
    fpr = np.clip(fpr, 0.0, 1.0)
    tpr = np.concatenate([[0.0], tpr, [1.0]])
    fpr = np.concatenate([[0.0], fpr, [1.0]])
    auc_value = float(np.trapz(tpr, fpr))
    bkg_eff_at_target = float(np.interp(target_sig_eff, tpr, fpr, left=fpr[0], right=fpr[-1]))
    return {"fpr": fpr, "tpr": tpr, "auc": auc_value, "bkg_eff_at_target": bkg_eff_at_target}


def _get_pt_val_weight(tree, pt_var, var, weight_var):
    pt_br = getattr(tree, pt_var)
    var_br = getattr(tree, var)
    try:
        n = len(pt_br)
    except TypeError:
        n = None
    if n is not None and n > 0:
        w_br = getattr(tree, weight_var) if weight_var else None
        try:
            weight_is_vector = w_br is not None and len(w_br) >= n
        except TypeError:
            weight_is_vector = False
        for j in range(n):
            pt = float(pt_br[j])
            val = float(var_br[j])
            if weight_is_vector:
                w = float(w_br[j])
            elif weight_var:
                w = float(w_br) if w_br is not None else 1.0
            else:
                w = 1.0
            yield pt, val, w
    else:
        pt = float(pt_br)
        val = float(var_br)
        w = float(getattr(tree, weight_var)) if weight_var else 1.0
        yield pt, val, w


print(f"Optimizing on: {var}")
pt_centers = []  # List of center pT bin for fitting a+b*pT
cuts = []        # List of the variable cut values for Target signal efficinecy in that pT bin
roc_curves = []  # list of (fpr_arr, tpr_arr, label)
dist_bins = []   # list of (ptmin, ptmax, sig_vals, bkg_vals, cut)

for (ptmin, ptmax) in pt_bins:
    sig_vals = []
    sig_wgts = []
    bkg_vals = []
    bkg_wgts = []

    for i in range(t_sig.GetEntries()):
        t_sig.GetEntry(i)
        for pt, val, w in _get_pt_val_weight(t_sig, pt_var, var, weight_var):
            if pt < ptmin or pt >= ptmax:
                continue
            sig_vals.append(val)
            sig_wgts.append(w)

    for i in range(t_bkg.GetEntries()):
        t_bkg.GetEntry(i)
        for pt, val, w in _get_pt_val_weight(t_bkg, pt_var, var, weight_var):
            if pt < ptmin or pt >= ptmax:
                continue
            bkg_vals.append(val)
            bkg_wgts.append(w)

    if len(sig_vals) < 10:
        continue

    sig_vals = np.array(sig_vals, dtype=float)
    sig_wgts = np.array(sig_wgts, dtype=float)
    bkg_vals = np.array(bkg_vals, dtype=float) if bkg_vals else np.array([], dtype=float)
    bkg_wgts = np.array(bkg_wgts, dtype=float) if bkg_wgts else np.array([], dtype=float)

    print(f"\n=== pT bin [{ptmin}, {ptmax}] diagnostics ===")
    summarize_values("signal", sig_vals, sig_wgts)
    summarize_values("background", bkg_vals, bkg_wgts)

    sig_vals_f, sig_wgts_f = filter_finite(sig_vals, sig_wgts)
    bkg_vals_f, bkg_wgts_f = filter_finite(bkg_vals, bkg_wgts)
    if var_min is not None:
        sig_vals_f, sig_wgts_f = filter_var_min(sig_vals_f, sig_wgts_f, var_min)
        bkg_vals_f, bkg_wgts_f = filter_var_min(bkg_vals_f, bkg_wgts_f, var_min)

    cut = find_cut(sig_vals_f, sig_wgts_f, target_eff)
    if len(bkg_vals_f) > 0 and np.sum(bkg_wgts_f) > 0:
        bkg_eff_at_cut = weighted_eff_below_cut(bkg_vals_f, bkg_wgts_f, cut)
        bkg_rejection = 1.0 - bkg_eff_at_cut
    else:
        bkg_rejection = np.nan

    pt_center = 0.5 * (ptmin + ptmax)
    pt_centers.append(pt_center)
    cuts.append(cut)
    dist_bins.append((ptmin, ptmax, np.array(sig_vals_f), np.array(bkg_vals_f), cut))

    roc = compute_roc(sig_vals_f, sig_wgts_f, bkg_vals_f, bkg_wgts_f, target_sig_eff=target_eff)
    if roc["fpr"].size > 0:
        roc_curves.append((roc["fpr"], roc["tpr"], f"pT [{ptmin},{ptmax}] GeV"))
    bkg_str = f"bkg_rejection = {bkg_rejection:.4f}" if not np.isnan(bkg_rejection) else "bkg_rejection = (no bkg)"
    roc_str = ""
    if roc["fpr"].size > 0 and not np.isnan(roc["auc"]):
        roc_str = f"  AUC = {roc['auc']:.4f}  bkg_eff@sig={target_eff:.2f} = {roc['bkg_eff_at_target']:.4g}"
    print(f"pT [{ptmin}, {ptmax}] → cut = {cut:.4f}  sig_eff = {target_eff:.2f}  {bkg_str}{roc_str}")

# Fit
graph = ROOT.TGraph(len(pt_centers))
print("\n=== pT center vs cut points ===")
out_txt = "pt_cut_points.txt"
with open(out_txt, "w") as f:
    f.write("# idx  pt_center_GeV  cut\n")
    for i in range(len(pt_centers)):
        graph.SetPoint(i, pt_centers[i], cuts[i])
        print(f"i={i:2d}  pt_center={pt_centers[i]:7.2f}  cut={cuts[i]:.6f}")
        f.write(f"{i:3d}  {pt_centers[i]:10.4f}  {cuts[i]:.8f}\n")
print(f"Saved points to {out_txt}")
fit = ROOT.TF1("fit", "[0] + [1]*x", 30, 200)
graph.Fit(fit, "Q")
a = fit.GetParameter(0)
b = fit.GetParameter(1)
print("\n=== FIT RESULT ===")
print(f"a = {a}")
print(f"b = {b}")
print("\n=== CMSSW PARAMETERS ===")
print(f"thrRegularEB = cms.vdouble({a})")
print(f"thrOverEE    = cms.vdouble({b})  (useEt = True)")

# ROC styled plot
if roc_curves:
    try:
        color_light_blue = ROOT.TColor.GetColor("#7ec8e3")
        color_dark_red = ROOT.TColor.GetColor("#a52a2a")
    except Exception:
        color_light_blue = ROOT.kAzure + 1
        color_dark_red = ROOT.kRed + 1
    colors = [color_light_blue, color_dark_red, ROOT.kGreen + 2, ROOT.kMagenta + 2, ROOT.kOrange + 2]
    markers = [21, 22, 23, 24, 25]

    canvas = ROOT.TCanvas("roc", "ROC: Background rejection vs Signal efficiency", 700, 600)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.06)
    canvas.SetBottomMargin(0.12)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridColor(ROOT.kGray + 1)
    ROOT.gStyle.SetTickLength(-0.02, "xyz")
    canvas.SetGridx(True)
    canvas.SetGridy(True)

    leg = ROOT.TLegend(0.18, 0.18, 0.52, 0.38)
    leg.SetBorderSize(1)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    graphs = []
    first = True
    for idx, (fpr, tpr, label) in enumerate(roc_curves):
        bkg_rej = 1.0 - np.asarray(fpr, dtype=float)
        sig_eff = np.asarray(tpr, dtype=float)
        g = ROOT.TGraph(len(sig_eff), sig_eff, bkg_rej)
        c = colors[idx % len(colors)]
        mk = markers[idx % len(markers)]
        g.SetLineColor(c)
        g.SetLineWidth(2)
        g.SetMarkerColor(c)
        g.SetMarkerStyle(mk)
        g.SetMarkerSize(0.8)
        graphs.append(g)
        if first:
            g.SetTitle(";Signal efficiency;Background rejection")
            g.Draw("APL")
            h = g.GetHistogram()
            h.SetMinimum(0.0)
            h.SetMaximum(1.0)
            h.GetXaxis().SetRangeUser(roc_xmin, roc_xmax)
            first = False
        else:
            g.Draw("PL same")
        leg.AddEntry(g, label, "lp")
    leg.Draw()
    pave = ROOT.TPaveText(0.18, 0.42, 0.55, 0.58, "NDC")
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    pave.SetTextAlign(12)
    pave.SetTextFont(42)
    pave.SetTextSize(0.032)
    pave.AddText(f"Variable: {var}")
    if var_min is not None:
        pave.AddText(f"var >= {var_min:g}")
    pave.AddText(f"Signal efficiency target: {target_eff}")
    pave.Draw()
    canvas._keepalive = [leg, pave] + graphs
    canvas.SaveAs("ROC_curves.png")
    print("\n   ROC curves saved to ROC_curves.png")
    canvas.Close()

# Distribution plots
if dist_bins:
    try:
        color_sig = ROOT.TColor.GetColor("#7ec8e3")
        color_bkg = ROOT.TColor.GetColor("#a52a2a")
    except Exception:
        color_sig = ROOT.kAzure + 1
        color_bkg = ROOT.kRed + 1
    n_dist_bins = 60
    for ptmin, ptmax, sig_vals_f, bkg_vals_f, cut in dist_bins:
        all_vals = np.concatenate([sig_vals_f, bkg_vals_f]) if (len(sig_vals_f) and len(bkg_vals_f)) else (sig_vals_f if len(sig_vals_f) else bkg_vals_f)
        if len(all_vals) < 2:
            continue
        vmin = float(np.min(all_vals))
        vmax = float(np.max(all_vals))
        if vmax <= vmin:
            vmax = vmin + 1.0
        c_dist = ROOT.TCanvas("dist", f"{var} pT [{ptmin},{ptmax}]", 700, 600)
        c_dist.SetGridx(True)
        c_dist.SetGridy(True)
        ROOT.gStyle.SetGridStyle(3)
        h_sig = ROOT.TH1F("h_sig", "", n_dist_bins, vmin, vmax)
        h_bkg = ROOT.TH1F("h_bkg", "", n_dist_bins, vmin, vmax)
        for x in sig_vals_f:
            h_sig.Fill(float(x))
        for x in bkg_vals_f:
            h_bkg.Fill(float(x))
        if h_sig.GetSumOfWeights() > 0:
            h_sig.Scale(1.0 / h_sig.GetSumOfWeights())
        if h_bkg.GetSumOfWeights() > 0:
            h_bkg.Scale(1.0 / h_bkg.GetSumOfWeights())
        h_sig.SetStats(0)
        h_bkg.SetStats(0)
        h_sig.SetLineColor(color_sig)
        h_sig.SetLineWidth(2)
        h_sig.SetFillStyle(0)
        h_bkg.SetLineColor(color_bkg)
        h_bkg.SetLineWidth(2)
        h_bkg.SetFillStyle(0)
        h_sig.SetTitle(f";{var};a.u.")
        y_max = max(h_sig.GetMaximum(), h_bkg.GetMaximum())
        h_sig.SetMaximum(1.10 * y_max if y_max > 0 else 1.0)
        h_sig.Draw("HIST")
        h_bkg.Draw("HIST SAME")
        line = ROOT.TLine(cut, 0, cut, max(h_sig.GetMaximum(), h_bkg.GetMaximum()) * 1.05)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.SetLineColor(ROOT.kBlack)
        line.Draw()
        leg_d = ROOT.TLegend(0.6, 0.72, 0.88, 0.88)
        leg_d.SetBorderSize(1)
        leg_d.SetFillStyle(0)
        leg_d.AddEntry(h_sig, "Signal", "l")
        leg_d.AddEntry(h_bkg, "Background", "l")
        leg_d.AddEntry(line, f"Cut (sig eff={target_eff})", "l")
        leg_d.Draw()
        pave_d = ROOT.TPaveText(0.18, 0.82, 0.5, 0.92, "NDC")
        pave_d.SetBorderSize(0)
        pave_d.SetFillStyle(0)
        pave_d.SetTextFont(42)
        pave_d.SetTextSize(0.032)
        pave_d.AddText(f"pT [{ptmin}, {ptmax}] GeV")
        pave_d.Draw()
        c_dist._keepalive = [h_sig, h_bkg, line, leg_d, pave_d]
        out_name = f"distributions_pt{ptmin:.0f}_{ptmax:.0f}.png"
        c_dist.SaveAs(out_name)
        print(f"   Distributions saved to {out_name}")
        c_dist.Close()