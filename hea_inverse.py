
import math
import random
import streamlit as st
import pandas as pd
import numpy as np

# ============================================================
#  Element Database
# ============================================================

ELEMENTS_DATA = {
    "Nb": {"r": 146, "VEC": 5, "ea": 1.51, "chi": 1.6,  "Tm": 2477+273.15, "B": 170.3, "alpha": 7.1e-6,
           "E_BCC": -10.0466, "E_FCC": -9.7232,  "E_HCP": -9.7551},
    "Mo": {"r": 139, "VEC": 6, "ea": 1.87, "chi": 2.16, "Tm": 2623+273.15, "B": 115.8, "alpha": 4.9e-6,
           "E_BCC": -10.7799, "E_FCC": -10.3784, "E_HCP": -10.3666},
    "Ta": {"r": 146, "VEC": 5, "ea": 1.50, "chi": 1.5,  "Tm": 3017+273.15, "B": 200.1, "alpha": 6.5e-6,
           "E_BCC": -11.7358, "E_FCC": -11.1896, "E_HCP": -11.4579},
    "W":  {"r": 139, "VEC": 6, "ea": 1.86, "chi": 2.36, "Tm": 3422+273.15, "B": 323.3, "alpha": 4.5e-6,
           "E_BCC": -12.7781, "E_FCC": -12.3115, "E_HCP": -12.2928},
    "V":  {"r": 135, "VEC": 5, "ea": 1.50, "chi": 1.63, "Tm": 1910+273.15, "B": 162.0, "alpha": 8.7e-6,
           "E_BCC": -8.9632,  "E_FCC": -8.7150,  "E_HCP": -8.7095},
    "Cr": {"r": 128, "VEC": 6, "ea": 1.87, "chi": 1.66, "Tm": 1907+273.15, "B": 116.7, "alpha": 6.1e-6,
           "E_BCC": -9.4655,  "E_FCC": -9.0845,  "E_HCP": -9.0751},
    "Ti": {"r": 147, "VEC": 4, "ea": 1.33, "chi": 1.54, "Tm": 1668+273.15, "B": 105.2, "alpha": 8.6e-6,
           "E_BCC": -2.2301,  "E_FCC": -2.2155,  "E_HCP": -2.2343},
    "Zr": {"r": 160, "VEC": 4, "ea": 1.33, "chi": 1.33, "Tm": 1855+273.15, "B":  83.35, "alpha": 5.5e-6,
           "E_BCC": -8.3598,  "E_FCC": -8.3972,  "E_HCP": -8.4354},
    "Hf": {"r": 159, "VEC": 4, "ea": 1.32, "chi": 1.30, "Tm": 2233+273.15, "B": 108.9, "alpha": 5.9e-6,
           "E_BCC": -9.6562,  "E_FCC": -9.7613,  "E_HCP": -9.8320},
}

DELTA_H_IJ = {
    ("Nb","Mo"): -6,  ("Mo","Nb"): -6,
    ("Nb","Ta"):  0,  ("Ta","Nb"):  0,
    ("Nb","W"):  -8,  ("W","Nb"):  -8,
    ("Nb","V"):  -1,  ("V","Nb"):  -1,
    ("Nb","Cr"): -7,  ("Cr","Nb"): -7,
    ("Nb","Ti"):  2,  ("Ti","Nb"):  2,
    ("Nb","Zr"):  4,  ("Zr","Nb"):  4,
    ("Nb","Hf"):  4,  ("Hf","Nb"):  4,
    ("Mo","Ta"): -5,  ("Ta","Mo"): -5,
    ("Mo","W"):   0,  ("W","Mo"):   0,
    ("Mo","V"):   0,  ("V","Mo"):   0,
    ("Mo","Cr"):  0,  ("Cr","Mo"):  0,
    ("Mo","Ti"): -4,  ("Ti","Mo"): -4,
    ("Mo","Zr"): -6,  ("Zr","Mo"): -6,
    ("Mo","Hf"): -4,  ("Hf","Mo"): -4,
    ("Ta","W"):  -7,  ("W","Ta"):  -7,
    ("Ta","V"):  -1,  ("V","Ta"):  -1,
    ("Ta","Cr"): -7,  ("Cr","Ta"): -7,
    ("Ta","Ti"):  1,  ("Ti","Ta"):  1,
    ("Ta","Zr"):  3,  ("Zr","Ta"):  3,
    ("Ta","Hf"):  3,  ("Hf","Ta"):  3,
    ("W","V"):   -1,  ("V","W"):   -1,
    ("W","Cr"):   1,  ("Cr","W"):   1,
    ("W","Ti"):  -6,  ("Ti","W"):  -6,
    ("W","Zr"):  -9,  ("Zr","W"):  -9,
    ("W","Hf"):  -6,  ("Hf","W"):  -6,
    ("V","Cr"):  -2,  ("Cr","V"):  -2,
    ("V","Ti"):  -2,  ("Ti","V"):  -2,
    ("V","Zr"):  -4,  ("Zr","V"):  -4,
    ("V","Hf"):  -2,  ("Hf","V"):  -2,
    ("Cr","Ti"): -7,  ("Ti","Cr"): -7,
    ("Cr","Zr"): -12, ("Zr","Cr"): -12,
    ("Cr","Hf"): -9,  ("Hf","Cr"): -9,
    ("Ti","Zr"):  0,  ("Zr","Ti"):  0,
    ("Ti","Hf"):  0,  ("Hf","Ti"):  0,
    ("Zr","Hf"):  0,  ("Hf","Zr"):  0,
}

R = 8.314

# ============================================================
#  Physics Functions
# ============================================================

def calc_r_bar(els, fracs):
    return sum(fracs[i] * ELEMENTS_DATA[els[i]]["r"] for i in range(len(els)))

def calc_delta_r(els, fracs):
    r_bar = calc_r_bar(els, fracs)
    return math.sqrt(sum(fracs[i] * (1 - ELEMENTS_DATA[els[i]]["r"] / r_bar)**2
                         for i in range(len(els)))) * 100

def calc_delta_S_mix(fracs):
    return -R * sum(c * math.log(c) for c in fracs if c > 0)

def calc_delta_H_mix(els, fracs):
    total = 0
    n = len(els)
    for i in range(n):
        for j in range(n):
            if i != j:
                total += 4 * DELTA_H_IJ.get((els[i], els[j]), 0) * fracs[i] * fracs[j]
    return total / 2

def calc_VEC(els, fracs):
    return sum(fracs[i] * ELEMENTS_DATA[els[i]]["VEC"] for i in range(len(els)))

def calc_E_DFT(els, fracs, structure="BCC"):
    return sum(fracs[i] * ELEMENTS_DATA[els[i]][f"E_{structure}"] for i in range(len(els)))

def calc_delta_H_elastic(els, fracs, T):
    data = ELEMENTS_DATA
    def v0(el): return (4/3) * math.pi * (data[el]["r"] * 1e-12)**3
    def v_T(el): return v0(el) * (1 + data[el]["alpha"] * (T - 298))
    v_vals = [v_T(el) for el in els]
    B_vals = [data[el]["B"] * 1e9 for el in els]
    v_bar = sum(fracs[i]*B_vals[i]*v_vals[i] for i in range(len(els))) / \
            sum(fracs[i]*B_vals[i] for i in range(len(els)))
    total = sum(fracs[i]*B_vals[i]*(v_vals[i]-v_bar)**2/(2*v_vals[i]) for i in range(len(els)))
    return total / 1.602e-19 / 6.022e23 * 1000

def calc_gamma(els, fracs):
    r_bar = calc_r_bar(els, fracs)
    r_vals = [ELEMENTS_DATA[el]["r"] for el in els]
    r_min, r_max = min(r_vals), max(r_vals)
    num = 1 - (r_min / (r_min + r_bar))**2
    den = 1 - (r_max / (r_max + r_bar))**2
    return float('inf') if den == 0 else num / den

def calc_Tm_mean(els, fracs):
    return sum(fracs[i] * ELEMENTS_DATA[els[i]]["Tm"] for i in range(len(els)))

def calc_Omega(Tm, dS, dH):
    return float('inf') if abs(dH) < 1e-10 else Tm * dS / (abs(dH) * 1000)

def calc_delta_chi(els, fracs):
    chi_vals = [ELEMENTS_DATA[el]["chi"] for el in els]
    chi_bar = sum(fracs[i] * chi_vals[i] for i in range(len(els)))
    return math.sqrt(sum(fracs[i] * (1 - chi_vals[i]/chi_bar)**2
                         for i in range(len(els)))) * 100

def calc_ea(els, fracs):
    return sum(fracs[i] * ELEMENTS_DATA[els[i]]["ea"] for i in range(len(els)))

# ============================================================
#  Criteria Checkers — return dict of {key: bool}
# ============================================================

def get_all_criteria(els, fracs, T):
    dH        = calc_delta_H_mix(els, fracs)
    dS        = calc_delta_S_mix(fracs)
    delta_r   = calc_delta_r(els, fracs)
    VEC       = calc_VEC(els, fracs)
    gamma     = calc_gamma(els, fracs)
    delta_chi = calc_delta_chi(els, fracs)
    ea        = calc_ea(els, fracs)
    Tm        = calc_Tm_mean(els, fracs)
    Omega     = calc_Omega(Tm, dS, dH)
    ratio     = T / Tm

    results = {}

    # --- SS criteria ---
    results["MC1 Zhang (2008)"]    = (0.5 < delta_r < 6.5) and (-17.5 < dH < 5)
    results["MC2 Yang (2012)"]     = (Omega >= 1.1) and (delta_r <= 6.6)
    results["MC3 Guo (2013)"]      = (-11.6 < dH < 3.2) and (delta_r <= 6.6)
    results["MC4 Wang (2014) γ"]   = gamma <= 1.175
    results["MC5 Poletti (2014)"]  = (1 < delta_r < 6) and (3 < delta_chi < 6)

    if ratio > 0.9:
        results["MC7 Wang (2014)"] = (-15 <= dH <= 5) and (delta_r <= 6.6)
    elif 0.5 <= ratio < 0.9:
        results["MC7 Wang (2014)"] = (dH >= -7.5) and (delta_r <= 3.3)
    else:
        results["MC7 Wang (2014)"] = False

    # --- BCC criteria ---
    if VEC < 6.87:
        results["LSS1 Guo (2011)"] = True
    elif VEC >= 8:
        results["LSS1 Guo (2011)"] = False
    else:
        results["LSS1 Guo (2011)"] = None   # mixed

    if VEC < 7.5 and (1.8 < ea < 2.3):
        results["LSS2 Poletti (2014)"] = True
    elif VEC > 7.5 and (1.6 < ea < 1.8):
        results["LSS2 Poletti (2014)"] = False
    else:
        results["LSS2 Poletti (2014)"] = None

    if ratio > 0.9:
        if VEC < 6.87:   results["LSS3 Wang (2014)"] = True
        elif VEC > 7.84: results["LSS3 Wang (2014)"] = False
        else:            results["LSS3 Wang (2014)"] = None
    elif 0.5 < ratio < 0.9:
        if VEC < 6:     results["LSS3 Wang (2014)"] = True
        elif VEC > 7.8: results["LSS3 Wang (2014)"] = False
        else:           results["LSS3 Wang (2014)"] = None
    else:
        results["LSS3 Wang (2014)"] = None

    return results

SS_KEYS  = ["MC1 Zhang (2008)", "MC2 Yang (2012)", "MC3 Guo (2013)",
            "MC4 Wang (2014) γ", "MC5 Poletti (2014)", "MC7 Wang (2014)"]
BCC_KEYS = ["LSS1 Guo (2011)", "LSS2 Poletti (2014)", "LSS3 Wang (2014)"]
ALL_KEYS = SS_KEYS + BCC_KEYS

# ============================================================
#  Monte Carlo Sampler
# ============================================================

def random_composition(n=5, lo=5, hi=35):
    """Generate a random composition summing to 100 with each value in [lo, hi]."""
    while True:
        vals = [random.uniform(lo, hi) for _ in range(n)]
        s = sum(vals)
        normed = [v / s * 100 for v in vals]
        if all(lo <= v <= hi for v in normed):
            return normed

def run_monte_carlo(elements, T, n_samples, progress_cb=None):
    """
    Returns a dict: {criterion_key: list_of_valid_compositions}
    Each composition is a list of fractions (%) matching elements order.
    Also returns 'ALL' key for compositions satisfying all criteria simultaneously.
    """
    buckets = {k: [] for k in ALL_KEYS}
    buckets["ALL"] = []

    n = len(elements)
    for idx in range(n_samples):
        if progress_cb and idx % max(1, n_samples // 200) == 0:
            progress_cb(idx / n_samples)

        pcts  = random_composition(n)
        fracs = [p / 100 for p in pcts]
        res   = get_all_criteria(elements, fracs, T)

        all_pass = True
        for key in ALL_KEYS:
            val = res.get(key)
            if val is True:
                buckets[key].append(pcts)
            elif val is None:
                buckets[key].append(pcts)   # mixed counts as partial pass
            else:
                all_pass = False

        if all_pass:
            buckets["ALL"].append(pcts)

    if progress_cb:
        progress_cb(1.0)
    return buckets

# ============================================================
#  Statistics helper
# ============================================================

def composition_stats(compositions, elements):
    """Return a DataFrame with min/max/mean/std for each element."""
    if not compositions:
        return pd.DataFrame()
    arr = np.array(compositions)   # shape (N, 5)
    rows = []
    for i, el in enumerate(elements):
        rows.append({
            "Element": el,
            "Min (%)":  round(float(arr[:, i].min()), 2),
            "Max (%)":  round(float(arr[:, i].max()), 2),
            "Mean (%)": round(float(arr[:, i].mean()), 2),
            "Std (%)":  round(float(arr[:, i].std()), 2),
        })
    return pd.DataFrame(rows)

# ============================================================
#  Streamlit App
# ============================================================

st.set_page_config(page_title="HEA Inverse Calculator", page_icon="🔁", layout="wide")

st.markdown("""
<style>
    .main { background-color: #0e1117; }
    .block-container { padding-top: 1.8rem; }
    .section-header {
        font-size: 1.05rem; font-weight: 700; color: #64ffda;
        border-bottom: 2px solid #64ffda33; padding-bottom: 8px;
        margin: 20px 0 14px 0; font-family: monospace;
        letter-spacing: 0.08em; text-transform: uppercase;
    }
    .info-card {
        background: linear-gradient(135deg, #1e2130, #252a3a);
        border: 1px solid #3a3f55; border-radius: 12px;
        padding: 16px 20px; margin-bottom: 10px;
    }
    .badge-yes {
        background: linear-gradient(135deg,#0d4f3c,#1a7a5e);
        border:1px solid #64ffda; color:#64ffda;
        padding:4px 12px; border-radius:16px;
        font-weight:700; font-size:0.8rem; display:inline-block;
    }
    .badge-no {
        background: linear-gradient(135deg,#4f0d0d,#7a1a1a);
        border:1px solid #ff6b6b; color:#ff6b6b;
        padding:4px 12px; border-radius:16px;
        font-weight:700; font-size:0.8rem; display:inline-block;
    }
    .count-chip {
        font-size:1.8rem; font-weight:900;
        font-family:monospace; color:#64ffda;
    }
    .el-bar-label {
        font-weight:700; font-family:monospace; font-size:0.9rem;
    }
    hr { border-color:#3a3f55; }
</style>
""", unsafe_allow_html=True)

# ── Header ───────────────────────────────────────────────────
st.markdown("""
<div style='text-align:center; padding:18px 0 8px 0;'>
    <div style='font-size:2.5rem; font-weight:900;
                background:linear-gradient(90deg,#ff79c6,#bd93f9,#64ffda);
                -webkit-background-clip:text; -webkit-text-fill-color:transparent;
                font-family:monospace; letter-spacing:0.1em;'>
        🔁 HEA INVERSE CALCULATOR
    </div>
    <div style='color:#8892b0; font-size:0.9rem; margin-top:6px;'>
        Find compositions that satisfy SS &amp; BCC criteria · Monte Carlo Search
    </div>
</div>
""", unsafe_allow_html=True)
st.markdown("---")

# ── Sidebar ──────────────────────────────────────────────────
with st.sidebar:
    st.markdown("<div style='text-align:center; font-size:1.15rem; font-weight:700; color:#ff79c6; font-family:monospace; padding:6px 0;'>🔬 SEARCH SETTINGS</div>",
                unsafe_allow_html=True)
    st.markdown("---")

    T = st.slider("🌡️ Temperature (K)", 300, 3000, 1000, 50)

    st.markdown("---")
    st.markdown("<div style='color:#8892b0; font-size:0.82rem; font-weight:600; letter-spacing:0.06em;'>SELECT 5 ELEMENTS</div>",
                unsafe_allow_html=True)
    selected = st.multiselect("", options=list(ELEMENTS_DATA.keys()),
                               default=["Nb","Mo","Ta","W","V"], max_selections=5)

    st.markdown("---")
    st.markdown("<div style='color:#8892b0; font-size:0.82rem; font-weight:600; letter-spacing:0.06em;'>COMPOSITION RANGE (%)</div>",
                unsafe_allow_html=True)
    lo = st.slider("Min per element (%)", 5, 20, 5, 1)
    hi = st.slider("Max per element (%)", 20, 40, 35, 1)
    if lo >= hi:
        st.error("Min must be less than Max")

    st.markdown("---")
    n_samples = st.select_slider(
        "🎲 Monte Carlo Samples",
        options=[10_000, 20_000, 50_000, 100_000, 200_000],
        value=50_000,
        help="More samples = more accurate ranges but slower."
    )

    st.markdown("---")
    run_btn = st.button("🚀 RUN SEARCH", use_container_width=True,
                        disabled=(len(selected) != 5 or lo >= hi))

# ── Welcome ──────────────────────────────────────────────────
if not run_btn:
    st.markdown("""
    <div style='text-align:center; padding:60px 20px;'>
        <div style='font-size:4rem;'>🔁</div>
        <div style='font-size:1.15rem; color:#8892b0; margin-top:16px; line-height:2.2;'>
            Select <b style='color:#ff79c6;'>5 elements</b>, set temperature and sample count,<br>
            then press <b style='color:#64ffda;'>🚀 RUN SEARCH</b><br>
            <span style='font-size:0.9rem;'>The engine will search for compositions that produce <b style='color:#bd93f9;'>SS + BCC</b></span>
        </div>
    </div>
    """, unsafe_allow_html=True)
    st.stop()

# ── Validation ───────────────────────────────────────────────
if len(selected) != 5:
    st.error("Please select exactly 5 elements.")
    st.stop()

# ── Run ──────────────────────────────────────────────────────
elements = selected

progress_bar = st.progress(0, text="🔄 Running Monte Carlo search...")

def update_progress(val):
    progress_bar.progress(val, text=f"🔄 Searching... {int(val*100)}%")

# Monkey-patch random_composition with user lo/hi
def random_composition_custom(n=5):
    _lo, _hi = lo, hi
    while True:
        vals = [random.uniform(_lo, _hi) for _ in range(n)]
        s = sum(vals)
        normed = [v / s * 100 for v in vals]
        if all(_lo <= v <= _hi for v in normed):
            return normed

import types
old_rc = random_composition

def run_mc_custom(elements, T, n_samples, progress_cb=None):
    buckets = {k: [] for k in ALL_KEYS}
    buckets["ALL"] = []
    n = len(elements)
    for idx in range(n_samples):
        if progress_cb and idx % max(1, n_samples // 200) == 0:
            progress_cb(idx / n_samples)
        pcts  = random_composition_custom(n)
        fracs = [p / 100 for p in pcts]
        res   = get_all_criteria(elements, fracs, T)
        all_pass = True
        for key in ALL_KEYS:
            val = res.get(key)
            if val is True:
                buckets[key].append(pcts)
            elif val is None:
                buckets[key].append(pcts)
            else:
                all_pass = False
        if all_pass:
            buckets["ALL"].append(pcts)
    if progress_cb:
        progress_cb(1.0)
    return buckets

buckets = run_mc_custom(elements, T, n_samples, update_progress)
progress_bar.empty()

# ── Results Overview ─────────────────────────────────────────
colors_el = ["#64ffda","#bd93f9","#ff79c6","#ffb86c","#8be9fd"]

st.markdown(f"""
<div style='background:linear-gradient(135deg,#1e2130,#252a3a);
            border:1px solid #3a3f55; border-radius:14px;
            padding:18px 26px; margin-bottom:20px;
            display:flex; align-items:center; gap:30px; flex-wrap:wrap;'>
    <div>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace; letter-spacing:0.1em;'>ELEMENTS</div>
        <div style='font-size:1.4rem; font-weight:800; color:#ccd6f6; font-family:monospace;'>
            {"  ·  ".join(elements)}
        </div>
    </div>
    <div>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace;'>TEMPERATURE</div>
        <div style='color:#64ffda; font-size:1.3rem; font-weight:700; font-family:monospace;'>{T} K</div>
    </div>
    <div>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace;'>SAMPLES</div>
        <div style='color:#ffb86c; font-size:1.3rem; font-weight:700; font-family:monospace;'>{n_samples:,}</div>
    </div>
    <div style='margin-left:auto; text-align:right;'>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace;'>ALL CRITERIA PASS</div>
        <div style='color:#ff79c6; font-size:2rem; font-weight:900; font-family:monospace;'>
            {len(buckets["ALL"]):,}
        </div>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace;'>
            ({len(buckets["ALL"])/n_samples*100:.1f}% of samples)
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Criteria count summary ────────────────────────────────────
st.markdown("<div class='section-header'>📊 Valid Compositions per Criterion</div>", unsafe_allow_html=True)

crit_cols = st.columns(3)
all_criteria_display = [
    ("MC1 Zhang (2008)",    "SS",  "#64ffda"),
    ("MC2 Yang (2012)",     "SS",  "#64ffda"),
    ("MC3 Guo (2013)",      "SS",  "#64ffda"),
    ("MC4 Wang (2014) γ",   "SS",  "#64ffda"),
    ("MC5 Poletti (2014)",  "SS",  "#64ffda"),
    ("MC7 Wang (2014)",     "SS",  "#64ffda"),
    ("LSS1 Guo (2011)",     "BCC", "#bd93f9"),
    ("LSS2 Poletti (2014)", "BCC", "#bd93f9"),
    ("LSS3 Wang (2014)",    "BCC", "#bd93f9"),
]

for idx, (key, kind, clr) in enumerate(all_criteria_display):
    cnt  = len(buckets[key])
    pct  = cnt / n_samples * 100
    with crit_cols[idx % 3]:
        st.markdown(f"""
        <div class='info-card' style='border-left:4px solid {clr};'>
            <div style='color:#8892b0; font-size:0.75rem; font-family:monospace;'>{kind} · {key}</div>
            <div style='font-size:1.7rem; font-weight:900; color:{clr}; font-family:monospace;'>{cnt:,}</div>
            <div style='color:#8892b0; font-size:0.78rem;'>{pct:.1f}% of samples</div>
        </div>
        """, unsafe_allow_html=True)

# ── Tabs for each criterion ───────────────────────────────────
st.markdown("<div class='section-header'>🔎 Composition Ranges per Criterion</div>", unsafe_allow_html=True)

tab_labels = [k for k,_,_ in all_criteria_display] + ["⭐ ALL CRITERIA"]
tabs = st.tabs(tab_labels)

def render_tab(tab, key, elements, colors_el):
    comps = buckets[key]
    with tab:
        if not comps:
            st.markdown("""
            <div style='text-align:center; padding:40px; color:#ff6b6b;
                        font-family:monospace; font-size:1rem;'>
                ❌ No valid compositions found for this criterion.<br>
                Try increasing sample count or changing elements.
            </div>
            """, unsafe_allow_html=True)
            return

        st.markdown(f"""
        <div style='color:#8892b0; font-size:0.85rem; margin-bottom:12px;'>
            <b style='color:#64ffda;'>{len(comps):,}</b> compositions satisfy this criterion
        </div>
        """, unsafe_allow_html=True)

        df = composition_stats(comps, elements)

        # ── Stats table ──
        st.dataframe(
            df.style
              .format({"Min (%)": "{:.2f}", "Max (%)": "{:.2f}",
                       "Mean (%)": "{:.2f}", "Std (%)": "{:.2f}"})
              .set_properties(**{"font-family": "monospace", "font-size": "14px"})
              .background_gradient(subset=["Mean (%)"], cmap="Blues"),
            use_container_width=True, hide_index=True
        )

        # ── Visual range bars ──
        st.markdown("<div style='margin-top:18px;'></div>", unsafe_allow_html=True)
        arr = np.array(comps)
        for i, el in enumerate(elements):
            mn  = float(arr[:, i].min())
            mx  = float(arr[:, i].max())
            avg = float(arr[:, i].mean())
            col_l, col_r = st.columns([2, 10])
            with col_l:
                st.markdown(f"<div style='font-weight:700; color:{colors_el[i]}; font-family:monospace; padding-top:4px; font-size:1rem;'>{el}</div>",
                            unsafe_allow_html=True)
            with col_r:
                # range bar: fill from min% to max% of [lo,hi] scale
                scale = hi - lo if hi > lo else 1
                start_pct = (mn - lo) / scale * 100
                width_pct  = (mx - mn) / scale * 100
                avg_pos    = (avg - lo) / scale * 100
                st.markdown(f"""
                <div style='position:relative; background:#1e2130; border-radius:8px;
                            height:26px; overflow:visible; margin-top:4px;'>
                    <div style='position:absolute; left:{start_pct:.1f}%; width:{width_pct:.1f}%;
                                height:100%; background:{colors_el[i]}33;
                                border:1px solid {colors_el[i]}88; border-radius:6px;'></div>
                    <div style='position:absolute; left:{avg_pos:.1f}%; width:3px;
                                height:100%; background:{colors_el[i]}; border-radius:2px;'></div>
                    <div style='position:absolute; right:8px; top:4px;
                                color:#8892b0; font-size:0.75rem; font-family:monospace;'>
                        {mn:.1f}% – {mx:.1f}%  (avg {avg:.1f}%)
                    </div>
                </div>
                """, unsafe_allow_html=True)
            st.markdown("<div style='margin-bottom:6px;'></div>", unsafe_allow_html=True)

        # ── Best 10 compositions ──
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown("<div style='color:#8892b0; font-size:0.82rem; font-weight:600; letter-spacing:0.06em; margin-bottom:8px;'>TOP 10 SAMPLE COMPOSITIONS</div>",
                    unsafe_allow_html=True)
        top10 = comps[:10]
        top_df = pd.DataFrame(top10, columns=elements)
        top_df = top_df.round(2)
        top_df.index = top_df.index + 1
        st.dataframe(top_df.style.format("{:.2f}%")
                              .set_properties(**{"font-family": "monospace"}),
                     use_container_width=True)

for idx, (key, _, _) in enumerate(all_criteria_display):
    render_tab(tabs[idx], key, elements, colors_el)

# ALL tab
render_tab(tabs[-1], "ALL", elements, colors_el)

# ── Summary Banner ────────────────────────────────────────────
st.markdown("---")
all_cnt = len(buckets["ALL"])
best_mc = max(SS_KEYS,  key=lambda k: len(buckets[k]))
best_ls = max(BCC_KEYS, key=lambda k: len(buckets[k]))

st.markdown(f"""
<div style='background:linear-gradient(135deg,#1e2130,#252a3a);
            border:1px solid #3a3f55; border-radius:14px;
            padding:22px 28px; display:flex; gap:40px;
            align-items:center; flex-wrap:wrap;'>
    <div style='font-size:1rem; font-weight:700; color:#8892b0; font-family:monospace;'>SUMMARY</div>
    <div>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace;'>ALL CRITERIA</div>
        <div style='color:#ff79c6; font-size:1.4rem; font-weight:900; font-family:monospace;'>
            {all_cnt:,} compositions
        </div>
    </div>
    <div>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace;'>BEST SS CRITERION</div>
        <div style='color:#64ffda; font-size:1rem; font-weight:700; font-family:monospace;'>
            {best_mc}  ({len(buckets[best_mc]):,})
        </div>
    </div>
    <div>
        <div style='color:#8892b0; font-size:0.72rem; font-family:monospace;'>BEST BCC CRITERION</div>
        <div style='color:#bd93f9; font-size:1rem; font-weight:700; font-family:monospace;'>
            {best_ls}  ({len(buckets[best_ls]):,})
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

st.markdown("<br><br>", unsafe_allow_html=True)
