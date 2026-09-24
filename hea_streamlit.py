
import math
import streamlit as st
import pandas as pd

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
#  Calculation Functions
# ============================================================

def calc_r_bar(elements, fractions):
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]]["r"] for i in range(len(elements)))

def calc_delta_r(elements, fractions):
    r_bar = calc_r_bar(elements, fractions)
    return math.sqrt(sum(fractions[i] * (1 - ELEMENTS_DATA[elements[i]]["r"] / r_bar)**2
                         for i in range(len(elements)))) * 100

def calc_delta_S_mix(fractions):
    return -R * sum(c * math.log(c) for c in fractions if c > 0)

def calc_delta_H_mix(elements, fractions):
    total = 0
    n = len(elements)
    for i in range(n):
        for j in range(n):
            if i != j:
                hij = DELTA_H_IJ.get((elements[i], elements[j]), 0)
                total += 4 * hij * fractions[i] * fractions[j]
    return total / 2

def calc_VEC(elements, fractions):
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]]["VEC"] for i in range(len(elements)))

def calc_E_DFT(elements, fractions, structure="BCC"):
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]][f"E_{structure}"] for i in range(len(elements)))

def calc_delta_H_elastic(elements, fractions, T):
    data = ELEMENTS_DATA
    def v0(el): return (4/3) * math.pi * (data[el]["r"] * 1e-12)**3
    def v_T(el): return v0(el) * (1 + data[el]["alpha"] * (T - 298))
    v_vals = [v_T(el) for el in elements]
    B_vals = [data[el]["B"] * 1e9 for el in elements]
    c = fractions
    v_bar = sum(c[i]*B_vals[i]*v_vals[i] for i in range(len(elements))) / \
            sum(c[i]*B_vals[i] for i in range(len(elements)))
    total = sum(c[i]*B_vals[i]*(v_vals[i]-v_bar)**2/(2*v_vals[i]) for i in range(len(elements)))
    return total / 1.602e-19 / 6.022e23 * 1000

def calc_gamma(elements, fractions):
    r_bar = calc_r_bar(elements, fractions)
    r_vals = [ELEMENTS_DATA[el]["r"] for el in elements]
    r_min, r_max = min(r_vals), max(r_vals)
    num = 1 - (r_min / (r_min + r_bar))**2
    den = 1 - (r_max / (r_max + r_bar))**2
    return float('inf') if den == 0 else num / den

def calc_Tm_mean(elements, fractions):
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]]["Tm"] for i in range(len(elements)))

def calc_Omega(Tm, dS, dH):
    return float('inf') if abs(dH) < 1e-10 else Tm * dS / (abs(dH) * 1000)

def calc_delta_chi(elements, fractions):
    chi_vals = [ELEMENTS_DATA[el]["chi"] for el in elements]
    chi_bar = sum(fractions[i] * chi_vals[i] for i in range(len(elements)))
    return math.sqrt(sum(fractions[i] * (1 - chi_vals[i]/chi_bar)**2
                         for i in range(len(elements)))) * 100

def calc_ea(elements, fractions):
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]]["ea"] for i in range(len(elements)))

# ============================================================
#  Criteria
# ============================================================

def check_SS(dH, delta_r, Omega, gamma, delta_chi, VEC, ea, T, Tm):
    results = {}
    results["MC1 — Zhang (2008)"] = ("SS", (0.5 < delta_r < 6.5) and (-17.5 < dH < 5))
    results["MC2 — Yang (2012)"]  = ("SS", (Omega >= 1.1) and (delta_r <= 6.6))
    c3 = (-11.6 < dH < 3.2) and (delta_r <= 6.6)
    results["MC3 — Guo (2013)"]   = ("SS", c3)
    results["MC4 — Wang (2014) γ"] = ("SS", gamma <= 1.175)
    results["MC5 — Poletti (2014)"] = ("SS", (1 < delta_r < 6) and (3 < delta_chi < 6))
    ratio = T / Tm
    if ratio > 0.9:
        results["MC7 — Wang (2014) [T/Tm>0.9]"] = ("SS", (-15 <= dH <= 5) and (delta_r <= 6.6))
    elif 0.5 <= ratio < 0.9:
        results["MC7 — Wang (2014) [0.5≤T/Tm<0.9]"] = ("SS", (dH >= -7.5) and (delta_r <= 3.3))
    else:
        results["MC7 — Wang (2014)"] = ("OUT_OF_RANGE", False)
    return results

def check_BCC(VEC, ea, T, Tm):
    results = {}
    if VEC < 6.87:
        results["LSS1 — Guo (2011)"] = ("BCC", True)
    elif VEC >= 8:
        results["LSS1 — Guo (2011)"] = ("FCC", False)
    else:
        results["LSS1 — Guo (2011)"] = ("FCC+BCC", None)

    if VEC < 7.5 and (1.8 < ea < 2.3):
        results["LSS2 — Poletti (2014)"] = ("BCC", True)
    elif VEC > 7.5 and (1.6 < ea < 1.8):
        results["LSS2 — Poletti (2014)"] = ("FCC", False)
    else:
        results["LSS2 — Poletti (2014)"] = ("Unresolved", None)

    ratio = T / Tm
    if ratio > 0.9:
        if VEC < 6.87:   lss3 = ("BCC", True)
        elif VEC > 7.84: lss3 = ("FCC", False)
        else:            lss3 = ("Mixed", None)
        results["LSS3 — Wang (2014) [T/Tm>0.9]"] = lss3
    elif 0.5 < ratio < 0.9:
        if VEC < 6:     lss3 = ("BCC", True)
        elif VEC > 7.8: lss3 = ("FCC", False)
        else:           lss3 = ("Mixed", None)
        results["LSS3 — Wang (2014) [0.5<T/Tm<0.9]"] = lss3
    else:
        results["LSS3 — Wang (2014)"] = ("Out of range", None)
    return results

# ============================================================
#  Streamlit App
# ============================================================

st.set_page_config(
    page_title="HEA Calculator",
    page_icon="⚗️",
    layout="wide",
)

# ── Custom CSS ──────────────────────────────────────────────
st.markdown("""
<style>
    .main { background-color: #0e1117; }
    .block-container { padding-top: 2rem; }
    .metric-card {
        background: linear-gradient(135deg, #1e2130, #252a3a);
        border: 1px solid #3a3f55;
        border-radius: 12px;
        padding: 18px 20px;
        margin-bottom: 10px;
        box-shadow: 0 4px 15px rgba(0,0,0,0.3);
    }
    .metric-label {
        font-size: 0.78rem;
        color: #8892b0;
        font-family: monospace;
        letter-spacing: 0.05em;
        text-transform: uppercase;
    }
    .metric-value {
        font-size: 1.5rem;
        font-weight: 700;
        color: #e6f1ff;
        font-family: monospace;
    }
    .metric-unit {
        font-size: 0.75rem;
        color: #64ffda;
        margin-left: 6px;
    }
    .badge-yes {
        background: linear-gradient(135deg, #0d4f3c, #1a7a5e);
        border: 1px solid #64ffda;
        color: #64ffda;
        padding: 5px 14px;
        border-radius: 20px;
        font-weight: 700;
        font-size: 0.82rem;
        display: inline-block;
    }
    .badge-no {
        background: linear-gradient(135deg, #4f0d0d, #7a1a1a);
        border: 1px solid #ff6b6b;
        color: #ff6b6b;
        padding: 5px 14px;
        border-radius: 20px;
        font-weight: 700;
        font-size: 0.82rem;
        display: inline-block;
    }
    .badge-mix {
        background: linear-gradient(135deg, #3a3000, #5a4a00);
        border: 1px solid #ffd700;
        color: #ffd700;
        padding: 5px 14px;
        border-radius: 20px;
        font-weight: 700;
        font-size: 0.82rem;
        display: inline-block;
    }
    .section-header {
        font-size: 1.1rem;
        font-weight: 700;
        color: #64ffda;
        border-bottom: 2px solid #64ffda33;
        padding-bottom: 8px;
        margin-bottom: 16px;
        font-family: monospace;
        letter-spacing: 0.08em;
        text-transform: uppercase;
    }
    .stable-badge {
        background: linear-gradient(135deg, #1a0533, #2d0a57);
        border: 2px solid #bd93f9;
        color: #bd93f9;
        padding: 8px 20px;
        border-radius: 25px;
        font-size: 1.1rem;
        font-weight: 700;
        display: inline-block;
        margin-top: 4px;
    }
    .alloy-title {
        font-size: 1.6rem;
        font-weight: 800;
        color: #ccd6f6;
        font-family: monospace;
        letter-spacing: 0.12em;
    }
    div[data-testid="stSlider"] label {
        color: #8892b0 !important;
        font-size: 0.85rem !important;
    }
    .stSelectbox label, .stMultiSelect label {
        color: #8892b0 !important;
    }
    hr { border-color: #3a3f55; }
</style>
""", unsafe_allow_html=True)

# ── Header ──────────────────────────────────────────────────
st.markdown("""
<div style='text-align:center; padding: 20px 0 10px 0;'>
    <div style='font-size:2.8rem; font-weight:900; 
                background: linear-gradient(90deg, #64ffda, #bd93f9, #ff79c6);
                -webkit-background-clip: text; -webkit-text-fill-color: transparent;
                font-family: monospace; letter-spacing: 0.1em;'>
        ⚗️ HEA PROPERTY CALCULATOR
    </div>
    <div style='color:#8892b0; font-size:0.95rem; margin-top:6px;'>
        High-Entropy Alloy · Refractory System · Nb · Mo · Ta · W · V · Cr · Ti · Zr · Hf
    </div>
</div>
""", unsafe_allow_html=True)

st.markdown("---")

# ── Sidebar ─────────────────────────────────────────────────
with st.sidebar:
    st.markdown("""
    <div style='text-align:center; padding:10px 0;'>
        <div style='font-size:1.3rem; font-weight:700; color:#64ffda; font-family:monospace;'>
            🔬 INPUT PARAMETERS
        </div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("---")

    T = st.slider("🌡️ Temperature (K)", min_value=300, max_value=3000,
                  value=1000, step=50,
                  help="Processing temperature used in elastic enthalpy, Omega, and structure criteria.")

    st.markdown("---")
    st.markdown("<div style='color:#8892b0; font-size:0.85rem; font-weight:600;'>SELECT 5 ELEMENTS</div>",
                unsafe_allow_html=True)

    selected_elements = st.multiselect(
        "",
        options=list(ELEMENTS_DATA.keys()),
        default=["Nb", "Mo", "Ta", "W", "Ti"],
        max_selections=5,
    )

    st.markdown("---")

    fractions_pct = {}
    total_pct = 0.0

    if len(selected_elements) == 5:
        st.markdown("<div style='color:#8892b0; font-size:0.85rem; font-weight:600;'>MOLAR FRACTIONS (%)</div>",
                    unsafe_allow_html=True)
        for el in selected_elements:
            fractions_pct[el] = st.slider(f"{el}", min_value=5, max_value=35,
                                           value=20, step=1, key=f"frac_{el}")
        total_pct = sum(fractions_pct.values())

        if abs(total_pct - 100) < 0.01:
            st.success(f"✅ Total = {total_pct}%")
        else:
            st.error(f"❌ Total = {total_pct}% (must be 100%)")
    else:
        st.warning("⚠️ Please select exactly 5 elements")

    st.markdown("---")
    calculate = st.button("⚡ CALCULATE", use_container_width=True,
                          disabled=(len(selected_elements) != 5 or abs(total_pct - 100) > 0.01))

# ── Main Content ─────────────────────────────────────────────
if not calculate:
    st.markdown("""
    <div style='text-align:center; padding: 60px 20px;'>
        <div style='font-size:4rem;'>🧪</div>
        <div style='font-size:1.3rem; color:#8892b0; margin-top:16px;'>
            Select 5 elements, set molar fractions (total = 100%),<br>
            then press <b style='color:#64ffda;'>CALCULATE</b>
        </div>
    </div>
    """, unsafe_allow_html=True)
    st.stop()

# ── Run Calculations ─────────────────────────────────────────
elements  = selected_elements
fractions = [fractions_pct[el] / 100 for el in elements]

delta_r   = calc_delta_r(elements, fractions)
dS        = calc_delta_S_mix(fractions)
dH        = calc_delta_H_mix(elements, fractions)
VEC       = calc_VEC(elements, fractions)
dH_el     = calc_delta_H_elastic(elements, fractions, T)
gamma     = calc_gamma(elements, fractions)
Tm        = calc_Tm_mean(elements, fractions)
Omega     = calc_Omega(Tm, dS, dH)
delta_chi = calc_delta_chi(elements, fractions)
ea        = calc_ea(elements, fractions)
dG        = dH * 1000 - T * dS
dG_MIDMA  = (dH + dH_el) * 1000 - T * dS
E_BCC     = calc_E_DFT(elements, fractions, "BCC")
E_FCC     = calc_E_DFT(elements, fractions, "FCC")
E_HCP     = calc_E_DFT(elements, fractions, "HCP")
stable    = min([("BCC", E_BCC), ("FCC", E_FCC), ("HCP", E_HCP)], key=lambda x: x[1])

# ── Alloy Name & Composition Bar ─────────────────────────────
alloy_name = "".join([f"{el}{int(fractions_pct[el])}" for el in elements])
st.markdown(f"""
<div style='background: linear-gradient(135deg,#1e2130,#252a3a);
            border:1px solid #3a3f55; border-radius:14px;
            padding:20px 28px; margin-bottom:24px;
            display:flex; align-items:center; gap:20px;'>
    <div>
        <div style='color:#8892b0; font-size:0.75rem; text-transform:uppercase;
                    letter-spacing:0.1em; font-family:monospace;'>Alloy Composition</div>
        <div class='alloy-title'>{alloy_name}</div>
    </div>
    <div style='margin-left:auto; text-align:right;'>
        <div style='color:#8892b0; font-size:0.75rem; font-family:monospace;'>Temperature</div>
        <div style='color:#64ffda; font-size:1.4rem; font-weight:700; font-family:monospace;'>{T} K</div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Composition Progress Bars ─────────────────────────────────
colors = ["#64ffda","#bd93f9","#ff79c6","#ffb86c","#8be9fd"]
bars_html = ""
for i, el in enumerate(elements):
    pct = fractions_pct[el]
    bars_html += f"""
    <div style='display:flex; align-items:center; gap:12px; margin-bottom:8px;'>
        <div style='width:32px; font-weight:700; color:{colors[i]}; font-family:monospace;'>{el}</div>
        <div style='flex:1; background:#1e2130; border-radius:8px; height:22px; overflow:hidden;'>
            <div style='width:{pct}%; background:{colors[i]}; height:100%;
                        border-radius:8px; display:flex; align-items:center;
                        padding-left:8px; font-size:0.75rem; font-weight:700; color:#0e1117;'>
                {pct}%
            </div>
        </div>
    </div>"""

st.markdown(f"""
<div style='background:#1a1e2e; border:1px solid #3a3f55; border-radius:12px; padding:18px 22px; margin-bottom:24px;'>
    <div class='section-header'>Composition</div>
    {bars_html}
</div>
""", unsafe_allow_html=True)

# ── Calculated Parameters ─────────────────────────────────────
st.markdown("<div class='section-header'>📊 Calculated Parameters</div>", unsafe_allow_html=True)

params = [
    ("δr",      f"{delta_r:.4f}", "%",          "Atomic size difference"),
    ("ΔS_mix",  f"{dS:.4f}",      "J/(mol·K)",  "Mixing entropy"),
    ("ΔH_mix",  f"{dH:.4f}",      "kJ/mol",     "Mixing enthalpy"),
    ("VEC",     f"{VEC:.4f}",     "",            "Avg valence electron concentration"),
    ("ΔH_el",   f"{dH_el:.4f}",   "kJ/mol",     "Elastic enthalpy"),
    ("γ",       f"{gamma:.4f}",   "",            "Gamma parameter"),
    ("Ω",       f"{Omega:.4f}",   "",            "Omega parameter"),
]

cols = st.columns(4)
for i, (label, val, unit, desc) in enumerate(params):
    with cols[i % 4]:
        st.markdown(f"""
        <div class='metric-card'>
            <div class='metric-label'>{desc}</div>
            <div class='metric-value'>{label} = {val}<span class='metric-unit'>{unit}</span></div>
        </div>
        """, unsafe_allow_html=True)

# ── DFT Energy & Stable Structure ────────────────────────────
st.markdown("<br>", unsafe_allow_html=True)
c1, c2, c3, c4 = st.columns(4)

dft_info = [
    (c1, "E_DFT  BCC", f"{E_BCC:.4f}", "eV/atom", stable[0]=="BCC"),
    (c2, "E_DFT  FCC", f"{E_FCC:.4f}", "eV/atom", stable[0]=="FCC"),
    (c3, "E_DFT  HCP", f"{E_HCP:.4f}", "eV/atom", stable[0]=="HCP"),
]
for col, label, val, unit, is_stable in dft_info:
    border = "#bd93f9" if is_stable else "#3a3f55"
    glow   = "box-shadow: 0 0 18px #bd93f944;" if is_stable else ""
    with col:
        st.markdown(f"""
        <div style='background:linear-gradient(135deg,#1e2130,#252a3a);
                    border:2px solid {border}; border-radius:12px;
                    padding:18px 20px; {glow}'>
            <div class='metric-label'>{label}</div>
            <div class='metric-value'>{val}<span class='metric-unit'>{unit}</span></div>
            {"<div class='stable-badge'>★ STABLE</div>" if is_stable else ""}
        </div>
        """, unsafe_allow_html=True)

with c4:
    st.markdown(f"""
    <div style='background:linear-gradient(135deg,#1a0533,#2d0a57);
                border:2px solid #bd93f9; border-radius:12px;
                padding:18px 20px; box-shadow:0 0 18px #bd93f944;'>
        <div class='metric-label'>Most Stable Structure</div>
        <div style='font-size:2rem; font-weight:900; color:#bd93f9; font-family:monospace;
                    margin-top:6px;'>{stable[0]}</div>
        <div style='color:#8892b0; font-size:0.78rem; font-family:monospace;'>{stable[1]:.4f} eV/atom</div>
    </div>
    """, unsafe_allow_html=True)

# ── Gibbs Free Energy ─────────────────────────────────────────
st.markdown("<br>", unsafe_allow_html=True)
st.markdown("<div class='section-header'>🔥 Gibbs Free Energy</div>", unsafe_allow_html=True)

g1, g2 = st.columns(2)
with g1:
    color_g  = "#64ffda" if dG  < 0 else "#ff6b6b"
    st.markdown(f"""
    <div class='metric-card' style='border-color:{color_g}44;'>
        <div class='metric-label'>ΔG = ΔH_mix − T·ΔS_mix</div>
        <div style='font-size:1.6rem; font-weight:700; color:{color_g}; font-family:monospace;'>
            {dG:.2f} <span style='font-size:0.8rem;'>J/mol</span>
        </div>
    </div>
    """, unsafe_allow_html=True)
with g2:
    color_gm = "#64ffda" if dG_MIDMA < 0 else "#ff6b6b"
    st.markdown(f"""
    <div class='metric-card' style='border-color:{color_gm}44;'>
        <div class='metric-label'>ΔG_MIDMA = (ΔH_mix + ΔH_el) − T·ΔS_mix</div>
        <div style='font-size:1.6rem; font-weight:700; color:{color_gm}; font-family:monospace;'>
            {dG_MIDMA:.2f} <span style='font-size:0.8rem;'>J/mol</span>
        </div>
    </div>
    """, unsafe_allow_html=True)

# ── SS Criteria ───────────────────────────────────────────────
st.markdown("<br>", unsafe_allow_html=True)
st.markdown("<div class='section-header'>🔍 Solid Solution (SS) Criteria</div>", unsafe_allow_html=True)

ss_results = check_SS(dH, delta_r, Omega, gamma, delta_chi, VEC, ea, T, Tm)
ss_yes = sum(1 for _, v in ss_results.items() if v[1] is True)

ss_rows = []
for model, (label, result) in ss_results.items():
    if result is True:
        badge = f"<span class='badge-yes'>✓ SS</span>"
    elif result is False:
        badge = f"<span class='badge-no'>✗ {label}</span>"
    else:
        badge = f"<span class='badge-mix'>~ OUT OF RANGE</span>"
    ss_rows.append((model, badge))

table_rows = "".join([f"""
<tr style='border-bottom:1px solid #2a2f42;'>
    <td style='padding:10px 16px; color:#ccd6f6; font-family:monospace; font-size:0.88rem;'>{m}</td>
    <td style='padding:10px 16px;'>{b}</td>
</tr>""" for m, b in ss_rows])

st.markdown(f"""
<div style='background:#1a1e2e; border:1px solid #3a3f55; border-radius:12px; overflow:hidden; margin-bottom:12px;'>
    <table style='width:100%; border-collapse:collapse;'>
        <thead>
            <tr style='background:#252a3a;'>
                <th style='padding:12px 16px; text-align:left; color:#64ffda;
                           font-family:monospace; font-size:0.82rem; letter-spacing:0.08em;'>MODEL</th>
                <th style='padding:12px 16px; text-align:left; color:#64ffda;
                           font-family:monospace; font-size:0.82rem; letter-spacing:0.08em;'>RESULT</th>
            </tr>
        </thead>
        <tbody>{table_rows}</tbody>
    </table>
</div>
<div style='color:#8892b0; font-size:0.88rem; padding:4px 2px;'>
    SS confirmed by <b style='color:#64ffda;'>{ss_yes}</b> out of
    <b style='color:#ccd6f6;'>{len(ss_results)}</b> models
</div>
""", unsafe_allow_html=True)

# ── BCC Criteria ──────────────────────────────────────────────
st.markdown("<br>", unsafe_allow_html=True)
st.markdown("<div class='section-header'>🔷 BCC Structure Criteria</div>", unsafe_allow_html=True)

bcc_results = check_BCC(VEC, ea, T, Tm)
bcc_yes = sum(1 for _, v in bcc_results.items() if v[1] is True)

bcc_rows = []
for model, (label, result) in bcc_results.items():
    if result is True:
        badge = f"<span class='badge-yes'>✓ BCC</span>"
    elif result is False:
        badge = f"<span class='badge-no'>✗ {label}</span>"
    else:
        badge = f"<span class='badge-mix'>~ {label}</span>"
    bcc_rows.append((model, badge))

table_rows_bcc = "".join([f"""
<tr style='border-bottom:1px solid #2a2f42;'>
    <td style='padding:10px 16px; color:#ccd6f6; font-family:monospace; font-size:0.88rem;'>{m}</td>
    <td style='padding:10px 16px;'>{b}</td>
</tr>""" for m, b in bcc_rows])

st.markdown(f"""
<div style='background:#1a1e2e; border:1px solid #3a3f55; border-radius:12px; overflow:hidden; margin-bottom:12px;'>
    <table style='width:100%; border-collapse:collapse;'>
        <thead>
            <tr style='background:#252a3a;'>
                <th style='padding:12px 16px; text-align:left; color:#64ffda;
                           font-family:monospace; font-size:0.82rem; letter-spacing:0.08em;'>MODEL</th>
                <th style='padding:12px 16px; text-align:left; color:#64ffda;
                           font-family:monospace; font-size:0.82rem; letter-spacing:0.08em;'>RESULT</th>
            </tr>
        </thead>
        <tbody>{table_rows_bcc}</tbody>
    </table>
</div>
<div style='color:#8892b0; font-size:0.88rem; padding:4px 2px;'>
    BCC confirmed by <b style='color:#64ffda;'>{bcc_yes}</b> out of
    <b style='color:#ccd6f6;'>{len(bcc_results)}</b> models
</div>
""", unsafe_allow_html=True)

# ── Summary Banner ────────────────────────────────────────────
st.markdown("<br>", unsafe_allow_html=True)
ss_color  = "#64ffda" if ss_yes  >= len(ss_results)//2+1  else "#ff6b6b"
bcc_color = "#64ffda" if bcc_yes >= len(bcc_results)//2+1 else "#ff6b6b"

st.markdown(f"""
<div style='background:linear-gradient(135deg,#1e2130,#252a3a);
            border:1px solid #3a3f55; border-radius:14px;
            padding:22px 28px; display:flex; gap:40px; align-items:center;'>
    <div style='font-size:1.05rem; font-weight:700; color:#8892b0; font-family:monospace;'>SUMMARY</div>
    <div>
        <div style='color:#8892b0; font-size:0.75rem; font-family:monospace;'>SOLID SOLUTION</div>
        <div style='color:{ss_color}; font-size:1.3rem; font-weight:800; font-family:monospace;'>
            {ss_yes}/{len(ss_results)} models confirm SS
        </div>
    </div>
    <div>
        <div style='color:#8892b0; font-size:0.75rem; font-family:monospace;'>BCC STRUCTURE</div>
        <div style='color:{bcc_color}; font-size:1.3rem; font-weight:800; font-family:monospace;'>
            {bcc_yes}/{len(bcc_results)} models confirm BCC
        </div>
    </div>
    <div style='margin-left:auto;'>
        <div style='color:#8892b0; font-size:0.75rem; font-family:monospace;'>STABLE PHASE</div>
        <div style='color:#bd93f9; font-size:1.3rem; font-weight:800; font-family:monospace;'>{stable[0]}</div>
    </div>
</div>
""", unsafe_allow_html=True)

st.markdown("<br><br>", unsafe_allow_html=True)
