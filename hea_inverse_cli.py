
"""
HEA Inverse Calculator — Command Line (No Streamlit)
=====================================================
Finds compositions of 5 refractory elements that satisfy
SS (Solid Solution) and BCC criteria via Monte Carlo search.

Usage:
    python hea_inverse_cli.py

Requirements:
    pip install pandas numpy
"""

import math
import random
import os
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
#  Colors for terminal (ANSI)
# ============================================================

class C:
    RESET  = "\033[0m"
    BOLD   = "\033[1m"
    CYAN   = "\033[96m"
    GREEN  = "\033[92m"
    YELLOW = "\033[93m"
    RED    = "\033[91m"
    PURPLE = "\033[95m"
    BLUE   = "\033[94m"
    GREY   = "\033[90m"
    WHITE  = "\033[97m"

def cprint(text, color=C.WHITE):
    print(color + text + C.RESET)

def header(title):
    w = 64
    print()
    cprint("=" * w, C.CYAN)
    cprint(f"  {title}", C.BOLD + C.CYAN)
    cprint("=" * w, C.CYAN)

def subheader(title):
    print()
    cprint(f"  ── {title} ──", C.YELLOW)
    cprint("  " + "─" * 58, C.GREY)

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
#  All Criteria
# ============================================================

SS_KEYS  = ["MC1 Zhang (2008)", "MC2 Yang (2012)", "MC3 Guo (2013)",
            "MC4 Wang γ (2014)", "MC5 Poletti (2014)", "MC7 Wang (2014)"]
BCC_KEYS = ["LSS1 Guo (2011)", "LSS2 Poletti (2014)", "LSS3 Wang (2014)"]
ALL_KEYS = SS_KEYS + BCC_KEYS

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
    results["MC1 Zhang (2008)"]    = (0.5 < delta_r < 6.5) and (-17.5 < dH < 5)
    results["MC2 Yang (2012)"]     = (Omega >= 1.1) and (delta_r <= 6.6)
    results["MC3 Guo (2013)"]      = (-11.6 < dH < 3.2) and (delta_r <= 6.6)
    results["MC4 Wang γ (2014)"]   = gamma <= 1.175
    results["MC5 Poletti (2014)"]  = (1 < delta_r < 6) and (3 < delta_chi < 6)

    if ratio > 0.9:
        results["MC7 Wang (2014)"] = (-15 <= dH <= 5) and (delta_r <= 6.6)
    elif 0.5 <= ratio < 0.9:
        results["MC7 Wang (2014)"] = (dH >= -7.5) and (delta_r <= 3.3)
    else:
        results["MC7 Wang (2014)"] = False

    # BCC
    if VEC < 6.87:
        results["LSS1 Guo (2011)"] = True
    elif VEC >= 8:
        results["LSS1 Guo (2011)"] = False
    else:
        results["LSS1 Guo (2011)"] = None

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

# ============================================================
#  Monte Carlo
# ============================================================

def random_composition(n, lo, hi):
    while True:
        vals = [random.uniform(lo, hi) for _ in range(n)]
        s = sum(vals)
        normed = [v / s * 100 for v in vals]
        if all(lo <= v <= hi for v in normed):
            return normed

def run_monte_carlo(elements, T, n_samples, lo, hi):
    buckets = {k: [] for k in ALL_KEYS}
    buckets["ALL"] = []
    n = len(elements)
    bar_width = 40

    print()
    for idx in range(n_samples):
        # progress bar
        if idx % max(1, n_samples // 100) == 0:
            done = int(idx / n_samples * bar_width)
            bar  = "█" * done + "░" * (bar_width - done)
            pct  = idx / n_samples * 100
            print(f"\r  {C.CYAN}[{bar}]{C.RESET} {C.YELLOW}{pct:5.1f}%{C.RESET}  "
                  f"{C.GREY}{idx:,}/{n_samples:,}{C.RESET}", end="", flush=True)

        pcts  = random_composition(n, lo, hi)
        fracs = [p / 100 for p in pcts]
        res   = get_all_criteria(elements, fracs, T)

        all_pass = True
        for key in ALL_KEYS:
            val = res.get(key)
            if val is True or val is None:
                buckets[key].append(pcts)
            else:
                all_pass = False

        if all_pass:
            buckets["ALL"].append(pcts)

    # final bar
    bar = "█" * bar_width
    print(f"\r  {C.GREEN}[{bar}]{C.RESET} {C.GREEN}100.0%{C.RESET}  "
          f"{C.GREY}{n_samples:,}/{n_samples:,}{C.RESET}  ✓")
    return buckets

# ============================================================
#  Display helpers
# ============================================================

def print_stats(compositions, elements, key):
    if not compositions:
        cprint(f"    ❌  No valid compositions found for: {key}", C.RED)
        return
    arr  = np.array(compositions)
    cnt  = len(compositions)
    cprint(f"    Valid compositions: {cnt:,}", C.GREEN)
    print()

    # header row
    el_w = 8
    print("  " + C.BOLD + C.WHITE +
          f"  {'Element':<{el_w}}  {'Min %':>8}  {'Max %':>8}  {'Mean %':>8}  {'Std %':>8}" +
          C.RESET)
    print("  " + C.GREY + "  " + "─"*50 + C.RESET)

    for i, el in enumerate(elements):
        mn  = arr[:, i].min()
        mx  = arr[:, i].max()
        avg = arr[:, i].mean()
        std = arr[:, i].std()

        # ASCII bar (range)
        bar_len   = 20
        scale     = 30  # total pct scale
        start_pos = int(mn / scale * bar_len)
        end_pos   = int(mx / scale * bar_len)
        avg_pos   = int(avg / scale * bar_len)
        bar = list(" " * bar_len)
        for b in range(start_pos, min(end_pos+1, bar_len)):
            bar[b] = "▓"
        if 0 <= avg_pos < bar_len:
            bar[avg_pos] = "█"
        bar_str = "".join(bar)

        print("  " +
              C.CYAN  + f"  {el:<{el_w}}" + C.RESET +
              C.WHITE + f"  {mn:>7.2f}%  {mx:>7.2f}%  {avg:>7.2f}%  {std:>7.2f}%" + C.RESET +
              C.GREY  + f"  [{bar_str}]" + C.RESET)

def print_top10(compositions, elements):
    if not compositions:
        return
    cprint("\n  Top 10 sample compositions:", C.YELLOW)
    header_row = "  " + C.BOLD + C.WHITE + f"  {'#':>3}  "
    for el in elements:
        header_row += f"{el:>8}"
    print(header_row + C.RESET)
    print("  " + C.GREY + "  " + "─" * (6 + 8*len(elements)) + C.RESET)

    for rank, comp in enumerate(compositions[:10], 1):
        row = "  " + C.CYAN + f"  {rank:>3}  " + C.RESET
        for val in comp:
            row += C.WHITE + f"{val:>7.2f}%" + C.RESET
        print(row)

def save_csv(buckets, elements, filename="hea_inverse_results.csv"):
    rows = []
    for key in ALL_KEYS + ["ALL"]:
        for comp in buckets[key]:
            row = {"Criterion": key}
            for i, el in enumerate(elements):
                row[el] = round(comp[i], 4)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(filename, index=False)
    return filename

# ============================================================
#  User Input Helpers
# ============================================================

def ask_elements():
    all_els = list(ELEMENTS_DATA.keys())
    cprint("\n  Available elements: " + "  ".join(all_els), C.CYAN)
    while True:
        raw = input(C.YELLOW + "\n  Enter 5 elements separated by spaces (e.g. Nb Mo Ta W V): " + C.RESET).strip()
        chosen = [e.strip().capitalize() for e in raw.split()]
        # fix capitalization: "nb" -> "Nb", "NB" -> "Nb"
        chosen = [e[0].upper() + e[1:].lower() if len(e) > 1 else e.upper() for e in chosen]
        if len(chosen) != 5:
            cprint("  ❌ Please enter exactly 5 elements.", C.RED)
            continue
        invalid = [e for e in chosen if e not in ELEMENTS_DATA]
        if invalid:
            cprint(f"  ❌ Unknown elements: {', '.join(invalid)}", C.RED)
            continue
        if len(set(chosen)) != 5:
            cprint("  ❌ Please choose 5 different elements.", C.RED)
            continue
        return chosen

def ask_int(prompt, default, lo=None, hi=None):
    while True:
        raw = input(C.YELLOW + f"  {prompt} [default={default}]: " + C.RESET).strip()
        if raw == "":
            return default
        try:
            val = int(raw)
            if lo is not None and val < lo:
                cprint(f"  ❌ Value must be >= {lo}", C.RED); continue
            if hi is not None and val > hi:
                cprint(f"  ❌ Value must be <= {hi}", C.RED); continue
            return val
        except ValueError:
            cprint("  ❌ Please enter a valid integer.", C.RED)

def ask_choice(prompt, options, default):
    opts_str = " / ".join(str(o) for o in options)
    while True:
        raw = input(C.YELLOW + f"  {prompt} [{opts_str}] [default={default}]: " + C.RESET).strip()
        if raw == "":
            return default
        try:
            val = int(raw)
            if val in options:
                return val
            cprint(f"  ❌ Choose from: {opts_str}", C.RED)
        except ValueError:
            cprint("  ❌ Please enter a number.", C.RED)

# ============================================================
#  Main
# ============================================================

def main():
    os.system("cls" if os.name == "nt" else "clear")

    header("HEA INVERSE CALCULATOR  —  CLI Edition")
    cprint("  Finds compositions satisfying SS + BCC criteria", C.GREY)
    cprint("  via Monte Carlo random search", C.GREY)

    # ── Inputs ──────────────────────────────────────────────
    subheader("STEP 1 — Select 5 Elements")
    elements = ask_elements()
    cprint(f"\n  Selected: {' · '.join(elements)}", C.GREEN)

    subheader("STEP 2 — Temperature")
    T = ask_int("Temperature in K", default=1000, lo=300, hi=3000)

    subheader("STEP 3 — Composition Range per Element (%)")
    lo = ask_int("Minimum % per element", default=5, lo=1, hi=20)
    hi = ask_int("Maximum % per element", default=35, lo=21, hi=50)

    subheader("STEP 4 — Monte Carlo Samples")
    n_samples = ask_choice(
        "Number of random samples",
        options=[10000, 20000, 50000, 100000, 200000],
        default=50000
    )

    subheader("STEP 5 — Save CSV?")
    save_yn = input(C.YELLOW + "  Save results to CSV file? [y/N]: " + C.RESET).strip().lower()
    do_save = save_yn == "y"

    # ── Run ─────────────────────────────────────────────────
    header("RUNNING MONTE CARLO SEARCH")
    cprint(f"  Elements  : {' · '.join(elements)}", C.WHITE)
    cprint(f"  Temperature: {T} K", C.WHITE)
    cprint(f"  Range     : {lo}% – {hi}% per element", C.WHITE)
    cprint(f"  Samples   : {n_samples:,}", C.WHITE)

    buckets = run_monte_carlo(elements, T, n_samples, lo, hi)

    # ── Results overview ─────────────────────────────────────
    header("RESULTS OVERVIEW")
    total_all = len(buckets["ALL"])
    cprint(f"  ⭐ Compositions passing ALL criteria: "
           f"{C.BOLD}{C.GREEN}{total_all:,}{C.RESET}  "
           f"({total_all/n_samples*100:.2f}% of samples)", C.WHITE)

    print()
    cprint("  SS Criteria:", C.CYAN)
    for key in SS_KEYS:
        cnt = len(buckets[key])
        bar_done = int(cnt / n_samples * 30)
        bar = "█"*bar_done + "░"*(30-bar_done)
        print(f"  {C.GREY}{key:<26}{C.RESET}  "
              f"{C.GREEN}{cnt:>7,}{C.RESET}  "
              f"{C.GREY}[{bar}] {cnt/n_samples*100:.1f}%{C.RESET}")

    print()
    cprint("  BCC Criteria:", C.PURPLE)
    for key in BCC_KEYS:
        cnt = len(buckets[key])
        bar_done = int(cnt / n_samples * 30)
        bar = "█"*bar_done + "░"*(30-bar_done)
        print(f"  {C.GREY}{key:<26}{C.RESET}  "
              f"{C.PURPLE}{cnt:>7,}{C.RESET}  "
              f"{C.GREY}[{bar}] {cnt/n_samples*100:.1f}%{C.RESET}")

    # ── Detailed stats per criterion ─────────────────────────
    all_crit_keys = ALL_KEYS + ["ALL"]
    all_crit_labels = {
        **{k: f"SS  · {k}" for k in SS_KEYS},
        **{k: f"BCC · {k}" for k in BCC_KEYS},
        "ALL": "⭐ ALL CRITERIA COMBINED",
    }

    while True:
        header("DETAILED VIEW")
        cprint("  Select a criterion to see composition ranges:\n", C.WHITE)
        for i, key in enumerate(all_crit_keys, 1):
            cnt = len(buckets[key])
            label = all_crit_labels[key]
            color = C.GREEN if "MC" in key else C.PURPLE if "LSS" in key else C.YELLOW
            print(f"  {C.GREY}{i:>2}.{C.RESET}  {color}{label:<40}{C.RESET}  "
                  f"{C.WHITE}{cnt:,} valid{C.RESET}")
        print(f"\n  {C.GREY} 0.  Exit{C.RESET}")

        choice_raw = input(C.YELLOW + "\n  Enter number: " + C.RESET).strip()
        try:
            choice = int(choice_raw)
        except ValueError:
            continue

        if choice == 0:
            break
        if 1 <= choice <= len(all_crit_keys):
            key = all_crit_keys[choice - 1]
            subheader(all_crit_labels[key])
            print_stats(buckets[key], elements, key)
            print_top10(buckets[key], elements)
            input(C.GREY + "\n  Press Enter to continue..." + C.RESET)

    # ── Save CSV ─────────────────────────────────────────────
    if do_save:
        fname = save_csv(buckets, elements)
        header("CSV SAVED")
        cprint(f"  ✅ Results saved to: {fname}", C.GREEN)
        cprint(f"     Rows: {sum(len(buckets[k]) for k in all_crit_keys):,}", C.GREY)

    header("DONE")
    cprint("  Thank you for using HEA Inverse Calculator!", C.CYAN)
    print()


if __name__ == "__main__":
    main()
