
import math

# ============================================================
#  Element Database
# ============================================================

ELEMENTS_DATA = {
    "Nb": {"r": 146, "VEC": 5, "ea": 1.51, "chi": 1.6,  "Tm": 2477 + 273.15, "B": 170.3, "alpha": 7.1e-6,
           "E_BCC": -10.0466, "E_FCC": -9.7232, "E_HCP": -9.7551},
    "Mo": {"r": 139, "VEC": 6, "ea": 1.87, "chi": 2.16, "Tm": 2623 + 273.15, "B": 115.8, "alpha": 4.9e-6,
           "E_BCC": -10.7799, "E_FCC": -10.3784, "E_HCP": -10.3666},
    "Ta": {"r": 146, "VEC": 5, "ea": 1.50, "chi": 1.5,  "Tm": 3017 + 273.15, "B": 200.1, "alpha": 6.5e-6,
           "E_BCC": -11.7358, "E_FCC": -11.1896, "E_HCP": -11.4579},
    "W":  {"r": 139, "VEC": 6, "ea": 1.86, "chi": 2.36, "Tm": 3422 + 273.15, "B": 323.3, "alpha": 4.5e-6,
           "E_BCC": -12.7781, "E_FCC": -12.3115, "E_HCP": -12.2928},
    "V":  {"r": 135, "VEC": 5, "ea": 1.50, "chi": 1.63, "Tm": 1910 + 273.15, "B": 162.0, "alpha": 8.7e-6,
           "E_BCC": -8.9632, "E_FCC": -8.7150, "E_HCP": -8.7095},
    "Cr": {"r": 128, "VEC": 6, "ea": 1.87, "chi": 1.66, "Tm": 1907 + 273.15, "B": 116.7, "alpha": 6.1e-6,
           "E_BCC": -9.4655, "E_FCC": -9.0845, "E_HCP": -9.0751},
    "Ti": {"r": 147, "VEC": 4, "ea": 1.33, "chi": 1.54, "Tm": 1668 + 273.15, "B": 105.2, "alpha": 8.6e-6,
           "E_BCC": -2.2301, "E_FCC": -2.2155, "E_HCP": -2.2343},
    "Zr": {"r": 160, "VEC": 4, "ea": 1.33, "chi": 1.33, "Tm": 1855 + 273.15, "B":  83.35, "alpha": 5.5e-6,
           "E_BCC": -8.3598, "E_FCC": -8.3972, "E_HCP": -8.4354},
    "Hf": {"r": 159, "VEC": 4, "ea": 1.32, "chi": 1.30, "Tm": 2233 + 273.15, "B": 108.9, "alpha": 5.9e-6,
           "E_BCC": -9.6562, "E_FCC": -9.7613, "E_HCP": -9.8320},
}

# Binary mixing enthalpy (kJ/mol)
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

R = 8.314  # J/(mol·K)

# ============================================================
#  Calculation Functions
# ============================================================

def calc_r_bar(elements, fractions):
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]]["r"] for i in range(len(elements)))

def calc_delta_r(elements, fractions):
    r_bar = calc_r_bar(elements, fractions)
    return math.sqrt(sum(fractions[i] * (1 - ELEMENTS_DATA[elements[i]]["r"] / r_bar) ** 2
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
    key = f"E_{structure}"
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]][key] for i in range(len(elements)))

def calc_delta_H_elastic(elements, fractions, T=1000):
    data = ELEMENTS_DATA
    def v0(el): 
        return (4/3) * math.pi * (data[el]["r"] * 1e-12) ** 3
    def v_T(el): 
        return v0(el) * (1 + data[el]["alpha"] * (T - 298))

    v_vals = [v_T(el) for el in elements]
    B_vals = [data[el]["B"] * 1e9 for el in elements]   # GPa → Pa
    c = fractions

    v_bar = sum(c[i] * B_vals[i] * v_vals[i] for i in range(len(elements))) / \
            sum(c[i] * B_vals[i] for i in range(len(elements)))

    total = sum(c[i] * B_vals[i] * (v_vals[i] - v_bar) ** 2 / (2 * v_vals[i])
                for i in range(len(elements)))

    N_A = 6.02214076e23
    return total * N_A / 1000   # J/atom → kJ/mol

def calc_gamma(elements, fractions):
    r_bar = calc_r_bar(elements, fractions)
    r_vals = [ELEMENTS_DATA[el]["r"] for el in elements]
    r_min = min(r_vals)
    r_max = max(r_vals)
    num = 1 - ((r_min / (r_min + r_bar)) ** 2)
    den = 1 - ((r_max / (r_max + r_bar)) ** 2)
    if den == 0:
        return float('inf')
    return num / den

def calc_Tm_mean(elements, fractions):
    # Weighted average melting temperature — used internally for Omega, MC7, LSS3
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]]["Tm"] for i in range(len(elements)))

def calc_Omega(Tm, dS, dH):
    if abs(dH) < 1e-10:
        return float('inf')
    return Tm * dS / (abs(dH) * 1000)

def calc_delta_chi(elements, fractions):
    chi_vals = [ELEMENTS_DATA[el]["chi"] for el in elements]
    chi_bar = sum(fractions[i] * chi_vals[i] for i in range(len(elements)))
    return math.sqrt(sum(fractions[i] * (1 - chi_vals[i] / chi_bar) ** 2
                         for i in range(len(elements)))) * 100

def calc_ea(elements, fractions):
    return sum(fractions[i] * ELEMENTS_DATA[elements[i]]["ea"] for i in range(len(elements)))

def calc_delta_G(dH, T, dS):
    return dH * 1000 - T * dS  # J/mol

def calc_delta_G_MIDMA(dH_mix, dH_el, T, dS):
    return (dH_mix + dH_el) * 1000 - T * dS  # J/mol

# ============================================================
#  Solid Solution (SS) Criteria
# ============================================================

def check_SS(dH, delta_r, Omega, gamma, delta_chi, VEC, ea, T, Tm):
    results = {}

    # MC1 Zhang (2008)
    c1 = (0.5 < delta_r < 6.5) and (-17.5 < dH < 5)
    results["MC1 Zhang (2008)"] = "YES - SS" if c1 else "NO  - IM/BMG"

    # MC2 Yang (2012)
    c2 = (Omega >= 1.1) and (delta_r <= 6.6)
    results["MC2 Yang (2012)"] = "YES - SS" if c2 else "NO  - IM"

    # MC3 Guo (2013)
    c3 = (-11.6 < dH < 3.2) and (delta_r <= 6.6)
    results["MC3 Guo (2013)"] = "YES - SS" if c3 else ("NO  - BMG" if (dH < -12 and delta_r > 6.4) else "NO  - IM/BMG")

    # MC4 Wang (2014) - gamma
    results["MC4 Wang (2014) [gamma]"] = "YES - SS" if (gamma <= 1.175) else "NO  - IM/BMG"

    # MC5 Poletti (2014)
    c5 = (1 < delta_r < 6) and (3 < delta_chi < 6)
    results["MC5 Poletti (2014)"] = "YES - SS" if c5 else "NO  - IM"

    # MC7 Wang (2014) - T/Tm
    ratio = T / Tm
    if ratio > 0.9:
        c7 = (-15 <= dH <= 5) and (delta_r <= 6.6)
        results["MC7 Wang (2014) [T/Tm>0.9]"] = "YES - SS" if c7 else "NO  - IM"
    elif 0.5 <= ratio < 0.9:
        c7 = (dH >= -7.5) and (delta_r <= 3.3)
        results["MC7 Wang (2014) [0.5<=T/Tm<0.9]"] = "YES - SS" if c7 else "NO  - IM"
    else:
        results["MC7 Wang (2014)"] = "WARN - T/Tm < 0.5 (out of range)"

    return results

# ============================================================
#  BCC Structure Criteria
# ============================================================

def check_BCC(VEC, ea, T, Tm):
    results = {}

    # LSS1 Guo (2011)
    if VEC < 6.87:
        results["LSS1 Guo (2011)"] = "YES - BCC"
    elif VEC >= 8:
        results["LSS1 Guo (2011)"] = "NO  - FCC"
    else:
        results["LSS1 Guo (2011)"] = "MIX - FCC+BCC"

    # LSS2 Poletti (2014)
    if VEC < 7.5 and (1.8 < ea < 2.3):
        results["LSS2 Poletti (2014)"] = "YES - BCC"
    elif VEC > 7.5 and (1.6 < ea < 1.8):
        results["LSS2 Poletti (2014)"] = "NO  - FCC"
    else:
        results["LSS2 Poletti (2014)"] = "UNRESOLVED"

    # LSS3 Wang (2014)
    ratio = T / Tm
    if ratio > 0.9:
        if VEC < 6.87:
            lss3 = "YES - BCC"
        elif VEC > 7.84:
            lss3 = "NO  - FCC"
        else:
            lss3 = "MIX - Mixed"
        results["LSS3 Wang (2014) [T/Tm>0.9]"] = lss3
    elif 0.5 < ratio < 0.9:
        if VEC < 6:
            lss3 = "YES - BCC"
        elif VEC > 7.8:
            lss3 = "NO  - FCC"
        else:
            lss3 = "MIX - Mixed"
        results["LSS3 Wang (2014) [0.5<T/Tm<0.9]"] = lss3
    else:
        results["LSS3 Wang (2014)"] = "WARN - T/Tm out of range"

    return results

# ============================================================
#  User Interface
# ============================================================

AVAILABLE_ELEMENTS = list(ELEMENTS_DATA.keys())

def get_temperature():
    print("\n" + "="*60)
    print("  Enter temperature (Kelvin)")
    print("  [Default: 1000 K — press Enter to confirm]")
    val = input("  T (K): ").strip()
    if val == "":
        return 1000.0
    try:
        return float(val)
    except ValueError:
        print("  Invalid input — using default 1000 K")
        return 1000.0

def select_elements():
    print("\n" + "="*60)
    print("  Available elements:")
    for i, el in enumerate(AVAILABLE_ELEMENTS, 1):
        print(f"   {i:2}. {el}")
    print("="*60)
    print("  Select exactly 5 elements (comma or space separated)")
    print("  Example:  Nb Mo Ta W V   or   1 2 3 4 5")

    while True:
        raw = input("\n  Your selection: ").strip()
        tokens = raw.replace(",", " ").split()
        selected = []
        valid = True
        for t in tokens:
            if t.isdigit():
                idx = int(t) - 1
                if 0 <= idx < len(AVAILABLE_ELEMENTS):
                    selected.append(AVAILABLE_ELEMENTS[idx])
                else:
                    print(f"  Invalid number: {t}")
                    valid = False
                    break
            elif t in AVAILABLE_ELEMENTS:
                selected.append(t)
            else:
                print(f"  Unknown element: '{t}'")
                valid = False
                break

        if not valid:
            continue
        if len(selected) != 5:
            print(f"  Exactly 5 elements required (you selected {len(selected)})")
            continue
        if len(set(selected)) != 5:
            print("  Duplicate elements detected")
            continue

        print(f"\n  Selected: {' | '.join(selected)}")
        return selected

def get_fractions(elements):
    print("\n" + "="*60)
    print("  Enter molar fraction for each element (%)")
    print("  Each value must be between 5 and 35 — total must equal 100")

    while True:
        fractions_pct = []
        for el in elements:
            while True:
                try:
                    val = float(input(f"  {el} (%): "))
                    if not (5 <= val <= 35):
                        print("  Value must be between 5 and 35")
                        continue
                    fractions_pct.append(val)
                    break
                except ValueError:
                    print("  Please enter a valid number")

        total = sum(fractions_pct)
        if abs(total - 100) > 0.01:
            print(f"\n  Total = {total:.2f}% — must equal 100. Please re-enter.\n")
            continue

        fractions = [p / 100 for p in fractions_pct]
        print(f"\n  Total = 100% OK")
        return fractions

def section(title):
    print("\n" + "="*60)
    print(f"  {title}")
    print("="*60)

def main():
    print("\n" + "="*60)
    print("  High-Entropy Alloy (HEA) Property Calculator")
    print("  Elements: Nb, Mo, Ta, W, V, Cr, Ti, Zr, Hf")
    print("="*60)

    T        = get_temperature()
    elements = select_elements()
    fractions = get_fractions(elements)

    # ── Calculations ──────────────────────────────────────────
    delta_r   = calc_delta_r(elements, fractions)
    dS        = calc_delta_S_mix(fractions)
    dH        = calc_delta_H_mix(elements, fractions)
    VEC       = calc_VEC(elements, fractions)
    dH_el     = calc_delta_H_elastic(elements, fractions, T)
    gamma     = calc_gamma(elements, fractions)
    Tm        = calc_Tm_mean(elements, fractions)   # internal use only
    Omega     = calc_Omega(Tm, dS, dH)
    delta_chi = calc_delta_chi(elements, fractions)  # internal use only
    ea        = calc_ea(elements, fractions)         # internal use only
    dG        = calc_delta_G(dH, T, dS)
    dG_MIDMA  = calc_delta_G_MIDMA(dH, dH_el, T, dS)

    # DFT energy for 3 structures
    E_BCC = calc_E_DFT(elements, fractions, "BCC")
    E_FCC = calc_E_DFT(elements, fractions, "FCC")
    E_HCP = calc_E_DFT(elements, fractions, "HCP")
    stable = min([("BCC", E_BCC), ("FCC", E_FCC), ("HCP", E_HCP)], key=lambda x: x[1])

    # ── Output ────────────────────────────────────────────────
    section("Calculated Parameters")
    print(f"   delta_r  (Atomic size difference)     : {delta_r:.4f} %")
    print(f"   dS_mix   (Mixing entropy)             : {dS:.4f} J/(mol·K)")
    print(f"   dH_mix   (Mixing enthalpy)            : {dH:.4f} kJ/mol")
    print(f"   VEC      (Avg valence electron conc.) : {VEC:.4f}")
    print(f"   E_DFT    BCC                          : {E_BCC:.4f} eV/atom")
    print(f"   E_DFT    FCC                          : {E_FCC:.4f} eV/atom")
    print(f"   E_DFT    HCP                          : {E_HCP:.4f} eV/atom")
    print(f"   Stable structure (lowest E_DFT)       : {stable[0]} ({stable[1]:.4f} eV/atom)")
    print(f"   dH_el    (Elastic enthalpy)           : {dH_el:.4f} kJ/mol")
    print(f"   gamma    (Gamma parameter)            : {gamma:.4f}")
    print(f"   Omega    (Omega parameter)            : {Omega:.4f}")

    section("Gibbs Free Energy")
    print(f"   dG       = dH_mix - T*dS_mix          : {dG:.2f} J/mol")
    print(f"   dG_MIDMA = (dH_mix + dH_el) - T*dS   : {dG_MIDMA:.2f} J/mol")

    section("Solid Solution (SS) Criteria")
    ss_results = check_SS(dH, delta_r, Omega, gamma, delta_chi, VEC, ea, T, Tm)
    ss_yes = 0
    for model, res in ss_results.items():
        print(f"   {model:<42}: {res}")
        if res.startswith("YES"):
            ss_yes += 1
    print(f"\n   SS confirmed by {ss_yes} out of {len(ss_results)} models")

    section("BCC Structure Criteria")
    bcc_results = check_BCC(VEC, ea, T, Tm)
    bcc_yes = 0
    for model, res in bcc_results.items():
        print(f"   {model:<42}: {res}")
        if res.startswith("YES"):
            bcc_yes += 1
    print(f"\n   BCC confirmed by {bcc_yes} out of {len(bcc_results)} models")

    print("\n" + "="*60)
    print("  Calculation complete.")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
