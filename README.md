# HEA Thermodynamic Design & Phase Prediction Suite

A semi-empirical thermodynamic calculation and inverse-design toolkit for **refractory high-entropy alloys (RHEAs)**, based on the **Melnick–Soolshenko MIDMA model** and the **HEAPS parameter/criteria framework**. The package predicts solid-solution formation, BCC lattice stability, and approximate DFT-based phase stability for alloys based on:

```text
Nb · Mo · Ta · W · V · Cr · Ti · Zr · Hf
```

The suite contains three coordinated Python programs:

| File | Role | Problem type | Interface |
|---|---|---|---|
| `hea_direct_cli.py` | Direct property and criteria calculator for a user-defined composition | Forward problem | Terminal / CLI |
| `hea_streamlit.py` | Interactive GUI version of the direct calculator | Forward problem | Streamlit |
| `hea_inverse.py` | Monte Carlo inverse search for compositions satisfying all SS + BCC criteria | Inverse problem | Streamlit |

This project is designed for **thermodynamic screening and alloy design**, not as a replacement for CALPHAD, full DFT, or experimental validation.

---

## Contents

1. [Theoretical Background](#1-theoretical-background)
2. [Element Database and Binary Mixing Enthalpy](#2-element-database-and-binary-mixing-enthalpy)
3. [Calculated Parameters](#3-calculated-parameters)
4. [Solid-Solution Formation Criteria](#4-solid-solution-formation-criteria)
5. [BCC Structure Criteria](#5-bcc-structure-criteria)
6. [DFT-Based Lattice Stability](#6-dft-based-lattice-stability)
7. [MIDMA Thermodynamic Model](#7-midma-thermodynamic-model)
8. [Program Structure](#8-program-structure)
9. [`hea_direct_cli.py`](#9-heap_direct_clipy)
10. [`hea_streamlit.py`](#10-hea_streamlitpy)
11. [`hea_inverse.py`](#11-hea_inversepy)
12. [Installation and Execution](#12-installation-and-execution)
13. [Inputs and Outputs](#13-inputs-and-outputs)
14. [Known Limitations and Items to Verify](#14-known-limitations-and-items-to-verify)
15. [Recommended Next Steps](#15-recommended-next-steps)
16. [References](#16-references)

---

# 1. Theoretical Background

## 1.1 High-entropy alloys and refractory HEAs

High-entropy alloys are usually defined as alloys containing at least five principal elements, each with an atomic fraction between 5 and 35 at.%. The high configurational entropy of mixing may stabilize disordered solid solutions over intermetallic compounds, but entropy alone is not sufficient for phase prediction.

Refractory high-entropy alloys are based on elements such as:

```text
W, Ta, Mo, Nb, V, Cr, Ti, Zr, Hf
```

These alloys are attractive for high-temperature structural applications because of their high melting points, high strength, and good creep resistance. However, they also have a strong tendency to form topologically close-packed phases such as Laves, sigma, and mu phases, which must be controlled during alloy design.

## 1.2 Thermodynamic parameters governing phase formation

The formation of solid solutions in HEAs is governed by several competing factors:

| Parameter | Symbol | Physical role |
|---|---|---|
| Atomic size mismatch | $\delta r$ | Lattice distortion and elastic strain energy |
| Mixing entropy | $\Delta S_{mix}$ | Driving force for random solid solution |
| Mixing enthalpy | $\Delta H_{mix}$ | Tendency toward intermetallics or segregation |
| Valence electron concentration | $VEC$ | Prediction of BCC vs FCC stability |
| Electronegativity mismatch | $\delta \chi$ | Bonding tendency and intermetallic formation |
| Omega parameter | $\Omega$ | Competition between entropy and enthalpy |
| Gamma parameter | $\gamma$ | Lattice distortion from largest/smallest atoms |
| Gibbs free energy | $\Delta G$ | Thermodynamic stability criterion |

## 1.3 Design philosophy

The present package follows the philosophy used in the HEAPS software: no single criterion is sufficiently accurate for phase prediction. Instead, several semi-empirical parameters and criteria are evaluated **simultaneously**, and an alloy is considered promising only when multiple criteria agree.

This is consistent with the strategy proposed by Martin et al. in the HEAPS tool, where the simultaneous use of several criteria improves the reliability of phase prediction.

---

# 2. Element Database and Binary Mixing Enthalpy

## 2.1 Element database

The package includes a built-in database of 9 refractory and related elements:

```text
Nb, Mo, Ta, W, V, Cr, Ti, Zr, Hf
```

For each element, the following properties are stored:

| Property | Code key | Unit | Description |
|---|---|---|---|
| Atomic radius | `r` | pm | Atomic radius |
| Valence electron concentration | `VEC` | — | Number of valence electrons |
| Electron-to-atom ratio | `ea` | — | e/a |
| Electronegativity | `chi` | — | Pauling scale |
| Melting temperature | `Tm` | K | Melting point |
| Bulk modulus | `B` | GPa | Used for elastic enthalpy |
| Linear thermal expansion coefficient | `alpha` | K⁻¹ | Used to correct atomic volume at T |
| DFT energy of BCC | `E_BCC` | eV/atom | Ab initio total energy |
| DFT energy of FCC | `E_FCC` | eV/atom | Ab initio total energy |
| DFT energy of HCP | `E_HCP` | eV/atom | Ab initio total energy |

Example entry:

```python
"Nb": {
    "r": 146,
    "VEC": 5,
    "ea": 1.51,
    "chi": 1.6,
    "Tm": 2477 + 273.15,
    "B": 170.3,
    "alpha": 7.1e-6,
    "E_BCC": -10.0466,
    "E_FCC": -9.7232,
    "E_HCP": -9.7551
}
```

## 2.2 Binary mixing enthalpy

The binary mixing enthalpies $\Delta H_{ij}$ are stored in `DELTA_H_IJ` in kJ/mol. These values are based on the semi-empirical Miedema model as tabulated by Takeuchi and Inoue.

Negative values indicate attractive interaction between the pair, while positive values indicate repulsion.

Example:

```python
("Nb", "W"): -8
("W", "Nb"): -8
("Cr", "Zr"): -12
("Ti", "Zr"): 0
("Nb", "Ti"): 2
```

---

# 3. Calculated Parameters

## 3.1 Average atomic radius

The average atomic radius is calculated by the rule of mixtures:

$$
\bar{r} = \sum_{i=1}^{n} c_i r_i
$$

where $c_i$ is the atomic fraction of element $i$ and $r_i$ is its atomic radius.

## 3.2 Atomic size mismatch $\delta r$

The atomic size mismatch is one of the most important parameters for predicting solid-solution stability:

$$
\delta r = \sqrt{\sum_{i=1}^{n} c_i \left(1 - \frac{r_i}{\bar{r}}\right)^2} \times 100
$$

In the code:

```python
def calc_delta_r(elements, fractions):
    r_bar = calc_r_bar(elements, fractions)
    return math.sqrt(
        sum(
            fractions[i] *
            (1 - ELEMENTS_DATA[elements[i]]["r"] / r_bar) ** 2
            for i in range(len(elements))
        )
    ) * 100
```

According to Melnick and Soolshenko, solid-solution formation is usually favored when:

$$
\delta < 6\%
$$

## 3.3 Mixing entropy $\Delta S_{mix}$

The ideal configurational mixing entropy is calculated using the Boltzmann expression:

$$
\Delta S_{mix} = -R \sum_{i=1}^{n} c_i \ln c_i
$$

where $R = 8.314 \, J/(mol \cdot K)$.

In the code:

```python
def calc_delta_S_mix(fractions):
    return -R * sum(c * math.log(c) for c in fractions if c > 0)
```

## 3.4 Mixing enthalpy $\Delta H_{mix}$

The mixing enthalpy is calculated from the binary mixing enthalpies:

$$
\Delta H_{mix} = \sum_{i,j=1, i \neq j}^{n} 4 \Delta H_{ij}^{mix} c_i c_j
$$

In the code, because each pair is counted twice, the sum is divided by 2:

```python
def calc_delta_H_mix(elements, fractions):
    total = 0
    n = len(elements)
    for i in range(n):
        for j in range(n):
            if i != j:
                hij = DELTA_H_IJ.get((elements[i], elements[j]), 0)
                total += 4 * hij * fractions[i] * fractions[j]
    return total / 2
```

According to Melnick and Soolshenko, solid-solution formation is usually favored when:

$$
\Delta H_{mix} > -10 \, kJ/mol
$$

## 3.5 Valence electron concentration VEC

The average valence electron concentration is:

$$
VEC = \sum_{i=1}^{n} c_i VEC_i
$$

According to Guo et al., the VEC parameter can be used to predict the lattice structure of the solid solution:

```text
VEC < 6.87  → BCC
VEC ≥ 8     → FCC
6.87 ≤ VEC < 8 → mixed FCC + BCC
```

## 3.6 Electronegativity mismatch $\delta \chi$

The electronegativity mismatch is calculated in the same form as $\delta r$:

$$
\delta \chi = \sqrt{\sum_{i=1}^{n} c_i \left(1 - \frac{\chi_i}{\bar{\chi}}\right)^2} \times 100
$$

In the code:

```python
def calc_delta_chi(elements, fractions):
    chi_vals = [ELEMENTS_DATA[el]["chi"] for el in elements]
    chi_bar = sum(fractions[i] * chi_vals[i] for i in range(len(elements)))
    return math.sqrt(
        sum(
            fractions[i] * (1 - chi_vals[i] / chi_bar) ** 2
            for i in range(len(elements))
        )
    ) * 100
```

## 3.7 Average melting temperature $T_m$

The average melting temperature is calculated by the rule of mixtures:

$$
T_m = \sum_{i=1}^{n} c_i T_{m,i}
$$

This is used internally for the $\Omega$ parameter and for temperature-dependent criteria.

## 3.8 Omega parameter $\Omega$

The $\Omega$ parameter was proposed by Yang and Zhang to quantify the competition between entropy and enthalpy:

$$
\Omega = \frac{T_m \Delta S_{mix}}{|\Delta H_{mix}|}
$$

In the code:

```python
def calc_Omega(Tm, dS, dH):
    if abs(dH) < 1e-10:
        return float('inf')
    return Tm * dS / (abs(dH) * 1000)
```

According to Yang and Zhang:

```text
Ω ≥ 1.1 and δr ≤ 6.6 → solid solution likely
```

## 3.9 Gamma parameter $\gamma$

The $\gamma$ parameter was proposed by Wang et al. to quantify lattice distortion by comparing the largest and smallest atoms with the average atomic size:

$$
\gamma =
\frac{
1 - \left(\frac{r_{min}}{r_{min} + \bar{r}}\right)^2
}{
1 - \left(\frac{r_{max}}{r_{max} + \bar{r}}\right)^2
}
$$

In the code:

```python
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
```

According to Wang et al.:

```text
γ ≤ 1.175 → solid solution likely
```

## 3.10 Gibbs free energy

The simple Gibbs free energy is:

$$
\Delta G = \Delta H_{mix} - T \Delta S_{mix}
$$

In the code:

```python
def calc_delta_G(dH, T, dS):
    return dH * 1000 - T * dS
```

The unit of the output is J/mol.

---

# 4. Solid-Solution Formation Criteria

The package implements several criteria for predicting solid-solution formation. These are based on the criteria summarized in the HEAPS software.

## 4.1 MC1 — Zhang 2008

This criterion uses $\delta r$ and $\Delta H_{mix}$.

Solid solution is predicted when:

$$
0.5 < \delta r < 6.5
$$

and

$$
-17.5 < \Delta H_{mix} < 5
$$

In the code:

```python
c1 = (0.5 < delta_r < 6.5) and (-17.5 < dH < 5)
```

## 4.2 MC2 — Yang 2012

This criterion uses $\Omega$ and $\delta r$.

Solid solution is predicted when:

$$
\Omega \ge 1.1
$$

and

$$
\delta r \le 6.6
$$

In the code:

```python
c2 = (Omega >= 1.1) and (delta_r <= 6.6)
```

## 4.3 MC3 — Guo 2013

This criterion also uses $\Delta H_{mix}$ and $\delta r$.

Solid solution is predicted when:

$$
-11.6 < \Delta H_{mix} < 3.2
$$

and

$$
\delta r \le 6.6
$$

In the code:

```python
c3 = (-11.6 < dH < 3.2) and (delta_r <= 6.6)
```

## 4.4 MC4 — Wang 2014 Gamma

This criterion uses only the $\gamma$ parameter.

Solid solution is predicted when:

$$
\gamma \le 1.175
$$

In the code:

```python
results["MC4 Wang (2014) [gamma]"] = gamma <= 1.175
```

## 4.5 MC5 — Poletti 2014

This criterion uses $\delta r$ and $\delta \chi$.

Solid solution is predicted when:

$$
1 < \delta r < 6
$$

and

$$
3 < \delta \chi < 6
$$

In the code:

```python
c5 = (1 < delta_r < 6) and (3 < delta_chi < 6)
```

## 4.6 MC7 — Wang 2014 with temperature dependence

This criterion considers the effect of evaluation temperature through the ratio $T/T_m$.

If:

$$
T/T_m > 0.9
$$

then solid solution is predicted when:

$$
-15 \le \Delta H_{mix} \le 5
$$

and

$$
\delta r \le 6.6
$$

If:

$$
0.5 \le T/T_m < 0.9
$$

then solid solution is predicted when:

$$
\Delta H_{mix} \ge -7.5
$$

and

$$
\delta r \le 3.3
$$

In the code:

```python
ratio = T / Tm

if ratio > 0.9:
    c7 = (-15 <= dH <= 5) and (delta_r <= 6.6)
elif 0.5 <= ratio < 0.9:
    c7 = (dH >= -7.5) and (delta_r <= 3.3)
else:
    c7 = False
```

## 4.7 Summary of solid-solution criteria

| Code | Reference | Parameters | Solid-solution condition |
|---|---|---|---|
| MC1 | Zhang 2008 | $\delta r, \Delta H_{mix}$ | $0.5 < \delta r < 6.5$ and $-17.5 < \Delta H_{mix} < 5$ |
| MC2 | Yang 2012 | $\Omega, \delta r$ | $\Omega \ge 1.1$ and $\delta r \le 6.6$ |
| MC3 | Guo 2013 | $\Delta H_{mix}, \delta r$ | $-11.6 < \Delta H_{mix} < 3.2$ and $\delta r \le 6.6$ |
| MC4 | Wang 2014 | $\gamma$ | $\gamma \le 1.175$ |
| MC5 | Poletti 2014 | $\delta r, \delta \chi$ | $1 < \delta r < 6$ and $3 < \delta \chi < 6$ |
| MC7 | Wang 2014 | $T/T_m, \Delta H_{mix}, \delta r$ | Depends on $T/T_m$ |

---

# 5. BCC Structure Criteria

These criteria are used to predict whether the solid solution will have a BCC lattice structure. For refractory HEAs, BCC is usually the desired structure.

## 5.1 LSS1 — Guo 2011

This criterion is based only on VEC.

```text
If VEC < 6.87 → BCC
If VEC ≥ 8 → FCC
If 6.87 ≤ VEC < 8 → mixed FCC + BCC
```

In the code:

```python
if VEC < 6.87:
    results["LSS1 Guo (2011)"] = True
elif VEC >= 8:
    results["LSS1 Guo (2011)"] = False
else:
    results["LSS1 Guo (2011)"] = None
```

The value `None` represents a mixed or unresolved case.

## 5.2 LSS2 — Poletti 2014

This criterion uses both VEC and e/a.

BCC is predicted when:

$$
VEC < 7.5
$$

and

$$
1.8 < e/a < 2.3
$$

FCC is predicted when:

$$
VEC > 7.5
$$

and

$$
1.6 < e/a < 1.8
$$

Otherwise the result is unresolved.

In the code:

```python
if VEC < 7.5 and (1.8 < ea < 2.3):
    results["LSS2 Poletti (2014)"] = True
elif VEC > 7.5 and (1.6 < ea < 1.8):
    results["LSS2 Poletti (2014)"] = False
else:
    results["LSS2 Poletti (2014)"] = None
```

## 5.3 LSS3 — Wang 2014 with temperature dependence

This criterion also considers the ratio $T/T_m$.

If:

$$
T/T_m > 0.9
$$

then:

```text
VEC < 6.87 → BCC
VEC > 7.84 → FCC
otherwise → mixed
```

If:

$$
0.5 < T/T_m < 0.9
$$

then:

```text
VEC < 6 → BCC
VEC > 7.8 → FCC
otherwise → mixed
```

## 5.4 Summary of BCC criteria

| Code | Reference | Parameters | BCC condition |
|---|---|---|---|
| LSS1 | Guo 2011 | VEC | $VEC < 6.87$ |
| LSS2 | Poletti 2014 | VEC, e/a | $VEC < 7.5$ and $1.8 < e/a < 2.3$ |
| LSS3 | Wang 2014 | VEC, $T/T_m$ | Depends on temperature |

---

# 6. DFT-Based Lattice Stability

## 6.1 Concept

The relative stability of BCC, FCC, and HCP structures for each pure element can be obtained from ab initio calculations. In the paper by Wang et al., total energies for 78 elements were calculated for BCC, FCC, and HCP structures and compared with CALPHAD/SGTE lattice stabilities.

In this package, the DFT energies of the pure elements for BCC, FCC, and HCP are stored in the element database.

## 6.2 Approximate alloy DFT energy

The approximate DFT energy of the alloy for a given structure is calculated as a composition-weighted average:

$$
E_{structure}^{alloy} = \sum_{i=1}^{n} c_i E_{structure}^{i}
$$

In the code:

```python
def calc_E_DFT(elements, fractions, structure="BCC"):
    key = f"E_{structure}"
    return sum(
        fractions[i] * ELEMENTS_DATA[elements[i]][key]
        for i in range(len(elements))
    )
```

## 6.3 Stable structure prediction

Three energies are calculated:

```text
E_DFT BCC
E_DFT FCC
E_DFT HCP
```

The structure with the lowest energy is selected as the most stable:

```python
stable = min(
    [("BCC", E_BCC), ("FCC", E_FCC), ("HCP", E_HCP)],
    key=lambda x: x[1]
)
```

## 6.4 Limitation of the weighted-average approach

The weighted-average DFT energy is a simple approximation. It does not include:

- atomic ordering effects
- explicit binary and higher-order interactions
- lattice distortion
- temperature effects
- configurational entropy beyond the ideal mixing term

Therefore this part is suitable for preliminary screening only, not for final phase prediction.

Wang et al. also showed that for some transition elements, there are significant discrepancies between ab initio lattice stability and CALPHAD/SGTE data, especially for unstable structures.

---

# 7. MIDMA Thermodynamic Model

## 7.1 Basis of the model

The thermodynamic model used in this package is based on the paper by Melnick and Soolshenko:

```text
Thermodynamic design of high-entropy refractory alloys
Journal of Alloys and Compounds 694, 223–227, 2017
```

In this model, the stability of a multicomponent substitutional solid solution is evaluated by considering three main contributions:

1. Mixing enthalpy $\Delta H_{mix}$
2. Elastic enthalpy $\Delta H_{el}$
3. Mixing entropy $\Delta S_{mix}$

The Gibbs free energy change is:

$$
\Delta G =
\Delta H_{mix} +
\Delta H_{el} -
T \Delta S_{mix}
$$

The composition with the lowest $\Delta G$ is considered the most thermodynamically stable.

## 7.2 Mixing enthalpy in regular solution approximation

The mixing enthalpy is written as:

$$
\Delta H_{mix} =
\sum_{i,j=1, i \neq j}^{n}
c_i c_j \Omega_{ij}
$$

where:

$$
\Omega_{ij} = 4 \Delta H_{ij}^{mix}
$$

The values of $\Delta H_{ij}^{mix}$ are taken from semi-empirical Miedema-based tables.

## 7.3 Atomic volume at temperature T

The atomic volume of each element at the reference temperature is approximated as a sphere:

$$
V_{0i} = \frac{4}{3} \pi r_i^3
$$

Then it is corrected to temperature T using the linear thermal expansion coefficient:

$$
V_i(T) = V_{0i} \left(1 + \alpha_i (T - 298)\right)
$$

In the code:

```python
def v0(el):
    return (4/3) * math.pi * (data[el]["r"] * 1e-12) ** 3

def v_T(el):
    return v0(el) * (1 + data[el]["alpha"] * (T - 298))
```

## 7.4 Average alloy volume

The average atomic volume of the alloy is obtained from mechanical equilibrium:

$$
\bar{V}(T) =
\frac{
\sum_{i=1}^{n} c_i B_i V_i(T)
}{
\sum_{i=1}^{n} c_i B_i
}
$$

In the code:

```python
v_bar = sum(c[i] * B_vals[i] * v_vals[i] for i in range(len(elements))) / \
        sum(c[i] * B_vals[i] for i in range(len(elements)))
```

## 7.5 Elastic enthalpy

The elastic distortion energy is:

$$
\Delta H_{el} =
\sum_{i=1}^{n}
c_i B_i
\frac{
\left(V_i(T) - \bar{V}(T)\right)^2
}{
2 V_i(T)
}
$$

In the code:

```python
total = sum(
    c[i] * B_vals[i] * (v_vals[i] - v_bar) ** 2 / (2 * v_vals[i])
    for i in range(len(elements))
)
```

In the terminal version, the unit conversion is:

```python
N_A = 6.02214076e23
return total * N_A / 1000
```

The output unit is kJ/mol.

## 7.6 Gibbs free energy with elastic term

The MIDMA Gibbs free energy is:

$$
\Delta G_{MIDMA} =
\Delta H_{mix} +
\Delta H_{el} -
T \Delta S_{mix}
$$

In the code:

```python
def calc_delta_G_MIDMA(dH_mix, dH_el, T, dS):
    return (dH_mix + dH_el) * 1000 - T * dS
```

The output unit is J/mol.

## 7.7 Monte Carlo design philosophy

In the Melnick and Soolshenko paper, the stable composition is found by minimizing $\Delta G$ using a Monte Carlo procedure:

1. Start from an initial composition.
2. Calculate $\Delta G$.
3. Randomly perturb the composition with small steps.
4. Calculate the new $\Delta G$.
5. If $\Delta G$ decreases, accept the new composition.
6. Otherwise keep the old composition.
7. Repeat until a local or global minimum is found.

The `hea_inverse.py` code also uses Monte Carlo, but instead of directly minimizing $\Delta G$, it searches for compositions that satisfy the full set of SS and BCC criteria.

---

# 8. Program Structure

The package contains three main programs, each covering a different level of usage.

```text
hea_direct_cli.py
hea_streamlit.py
hea_inverse.py
```

## 8.1 Relationship between the programs

```text
                    +----------------------+
                    |   Element Database   |
                    |   DELTA_H_IJ         |
                    +----------+-----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
+---------v--------+  +--------v---------+  +-------v----------+
| hea_direct_cli   |  | hea_streamlit    |  | hea_inverse      |
| Forward terminal |  | Forward GUI      |  | Reverse search   |
+------------------+  +------------------+  +------------------+
```

| File | Input | Processing | Output |
|---|---|---|---|
| `hea_direct_cli.py` | Temperature, elements, fractions | Parameter and criteria calculation | Text report |
| `hea_streamlit.py` | Temperature, elements, fractions | Parameter and criteria calculation | GUI dashboard |
| `hea_inverse.py` | Temperature, elements, composition range, sample count | Monte Carlo search | Valid compositions |

---

# 9. `hea_direct_cli.py`

## 9.1 Renamed file

The original name of this file was:

```text
english.py
```

It has been renamed to:

```text
hea_direct_cli.py
```

Reason:

- `hea` stands for High-Entropy Alloy.
- `direct` indicates that this code solves the forward problem: the user gives a composition and the code calculates the parameters and criteria.
- `cli` indicates that it is a command-line interface.
- The name is consistent with `hea_inverse.py`, which solves the inverse problem.

## 9.2 Role

This code is a terminal-based calculator for analyzing a user-defined composition. The user:

1. Enters the temperature.
2. Selects exactly 5 elements.
3. Enters the molar fraction of each element in percent.
4. The program calculates the parameters and evaluates the criteria.

## 9.3 Main functions

| Function | Purpose |
|---|---|
| `calc_r_bar` | Average atomic radius |
| `calc_delta_r` | Atomic size mismatch |
| `calc_delta_S_mix` | Mixing entropy |
| `calc_delta_H_mix` | Mixing enthalpy |
| `calc_VEC` | Valence electron concentration |
| `calc_E_DFT` | DFT energy for BCC/FCC/HCP |
| `calc_delta_H_elastic` | Elastic enthalpy |
| `calc_gamma` | Gamma parameter |
| `calc_Tm_mean` | Average melting temperature |
| `calc_Omega` | Omega parameter |
| `calc_delta_chi` | Electronegativity mismatch |
| `calc_ea` | Average e/a |
| `calc_delta_G` | Simple Gibbs free energy |
| `calc_delta_G_MIDMA` | Gibbs free energy with elastic term |
| `check_SS` | Solid-solution criteria |
| `check_BCC` | BCC structure criteria |

## 9.4 Terminal user interface

The code has four main input stages:

### 9.4.1 Temperature

```text
Enter temperature (Kelvin)
[Default: 1000 K — press Enter to confirm]
```

If the user presses Enter, the default value is:

```text
T = 1000 K
```

### 9.4.2 Element selection

The user must select exactly 5 elements.

Example:

```text
Nb Mo Ta W V
```

or by number:

```text
1 2 3 4 5
```

The code checks that:

- exactly 5 elements are selected
- no duplicate elements are selected
- all element names are valid

### 9.4.3 Molar fractions

For each element, the user enters a molar percentage between 5 and 35.

Example:

```text
Nb (%): 20
Mo (%): 20
Ta (%): 20
W (%): 20
V (%): 20
```

The sum must be 100.

## 9.5 Output

The output contains:

### 9.5.1 Calculated parameters

```text
delta_r
dS_mix
dH_mix
VEC
E_DFT BCC
E_DFT FCC
E_DFT HCP
Stable structure
dH_el
gamma
Omega
```

### 9.5.2 Gibbs free energy

```text
dG       = dH_mix - T*dS_mix
dG_MIDMA = (dH_mix + dH_el) - T*dS
```

### 9.5.3 Solid-solution criteria

```text
MC1 Zhang (2008)
MC2 Yang (2012)
MC3 Guo (2013)
MC4 Wang (2014) [gamma]
MC5 Poletti (2014)
MC7 Wang (2014)
```

### 9.5.4 BCC criteria

```text
LSS1 Guo (2011)
LSS2 Poletti (2014)
LSS3 Wang (2014)
```

---

# 10. `hea_streamlit.py`

## 10.1 Role

This code implements the same calculation logic as the terminal code, but with an interactive Streamlit GUI.

The user can:

- select the temperature with a slider
- select 5 elements from a multiselect menu
- adjust the molar fraction of each element with sliders
- see the results as cards, tables, badges, and colored indicators

## 10.2 Inputs

| Input | Type | Range |
|---|---|---|
| Temperature | Slider | 300 to 3000 K |
| Elements | Multi-select | Select 5 elements |
| Molar fraction per element | Slider | 5 to 35 % |

The sum of the molar fractions must be 100.

## 10.3 Output sections

### 10.3.1 Alloy composition

The program automatically builds the alloy name.

Example:

```text
Nb20Mo20Ta20W20V20
```

### 10.3.2 Thermodynamic parameters

The following parameters are shown as cards:

```text
δr
ΔS_mix
ΔH_mix
VEC
ΔH_el
γ
Ω
```

### 10.3.3 DFT energies

Three structures are compared:

```text
BCC
FCC
HCP
```

The structure with the lowest energy is marked with:

```text
★ STABLE
```

### 10.3.4 Gibbs free energy

Two quantities are displayed:

```text
ΔG = ΔH_mix − T·ΔS_mix
ΔG_MIDMA = (ΔH_mix + ΔH_el) − T·ΔS_mix
```

If the value is negative, it is thermodynamically favorable.

### 10.3.5 Solid-solution criteria

The results of each criterion are displayed as green, red, or yellow badges:

```text
✓ SS
✗ IM/BMG
~ OUT OF RANGE
```

### 10.3.6 BCC criteria

The results of LSS1, LSS2, and LSS3 are displayed as:

```text
✓ BCC
✗ FCC
~ Mixed
```

## 10.4 Difference from the terminal code

| Feature | `hea_direct_cli.py` | `hea_streamlit.py` |
|---|---|---|
| Interface | Terminal | Browser |
| Input | Manual typing | Sliders and menus |
| Output | Plain text | Cards, tables, colors |
| Extra libraries | None | Streamlit, pandas |
| Best use | Quick calculation | Presentation and visual analysis |

---

# 11. `hea_inverse.py`

## 11.1 Role

This code is an inverse-design tool. Instead of the user giving a composition and checking whether it is good, the program generates thousands of random compositions and finds those that satisfy the desired criteria.

This code is philosophically similar to the Monte Carlo approach used by Melnick and Soolshenko, but its objective is to find compositions satisfying a set of SS and BCC criteria rather than directly minimizing $\Delta G$.

## 11.2 Inputs

| Input | Type | Description |
|---|---|---|
| Temperature | Slider | Evaluation temperature |
| Elements | Multi-select | Select exactly 5 elements |
| Minimum composition | Slider | Minimum percent per element |
| Maximum composition | Slider | Maximum percent per element |
| Monte Carlo samples | Select slider | 10,000 to 200,000 |

## 11.3 Random composition generation

The main function is:

```python
def random_composition(n=5, lo=5, hi=35):
    while True:
        vals = [random.uniform(lo, hi) for _ in range(n)]
        s = sum(vals)
        normed = [v / s * 100 for v in vals]
        if all(lo <= v <= hi for v in normed):
            return normed
```

This function:

1. Generates five random numbers between the minimum and maximum.
2. Normalizes them so that their sum is 100.
3. Checks that all values are still within the allowed range.
4. If not, it repeats the sampling.

## 11.4 Criteria evaluation

For each random composition, the following function is called:

```python
res = get_all_criteria(elements, fracs, T)
```

This function returns one of the following for each criterion:

| Value | Meaning |
|---|---|
| `True` | Criterion passed |
| `False` | Criterion failed |
| `None` | Mixed or unresolved case |

## 11.5 Result storage

The results are stored in a dictionary called `buckets`:

```python
buckets = {
    "MC1 Zhang (2008)": [],
    "MC2 Yang (2012)": [],
    "MC3 Guo (2013)": [],
    "MC4 Wang (2014) γ": [],
    "MC5 Poletti (2014)": [],
    "MC7 Wang (2014)": [],
    "LSS1 Guo (2011)": [],
    "LSS2 Poletti (2014)": [],
    "LSS3 Wang (2014)": [],
    "ALL": []
}
```

If a composition satisfies all criteria, or only has `True` and `None` results, it is stored in the `"ALL"` list.

## 11.6 Outputs

### 11.6.1 Number of valid compositions per criterion

For each criterion, the number of valid compositions is displayed.

Example:

```text
MC1 Zhang (2008): 35,000 compositions
MC2 Yang (2012): 42,000 compositions
LSS1 Guo (2011): 39,000 compositions
ALL: 18,000 compositions
```

### 11.6.2 Composition statistics

For each element, the following statistics are calculated:

```text
Min (%)
Max (%)
Mean (%)
Std (%)
```

### 11.6.3 Composition range bars

For each element, a bar is displayed showing the range of that element in the valid compositions.

### 11.6.4 Ten sample compositions

In each tab, the first ten compositions are shown as examples.

## 11.7 Difference from the other two codes

| Feature | `hea_direct_cli.py` | `hea_streamlit.py` | `hea_inverse.py` |
|---|---|---|---|
| Problem type | Forward | Forward | Inverse |
| Input | Fixed composition | Fixed composition | Composition range |
| Method | Single calculation | Single calculation | Monte Carlo |
| Output | Results for one alloy | Results for one alloy | Set of valid compositions |
| Use | Quick analysis | Presentation and visual check | Design and screening |

---

# 12. Installation and Execution

## 12.1 Requirements

### For the terminal code

```bash
python
```

No external libraries are needed besides the standard `math` module.

### For the Streamlit codes

```bash
python
streamlit
pandas
numpy
```

Install the libraries:

```bash
pip install streamlit pandas numpy
```

or:

```bash
python -m pip install streamlit pandas numpy
```

## 12.2 Running the terminal code

After renaming the file to `hea_direct_cli.py`:

```bash
python hea_direct_cli.py
```

or on systems using Python 3:

```bash
python3 hea_direct_cli.py
```

## 12.3 Running the direct GUI code

```bash
streamlit run hea_streamlit.py
```

Then open the browser at:

```text
http://localhost:8501
```

If the `streamlit` command is not found:

```bash
python -m streamlit run hea_streamlit.py
```

## 12.4 Running the inverse design code

```bash
streamlit run hea_inverse.py
```

Then open the browser at:

```text
http://localhost:8501
```

If the `streamlit` command is not found:

```bash
python -m streamlit run hea_inverse.py
```

---

# 13. Inputs and Outputs

## 13.1 Common inputs

| Input | Unit | Range | Description |
|---|---|---|---|
| Temperature | K | Usually 300 to 3000 | Evaluation or processing temperature |
| Elements | — | Select 5 elements | Principal alloy elements |
| Molar fraction | % | 5 to 35 per element | Must sum to 100 |

## 13.2 Output parameters

| Parameter | Unit | Description |
|---|---|---|
| $\delta r$ | % | Atomic size mismatch |
| $\Delta S_{mix}$ | J/(mol·K) | Mixing entropy |
| $\Delta H_{mix}$ | kJ/mol | Mixing enthalpy |
| VEC | — | Valence electron concentration |
| e/a | — | Electron-to-atom ratio |
| $\Delta H_{el}$ | kJ/mol | Elastic enthalpy |
| $\gamma$ | — | Gamma parameter |
| $\Omega$ | — | Omega parameter |
| $\Delta G$ | J/mol | Simple Gibbs free energy |
| $\Delta G_{MIDMA}$ | J/mol | Gibbs free energy with elastic term |
| $E_{DFT}^{BCC}$ | eV/atom | DFT energy of BCC |
| $E_{DFT}^{FCC}$ | eV/atom | DFT energy of FCC |
| $E_{DFT}^{HCP}$ | eV/atom | DFT energy of HCP |
| Stable structure | BCC/FCC/HCP | Structure with lowest DFT energy |

## 13.3 Criteria outputs

| Category | Outputs |
|---|---|
| Solid-solution criteria | MC1, MC2, MC3, MC4, MC5, MC7 |
| BCC criteria | LSS1, LSS2, LSS3 |
| Final output | Number of criteria confirming SS or BCC |

---

# 14. Known Limitations and Items to Verify

## 14.1 Simplification of DFT energy by weighted average

The alloy DFT energy is calculated as a composition-weighted average:

$$
E_{alloy} = \sum c_i E_i
$$

This is a simple approximation and does not include:

- atomic ordering effects
- explicit binary and higher-order interactions
- lattice distortion
- temperature effects
- disordered vs ordered phases
- volume and stress effects

Therefore this part is suitable for preliminary screening only, not for final prediction.

## 14.2 DFT data mismatch in the Streamlit version

In `hea_streamlit.py`, some `E_BCC` and `E_FCC` values differ from those in `hea_direct_cli.py` and `hea_inverse.py`.

For example:

| Element | `E_BCC` in terminal version | `E_BCC` in Streamlit version |
|---|---:|---:|
| Nb | -10.0466 | -4.1302 |
| Mo | -10.7799 | -5.2954 |
| Ta | -11.7358 | -5.7318 |
| W | -12.7781 | -11.9107 |
| V | -8.9632 | -5.2954 |

This difference may be due to a different source, data entry error, or a change of scale. To avoid inconsistent results, it is recommended to unify the DFT values across all files.

The main reference for these energies is the paper by Wang et al.

## 14.3 Difference in elastic enthalpy unit conversion

In the terminal code:

```python
return total * N_A / 1000
```

In the Streamlit and inverse codes:

```python
return total / 1.602e-19 / 6.022e23 * 1000
```

These two expressions may give different results. The physically appropriate conversion from J/atom to kJ/mol is:

$$
\Delta H_{el} \left(\frac{kJ}{mol}\right)
=
\Delta H_{el} \left(\frac{J}{atom}\right)
\times
N_A
\times
\frac{1}{1000}
$$

Therefore it is recommended to unify the unit conversion across all files.

## 14.4 Possible execution error in the terminal code

If the terminal code contains:

```python
if name == "main":
    main()
```

it should be corrected to:

```python
if __name__ == "__main__":
    main()
```

Otherwise the program may not run.

## 14.5 Accuracy of semi-empirical criteria

The criteria used are based on limited experimental data and simple approximations. Therefore:

- no single criterion is definitive
- simultaneous use of multiple criteria is recommended
- final design should be verified with CALPHAD, DFT, or experiments

## 14.6 No prediction of complex phases

The current codes do not independently predict the following phases:

```text
Laves
Sigma
B2
TCP
FCC + L12
HCP
BMG
```

However, criteria for these phases can be added from other papers.

## 14.7 Additional criteria from other sources

### BCC + Laves criteria from Zareipour et al.

For predicting Laves phase formation in refractory HEAs:

$$
-20 \le \Delta H_{mix} \le -3 \, kJ/mol
$$

and

$$
4 < \delta < 10
$$

Also, it was reported that:

$$
VEC \le 6.78
$$

is usually associated with BCC stability.

### FCC + L1₂ criteria from Tao et al.

For designing superalloys with FCC + L1₂ microstructure:

$$
VEC > 8
$$

$$
-16.0 < \Delta H_{mix} < -9.7 \, kJ/mol
$$

$$
1671 < T_m < 1822 \, K
$$

These criteria can be added as separate modules to the package.

---

# 15. Recommended Next Steps

## 15.1 Unify the database

It is recommended that all files use a common database.

Best structure:

```text
hea_database.py
hea_calculations.py
hea_criteria.py
hea_direct_cli.py
hea_streamlit.py
hea_inverse.py
```

## 15.2 Add output saving

In the Streamlit versions, it would be useful to add export options:

```text
CSV
Excel
JSON
```

## 15.3 Add more criteria

The following criteria can be added:

```text
MC6 Singh
MC8 Ye
MC9 Troparevsky
MC10 Senkov
LSS4 Ye
LSS5 Zeng
ImF1 Dong
ImF2 Tsai
ImF3 Yurchenko
FMP1 Sheikh
Laves criteria
FCC + L12 criteria
```

## 15.4 Add Gibbs free energy minimization

An optimization module can be added to find the composition with the lowest:

$$
\Delta G_{MIDMA}
$$

similar to the Monte Carlo method used by Melnick and Soolshenko.

## 15.5 Add CALPHAD or Thermo-Calc comparison

To increase accuracy, the semi-empirical results can be compared with CALPHAD predictions.

---

# 16. References

1. A.B. Melnick, V.K. Soolshenko, *Thermodynamic design of high-entropy refractory alloys*, Journal of Alloys and Compounds, 694, 223–227, 2017.

2. P. Martin, C.E. Madrid-Cortes, C. Cáceres, N. Araya, C. Aguilar, J.M. Cabrera, *HEAPS: A user-friendly tool for the design and exploration of high-entropy alloys based on semi-empirical parameters*, Computer Physics Communications, 278, 108398, 2022.

3. Y. Wang, S. Curtarolo, C. Jiang, R. Arroyave, T. Wang, G. Ceder, L.-Q. Chen, Z.-K. Liu, *Ab initio lattice stability in comparison with CALPHAD lattice stability*, Computer Coupling of Phase Diagrams and Thermochemistry, 28, 79–90, 2004.

4. F. Zareipour, H. Shahmir, Y. Huang, *Formation and significance of topologically close-packed Laves phases in refractory high-entropy alloys*, Journal of Alloys and Compounds, 986, 174148, 2024.

5. Q. Tao, X. Yang, L. Bao, Y. Zhou, T. Yang, Y. Zhao, R. Shi, Z. Yao, X. Liu, *Transforming machine learning model knowledge into material insights for multi-principal-element superalloy phase design*, npj Computational Materials, 11, 99, 2025.

6. Y. Zhang, Y.J. Zhou, J.P. Lin, G.L. Chen, P.K. Liaw, *Solid-solution phase formation rules for multi-component alloys*, Advanced Engineering Materials, 10, 534–538, 2008.

7. X. Yang, Y. Zhang, *Prediction of high-entropy stabilized solid-solution in multi-component alloys*, Materials Chemistry and Physics, 132, 233–238, 2012.

8. S. Guo, Q. Hu, C. Ng, C.T. Liu, *More than entropy in high-entropy alloys: forming solid solutions or amorphous phase*, Intermetallics, 41, 96–103, 2013.

9. Z. Wang, Y. Huang, Y. Yang, J. Wang, C.T. Liu, *Atomic-size effect and solid solubility of multicomponent alloys*, Scripta Materialia, 94, 28–31, 2014.

10. M.G. Poletti, L. Battezzati, *Electronic and thermodynamic criteria for the occurrence of high entropy alloys in metallic systems*, Acta Materialia, 75, 297–306, 2014.

11. S. Guo, C. Ng, J. Lu, C.T. Liu, *Effect of valence electron concentration on stability of fcc or bcc phase in high entropy alloys*, Journal of Applied Physics, 109, 103505, 2011.

12. A. Takeuchi, A. Inoue, *Classification of bulk metallic glasses by atomic size difference, heat of mixing and period of constituent elements*, Materials Transactions, 46, 2817–2829, 2005.