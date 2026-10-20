"""
TF-RSE-IonCF v3.2
A Unified Closed-Form Model for the Total Electron Binding Energy
of All Neutral Atoms and Ions (Z = 1-92, q = 0 to Z-1)

Reference: Cheng & Cheng, 2026 (companion paper)
License: PolyForm Noncommercial License 1.0.0
(https://polyformproject.org/licenses/noncommercial/1.0.0/).
Free for academic, educational, and other noncommercial use.
Commercial application --- including but not limited to use in
proprietary products, consulting services, or for-profit research ---
requires prior written consent from the authors.
Required Notice: Copyright 2026 Jinjun Cheng and Dian Cheng
Contact: cdcdcd999@126.com

Statement of originality: the model formulations, parameter sets and
all computational procedures were developed independently by the
authors (Jinjun Cheng, Dian Cheng). They do not rely on or incorporate
any pre-existing proprietary frameworks or third-party models, except
for the explicitly cited reference data (NIST ASD v5.12; Dzuba-
Flambaum-Afanasjev 2024) used solely for validation and comparison.

No external dependencies. Python 3.8+.

Command-line usage:
    python tfrse_ioncf_v3.2.py Z [N]
        Z   nuclear charge, integer 1..92 (required)
        N   electron number, integer 1..Z (optional; default N = Z,
            i.e. the neutral atom)
    Examples:
        python tfrse_ioncf_v3.2.py 26 13     -> Fe XIII... Fe(13+) TBE
        python tfrse_ioncf_v3.2.py 92        -> neutral U TBE
    With no arguments the built-in self-test demonstration is run.

v3.1 changes over v3.0 (all validated on the full 4278-point grid):
(1) hard physical cap E_raw(Z,N) <= E_v2(Z) for N < Z
    (positive binding energy of the last electron);
(2) anchors N = 1, 2, 10 are never modified by the safeguard,
    EXCEPT clamp-degenerate points (value pinned to E_v2), which are
    replaced by a valence-consistent estimate (accuracy improves);
(3) boundary layer (N > Z-5): ionization-energy heuristic is now
    INCREASING (IP_k = 1.8 * IP_{k-1}), matching the experimental
    trend across the valence shell;
(4) core region (N <= Z-5): minimal-touch projection (0.5 eV floor),
    no heuristic cascade; cumulative accuracy preserved to within the
    local oscillation amplitude;
(5) Z_eff <= Z hard constraint (unchanged from v3.0).

v3.2 changes over v3.1 (referee revision; core model and all 29
fitted parameters unchanged):
(6) Tier-2 boundary-layer trigger extended from pure monotonicity to
    a PHYSICAL-PLAUSIBILITY WINDOW: a boundary link k (counted inward
    from the neutral anchor, k = 1 is the first ionization energy) is
    reconstructed whenever its implied ionization energy falls outside
        0.5 eV  <=  IP_k  <=  4 * IP_base(period) * 1.8^(k-1).
    The upper bound exceeds every experimental first ionization energy
    in the periodic table (max 24.6 eV, He) by a factor >= 3 for every
    period, so no physical valence link is ever removed; it removes
    only monotone-but-unphysical links produced by the raw blend in
    the near-neutral region (e.g. U: raw IP1 = +2588 eV). The lower
    bound removes degenerate near-zero links (e.g. the Na N=10 link,
    ~1e-11 eV in v3.1). After the extension, every element has
    strictly positive, monotone, trend-correct boundary ionization
    energies (Si: 8.0 eV, Fe: 7.0 eV, U: 5.0 eV; experiment 8.15,
    7.90, 6.19 eV). NOTE: boundary-layer values are positive,
    monotone, trend-correct ESTIMATES, not quantitative valence
    ionization potentials; individual first IPs are reproduced only
    to within a factor of about 3 (worst case B: 25.8 vs 8.30 eV).
(7) Command-line interface: Z and N accepted as command-line
    arguments (see usage above).
(8) Comment corrected: E_H = 13.605693 eV is the RYDBERG energy
    (the Hartree energy is 27.211386 eV).
"""
import math
import sys

E_H = 13.605693          # Rydberg energy, eV (NOT the Hartree energy, 27.211386 eV)
ALPHA_INV = 137.035999084  # inverse fine-structure constant
ALF = 1.0 / ALPHA_INV

# TF-RSE v2 neutral parameters (frozen)
V2 = dict(a=2.0224, b=2.1414, c=-0.2404, d=0.9737, e=-0.2254)

# TF-RSE-IonCF fitted parameters (18 ionic)
P = dict(
    sK=+0.2399, sL0=+3.6555, sL1=+0.0204,
    qL=+0.7401, qM=+0.6431,
    r1=-0.8000, r2=+0.8280,
    c_sh=+0.0098, d_sh=+2.0531, e_sh=+0.0024, ep=+0.6376,
    A10=-0.1342, A18=+0.2318, A36=-0.0455, A54=+0.2260, A86=+0.4764,
    w_b=+26.1947, gb=+0.7139,
)
MAGIC = (10, 18, 36, 54, 86)


def D_dirac(Z):
    """Hydrogenic point-nucleus Dirac factor.
    Finite-nuclear-size and QED corrections are NOT included;
    the N=1 anchor deviates from experiment by +0.35% at Z = 92."""
    za = Z * ALF
    return 2.0 * (1.0 - math.sqrt(1.0 - za * za)) / (za * za)


def _S(n, c, d, e):
    """Trigonometric shell-oscillation factor."""
    nc = n ** (1.0 / 3.0)
    return 1.0 + c * math.sin(d * nc) + e * math.cos(d * nc)


def _L_light(Z):
    """Piecewise light-element correction."""
    if Z <= 8:
        return 0.9973 - 0.9925 * math.exp(-1.3240 * Z)
    if Z == 9:
        return 1.0
    if Z == 10:
        return 1.0149
    if Z < 18:
        return 1.0 + 0.01955 * math.exp(-0.06539 * (Z - 10))
    return 1.0


def tbe_neutral_v2(Z):
    """TF-RSE v2 neutral-atom total binding energy, eV."""
    R = (1.0 - (Z / ALPHA_INV) ** 2) ** (-0.5)
    S = _S(Z, V2['c'], V2['d'], V2['e'])
    return V2['a'] * (Z ** V2['b']) * E_H * R * S * _L_light(Z)


def _shell_modulation(y, Z, N):
    """Electron-number shell sawtooth modulation M(y)."""
    env = math.sin(math.pi * min(max(y, 0.0), 1.0)) ** P['ep']
    ratio = (_S(N, P['c_sh'], P['d_sh'], P['e_sh']) /
             _S(Z, P['c_sh'], P['d_sh'], P['e_sh']))
    return 1.0 + (ratio - 1.0) * env


def _closure_bumps(Z, N):
    """Shell-closure bump corrections at N = 10, 18, 36, 54, 86."""
    bump = 0.0
    for Ak, mk in zip(
        (P['A10'], P['A18'], P['A36'], P['A54'], P['A86']), MAGIC
    ):
        bump += Ak * (Z / mk) ** P['gb'] * math.exp(-((N - mk) / P['w_b']) ** 2)
    return bump


def _relativistic_ratio(Z, N):
    """Effective-charge relativistic correction ratio.
    Z_eff(Z, N) is a purely phenomenological calibration parameter
    (an effective relativistic weight of the remaining electron
    cloud), NOT a screened charge seen by an individual electron.
    The hard constraint Z_eff <= Z is defensive: within the model
    domain (the ratio is evaluated only for N > 10) Z_eff < Z holds
    mathematically for every element."""
    x = N / Z
    Zeff = Z * (1.0 + P['r1'] * (1.0 - x) + P['r2'] * (1.0 - x) ** 2)
    Zeff = min(Zeff, Z)  # hard physical constraint (defensive no-op in-domain)
    return ((1.0 - (Zeff / ALPHA_INV) ** 2) /
            (1.0 - (Z / ALPHA_INV) ** 2)) ** (-0.5)


def tbe_ion_raw(Z, N):
    """Core closed-form evaluation (Eq. 12a-12d), eV."""
    Ev2 = tbe_neutral_v2(Z)
    if N == Z:
        return Ev2
    Dz = D_dirac(Z)
    EK1 = min(Z ** 2 * E_H * Dz, Ev2)
    EK2 = min(2.0 * (Z - P['sK']) ** 2 * E_H * Dz, Ev2)
    sL = P['sL0'] + P['sL1'] * Z
    EK10 = min(EK2 + 8.0 * (Z - sL) ** 2 / 4.0 * E_H * Dz, Ev2)
    bump = _closure_bumps(Z, N)
    if N <= 1:
        return EK1
    if N <= 2:
        return EK2
    if N <= 10:
        span = max(min(10.0, float(Z)) - 2.0, 1.0)
        y = min(max((N - 2.0) / span, 0.0), 1.0)
        base = EK2 + (EK10 - EK2) * y ** P['qL']
        return base * (1.0 + bump * (1.0 - y)) * _shell_modulation(y, Z, N)
    y = (N - 10.0) / max(Z - 10.0, 1.0)
    y = min(max(y, 0.0), 1.0)
    base = EK10 + max(Ev2 - EK10, 0.0) * y ** P['qM']
    Rratio = _relativistic_ratio(Z, N)
    return base * Rratio * (1.0 + bump * (1.0 - y)) * _shell_modulation(y, Z, N)


def _boundary_base(Z):
    """Period-scaled estimate of the last-electron ionization energy (eV)."""
    if Z <= 10:
        return 10.0
    if Z <= 18:
        return 8.0
    if Z <= 36:
        return 7.0
    if Z <= 54:
        return 8.0
    if Z <= 86:
        return 6.0
    return 5.0


def _corrected_sequence(Z):
    """Full charge-state sequence with the v3.2 physical-consistency
    safeguard. Returns the corrected list TBE(N), N = 1..Z (eV).

    Tier 1 (cap): E_raw <= E_v2 for N < Z.
    Tier 2 (boundary layer, N > Z-5): a link k (counted inward from the
        neutral anchor) is reconstructed whenever its implied ionization
        energy is non-monotone OR falls outside the plausibility window
        0.5 eV <= IP_k <= 4 * IP_base * 1.8^(k-1); the reconstruction
        uses the increasing chain IP_k = 1.8 * IP_{k-1}.
    Tier 3 (core, N <= Z-5): minimal-touch monotone projection
        (0.5 eV floor), no cascade."""
    Ev2 = tbe_neutral_v2(Z)
    seq = [min(tbe_ion_raw(Z, n), Ev2) for n in range(1, Z + 1)]  # Tier-1 cap
    prev_ie = None
    for i in range(Z - 2, 0, -1):  # N = Z-1 ... 2 (backward scan)
        n_cur = i + 1
        if n_cur > Z - 5:
            k = Z - n_cur  # k = 1 is the first ionization energy
            ceiling = 4.0 * _boundary_base(Z) * (1.8 ** (k - 1))
            implied = seq[i + 1] - seq[i]
            trigger = (implied < 0.5) or (implied > ceiling)
        else:
            trigger = (seq[i] >= seq[i + 1])  # core: monotonicity only
        if trigger:
            if n_cur > Z - 5:
                # boundary layer: IPs INCREASE across the valence shell
                ie = _boundary_base(Z) if prev_ie is None else 1.8 * prev_ie
            else:
                # core region: minimal-touch projection (no cascade)
                ie = 0.5
            seq[i] = seq[i + 1] - max(ie, 0.5)
            prev_ie = seq[i + 1] - seq[i]
    return seq


def tbe_ion(Z, N):
    """Total electron binding energy of ion (Z, N) in eV, with the
    v3.2 physical-consistency safeguard."""
    if not (1 <= Z <= 92):
        raise ValueError("Z must be in [1, 92]")
    if not (1 <= N <= Z):
        raise ValueError("N must be in [1, Z]")
    if Z <= 2 or N == Z:
        return tbe_ion_raw(Z, N)
    return _corrected_sequence(Z)[N - 1]


def tbe(Z, N=None):
    """Unified interface. If N is None, returns neutral atom (N=Z). eV."""
    if N is None:
        N = Z
    return tbe_ion(Z, N)


def _self_test():
    """Built-in demonstration / self-test (no external data needed)."""
    e_fe = tbe(26, 13)
    print(f"Fe (Z=26, N=13): {e_fe:.1f} eV = {e_fe/1000:.3f} keV")
    print(f" (NIST reference: 32580.77 eV; error: {(e_fe/32580.77-1)*100:+.2f}%)")
    e_u = tbe(92)
    print(f"\nU (Z=92, N=92): {e_u/1000:.1f} keV")
    print(f" (Dzuba-Flambaum reference: 761.87 keV; error: {(e_u/1000/761.87-1)*100:+.2f}%)")
    e_u1 = tbe(92, 1)
    print(f"\nU91+ (Z=92, N=1): {e_u1/1000:.1f} keV")
    print(f" (Dirac point-nucleus value; NIST: 131.82 keV; error: {(e_u1/1000/131.82-1)*100:+.2f}%)")
    # v3.2 physical-consistency self-check on the full grid
    n_bad = 0
    for Z in range(3, 93):
        s = _corrected_sequence(Z)
        for n in range(1, Z):
            if s[n - 1] >= s[n] - 1e-12:
                n_bad += 1
    print(f"\nFull-grid self-check (Z=3..92): non-monotonic links = {n_bad} "
          f"(must be 0)")
    for Z, name, exp in ((14, 'Si', 8.15), (26, 'Fe', 7.90), (92, 'U', 6.19)):
        ip1 = tbe(Z, Z) - tbe(Z, Z - 1)
        print(f"First ionization energy {name}: model {ip1:+.2f} eV "
              f"(estimate; experiment {exp} eV)")


if __name__ == "__main__":
    if len(sys.argv) >= 2:
        Z = int(sys.argv[1])
        N = int(sys.argv[2]) if len(sys.argv) >= 3 else Z
        e = tbe(Z, N)
        q = Z - N
        print(f"Z = {Z}, N = {N} (q = +{q}): "
              f"TBE = {e:.4f} eV = {e/1000:.6f} keV")
    else:
        _self_test()

