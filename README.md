# TF-RSE-IonCF v3.2

[![DOI](https://zenodo.org/badge/1426650227.svg)](https://doi.org/10.5281/zenodo.23524783)

**A Unified Closed-Form Model for the Total Electron Binding Energy of All Neutral Atoms and Ions (Z = 1–92, q = 0 to Z−1)**

Companion code for the paper:

> J. Cheng and D. Cheng, *"TF-RSE-IonCF: A Unified Closed-Form Model for the Total Electron Binding Energy of All Neutral Atoms and Ions (Z = 1–92, q = 0 to Z−1)"*, **Atomic Data and Nuclear Data Tables** (2026), accepted for publication, manuscript ADNDT-D-26-00055R2.

## Overview

TF-RSE-IonCF returns the total electron binding energy (TBE) — the energy required to remove **all** electrons — of any neutral atom or ion with nuclear charge Z = 1–92 and electron number N = 1–Z, from a single closed-form expression:

- **29 fitted parameters**, no external dependencies, pure Python 3.8+
- Neutral anchor layer (TF-RSE v2): mean error **0.40%** over all 92 neutral atoms
- Full grid (4278 charge states): mean |error| **0.260%**, median **0.161%**, 97.8% of points within 1%, 100% within 5%
- Validated against NIST ASD v5.12 and Dzuba–Flambaum–Afanasjev (2024) reference data

## Installation

No dependencies. Copy `tfrse_ioncf_v3.2.py` into your project, or run it directly.

```bash
python tfrse_ioncf_v3.2.py            # built-in self-test
python tfrse_ioncf_v3.2.py 26 13      # Fe XIII (Z=26, N=13)
python tfrse_ioncf_v3.2.py 92         # neutral uranium
```

## Python API

```python
from tfrse_ioncf_v3_2 import tbe

e = tbe(26, 13)   # Fe 13+ total electron binding energy, in eV
e = tbe(92)       # neutral atom: N defaults to Z
```

- `tbe(Z, N=None)` — unified interface; N = electron number, defaults to the neutral atom
- `tbe_ion(Z, N)` — ionic TBE with the v3.2 physical-consistency safeguard
- `tbe_neutral_v2(Z)` — neutral-anchor layer alone

## Validation data

`validation_full.xlsx` contains the complete 4278-point validation grid:

| Sheet | Contents |
|---|---|
| `validation_full` | Z, element, N, q, TBE_ref_eV, TBE_model_eV, error_pct for every charge state |
| `summary` | Aggregate error statistics (full grid, neutral atoms, ions) |

## License

**PolyForm Noncommercial License 1.0.0** — free for academic, educational, and other noncommercial use (including use by universities, public research organizations, and government institutions). Commercial application — including use in proprietary products, consulting services, or for-profit research — requires **prior written consent from the authors**. See `LICENSE`.

## Citation

If you use this code or the accompanying data, please cite the paper (see `CITATION.cff`):

```bibtex
@article{Cheng2026TFRSEIonCF,
  author  = {Cheng, Jinjun and Cheng, Dian},
  title   = {TF-RSE-IonCF: A Unified Closed-Form Model for the Total Electron
             Binding Energy of All Neutral Atoms and Ions (Z = 1--92, q = 0 to Z-1)},
  journal = {Atomic Data and Nuclear Data Tables},
  year    = {2026},
  note    = {Accepted for publication}
}
```

## Contact

Jinjun Cheng — cdcdcd999@126.com
