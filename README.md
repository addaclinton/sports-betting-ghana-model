# Sports betting dynamics in Ghana: model code

Code and output data for the manuscript *A Compartmental Mathematical Model for Sports Betting Dynamics in Ghana* (under review; citation to be added).

The code integrates a deterministic five-compartment model of sports betting among Ghanaian adults, computes the basic reproduction number and the endemic equilibrium, and reproduces every number, table and figure dataset in the manuscript: the baseline simulation, the sensitivity indices, the Latin hypercube uncertainty analysis (20,000 draws) and the intervention comparison.

## Model

Compartments: Susceptible `S`, Experimenters `E`, Regular bettors `R`, Problem bettors `P`, Recovered `Q`. Time is in years; `N = S+E+R+P+Q` is constant.

```
dS/dt = Lam + delta*Q - beta*S*(E+R+P)/N - mu*S
dE/dt = beta*S*(E+R+P)/N - (a1+g1+mu)*E
dR/dt = a1*E - (a2+g2+mu)*R
dP/dt = a2*R - (g3+mu)*P
dQ/dt = g1*E + g2*R + g3*P - (delta+mu)*Q
```

with `a = a1+g1+mu`, `b = a2+g2+mu`, `c = g3+mu` and

```
R0 = beta * ( 1/a + a1/(a*b) + a1*a2/(a*b*c) )
```

The closed-form endemic equilibrium (exists iff `R0 > 1`) is implemented in `endemic()` in `src/model.py`.

## Parameters and their sources

| Class | Parameter | Value | Source |
|---|---|---|---|
| A | `N0` adults (18+) | 17,931,673 | Ghana Statistical Service, 2021 Population and Housing Census |
| A | `mu` exit rate | 1/(65.69 - 18) = 0.02097 | World Bank life expectancy at birth, 2024 |
| A | `Lam` recruitment | `mu * N0` | stationary adult population |
| B | `a1`, `a2` | 0.6351, 0.1100 | calibrated so that the equilibrium composition E:R:P is 16:44:40 (Manu et al., 2024, *BMC Public Health*; 403 bettors, Volta Region) |
| C | `R0` (hence `beta`) | 1.30; range 0.8-6.0 | illustrative; upper end covers the value (5.49) that reproduces 42% past-year betting (TGM Research, 2022) |
| C | `g1` | 0.50; range 0.125-2.0 | no published value; factor-4 range |
| C | `g2` | 0.10; range 0.025-0.40 | no published value; factor-4 range |
| C | `g3` | 0.10; range 0.025-0.54 | upper bound from 12-month abstinence of 41.6% (Müller et al., 2017, *European Psychiatry*) |
| C | `delta` | 0.15; range 0.034-0.88 | bounds from Bellringer et al. (2025, New Zealand) and Müller et al. (2017) |

**Class C values are not estimates.** No Ghanaian data exist for transmission, quitting or relapse rates (Badu et al., 2023, *Health Promotion International*). They are scenario values and ranges, and the uncertainty analysis shows how much they matter. See Table 1 and Section 4.1 of the manuscript.

## Repository layout

```
src/model.py              model, R0, RK4 integrator, endemic equilibrium
scripts/run_analysis.py   reproduces all results -> data/
scripts/verify.py         checks the numbers in data/summary.json (and, optionally, a manuscript .tex)
tests/test_model.py       unit tests (R0 vs next-generation matrix, equilibrium, conservation, convergence, ...)
data/                     output data used in the manuscript
```

## Quick start

Requires Python 3.9+ with NumPy and SciPy.

```bash
pip install -r requirements.txt
python -m unittest discover -s tests -v      # 9 tests, about 3 s
python scripts/run_analysis.py               # regenerates data/ (about 20 s)
python scripts/verify.py                     # numerical self-checks
```

To write the results somewhere else (for example to compare with the shipped files): `python scripts/run_analysis.py --outdir /path/to/folder`. `make test`, `make reproduce` and `make verify` do the same on systems with `make`.

## Outputs and where they appear in the manuscript

| File | Content | Manuscript |
|---|---|---|
| `data/fig2.csv` | baseline time evolution, 30 years (millions of adults) | Figure 2 |
| `data/fig4.csv` | equilibrium share of adults betting vs `R0` | Figure 3 |
| `data/fig3.csv` | problem bettors under the two interventions (thousands) | Figure 4 |
| `data/summary.json` | R0, equilibrium, 30-year state, sensitivity indices, intervention results, survey-matched scenario, uncertainty analysis (medians, percentiles, PRCC) | Tables 2-4 and numbers in the text |

Note that the file names `fig3.csv` and `fig4.csv` are in the reverse order of the figure numbers in the manuscript.

## Reproducibility

- Tested with Python 3.13.17, NumPy 2.5.3 and SciPy 1.18.1. Regenerated outputs were identical to the shipped files in this environment.
- The Latin hypercube sampling uses a fixed seed (`20261027`). Results should be reproducible to numerical precision, but different NumPy/SciPy versions may change the sampled values slightly.
- The integrator is classical RK4 (step 0.005 years for single runs, 0.05 for the 20,000-run uncertainty analysis). The vectorized and single-run integrators agree to better than 1e-10 at baseline, and RK4 agrees with SciPy's DOP853 to better than 1e-9 (`tests/test_model.py`).

## AI assistance

I confirm that all derivations were carried out and verified by the authors; Claude was used only as a checking aid and did not generate any derivation. All text was written by the authors, with Claude used solely for word- and sentence-level language correction per Section 11 of the MWSJ Publication Ethics policy.

## License

MIT (see `LICENSE`).
