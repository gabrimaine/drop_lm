# drop_lm — fitting lipid-droplet shapes with the Deslandes–Thiam–Forêt model

This repository fits a measured lipid-droplet (LD) profile `z(r)` with the
membrane-elasticity model of **Deslandes, Thiam & Forêt, *Biophys. J.* 113:15
(2017)** (Supporting Material, Sect. 3–5), and extracts:

- the droplet **volume** `V`,
- the tension ratio **γ_b / γ_m** (bilayer / monolayer),
- the **monolayer tension γ_m** and **bilayer tension γ_b**,

at a fixed monolayer bending rigidity **κ_m = 6 k_BT**.

The paper gives a free-energy functional, not a formula for `z(r)`. The code
minimises it by solving the Euler–Lagrange shape equations numerically and fits
the resulting shape to the data. Everything (code, data, results, and the
record of how each result was made) is in this repository.

> **New to Lightcone?** You do not need it to *look at* the results. They are
> plain files (JSON, CSV, PNG) once downloaded. Lightcone (`lc`) is only needed
> to *recompute* them, or to browse them in its viewer. Sections 2–3 are enough
> to explore; Sections 5–6 explain how to rerun.

## Quick start: read the report with MyST

The report is `index.md`, rendered with [MyST](https://mystmd.org). Every number
and figure in it is read from `results/` when the page is built, so the result
files must be downloaded first.

**Step 1: install the tools (once)**
- Node.js 20 or newer: `brew install node` (macOS), or from https://nodejs.org
- MyST CLI: `npm install -g mystmd`
- git-annex: `brew install git-annex` (macOS) or `apt install git-annex`

**Step 2: clone and download the results**
```bash
git clone git@github.com:gabrimaine/drop_lm.git
cd drop_lm
git annex init                      # once per clone
git annex enableremote github-lfs   # connect to the file storage on GitHub
git annex get results/ data/        # without this the figures are broken images
```

**Step 3: start the preview**
```bash
myst start
```
The first run downloads the pinned plugin (MySTRA v0.0.8) and theme
(astra-article-theme v0.0.16), which takes a minute. When it is ready it prints
`👉 http://localhost:3000 👈` (or the next free port).

**Step 4: open the printed URL in your browser.**

*Running on a remote machine (e.g. S3DF)?* Forward **both** ports: the page is
served on the printed port (3000) and the images on 3100.
```bash
ssh -L 3000:localhost:3000 -L 3100:localhost:3100 <remote-host>
```

**Step 5: switch between analysis choices.** Go to the section **"Changing
the analysis choices"** and click the **tabs** (Baseline, 298 K, Free
baseline, r ≤ 40 nm, Bootstrap). The decision boxes higher up (e.g.
"Temperature: 310 K / 298.15 K") are display-only; clicking them does nothing.

**Step 6: look at the code.** The last section of the report, **"Code"**, lists
every script. Click a file name to open it as plain text in a new tab (and save
it from there), or expand its box to read it in the page. On GitHub, see
[The code](#the-code) below.

**Static version (optional)**
```bash
myst build --html
cd _build/html && python3 -m http.server 8000   # then open http://localhost:8000
```
Don't open `_build/html/index.html` by double-clicking it. The page needs a
web server: from a `file://` URL nothing loads and the tabs do not work.

---

## Contents

1. [What is in here](#1-what-is-in-here)
2. [Get the repository and its files](#2-get-the-repository-and-its-files)
3. [Explore the results (no Lightcone needed)](#3-explore-the-results-no-lightcone-needed)
4. [The analysis choices ("decisions") and "universes"](#4-the-analysis-choices-decisions-and-universes)
5. [Browse interactively: the Lightcone viewer and the report](#5-browse-interactively-the-lightcone-viewer-and-the-report)
6. [Recompute the results](#6-recompute-the-results)
7. [Fit your own curve](#7-fit-your-own-curve)
8. [The physics in one page](#8-the-physics-in-one-page)
9. [Assumptions and caveats](#9-assumptions-and-caveats)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. What is in here

```
astra.yaml          the analysis specification: inputs, outputs, decisions,
                    and the exact command that produces each output
universes/          one YAML file per combination of decisions ("universe")
results/<universe>/ the outputs of each universe (see table below)
src/                the Python code (solver, fit, plots)
data/raw/           the measured profile (td_reset2.xvg), the paper (PDF), and
                    compare_fit.py, an independent fit used for cross-checking
pyproject.toml,     the Python environment (numpy, scipy, matplotlib),
uv.lock             pinned exactly
index.md, myst.yml  the report, written as a short paper (MyST); see Section 5
references.bib      its bibliography
AGENTS.md           working notes for the AI assistant that built this
```

Each universe directory contains the same seven outputs:

| File | What it is |
|---|---|
| `fit_params.json` | **The answer.** Best-fit V, γ_b/γ_m, γ_m, γ_b and derived quantities (contact radius R, thickness h, λ_m, pressure, contact angles), their uncertainties and correlation matrix, and fit-quality numbers. |
| `fit_summary.csv` | The same numbers as a flat table, one row per quantity: `quantity, label, value, value_std, unit`. The report reads its numbers from here. |
| `fit_plot.png` | Data vs best-fit shape, residuals, and the fitted values. |
| `fit_profile.csv` | The best-fit model contour: columns `s_nm, r_nm, z_nm, psi_deg, dpsi_ds_inv_nm, is_bilayer` (arclength, position, tangent angle, curvature; `is_bilayer = 1` marks the flat bilayer beyond the contact line). |
| `identifiability.png` | χ² map over (γ_m, γ_b/γ_m): shows that this profile constrains *both* tensions, not just a combination of them. |
| `solver_validation.png` | Checks of the shape solver against the paper's analytic small- and large-volume limits. It doesn't depend on the data, so it's identical in every universe. |
| `cross_check.csv` | The profile refitted by the independent code `data/raw/compare_fit.py`, next to our fit: for each quantity both values, both 1σ errors, and the difference in units of our σ. |

Next to each output there is a hidden `.<output>.manifest.json`: the run
record (exact command, decisions, git commit, input checksums, timestamps).

### The code

[![Browse the code](https://img.shields.io/badge/code-src%2F-blue)](src/)
[![Download ZIP](https://img.shields.io/badge/download-ZIP-green)](https://github.com/gabrimaine/drop_lm/archive/refs/heads/main.zip)

Click a file to read it on GitHub. Its page has **Raw** and **Download raw
file** buttons (top right of the file).

| File | What it does | Makes |
|---|---|---|
| [`src/lens.py`](src/lens.py) | Shape solver: shape equations, axis series start, contact line, lens-branch selection | used by every output |
| [`src/fitting.py`](src/fitting.py) | Least-squares model, initial grid search, derived physical quantities | used by `fit_params`, `identifiability` |
| [`src/fit.py`](src/fit.py) | Runs the fit and its uncertainty estimate | `fit_params` |
| [`src/fit_profile.py`](src/fit_profile.py) | Writes the best-fit contour | `fit_profile` |
| [`src/plot_fit.py`](src/plot_fit.py) | Data vs model figure | `fit_plot` |
| [`src/identifiability.py`](src/identifiability.py) | Profile Δχ² map over the two tensions | `identifiability` |
| [`src/validate_solver.py`](src/validate_solver.py) | Solver against the SI's analytic limits | `solver_validation` |
| [`src/summarize.py`](src/summarize.py) | Flattens the fit into one row per quantity | `fit_summary` |
| [`src/cross_check.py`](src/cross_check.py) | Refits with the independent code and tabulates both | `cross_check` |
| [`data/raw/compare_fit.py`](data/raw/compare_fit.py) | Independent implementation of the same fit (standalone; `python data/raw/compare_fit.py --test` runs its self-tests) | used by `cross_check` |

The exact command that runs each script is in [`astra.yaml`](astra.yaml).
Two caveats on the ZIP: the repository is private, so only collaborators can
download it; and it holds the scripts in full but only one-line annex pointers
for `data/` and `results/` (Section 2.2).

---

## 2. Get the repository and its files

### 2.1 Prerequisites

| Tool | Why | Install |
|---|---|---|
| `git` | the repository | your OS package manager |
| `git-annex` | downloads the data and result files (see below) | bundled with Lightcone on Linux (below), or `brew install git-annex` / `apt install git-annex` / `conda install -c conda-forge git-annex`. Version 8 or newer; 10.x recommended. |
| `uv` | Python environment, only to run code | `curl -LsSf https://astral.sh/uv/install.sh \| sh` ([docs](https://docs.astral.sh/uv/)) |
| `lc` (Lightcone CLI) | only for the viewer and recomputation | `uv tool install lightcone-cli==0.5.1` (on Linux this also installs `git-annex`). |

The repository on GitHub is **private**: ask the owner to add you as a
collaborator, and make sure your SSH key is registered with GitHub
(`ssh -T git@github.com` should greet you by name).

### 2.2 Why there is an extra step: git-annex

Git is bad at storing binary files, so data and results (everything under
`data/` and `results/`) are managed by **git-annex**. Git only stores a small
*pointer* for each such file, and the contents are stored in **GitHub LFS** on
the same repository, through a git-annex "special remote" named `github-lfs`.

Practical consequence: after a plain `git clone`, files like
`results/baseline/fit_plot.png` are one-line text pointers such as
`/annex/objects/SHA256E-s121856--e42a….png`. GitHub's web page shows the same
pointer. You need one more command to download the real contents.

### 2.3 Clone and download

```bash
git clone git@github.com:gabrimaine/drop_lm.git
cd drop_lm
git annex init                      # register this clone with git-annex
git annex enableremote github-lfs   # connect to the LFS storage on GitHub
git annex get .                     # download everything (~4 MB)
```

You can also download selectively, e.g. `git annex get results/baseline` or
`git annex get results/*/fit_params.json`.

Useful git-annex commands:

```bash
git annex whereis results/baseline/fit_plot.png   # where copies of a file live
git annex find --not --in=here                     # files not downloaded yet
git annex drop results/bootstrap                   # free local space (safe: refuses if no other copy exists)
```

Files are checksummed (SHA-256), and git-annex verifies every download.

---

## 3. Explore the results (no Lightcone needed)

After `git annex get .`, everything is ordinary files.

**Look at the figures:**

```bash
open results/baseline/fit_plot.png        # macOS
xdg-open results/baseline/fit_plot.png    # Linux
```

**Print the main numbers for every universe:**

```bash
python3 - <<'PY'
import json, glob, os
print(f"{'universe':10s} {'V [nm^3]':>14s} {'gb/gm':>8s} {'gm [mN/m]':>14s} {'gb [mN/m]':>14s}")
for f in sorted(glob.glob("results/*/fit_params.json")):
    d = json.load(open(f)); b = d["best_fit"]; s = d["uncertainty"]["sigma"]
    u = os.path.basename(os.path.dirname(f))
    print(f"{u:10s} {b['V_nm3']:8.0f} ± {s['V_nm3']:<3.0f} {b['gamma_ratio']:8.4f} "
          f"{b['gamma_m_mN_per_m']:7.3f} ± {s['gamma_m_mN_per_m']:.3f} "
          f"{b['gamma_b_mN_per_m']:7.3f} ± {s['gamma_b_mN_per_m']:.3f}")
PY
```

At the time of writing (commit `169a0c9`) this gives:

| Universe | What differs from baseline | V (nm³) | γ_b/γ_m | γ_m (mN/m) | γ_b (mN/m) |
|---|---|---|---|---|---|
| `baseline` | — | 11483 ± 3 | 0.8742 | 3.373 ± 0.027 | 2.949 ± 0.027 |
| `temp-298` | T = 298.15 K instead of 310 K | 11483 | 0.8742 | 3.244 | 2.836 |
| `z0-free` | free vertical offset of the bilayer plane | 11491 | 0.8737 | 3.359 | 2.935 |
| `r-le-40` | only r ≤ 40 nm fitted | 11483 | 0.8742 | 3.373 | 2.949 |
| `bootstrap` | bootstrap instead of Jacobian errors | 11483 | 0.8742 | 3.373 | 2.949 |

The JSON files are the source of truth; if this table and the files ever
disagree, the files win.

**Key fields in `fit_params.json`:**

- `best_fit.V_nm3`: total droplet volume (both halves of the lens), nm³
- `best_fit.gamma_ratio`: γ_b/γ_m
- `best_fit.gamma_m_mN_per_m`, `best_fit.gamma_b_mN_per_m`: tensions, mN/m
- `best_fit.lambda_nm`: λ_m = √(κ_m/γ_m), the bending length
- `best_fit.R_nm`, `best_fit.h_nm`: contact-line radius and total thickness
- `uncertainty.sigma`: 1σ errors; `uncertainty.correlation`: correlation matrix
- `fit_quality.rms_nm`: rms residual (≈ 8 pm for a droplet 11 nm high)
- `decisions`: which choices produced this file

**Load the model contour, e.g. to overlay it in your own plots:**

```python
import numpy as np
m = np.genfromtxt("results/baseline/fit_profile.csv", delimiter=",", names=True)
print(m.dtype.names)        # column names
```

---

## 4. The analysis choices ("decisions") and "universes"

Any analysis involves choices that are defensible but not unique. Here they are
written down explicitly in `astra.yaml` (the [ASTRA](https://github.com/LightconeResearch)
format) as **decisions**, each with options, a default, and a rationale:

| Decision | Options (default first) | Why it matters |
|---|---|---|
| `temperature` | `t310` (310 K), `t298` (298.15 K) | Converts κ_m = 6 k_BT to joules. Both tensions scale ∝ T; the lengths and the volume don't change. |
| `baseline` | `fixed` (bilayer at z = 0), `free` (fit an offset z0) | Tests whether a reference shift biases the fit. |
| `fit_range` | `all` (r ≤ 54 nm), `r_le_40` | The far tail drifts slightly below 0, which is not droplet shape. |
| `uncertainty` | `jacobian`, `residual_bootstrap` (200 refits) | Two error estimates. Both treat residuals as independent, so they are lower bounds. |

A **universe** is one complete choice of options, stored as a small file in
`universes/`:

```yaml
# universes/temp-298.yaml
id: temp-298
description: Baseline with kappa_m = 6 kBT converted at 298.15 K instead of 310 K
decisions:
  temperature: t298
  baseline: fixed
  fit_range: all
  uncertainty: jacobian
```

Every universe gets its own `results/<universe>/` directory. The four variants
here each change **one** decision from `baseline`, so you can see the effect of
each choice in isolation. Combinations of two or more changes are not
precomputed; Section 6 shows how to add one.

The `recipe:` of each output in `astra.yaml` is the exact command that made it,
with `{decisions.x}` filled in from the universe. That is how a decision reaches
the code.

---

## 5. Browse interactively: the Lightcone viewer and the report

```bash
uv tool install lightcone-cli==0.5.1   # once
lc view                                # prints a URL (http://localhost:<port>/lightcone?token=…)
```

Open the URL in a browser. The page renders `astra.yaml`: the inputs, the
outputs with their run records, the findings and the papers behind them, and
each decision with its options and rationale. It always shows the **baseline**
universe (lightcone-cli 0.5.1 fixes it, and the option lists cannot be
clicked to switch). To see another universe, open `results/<universe>/`, or
use the tabs in the report (below).
`lc view` returns immediately, and the server stops itself after 30 minutes
of inactivity.

**If the repository is on a remote machine** (e.g. a cluster login node), run
`lc view` there and forward the port from your laptop:

```bash
ssh -L 8889:localhost:8889 <remote-host>   # use the port and host that lc view printed
```

then open the printed URL locally.

### Read the report

`index.md` is a short paper about the analysis: model, method, results and
discussion, with a bibliography (`references.bib`). It contains no typed-in
numbers. Every value, figure and decision is pulled from `astra.yaml` and
`results/` when the page is built, so it stays in step with the analysis.
To view it you need Node.js 20 or newer and the MyST CLI:

```bash
npm install -g mystmd        # once
myst start                   # live preview; prints the URL (port 3000 or the next free one)
myst build --html            # or a static site in _build/html/ (serve it, see below)
```

A static build has to be served, e.g. `cd _build/html && python3 -m http.server 8000`;
opened from a `file://` URL the page does not load. On a remote machine,
forward the page port *and* port 3100, which serves the images (see the Quick
start).

The plugin (MySTRA v0.0.8) and the theme (astra-article-theme v0.0.16) are
pinned in `myst.yml` and downloaded on the first build. Run `myst build`
(plain) to see reference problems. They are printed as lines starting with
`⛔️`, and the build still exits 0.

The report shows the first universe in `universes/` (`baseline`). The section
"Changing the analysis choices" has one tab per universe with that
universe's own fit plot. Neither the report nor `lc view` can switch
decisions in the page. The section "Independent cross-check" compares the fit
with `data/raw/compare_fit.py`, and the last section, "Code", shows and links
every script. For any other combination, add a universe (Section 6)
and look in `results/<universe>/`.

---

## 6. Recompute the results

### 6.1 With Lightcone (recommended: keeps the run records)

Lightcone is tested here on Linux; it runs each recipe in a sandbox with only
the declared inputs visible.

```bash
uv tool install lightcone-cli==0.5.1
lc compute launch --wait        # start a local "cluster" named `local` (this machine)
lc status                       # current / stale / never-made, per universe and output
lc materialize local            # make whatever is stale, committing each output as it lands
lc materialize --check          # exit 0 only if everything is up to date
```

`lc materialize` needs a clean git tree (commit first). It remakes an output
only when its recipe, its decisions, or the content of one of its inputs
changed. Editing `src/fit.py` therefore remakes `fit_params` and everything
downstream, and nothing else.

You can target one output or one universe: `lc materialize local fit_params`
(all universes) or `lc materialize local temp-298/fit_plot`.

**Add a universe** (e.g. 298 K *and* a free offset):

```bash
cat > universes/t298-z0free.yaml <<'YAML'
id: t298-z0free
description: 298.15 K and a free bilayer offset
decisions:
  temperature: t298
  baseline: free
  fit_range: all
  uncertainty: jacobian
YAML
uvx astra-tools@0.2.18 validate      # check the spec and universes
git add universes/t298-z0free.yaml && git commit -m "Add t298-z0free universe"
lc materialize local                 # computes results/t298-z0free/
```

Keep `baseline.yaml` alphabetically first in `universes/`: the MyST report
renders the first file.

**Change or add a decision:** edit `astra.yaml`, make the code accept the new
option (scripts take each decision as a command-line flag), add the selection
to *every* universe file, validate, commit, materialize.

**Share new results:** `git annex copy --to github-lfs && git push origin main git-annex`.

Cost on 16 cores: `fit_params` ≈ 1 min (bootstrap ≈ 2.5 min), `identifiability`
≈ 2.5 min, everything else seconds. A universe takes about 4 minutes.

### 6.2 Without Lightcone (plain Python)

The scripts are ordinary command-line programs; `uv` creates the exact
environment from `uv.lock`:

```bash
uv sync
OUT=$(mktemp -d)        # outside the repository, so the git tree stays clean
uv run python src/fit.py --profile data/raw/td_reset2.xvg --kappa-kT 6 \
    --temperature t310 --baseline fixed --fit-range all --uncertainty jacobian \
    --procs 4 --output $OUT/fit_params.json
uv run python src/fit_profile.py --fit-params $OUT/fit_params.json --output $OUT/fit_profile.csv
uv run python src/plot_fit.py --profile data/raw/td_reset2.xvg \
    --fit-profile $OUT/fit_profile.csv --fit-params $OUT/fit_params.json \
    --output $OUT/fit_plot.png
uv run python src/summarize.py --fit-params $OUT/fit_params.json --output $OUT/fit_summary.csv
echo $OUT
```

The commands for the other outputs are the `recipe:` entries in `astra.yaml`.
Write outside the repository, not into `results/`: Lightcone treats files in
`results/` that it did not make as foreign and will remake them, and any
untracked file makes `lc materialize` refuse to start until you commit or remove it.

---

## 7. Fit your own curve

Input format: two whitespace-separated columns, `r` and `z` in **nm**; lines
starting with `#` or `@` (GROMACS `.xvg` headers) are ignored. `z` is the
height of the oil/monolayer interface above the flat bilayer plane, and `r` is
the distance from the droplet axis (see the assumptions in Section 9).

Quick, outside Lightcone:

```bash
uv run python src/fit.py --profile my_curve.xvg --kappa-kT 6 \
    --temperature t310 --baseline fixed --fit-range all --uncertainty jacobian \
    --output my_fit.json      # e.g. outside the repository
```

To make it a tracked part of the analysis, put the file in `data/raw/`, declare
it as an input in `astra.yaml` (copy the `profile` entry), add outputs for it,
and `lc materialize`.

`--kappa-kT` sets κ_m in units of k_BT. The tensions scale with κ_m (γ_m =
κ_m/λ_m²); the shape, and hence V and γ_b/γ_m, do not.

---

## 8. The physics in one page

The droplet is a lens of oil between the two leaflets of a bilayer. In the
**symmetric** model (SI Sect. 5, identical leaflets) the bilayer stays flat
and only one monolayer shape has to be found, by axisymmetric parametrisation
in arclength `s`: `r(s)`, `z(s)`, tangent angle `ψ(s)`, curvature
`u = dψ/ds`, and a Lagrange multiplier `t(s)`. Minimising the free energy at
fixed volume gives the shape equations (SI Eq. 18):

```
r' = cos ψ,   z' = −sin ψ,   ψ' = u,
u' = −u cos ψ / r + sin ψ cos ψ / r² − p r cos ψ / (2κ) + t sin ψ / (κ r),
t' = κ/2 (u² − sin²ψ / r²) + γ_m − p r sin ψ
```

with `p` the oil pressure (the Lagrange multiplier for the volume).

**Boundary conditions:**
- **Axis:** `r = ψ = 0`. The solver starts from the series solution ψ ≈ a s + (a − P/2) s³/8, where `a` is the curvature at the top.
- **Contact line:** `ψ` returns to 0, where the monolayer meets the flat bilayer.
- **Force balance at the contact line:** κ_m u(R)² = 2γ_m − γ_b. This is the condition that selects the physical shape.

**Numerics** (`src/lens.py`): in units of λ_m = √(κ_m/γ_m), a shape depends
only on (γ_b/γ_m, a). For each pair, the pressure is found by shooting from the
axis until the contact condition holds. Usually two solutions exist, and the
code keeps the lens (ψ < 90°) rather than the overhanging near-sphere. The
solver conserves the Hamiltonian (SI Eq. 19) to ~10⁻⁹ and reproduces the
paper's analytic limits (`solver_validation.png`).

**Fit** (`src/fitting.py`): least squares over (ln λ_m, γ_b/γ_m, ln a), plus z0
if `baseline: free`, after a coarse grid search for the starting point. Then
γ_m = κ_m/λ_m², γ_b = (γ_b/γ_m)·γ_m, and V is the integral of the fitted shape.

**Why both tensions can be measured:** for small droplets (V ≪ λ_m³) the shape
only constrains 2γ_m − γ_b (SI Eqs. 31–34). This droplet has V/λ_m³ ≈ 550, so
both tensions are constrained separately; `identifiability.png` shows a closed
χ² minimum.

---

## 9. Assumptions and caveats

- **Geometry:** `z(r)` is taken as the upper half of an up-down **symmetric**
  lens, measured from the flat bilayer plane (total thickness 2 z(0)). If the
  two leaflets differ, the full asymmetric contact conditions of the SI are
  needed and the bilayer bending rigidity enters.
- **Temperature:** the simulation temperature was not specified; 310 K is the
  default and 298.15 K a variant. Use the real value if you know it.
- **A factor of ½:** the SI text writes the contact-line tension as
  t = γR − κRu². Its own Eq. 19 and small-volume results are consistent only
  with t = γR − κRu²/2, which is what the code uses.
- **Error bars are lower bounds:** both methods treat residuals as independent,
  but neighbouring points along a profile are correlated. Read the ±0.03 mN/m
  on the tensions as a few-percent uncertainty at best. The spread across
  universes (≈ 0.5% for a free offset) is a better guide to systematics, and the
  temperature choice dominates the absolute tensions.
- κ_m = 6 k_BT is an input. Tensions scale linearly with it.

---

## 10. Troubleshooting

- **A result file contains `/annex/objects/SHA256E-…` instead of data:** the
  content isn't downloaded yet. Run `git annex get <file>`.
- **`git annex get` says "not available":** run `git annex enableremote github-lfs`
  first, and check `ssh -T git@github.com` works.
- **`Permission denied (publickey)`:** your SSH key isn't on GitHub, or you
  don't have access to the private repository.
- **`lc materialize` refuses to start:** the git tree isn't clean; commit or stash.
- **An output shows as `stale` after you wrote into `results/` by hand:**
  expected (a "foreign write"). `lc materialize` will remake it.
- **The report page is blank or unstyled:** you opened `_build/html/index.html`
  from disk. Serve it (`python3 -m http.server`) or use `myst start`.
- **The report's tabs switch but the figures are empty:** on a remote machine,
  forward port 3100 as well as the page port. Locally, the PNGs are probably
  not downloaded: `git annex get results/`.
- **Lots of `lmod/init/bash … Permission denied` lines in `lc` logs on HPC
  systems:** harmless. The sandbox blocks a few tools that the cluster's
  shell start-up scripts call.

## Reference

F. Deslandes, A. R. Thiam, L. Forêt, *Lipid Droplets Can Spontaneously Bud
Off from a Symmetric Bilayer*, Biophys. J. 113:15 (2017), and its Supporting
Material. PDF in `data/raw/deslandes_thiam_foret_2017_bpj.pdf`.
