# Cancer Research Tool

<p align="center">
  <img src="docs/demo.gif" alt="Cancer Protein Explorer — AlphaFold structure, dynamics, and ranked druggable targets" width="800">
</p>

> A research & education tool that leverages Google DeepMind's **AlphaFold** to resolve protein folding dynamics and rank druggable cancer-driver residues into a transparent, explainable shortlist — on a laptop.

<p align="center">
  <img src="https://github.com/Felixsavedra-1/Cancer-Research-tool/actions/workflows/ci.yml/badge.svg" alt="Tests">
  <img src="https://img.shields.io/badge/Python-3.11-3776AB?logo=python&logoColor=white" alt="Python 3.11">
  <img src="https://img.shields.io/badge/Streamlit-app-FF4B4B?logo=streamlit&logoColor=white" alt="Streamlit">
  <img src="https://img.shields.io/badge/AlphaFold-powered-0B3D91" alt="AlphaFold-powered">
  <img src="https://img.shields.io/badge/AlphaMissense-pathogenicity-0B3D91" alt="AlphaMissense">
  <img src="https://img.shields.io/badge/use-research%20%2F%20education-informational" alt="Research / Education use">
</p>

Type a gene — say `TP53`, `KRAS`, or `BRAF`. The tool then:

1. fetches its real **AlphaFold** structure,
2. computes the protein's intrinsic **folding dynamics**,
3. scores every variant with DeepMind's **AlphaMissense**,
4. finds its **druggable pockets**, and
5. fuses it all into a ranked, explainable shortlist of **candidate cancer-driver residues** —
   each with a plain-English rationale.

No GPU required.

> **Does the science hold up?** As a sanity check, TP53's top three ranked residues come out
> **R273, R248, R175** — the textbook p53 drivers. KRAS tops at **G12**, BRAF at **V600**, EGFR
> at the **exon-19 deletion / L858 / T790** sites. Rediscovering known biology is reassuring, but
> it's a sanity check, not validation — so it's also benchmarked.
>
> **The benchmark, honestly.** Against a curated gold standard of known oncogenic residues across
> a **27-gene panel**, the score separates true drivers from merely-recurrent residues at
> **composite AUPRC ≈ 0.68** — about 3× the 0.22 random baseline. Two caveats stated up front:
> recurrence *alone* is a strong baseline (AUPRC ≈ 0.70) that the composite doesn't quite beat on
> this recurrence-defined set, and the ENM structural axis carried no ranking signal, so it was
> **dropped from the score** and kept only as a visualisation. Full harness and ablation live in
> [`benchmark/REPORT.md`](benchmark/REPORT.md) (regenerate with `python scripts/benchmark.py`).

## Highlights

- **Folding dynamics, no simulation.** Elastic Network Model normal-mode analysis (ProDy)
  extracts per-residue flexibility, hinge/domain pivots, and collective motion from a single
  AlphaFold structure — in milliseconds on a CPU. It's the laptop-grade route to the motions a
  GPU molecular-dynamics run would reveal.
- **Variant effect from DeepMind AlphaMissense.** Every residue carries a predicted
  pathogenicity score and class (likely benign → pathogenic), pulled from the same AlphaFold DB
  endpoint as the structure.
- **Druggable-pocket detection.** A LIGSITE-style geometric scan finds enclosed cavities and the
  residues that line them — the "is this targetable?" axis of the score.
- **Target Priority Score — an explainable ranking.** One transparent composite:
  `0.35·recurrence + 0.41·pathogenicity + 0.24·druggability`. It's enriched with **Open Targets**
  tractability and disease links. (Folding dynamics is shown for context but is *not* a scored
  term — the benchmark found it added no ranking signal.)
- **Real AlphaFold structures, zero GPU.** It pulls *finished* models from DeepMind's
  ~200M-protein set, so all the heavy lifting is a download — nothing runs AlphaFold locally.
- **Confidence-aware 3D, three colour modes.** `py3Dmol` rendering you can paint by pLDDT
  confidence, NMA flexibility, or the Target Priority Score.
- **Pure, unit-tested scientific core.** The dynamics, pathogenicity parsing, and scoring formula
  have no UI or network dependency, and are covered by a network-free test suite.

**Tech stack:** Python 3.11 · Streamlit · py3Dmol · ProDy · NumPy/SciPy · pandas · requests · pytest

## Quickstart

The dashboard is a single self-contained file. Hosted on any static URL (e.g. GitHub
Pages) it gives the full experience — 3D structure, folding-dynamics colouring, and the
precomputed **Target Priority** ranking. Locally, serve the folder so it can load the
precomputed `data/` files:

```bash
python3 -m http.server      # then open http://localhost:8000/cancer-explorer.html
```

You can still just `open cancer-explorer.html` directly — it loads TP53 with the live
structure + hotspots, but browsers block `file://` from reading the precomputed `data/`
JSON, so dynamics/priority colouring need a served copy (or GitHub Pages).

<details>
<summary><b>Developer / Python alternative (Streamlit)</b></summary>

A Streamlit version of the same tool is also included for local Python development:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

</details>

<details>
<summary><b>Live-analysis API — any gene, not just the featured set (optional)</b></summary>

The precomputed `data/` covers the 27 featured genes (the benchmark panel). To analyse **any** human gene from
the standalone HTML, run the live-analysis service — it executes the same Python engine and
returns the same JSON schema, computing once and disk-caching under `data/cache/`.

```bash
pip install ".[api]"
uvicorn api.main:app --port 8100          # then GET http://localhost:8100/analyze/PIK3CA
# or, containerised:
docker compose up                          # builds the image and serves on :8100
```

Point the HTML at it with a `?api=` query param (or set `window.EXPLORER_API_BASE`):

```
http://localhost:8000/cancer-explorer.html?api=http://localhost:8100
```

The dashboard still tries the precomputed static file first and only falls back to the API for
other genes, so it degrades gracefully to zero-install behaviour when no API is configured.

</details>

---

<details>
<summary><b>How it runs on a laptop</b></summary>

Predicting protein structures with AlphaFold needs big GPUs. But DeepMind already
ran AlphaFold on ~200 million proteins and made the results free. This tool
**downloads finished structures** instead of computing them, so all the heavy
work is done by hosted services and the browser — no GPU required.

</details>

<details>
<summary><b>Data sources</b> (all free, public)</summary>

| Data | Source |
|------|--------|
| Gene → protein + sequence | [UniProt](https://www.uniprot.org) |
| 3D structures (with pLDDT confidence) | [AlphaFold DB](https://alphafold.ebi.ac.uk) |
| Per-variant pathogenicity | [AlphaMissense](https://alphafold.ebi.ac.uk) (DeepMind) |
| Recurrent cancer mutation hotspots | [cancerhotspots.org](https://www.cancerhotspots.org) |
| Target tractability & disease links | [Open Targets](https://platform.opentargets.org) |
| Folding dynamics (ENM/NMA) | computed locally with [ProDy](http://prody.csb.pitt.edu) |
| Druggable pockets | computed locally (LIGSITE-style geometric scan) |

</details>

<details>
<summary><b>Using the app</b></summary>

In the sidebar: pick or type a gene, enter mutations like `R175H` (separate several
with commas or spaces), and choose how to colour the structure — **pLDDT confidence**,
**NMA flexibility**, or the **Target Priority Score**. The ranked priority table is the
main deliverable; click a residue's row to drop it onto the 3D structure. Export the full
ranking (every sub-score + rationale) as **CSV**, or the whole analysis payload as **JSON**,
from the buttons on the Target Priority panel — for your own downstream analysis.

</details>

<details>
<summary><b>Architecture & project layout</b></summary>

```
app.py                 Streamlit UI — wires the modules + priority panel together
cancer_tool/
  uniprot.py           gene symbol -> UniProt accession + sequence
  structures.py        fetch AlphaFold structure (PDB text)
  mutations.py         parse/validate mutations + fetch hotspots
  dynamics.py          ENM/NMA folding dynamics (flexibility, hinges, modes)
  pathogenicity.py     DeepMind AlphaMissense per-variant pathogenicity
  pockets.py           druggable-pocket detection (LIGSITE-style / fpocket)
  scoring.py           Target Priority Score (pure, unit-tested composite)
  metrics.py           pure-NumPy ranking metrics (AUROC/AUPRC/P@k) for the benchmark
  targets.py           Open Targets tractability + disease links
  viewer.py            3D py3Dmol view (pLDDT / flexibility / priority colouring)
  pipeline.py          analyze_gene() — shared engine for precompute + the API
api/main.py            FastAPI live-analysis service (/health, /analyze/{gene})
Dockerfile             container image for the live-analysis API
scripts/precompute.py  precompute engine -> data/{GENE}.json for the HTML
scripts/benchmark.py   quantitative benchmark vs a gold standard -> benchmark/
data/                  precomputed per-gene results (powers cancer-explorer.html)
benchmark/             gold_standard.json + committed results.json / REPORT.md
docs/METHODS.md        scoring formula, ENM parameters, pocket caveats, benchmark
requirements.lock      exact-pinned scientific stack for reproducibility
tests/                 network-free unit tests (scoring, dynamics, metrics, benchmark, ...)
```

To refresh the standalone HTML's precomputed data:

```bash
python scripts/precompute.py            # all featured genes
python scripts/precompute.py TP53 KRAS  # specific genes
```

</details>

<details>
<summary><b>Testing</b></summary>

```bash
pytest
```

The suite needs no network — the scoring formula, ENM/NMA dynamics (on a synthetic
fixture), AlphaMissense parsing, the ranking metrics, and a full-pipeline reproducibility
check (KRAS → G12 replayed from captured API responses) all run with no Streamlit/UI or HTTP
dependency. CI also lints with `ruff` and enforces a coverage floor across Python 3.11/3.12.

To (re)run the quantitative benchmark against live data:

```bash
pip install ".[benchmark]"          # optional: enables the data-fit weight comparison
python scripts/benchmark.py         # -> benchmark/results.json + REPORT.md
```

</details>

<details>
<summary><b>Limitations</b></summary>

Honest scope, so the rankings aren't over-read:

- **Single-chain only.** Analysis runs on the AlphaFold monomer; assembly/interface
  dynamics of multimers (e.g. the p53 tetramer) are not modelled.
- **Pockets are geometric.** A LIGSITE-style scan finds enclosed cavities by shape; it
  does not assess lipophilicity, desolvation, or true ligandability.
- **Confidence-aware, not confidence-proof.** Structural criticality is down-weighted by
  AlphaFold pLDDT, but low-confidence regions still carry more uncertainty than the score
  conveys.
- **Not clinically validated.** The ranking rediscovers known driver biology and is
  benchmarked against a curated gold standard (see [`benchmark/REPORT.md`](benchmark/REPORT.md)),
  but the composite weights are expert-set, not trained on outcome data — and on a
  recurrence-defined benchmark, recurrence alone is a strong baseline. Treat the score as a
  research/education heuristic, meaningful *within* one protein, not a validated clinical metric.

Full method — the scoring formula and weights, ENM/NMA parameters, and pocket-detection
caveats — is written up in [`docs/METHODS.md`](docs/METHODS.md).

</details>

<details>
<summary><b>Disclaimer</b></summary>

⚠️ **Research and education use only.** This is not a medical device. It must
not be used for diagnosis or treatment decisions.

</details>

---

<table align="center" width="100%">
<tr>
<td width="50%" align="center" valign="middle">
  <img src="VRcompany.png" width="92%" alt="Vedra Research" />
</td>
<td width="50%" align="center" valign="middle">
  <img src="vr01-helix.gif" width="92%" alt="VR-01 — animated double-helix visual: crimson and bone strands winding with base-pair rungs, evoking the protein-structure science at the core of the explorer" />
</td>
</tr>
</table>
