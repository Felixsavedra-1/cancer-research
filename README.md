# Cancer Research Tool

<p align="center">
  <img src="docs/demo.gif" alt="Cancer Protein Explorer — AlphaFold structure, dynamics, and ranked druggable targets" width="800">
</p>

> Advanced research tool that leverages Google DeepMind's **AlphaFold** to resolve protein folding dynamics and rank druggable cancer-driver residues — accelerating target discovery, on a laptop.

<p align="center">
  <img src="https://github.com/Felixsavedra-1/Cancer-Research-tool/actions/workflows/ci.yml/badge.svg" alt="Tests">
  <img src="https://img.shields.io/badge/Python-3.11-3776AB?logo=python&logoColor=white" alt="Python 3.11">
  <img src="https://img.shields.io/badge/Streamlit-app-FF4B4B?logo=streamlit&logoColor=white" alt="Streamlit">
  <img src="https://img.shields.io/badge/AlphaFold-powered-0B3D91" alt="AlphaFold-powered">
  <img src="https://img.shields.io/badge/AlphaMissense-pathogenicity-0B3D91" alt="AlphaMissense">
  <img src="https://img.shields.io/badge/use-research%20%2F%20education-informational" alt="Research / Education use">
</p>

Type a gene (e.g. `TP53`, `KRAS`, `BRAF`) and the tool fetches its real **AlphaFold** structure, computes its intrinsic **folding dynamics**, scores every variant with DeepMind's **AlphaMissense**, finds **druggable pockets**, and fuses it all into a ranked shortlist of the **most actionable cancer-driver residues** — each with a plain-English rationale. No GPU required.

> **Does the science hold up?** On TP53 the top three ranked residues come out **R273, R248, R175** — the textbook p53 drivers. KRAS tops at **G12**, BRAF at **V600**, EGFR at the **exon-19 deletion / L858 / T790** sites. The ranking rediscovers known biology from first principles.

## Highlights

- **Folding dynamics, no simulation.** Elastic Network Model normal-mode analysis (ProDy) extracts per-residue flexibility, hinge/domain pivots, and collective motion from a single AlphaFold structure in milliseconds — the laptop-grade route to the motions a GPU molecular-dynamics run would reveal.
- **Variant effect from DeepMind AlphaMissense.** Every residue carries a predicted pathogenicity score + class (likely benign → pathogenic), pulled from the same AlphaFold DB endpoint as the structure.
- **Druggable-pocket detection.** A LIGSITE-style geometric scan finds enclosed, ligandable cavities and which residues line them — the actionability axis of the score.
- **Target Priority Score — the discovery engine.** One explainable composite: `0.30·recurrence + 0.35·pathogenicity + 0.20·druggability + 0.15·structural criticality`, ranking residues as druggable cancer drivers. Enriched with **Open Targets** tractability and disease links.
- **Real AlphaFold structures, zero GPU.** Pulls finished models from DeepMind's ~200M-protein set; all heavy lifting runs on hosted services or in milliseconds on CPU.
- **Confidence-aware 3D, three colour modes.** `py3Dmol` rendering you can paint by pLDDT confidence, NMA flexibility, or the Target Priority Score.
- **Pure, unit-tested scientific core.** Dynamics, pathogenicity parsing, and the scoring formula carry no UI/network dependency and are covered by a network-free test suite.

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
main deliverable; click a residue's row to drop it onto the 3D structure.

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
  targets.py           Open Targets tractability + disease links
  viewer.py            3D py3Dmol view (pLDDT / flexibility / priority colouring)
scripts/precompute.py  precompute engine -> data/{GENE}.json for the HTML
data/                  precomputed per-gene results (powers cancer-explorer.html)
tests/                 network-free unit tests (scoring, dynamics, pathogenicity, ...)
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
fixture), and AlphaMissense parsing all run with no Streamlit/UI or HTTP dependency.

</details>

<details>
<summary><b>Roadmap</b></summary>

- Cohort frequencies & cancer-type filter via cBioPortal
- Known-drug overlay (which approved drugs hit each pocket)
- All-atom molecular dynamics (OpenMM) for the shortlisted residues
- AlphaMissense heatmap per residue (all 19 substitutions)

</details>

<details>
<summary><b>Disclaimer</b></summary>

⚠️ **Research and education use only.** This is not a medical device. It must
not be used for diagnosis or treatment decisions.

</details>

---

<p align="center">
  <img src="docs/VRcompany.png" alt="Vedra Research" width="50%">
</p>
