"""Cancer Protein Explorer — Streamlit app. Run with: streamlit run app.py"""

from __future__ import annotations

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components

from cancer_tool import (
    __version__,
    dynamics,
    mutations as mut,
    pathogenicity,
    pockets,
    scoring,
    structures,
    targets,
    uniprot,
    viewer,
)

EXAMPLE_GENES = ["TP53", "KRAS", "BRAF", "KIT", "FLT3", "IDH2", "PTEN", "EGFR"]

_BRAND_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,300;0,400;0,500;1,300;1,400&family=Hanken+Grotesk:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap');

:root {
  --black: #0A0A0B; --black-2: #111113; --panel: #141417;
  --line: #26262B; --line-soft: #19191D;
  --bone: #E8E3D6; --grey: #8C8C93; --grey-dim: #5C5C63;
  --gold: #C6A15B; --gold-dim: #6e5c34; --crimson: #B23A3A;
  --display: 'Spectral', Georgia, serif;
  --body: 'Hanken Grotesk', system-ui, -apple-system, sans-serif;
  --mono: 'IBM Plex Mono', ui-monospace, monospace;
}

/* Background grid overlay */
.stApp::before {
  content: ""; position: fixed; inset: 0; z-index: 0; pointer-events: none;
  background-image:
    linear-gradient(var(--line-soft) 1px, transparent 1px),
    linear-gradient(90deg, var(--line-soft) 1px, transparent 1px);
  background-size: 64px 64px; opacity: .45;
  -webkit-mask-image: radial-gradient(120% 120% at 70% 0%, #000 0%, transparent 75%);
          mask-image: radial-gradient(120% 120% at 70% 0%, #000 0%, transparent 75%);
}

html, body, .stApp, [data-testid="stSidebar"] { font-family: var(--body); }
h1, h2, h3 {
  font-family: var(--display) !important; font-weight: 300 !important;
  letter-spacing: -.012em;
}
[data-testid="stCaptionContainer"], .stCaption { font-family: var(--mono); color: var(--grey); }

/* Brand masthead */
.vr-masthead { margin: 0 0 4px; }
.vr-designation {
  font-family: var(--mono); font-size: .72rem; letter-spacing: .18em;
  text-transform: uppercase; color: var(--gold);
}
.vr-mark {
  font-family: var(--mono); font-size: .78rem; letter-spacing: .14em;
  color: var(--grey); float: right;
}
.vr-mark a { color: var(--grey); text-decoration: none; }
.vr-mark a:hover { color: var(--bone); }

/* Buttons */
.stButton button, .stDownloadButton button {
  background: var(--gold); color: var(--black);
  border: 1px solid var(--gold); border-radius: 2px;
  font-family: var(--mono); font-weight: 500; font-size: .8rem;
  letter-spacing: .12em; text-transform: uppercase;
}
.stButton button:hover, .stDownloadButton button:hover {
  background: var(--bone); border-color: var(--bone); color: var(--black);
}

/* Inputs */
[data-testid="stTextInput"] input,
[data-testid="stTextArea"] textarea,
[data-baseweb="select"] > div {
  background: var(--black-2) !important; color: var(--bone) !important;
  border: 1px solid var(--line) !important; border-radius: 2px !important;
}
[data-testid="stTextInput"] input:focus,
[data-testid="stTextArea"] textarea:focus,
[data-baseweb="select"] > div:focus-within {
  border-color: var(--gold) !important; box-shadow: none !important;
}

/* Alerts / disclaimer */
[data-testid="stAlert"] {
  background: rgba(178,58,58,.07) !important;
  border: 1px solid var(--crimson) !important; border-radius: 2px !important;
  color: #d98d8d !important;
}

/* Tables & links */
[data-testid="stTable"] th, [data-testid="stDataFrame"] th {
  font-family: var(--mono); text-transform: uppercase; letter-spacing: .12em;
  font-size: .66rem; color: var(--grey);
}
a, a:visited { color: var(--gold); }
</style>
"""


def inject_brand_css() -> None:
    """Apply the Vedra Research brand fonts and component styling over the dark theme."""
    st.markdown(_BRAND_CSS, unsafe_allow_html=True)
    st.markdown(
        '<div class="vr-masthead">'
        '<span class="vr-designation">VR-05</span>'
        '<span class="vr-mark">'
        '<a href="https://vedraresearch.github.io/Vedra-Research/">Vedra Research ↗</a>'
        "</span></div>",
        unsafe_allow_html=True,
    )


@st.cache_resource
def http_session() -> requests.Session:
    session = requests.Session()
    session.headers.update({"User-Agent": f"cancer-protein-explorer/{__version__}"})
    return session


@st.cache_data(show_spinner=False)
def load_protein(gene: str):
    return uniprot.get_protein(gene, session=http_session())


@st.cache_data(show_spinner=False)
def load_structure(accession: str) -> str:
    return structures.fetch_alphafold_pdb(accession, session=http_session())


@st.cache_data(show_spinner=False)
def load_hotspots(gene: str):
    return mut.fetch_hotspots(gene, session=http_session())


@st.cache_data(show_spinner=False)
def load_pathogenicity(accession: str):
    return pathogenicity.fetch_alphamissense(accession, session=http_session())


@st.cache_data(show_spinner=False)
def load_dynamics(pdb_text: str):
    try:
        return dynamics.compute_dynamics(pdb_text)
    except dynamics.DynamicsError:
        return None


@st.cache_data(show_spinner=False)
def load_pockets(pdb_text: str):
    try:
        return pockets.detect_pockets(pdb_text)
    except Exception:
        return []


@st.cache_data(show_spinner=False)
def load_targets(gene: str):
    return targets.fetch_target_context(gene, session=http_session())


st.set_page_config(page_title="Cancer Protein Explorer", layout="wide")
inject_brand_css()

st.title("🧬 Cancer Protein Explorer")
st.caption(
    "AlphaFold structures + protein dynamics + AlphaMissense pathogenicity, fused into "
    "a ranked shortlist of druggable cancer-driver residues."
)
st.warning(
    "**Research and education use only.** This is not a medical device and must "
    "not be used for diagnosis or treatment decisions.",
    icon="⚠️",
)

with st.sidebar:
    st.header("Protein")
    gene = st.selectbox(
        "Gene symbol",
        options=EXAMPLE_GENES,
        index=0,
        accept_new_options=True,
        help="Pick an example or type any human gene symbol (e.g. PTEN).",
    )
    st.header("Mutations")
    mutation_text = st.text_area(
        "Mutations to highlight",
        value="R175H",
        help="One-letter form like R175H. Separate multiple with commas or spaces.",
    )
    st.header("Display")
    color_mode = st.radio(
        "Colour the structure by",
        options=["plddt", "flexibility", "priority"],
        format_func={
            "plddt": "AlphaFold confidence (pLDDT)",
            "flexibility": "Dynamics — flexibility (NMA)",
            "priority": "Target Priority Score",
        }.get,
        help=(
            "pLDDT: model confidence (🔴 low → 🔵 high). "
            "Flexibility: ENM/NMA intrinsic motion (🔵 rigid core → 🔴 mobile loops). "
            "Priority: composite druggable-driver score painted on scored residues."
        ),
    )
    show_pocket = st.toggle(
        "Show top druggable pocket",
        value=False,
        help="Overlay a translucent gold surface on the most druggable detected pocket.",
    )

if not gene:
    st.info("Choose or type a gene symbol in the sidebar to begin.")
    st.stop()

gene = gene.strip().upper()

with st.spinner(f"Looking up {gene}…"):
    try:
        protein = load_protein(gene)
    except requests.RequestException as exc:
        st.error(f"Could not reach UniProt: {exc}")
        st.stop()

if not protein:
    st.error(f"No reviewed human protein found for gene **{gene}**. Check the symbol.")
    st.stop()

with st.spinner("Fetching AlphaFold structure…"):
    try:
        pdb_text = load_structure(protein["accession"])
    except requests.RequestException as exc:
        st.error(
            f"No AlphaFold structure available for {protein['accession']} ({exc})."
        )
        st.stop()

parsed, unparseable = mut.parse_mutations(mutation_text)
highlights, invalid = [], []
for mutation in parsed:
    ok, message = mut.validate_mutation(mutation, protein["sequence"])
    if ok:
        highlights.append({"position": mutation["position"], "label": mutation["label"]})
    else:
        invalid.append(f"{mutation['label']}: {message}")

# --- Run the analysis engine (each piece cached, each fails soft) -------------
with st.spinner("Analysing structure, dynamics and variant effects…"):
    try:
        hotspots = load_hotspots(gene)
    except requests.RequestException:
        hotspots = []
    pathos = load_pathogenicity(protein["accession"])
    dyn = load_dynamics(pdb_text)
    pocket_list = load_pockets(pdb_text)
    priority_rows = scoring.score_residues(
        hotspots, pathos, dyn, pocket_list, protein["sequence"]
    )

flex_map = dynamics.flexibility_by_position(dyn) if dyn else {}
max_score = max((r["score"] for r in priority_rows), default=0.0) or 1.0
priority_map = {r["position"]: r["score"] / max_score for r in priority_rows}
pocket_residues = pocket_list[0]["residues"] if (show_pocket and pocket_list) else None

viewer_col, info_col = st.columns([3, 2], gap="large")

with viewer_col:
    st.subheader(f"{protein['name']} ({gene})")
    html = viewer.render_structure(
        pdb_text,
        highlights=highlights,
        color_mode=color_mode,
        flexibility=flex_map,
        priority=priority_map,
        pocket_residues=pocket_residues,
    )
    components.html(html, height=620, scrolling=False)
    _legend = {
        "plddt": "Coloured by pLDDT: 🔴 low confidence → 🔵 high confidence.",
        "flexibility": "Coloured by ENM/NMA flexibility: 🔵 rigid core → 🔴 mobile.",
        "priority": "Coloured by Target Priority Score: dark → 🟡 gold → 🔴 highest.",
    }
    st.caption(_legend[color_mode])

with info_col:
    st.subheader("Protein")
    st.markdown(
        f"- **Gene:** {gene}\n"
        f"- **UniProt:** [{protein['accession']}]"
        f"(https://www.uniprot.org/uniprotkb/{protein['accession']})\n"
        f"- **Length:** {protein['length']} residues\n"
        f"- **Structure:** [AlphaFold DB]"
        f"(https://alphafold.ebi.ac.uk/entry/{protein['accession']})"
    )

    if unparseable:
        st.warning("Could not parse: " + ", ".join(unparseable))
    if invalid:
        for problem in invalid:
            st.error(problem)
    if highlights:
        st.success("Highlighted: " + ", ".join(h["label"] for h in highlights))

    st.subheader("Folding dynamics")
    if dyn:
        rank = sorted(
            zip(dyn["residue_numbers"], dyn["flexibility"]), key=lambda x: x[1]
        )
        rigid = ", ".join(str(r) for r, _ in rank[:5])
        mobile = ", ".join(str(r) for r, _ in rank[-5:][::-1])
        hinge_txt = ", ".join(str(h) for h in dyn["hinges"][:8]) or "none detected"
        st.markdown(
            f"- **Most rigid (core):** {rigid}\n"
            f"- **Most mobile:** {mobile}\n"
            f"- **Hinge / domain pivots:** {hinge_txt}"
        )
        st.caption(
            f"Elastic Network Model normal-mode analysis ({dyn['n_modes']} modes) — "
            "intrinsic motion inferred from one AlphaFold structure, no simulation."
        )
    else:
        st.caption("Dynamics analysis unavailable for this structure.")

    st.subheader("Druggability")
    if pocket_list:
        top = pocket_list[0]
        st.markdown(
            f"- **Pockets detected:** {len(pocket_list)} "
            f"({top['source']})\n"
            f"- **Top pocket:** {top['volume']} Å³, druggability "
            f"{top['druggability']:.2f}, {len(top['residues'])} lining residues"
        )
    else:
        st.caption("No enclosed pocket detected.")

    context = load_targets(gene)
    if context:
        modalities = sorted({t["modality"] for t in context["tractability"]})
        if modalities:
            st.markdown("**Tractable modalities:** " + ", ".join(modalities))
        if context["top_diseases"]:
            diseases = " · ".join(
                f"{d['name']} ({d['score']:.2f})" for d in context["top_diseases"][:3]
            )
            st.markdown(f"**Top disease links:** {diseases}")
        st.caption("Target context from [Open Targets](https://platform.opentargets.org).")

# --- Target Priority panel (the discovery deliverable) -----------------------
st.divider()
st.subheader("🎯 Target Priority — ranked druggable cancer-driver residues")
if priority_rows:
    df = pd.DataFrame(priority_rows)[
        ["residue", "score", "recurrence", "pathogenicity", "druggability",
         "criticality", "tumours", "rationale"]
    ]
    df.columns = [
        "Residue", "Priority", "Recurrence", "Pathogenicity", "Druggability",
        "Criticality", "Tumours", "Why it ranks here",
    ]
    st.dataframe(
        df,
        hide_index=True,
        use_container_width=True,
        height=420,
        column_config={
            "Priority": st.column_config.ProgressColumn(
                "Priority", min_value=0, max_value=100, format="%.0f"
            ),
            "Recurrence": st.column_config.ProgressColumn(
                "Recurrence", min_value=0, max_value=1, format="%.2f"
            ),
            "Pathogenicity": st.column_config.ProgressColumn(
                "Pathogenicity", min_value=0, max_value=1, format="%.2f"
            ),
            "Druggability": st.column_config.ProgressColumn(
                "Druggability", min_value=0, max_value=1, format="%.2f"
            ),
            "Criticality": st.column_config.ProgressColumn(
                "Criticality", min_value=0, max_value=1, format="%.2f"
            ),
        },
    )
    st.caption(
        "Priority = 0.30·recurrence + 0.35·pathogenicity + 0.20·druggability + "
        "0.15·structural criticality. Recurrence: "
        "[cancerhotspots.org](https://www.cancerhotspots.org). Pathogenicity: "
        "DeepMind [AlphaMissense](https://alphafold.ebi.ac.uk). Druggability: "
        "geometric pocket detection. Criticality: ENM/NMA. "
        "Copy a residue (e.g. R175) into the sidebar to highlight it in 3D."
    )
else:
    st.caption("No recurrent cancer-mutation sites on record for this gene to rank.")
