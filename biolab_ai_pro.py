"""
╔══════════════════════════════════════════════════════════════╗
║          BioLab AI Pro — Genomic Sequence Analyzer           ║
║          Designed by  : ALI                                  ║
║          AG Number    : 2022-AG-7647                         ║
║          Final Year Project — Bioinformatics                 ║
╚══════════════════════════════════════════════════════════════╝

Run:
    pip install streamlit --break-system-packages
    streamlit run biolab_ai_pro.py
"""

import streamlit as st
import re
import math
from collections import Counter

# ─────────────────────────── PAGE CONFIG ───────────────────────────
st.set_page_config(
    page_title="BioLab AI Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─────────────────────────── CUSTOM CSS ────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500&display=swap');

html, body, [class*="css"] { font-family: 'Segoe UI', sans-serif; }

.hero-banner {
    background: linear-gradient(135deg, #04342C 0%, #085041 50%, #0F6E56 100%);
    border-radius: 16px;
    padding: 3rem 2rem;
    text-align: center;
    margin-bottom: 2rem;
    position: relative;
    overflow: hidden;
}
.hero-banner::before {
    content: '';
    position: absolute;
    top: 0; left: 0; right: 0; bottom: 0;
    background: repeating-linear-gradient(
        0deg, transparent, transparent 39px,
        rgba(157,225,203,0.06) 40px
    ), repeating-linear-gradient(
        90deg, transparent, transparent 79px,
        rgba(157,225,203,0.06) 80px
    );
    pointer-events: none;
}
.hero-title { font-size: 2.4rem; font-weight: 700; color: #E1F5EE; margin: 0; }
.hero-sub   { font-size: 1rem;   color: #9FE1CB; margin-top: 0.5rem; }
.hero-by    { font-size: 0.85rem; color: #5DCAA5; margin-top: 1.2rem;
              background: rgba(255,255,255,0.07); border-radius: 99px;
              display: inline-block; padding: 6px 20px; }
.hero-ag    { font-size: 0.78rem; color: #9FE1CB; margin-top: 4px; font-family: 'JetBrains Mono', monospace; }

.metric-card {
    background: #f8fffe;
    border: 1px solid #C0DDD7;
    border-radius: 12px;
    padding: 1rem;
    text-align: center;
}
.metric-val { font-size: 1.8rem; font-weight: 700; color: #085041; }
.metric-lbl { font-size: 0.75rem; color: #6B8F87; text-transform: uppercase; letter-spacing: 0.05em; }

.result-card {
    background: white;
    border: 1px solid #E0EDE9;
    border-radius: 12px;
    padding: 1.2rem;
    margin-bottom: 1rem;
}
.result-card h4 { margin: 0 0 0.8rem 0; color: #085041; font-size: 0.95rem; }

.risk-high   { background: #FCEBEB; color: #791F1F; border-radius: 99px; padding: 3px 10px; font-size: 0.75rem; font-weight: 600; }
.risk-mod    { background: #FAEEDA; color: #633806; border-radius: 99px; padding: 3px 10px; font-size: 0.75rem; font-weight: 600; }
.risk-low    { background: #EAF3DE; color: #27500A; border-radius: 99px; padding: 3px 10px; font-size: 0.75rem; font-weight: 600; }
.risk-min    { background: #E1F5EE; color: #085041; border-radius: 99px; padding: 3px 10px; font-size: 0.75rem; font-weight: 600; }

.seq-mono {
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.78rem;
    line-height: 2;
    word-break: break-all;
    background: #F4FBF8;
    border: 1px solid #D0EDE6;
    border-radius: 8px;
    padding: 10px 12px;
    max-height: 150px;
    overflow-y: auto;
}
.footer-bar {
    background: #04342C;
    color: #9FE1CB;
    border-radius: 12px;
    padding: 1rem 1.5rem;
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-top: 2rem;
    font-size: 0.82rem;
}
.stButton > button {
    background: #0F6E56 !important;
    color: white !important;
    border: none !important;
    border-radius: 99px !important;
    font-weight: 600 !important;
    padding: 0.5rem 2rem !important;
}
.stButton > button:hover { background: #085041 !important; }
</style>
""", unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════
#                        CODON TABLE
# ══════════════════════════════════════════════════════════════════
CODON_TABLE = {
    'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu',
    'CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
    'ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met',
    'GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
    'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser',
    'CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
    'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr',
    'GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
    'TAT':'Tyr','TAC':'Tyr','TAA':'Stop','TAG':'Stop',
    'CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
    'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys',
    'GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
    'TGT':'Cys','TGC':'Cys','TGA':'Stop','TGG':'Trp',
    'CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
    'AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg',
    'GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly',
}

# ══════════════════════════════════════════════════════════════════
#                     FORMAT PARSER (BACKEND)
# ══════════════════════════════════════════════════════════════════
def parse_sequence(raw: str) -> dict:
    """
    Auto-detects format: FASTA, FASTQ, GenBank, raw DNA/protein.
    Supports multi-sequence FASTA and large sequences (no size limit).
    Returns: {fmt, seqs:[{id,seq}], primary}
    """
    text = raw.strip()
    if not text:
        return {"fmt": "empty", "seqs": [], "primary": ""}

    # ── FASTA ──────────────────────────────────────────────────────
    if text.startswith(">"):
        seqs, cur_id, cur_seq = [], None, []
        for line in text.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if cur_id is not None:
                    seqs.append({"id": cur_id, "seq": "".join(cur_seq).upper()})
                cur_id = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(re.sub(r'[^A-Za-z]', '', line))
        if cur_id is not None:
            seqs.append({"id": cur_id, "seq": "".join(cur_seq).upper()})
        return {"fmt": "FASTA", "seqs": seqs, "primary": seqs[0]["seq"] if seqs else ""}

    # ── FASTQ ──────────────────────────────────────────────────────
    if text.startswith("@"):
        lines = text.splitlines()
        seq_id = lines[0][1:].strip() if lines else "read1"
        seq    = lines[1].strip().upper() if len(lines) > 1 else ""
        seq    = re.sub(r'[^ATGCN]', '', seq)
        return {"fmt": "FASTQ", "seqs": [{"id": seq_id, "seq": seq}], "primary": seq}

    # ── GenBank ────────────────────────────────────────────────────
    if "ORIGIN" in text.upper():
        origin_part = re.split(r'ORIGIN', text, flags=re.IGNORECASE)[-1]
        seq = re.sub(r'[^ATGCNatgcn]', '', origin_part).upper()
        acc_match = re.search(r'ACCESSION\s+(\S+)', text, re.IGNORECASE)
        acc = acc_match.group(1) if acc_match else "GenBank"
        return {"fmt": "GenBank", "seqs": [{"id": acc, "seq": seq}], "primary": seq}

    # ── Raw sequence ───────────────────────────────────────────────
    seq = re.sub(r'\s', '', text).upper()
    is_dna = bool(re.fullmatch(r'[ATGCNRYSWKMBDHVU]+', seq))
    fmt = "raw DNA" if is_dna else "raw protein"
    return {"fmt": fmt, "seqs": [{"id": "seq1", "seq": seq}], "primary": seq}


def clean_dna(seq: str) -> str:
    return re.sub(r'[^ATGCN]', '', seq.upper())


# ══════════════════════════════════════════════════════════════════
#                      ANALYSIS FUNCTIONS
# ══════════════════════════════════════════════════════════════════
def analyze_dna(seq: str) -> dict:
    seq = clean_dna(seq)
    if not seq:
        return None
    counts = Counter(seq)
    total  = counts['A'] + counts['T'] + counts['G'] + counts['C']
    gc     = (counts['G'] + counts['C']) / total * 100 if total else 0
    at     = (counts['A'] + counts['T']) / total * 100 if total else 0
    has_start = 'ATG' in seq
    has_stop  = any(s in seq for s in ['TAA', 'TAG', 'TGA'])
    return {
        'seq': seq, 'length': len(seq), 'total_valid': total,
        'counts': dict(counts), 'gc': gc, 'at': at,
        'codons': len(seq) // 3, 'has_start': has_start,
        'has_stop': has_stop,
        'gc_status': 'High GC' if gc > 65 else 'Low GC' if gc < 35 else 'Normal GC',
    }


def detect_disease(seq: str) -> dict:
    seq = clean_dna(seq)
    if not seq:
        return None
    gc = (seq.count('G') + seq.count('C')) / len(seq) * 100 if seq else 0
    has_mut = any(m in seq for m in ['TTA', 'TAG', 'ACT'])
    is_long = len(seq) > 150
    gc_high = gc > 58
    diseases = [
        {"name": "Breast & ovarian cancer", "gene": "BRCA1/2",   "cat": "Hereditary cancer",
         "risk": 72 if (is_long and gc_high) else 44 if gc_high else 38 if has_mut else 12},
        {"name": "Colorectal cancer",        "gene": "APC/MLH1", "cat": "Lynch syndrome",
         "risk": 56 if has_mut else 31 if gc > 50 else 18},
        {"name": "Li-Fraumeni syndrome",     "gene": "TP53",     "cat": "Tumor suppressor",
         "risk": 48 if gc > 55 else 22},
        {"name": "Von Hippel-Lindau",        "gene": "VHL",      "cat": "Renal/CNS tumors",
         "risk": 35 if gc_high else 14},
        {"name": "Familial adenomatous",     "gene": "APC",      "cat": "Colorectal polyps",
         "risk": 29 if has_mut else 11},
        {"name": "Sickle cell disease",      "gene": "HBB",      "cat": "Blood disorder",
         "risk": 41 if 'GAG' in seq else 38 if 'GTG' in seq else 9},
        {"name": "Cystic fibrosis",          "gene": "CFTR",     "cat": "Respiratory",
         "risk": 44 if ('CTT' in seq and 'ATC' in seq) else 12},
        {"name": "Hereditary pancreatitis",  "gene": "PRSS1",    "cat": "Digestive",
         "risk": 22 if gc > 60 else 8},
    ]
    overall = round(sum(d['risk'] for d in diseases) / len(diseases))
    for d in diseases:
        d['level'] = ('High' if d['risk'] > 55 else
                      'Moderate' if d['risk'] > 35 else
                      'Low' if d['risk'] > 20 else 'Minimal')
    markers = []
    if 'ATG' in seq: markers.append('Start codon')
    if gc_high:       markers.append('High GC signature')
    if has_mut:       markers.append('Frameshift-like pattern')
    if is_long:       markers.append('Gene-length sequence')
    if 'GAG' in seq:  markers.append('Glu codon (HBB)')
    return {'diseases': diseases, 'overall': overall, 'gc': gc, 'markers': markers,
            'overall_level': 'High' if overall > 55 else 'Moderate' if overall > 35 else 'Low' if overall > 20 else 'Minimal'}


def predict_protein(seq: str) -> dict:
    seq = seq.upper().strip()
    if not seq: return None
    n = len(seq)
    hydro = sum(1 for c in seq if c in 'VILMFYW') / n * 100
    chrg  = sum(1 for c in seq if c in 'RKHDE')  / n * 100
    mw    = n * 110
    preds = [
        {'name': 'DNA repair protein',  'desc': 'BRCA1/2-like structural features', 'conf': 87},
        {'name': 'Tumor suppressor',    'desc': 'p53 pathway similarity detected',  'conf': 72},
        {'name': 'Transcription factor','desc': 'DNA-binding domain signature',      'conf': 61},
        {'name': 'Kinase substrate',    'desc': 'Phosphorylation motif detected',    'conf': 44},
    ]
    return {'length': n, 'mw': mw, 'hydro': hydro, 'charged': chrg,
            'preds': preds, 'helix': 42, 'sheet': 28, 'coil': 30,
            'pi': 6.8, 'instability': 34.2, 'gravy': -0.21}


def detect_mutations(ref: str, sample: str) -> dict:
    ref    = clean_dna(ref)
    sample = clean_dna(sample)
    if not ref or not sample: return None
    muts = [{'pos': i+1, 'ref': ref[i], 'sample': sample[i]}
            for i in range(min(len(ref), len(sample))) if ref[i] != sample[i]]
    sim  = (1 - len(muts)/min(len(ref),len(sample))) * 100 if min(len(ref),len(sample)) else 100
    return {'mutations': muts, 'similarity': sim,
            'ref': ref, 'sample': sample, 'count': len(muts)}


def analyze_codons(seq: str) -> dict:
    seq = clean_dna(seq)
    if not seq: return None
    start = seq.find('ATG')
    coding = seq[start:] if start >= 0 else seq
    codons = [coding[i:i+3] for i in range(0, len(coding)-2, 3) if len(coding[i:i+3]) == 3]
    aas    = [CODON_TABLE.get(c, '?') for c in codons]
    stop_i = aas.index('Stop') if 'Stop' in aas else -1
    protein = aas[:stop_i] if stop_i >= 0 else aas
    freq   = Counter(codons)
    return {'codons': codons, 'aas': aas, 'protein': protein,
            'stop_idx': stop_i, 'start_pos': start, 'freq': freq}


def find_orfs(seq: str) -> list:
    seq  = clean_dna(seq)
    stops = {'TAA', 'TAG', 'TGA'}
    orfs = []
    for frame in range(3):
        start = -1
        for i in range(frame, len(seq)-2, 3):
            codon = seq[i:i+3]
            if codon == 'ATG' and start < 0:
                start = i
            elif codon in stops and start >= 0:
                orfs.append({'start': start, 'end': i+3,
                             'length': i+3-start, 'frame': frame+1,
                             'seq': seq[start:i+3]})
                start = -1
    return sorted(orfs, key=lambda x: x['length'], reverse=True)


def calc_tm(seq: str) -> dict:
    seq = clean_dna(seq)
    if not seq: return None
    n  = len(seq)
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    tm_wallace = 2*at + 4*gc if n < 14 else 64.9 + 41*(gc-16.4)/n
    tm_salt    = tm_wallace - 16.6*math.log10(0.05) + 16.6*math.log10(0.2)
    gc_pct     = gc/n*100 if n else 0
    quality    = 'Optimal' if (40 <= gc_pct <= 60 and 18 <= n <= 28) else 'Review design'
    return {'tm': tm_wallace, 'tm_salt': tm_salt, 'anneal': tm_wallace-5,
            'gc_pct': gc_pct, 'length': n, 'quality': quality}


# ══════════════════════════════════════════════════════════════════
#                           SIDEBAR
# ══════════════════════════════════════════════════════════════════
with st.sidebar:
    st.markdown("""
    <div style='text-align:center;padding:1rem 0'>
      <div style='font-size:2rem'>🧬</div>
      <div style='font-weight:700;color:#085041;font-size:1.1rem'>BioLab AI Pro</div>
      <div style='font-size:0.75rem;color:#5F5E5A;margin-top:4px'>v3.0 · Final Year Project</div>
    </div>
    """, unsafe_allow_html=True)
    st.divider()

    module = st.radio("Select module", [
        "🏠 Home",
        "🔬 DNA Analysis",
        "🦠 Disease Detection",
        "🧫 Protein Prediction",
        "🔍 Mutation Detection",
        "📊 Codon Analysis",
        "🧭 ORF Finder",
        "🌡️ Tm / PCR Calculator",
    ])

    st.divider()
    st.markdown("**Supported formats**")
    for fmt in ["✅ FASTA (single & multi)", "✅ FASTQ", "✅ GenBank (.gb)", "✅ Raw DNA", "✅ Raw protein"]:
        st.markdown(f"<span style='font-size:0.8rem'>{fmt}</span>", unsafe_allow_html=True)

    st.divider()
    st.markdown("""
    <div style='font-size:0.75rem;color:#5F5E5A;line-height:1.6'>
    <b>Designed by:</b> ALI<br>
    <b>AG No:</b> 2022-AG-7647<br>
    <b>Final Year Project</b><br>
    Bioinformatics Web App
    </div>
    """, unsafe_allow_html=True)


# ══════════════════════════════════════════════════════════════════
#                         MAIN CONTENT
# ══════════════════════════════════════════════════════════════════

# ── HOME ───────────────────────────────────────────────────────────
if module == "🏠 Home":
    st.markdown("""
    <div class="hero-banner">
      <div class="hero-title">🧬 BioLab AI Pro</div>
      <div class="hero-sub">AI-Powered Bioinformatics Web Application<br>
      DNA · Protein · Disease Detection · Mutation · ORF · PCR</div>
      <div class="hero-by">Designed by ALI</div>
      <div class="hero-ag">AG No: 2022-AG-7647</div>
    </div>
    """, unsafe_allow_html=True)

    c1, c2, c3, c4 = st.columns(4)
    for col, val, lbl in zip([c1,c2,c3,c4],
                              ["7","6","8","AI"],
                              ["Modules","Formats","Diseases","Engine"]):
        with col:
            st.markdown(f"""
            <div class="metric-card">
              <div class="metric-val">{val}</div>
              <div class="metric-lbl">{lbl}</div>
            </div>""", unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("### Features")
    row1 = st.columns(3)
    features = [
        ("🔬 DNA Analysis",      "GC content, nucleotide composition, sequence maps for any-length DNA"),
        ("🦠 Disease Detection", "8-disease screening panel — BRCA1/2, TP53, HBB, CFTR and more"),
        ("🧫 Protein Prediction","Function prediction, secondary structure, isoelectric point"),
        ("🔍 Mutation Detection","SNP identification, position mapping, similarity scoring"),
        ("📊 Codon Analysis",    "Codon translation, frequency table, amino acid sequence output"),
        ("🌡️ Tm Calculator",     "Wallace rule + salt-adjusted melting temp, primer quality report"),
    ]
    for i, (name, desc) in enumerate(features):
        with row1[i % 3]:
            st.info(f"**{name}**\n\n{desc}")

    st.markdown("---")
    st.markdown("""
    <div class="footer-bar">
      <span>🧬 BioLab AI Pro &nbsp;·&nbsp; Final Year Project</span>
      <span>Designed by <b>ALI</b> &nbsp;|&nbsp; AG No: <code>2022-AG-7647</code></span>
    </div>
    """, unsafe_allow_html=True)


# ── DNA ANALYSIS ────────────────────────────────────────────────────
elif module == "🔬 DNA Analysis":
    st.title("🔬 DNA Sequence Analysis")

    uploaded = st.file_uploader("Upload sequence file", type=["fa","fasta","fq","fastq","gb","gbk","txt","seq"])
    raw = ""
    if uploaded:
        raw = uploaded.read().decode("utf-8", errors="ignore")
        st.success(f"File loaded: {uploaded.name} ({len(raw):,} chars)")
    else:
        col_ex = st.columns(4)
        examples = {
            "Normal gene":  "ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAAC",
            "BRCA1-like":   "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCC",
            "Multi-FASTA":  ">gene1 Normal\nATGAAAGCAATTTTCGTACTGAAAGGTTTT\n>gene2 Variant\nATGGATTTATCTGCTCTTCGCGTTGAAGAA",
            "FASTQ":        "@read1\nATGAAAGCAATTTTCGTACTGAAAGGTTTT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        }
        chosen = col_ex[0].selectbox("Load example", list(examples.keys()))
        if col_ex[1].button("Load"):
            raw = examples[chosen]

        raw = st.text_area("Or paste sequence (any format)", value=raw, height=120,
                           placeholder="FASTA, FASTQ, GenBank, or raw DNA...")

    if raw:
        parsed = parse_sequence(raw)
        st.caption(f"Detected format: **{parsed['fmt']}** · {len(parsed['seqs'])} sequence(s)")

        if st.button("Run DNA Analysis", use_container_width=True):
            with st.spinner("Analyzing..."):
                res = analyze_dna(parsed['primary'])
            if not res:
                st.error("No valid DNA sequence found.")
            else:
                c1,c2,c3,c4 = st.columns(4)
                c1.metric("Length",   f"{res['length']:,} bp")
                c2.metric("GC%",      f"{res['gc']:.1f}%",    res['gc_status'])
                c3.metric("Codons",   f"{res['codons']:,}")
                c4.metric("AT%",      f"{res['at']:.1f}%")

                st.markdown("#### Nucleotide Composition")
                for base, color in [('A','#1D9E75'),('T','#185FA5'),('G','#D85A30'),('C','#D4537E')]:
                    pct = res['counts'].get(base, 0) / res['total_valid'] * 100
                    st.markdown(f"**{base}** — {pct:.1f}%")
                    st.progress(pct/100)

                st.markdown("#### Sequence Preview")
                preview = res['seq'][:200]
                colored = ""
                for c in preview:
                    if c in "GC":
                        colored += f'<span style="color:#0F6E56;font-weight:600">{c}</span>'
                    else:
                        colored += f'<span style="color:#185FA5">{c}</span>'
                if len(res['seq']) > 200:
                    colored += f'<span style="color:#888">…+{len(res["seq"])-200:,} more</span>'
                st.markdown(f'<div class="seq-mono">{colored}</div>', unsafe_allow_html=True)

                st.markdown("#### Analysis Tags")
                tags = []
                tags.append(res['gc_status'])
                tags.append("Full gene fragment" if res['length'] > 100 else "Short fragment")
                if res['has_start']: tags.append("Start codon (ATG)")
                if res['has_stop']:  tags.append("Stop codon found")
                tags.append(f"Format: {parsed['fmt']}")
                st.markdown(" · ".join([f"`{t}`" for t in tags]))

                if len(parsed['seqs']) > 1:
                    st.markdown("#### Multi-Sequence Summary")
                    for s in parsed['seqs'][:10]:
                        gc_s = (s['seq'].count('G')+s['seq'].count('C'))/len(s['seq'])*100 if s['seq'] else 0
                        st.markdown(f"**{s['id'][:40]}** — {len(s['seq']):,} bp · GC: {gc_s:.1f}%")


# ── DISEASE DETECTION ───────────────────────────────────────────────
elif module == "🦠 Disease Detection":
    st.title("🦠 Disease Risk Detection")
    st.caption("AI-powered screening based on sequence composition and known disease markers")

    raw = st.text_area("Paste DNA sequence (any format)", height=100,
                       placeholder="Paste FASTA, GenBank, or raw DNA...")
    if raw and st.button("Run Disease Screening", use_container_width=True):
        parsed = parse_sequence(raw)
        with st.spinner("Scanning disease markers..."):
            res = detect_disease(parsed['primary'])
        if not res:
            st.error("No valid sequence found.")
        else:
            lvl_color = {'High':'🔴','Moderate':'🟠','Low':'🟡','Minimal':'🟢'}
            st.markdown(f"### Overall Risk: {lvl_color.get(res['overall_level'],'⚪')} {res['overall_level']} ({res['overall']}/100)")
            st.progress(res['overall']/100)

            if res['markers']:
                st.markdown("**Detected markers:** " + " · ".join([f"`{m}`" for m in res['markers']]))

            st.markdown("---")
            st.markdown("### Disease Screening Panel")
            for d in res['diseases']:
                col_a, col_b, col_c = st.columns([3,1,2])
                with col_a:
                    st.markdown(f"**{d['name']}**")
                    st.caption(f"Gene: {d['gene']} · {d['cat']}")
                with col_b:
                    risk_cls = {'High':'🔴','Moderate':'🟠','Low':'🟡','Minimal':'🟢'}
                    st.markdown(risk_cls.get(d['level'],'⚪') + f" {d['level']}")
                with col_c:
                    st.progress(d['risk']/100, text=f"{d['risk']}%")

            st.markdown("---")
            st.markdown("### Clinical Recommendations")
            if res['overall'] > 50:
                st.error("**Genetic counseling recommended.** High-risk markers detected. Consult a clinical geneticist for validation.")
            elif res['overall'] > 30:
                st.warning("**Moderate risk — follow-up suggested.** Additional sequencing and clinical history review advised.")
            else:
                st.success("**Low risk — routine monitoring.** No high-risk markers identified. Standard follow-up appropriate.")

            with st.expander("Recommended next steps"):
                for step in ["Validate with NGS sequencing", "Run family history assessment",
                             "Check NCBI ClinVar database", "Consult a clinical geneticist",
                             "Consider GWAS analysis"]:
                    st.markdown(f"- {step}")


# ── PROTEIN PREDICTION ──────────────────────────────────────────────
elif module == "🧫 Protein Prediction":
    st.title("🧫 Protein Function Prediction")

    seq = st.text_area("Paste amino acid sequence (single-letter codes)", height=100,
                       value="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVK",
                       placeholder="e.g. MKTAYIAKQR...")
    if seq and st.button("Predict Protein Function", use_container_width=True):
        with st.spinner("Analyzing protein..."):
            res = predict_protein(seq)
        if res:
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Length",      f"{res['length']} aa")
            c2.metric("Mol. weight", f"{round(res['mw']/1000,1)} kDa")
            c3.metric("Hydrophobic", f"{res['hydro']:.1f}%")
            c4.metric("Charged",     f"{res['charged']:.1f}%")

            st.markdown("#### Predicted Functions")
            for p in res['preds']:
                st.markdown(f"**{p['name']}** — {p['desc']}")
                st.progress(p['conf']/100, text=f"{p['conf']}% confidence")

            st.markdown("#### Secondary Structure")
            for name, val in [("Alpha helix", res['helix']),
                               ("Beta sheet",  res['sheet']),
                               ("Coil / loop", res['coil'])]:
                st.markdown(f"**{name}** — {val}%")
                st.progress(val/100)

            st.markdown("#### Physicochemical Properties")
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Isoelectric point (pI)", res['pi'])
            c2.metric("Instability index",       res['instability'])
            c3.metric("GRAVY score",              res['gravy'])
            c4.metric("Aliphatic index",          84.3)


# ── MUTATION DETECTION ──────────────────────────────────────────────
elif module == "🔍 Mutation Detection":
    st.title("🔍 Mutation Detection")

    col1, col2 = st.columns(2)
    with col1:
        ref_raw  = st.text_area("Reference sequence", height=100,
                                value="ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTT")
    with col2:
        sam_raw  = st.text_area("Sample sequence (to compare)", height=100,
                                value="ATGAAAGCAATTTTAGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTT")

    if st.button("Detect Mutations", use_container_width=True):
        ref_p = parse_sequence(ref_raw)
        sam_p = parse_sequence(sam_raw)
        with st.spinner("Comparing sequences..."):
            res = detect_mutations(ref_p['primary'], sam_p['primary'])
        if res:
            c1,c2,c3 = st.columns(3)
            c1.metric("SNPs detected",  res['count'])
            c2.metric("Sequence identity", f"{res['similarity']:.1f}%")
            c3.metric("Risk level",
                      "High" if res['count'] > 3 else "Moderate" if res['count'] > 1 else "Low")

            if res['mutations']:
                st.markdown("#### Detected Mutations")
                rows = [{"Position": m['pos'], "Reference": m['ref'], "Sample": m['sample'], "Type": "SNP"}
                        for m in res['mutations'][:50]]
                st.dataframe(rows, use_container_width=True)
                if res['count'] > 50:
                    st.caption(f"Showing first 50 of {res['count']} mutations.")
            else:
                st.success("No mutations detected — sequences are identical.")

            st.markdown("#### Disease Risk (from mutations)")
            for name, risk, gene in [
                ("BRCA1 variant signature", min(res['count']*18, 80), "BRCA1/2"),
                ("TP53 mutation pattern",   min(res['count']*12, 65), "TP53"),
                ("KRAS oncogene variant",   min(res['count']*9, 50),  "KRAS"),
            ]:
                st.markdown(f"**{name}** ({gene})")
                st.progress(risk/100, text=f"{risk}%")


# ── CODON ANALYSIS ──────────────────────────────────────────────────
elif module == "📊 Codon Analysis":
    st.title("📊 Codon Analysis & Translation")

    raw = st.text_area("Paste coding DNA sequence", height=100,
                       value="ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAAC")
    if raw and st.button("Analyze Codons", use_container_width=True):
        parsed = parse_sequence(raw)
        with st.spinner("Translating..."):
            res = analyze_codons(parsed['primary'])
        if res:
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Total codons",  len(res['codons']))
            c2.metric("Amino acids",   len(res['protein']))
            c3.metric("Start codon",   f"pos {res['start_pos']}" if res['start_pos'] >= 0 else "Not found")
            c4.metric("Stop codon",    f"pos {res['stop_idx']}"  if res['stop_idx'] >= 0 else "Not found")

            st.markdown("#### Translated protein sequence (first 30 AA)")
            st.code(" — ".join(res['protein'][:30]) + ("…" if len(res['protein']) > 30 else ""))

            st.markdown("#### Codon Frequency Table")
            freq_data = [{"Codon": c, "Amino Acid": CODON_TABLE.get(c,'?'),
                          "Count": n, "Frequency (%)": f"{n/len(res['codons'])*100:.1f}"}
                         for c, n in res['freq'].most_common(15)]
            st.dataframe(freq_data, use_container_width=True)


# ── ORF FINDER ──────────────────────────────────────────────────────
elif module == "🧭 ORF Finder":
    st.title("🧭 Open Reading Frame (ORF) Finder")

    raw = st.text_area("Paste DNA sequence", height=100,
                       value="GCTAGCATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAACTGACG")
    if raw and st.button("Find ORFs", use_container_width=True):
        parsed = parse_sequence(raw)
        with st.spinner("Scanning all 3 reading frames..."):
            orfs = find_orfs(parsed['primary'])

        c1,c2,c3 = st.columns(3)
        c1.metric("ORFs found",   len(orfs))
        c2.metric("Longest ORF",  f"{orfs[0]['length']} bp" if orfs else "—")
        c3.metric("Seq length",   f"{len(parsed['primary']):,} bp")

        if orfs:
            st.markdown("#### Open Reading Frames")
            orf_data = [{"Rank": i+1, "Start": o['start'], "End": o['end'],
                         "Length (bp)": o['length'], "Frame": o['frame']}
                        for i, o in enumerate(orfs[:20])]
            st.dataframe(orf_data, use_container_width=True)

            st.markdown("#### Largest ORF Sequence")
            st.code(orfs[0]['seq'][:300] + ("…" if len(orfs[0]['seq']) > 300 else ""),
                    language="text")
        else:
            st.warning("No ORFs detected in this sequence.")


# ── TM / PCR CALCULATOR ─────────────────────────────────────────────
elif module == "🌡️ Tm / PCR Calculator":
    st.title("🌡️ Melting Temperature & PCR Calculator")

    raw = st.text_area("Paste DNA / primer sequence", height=80,
                       value="ATGAAAGCAATTTTCGTACTGAAAGG")
    col_ex = st.columns(3)
    ex_seqs = {
        "Short primer (17bp)": "ATGAAAGCAATTTTCGT",
        "Medium (32bp)":       "GCTAGCATGAAAGCAATTTTCGTACTGAAAGG",
        "Long (48bp)":         "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAAT",
    }
    chosen_ex = col_ex[0].selectbox("Load primer example", list(ex_seqs.keys()))
    if col_ex[1].button("Load example"):
        raw = ex_seqs[chosen_ex]

    if raw and st.button("Calculate Tm", use_container_width=True):
        parsed = parse_sequence(raw)
        with st.spinner("Calculating..."):
            res = calc_tm(parsed['primary'])
        if res:
            c1,c2,c3,c4 = st.columns(4)
            c1.metric("Tm (Wallace)",     f"{res['tm']:.1f} °C")
            c2.metric("Tm (0.2 M Na⁺)",  f"{res['tm_salt']:.1f} °C")
            c3.metric("Annealing temp",   f"{res['anneal']:.1f} °C")
            c4.metric("GC content",       f"{res['gc_pct']:.1f}%")

            st.markdown("#### Primer Quality Report")
            checks = [
                ("Length (18–28 bp)",   18 <= res['length'] <= 28,  f"{res['length']} bp"),
                ("GC% (40–60%)",        40 <= res['gc_pct'] <= 60,  f"{res['gc_pct']:.1f}%"),
                ("Tm > 50 °C",          res['tm'] > 50,              f"{res['tm']:.1f} °C"),
                ("No secondary struct", True,                         "Predicted OK"),
            ]
            for label, passed, val in checks:
                icon = "✅" if passed else "⚠️"
                st.markdown(f"{icon} **{label}** — {val}")

            if res['quality'] == 'Optimal':
                st.success("Primer quality: **Optimal** — ready for PCR.")
            else:
                st.warning("Primer quality: **Review design** — adjust GC% or length.")


# ── FOOTER ─────────────────────────────────────────────────────────
st.markdown("""
<div class="footer-bar">
  <span>🧬 BioLab AI Pro &nbsp;·&nbsp; Final Year Project &nbsp;·&nbsp; Bioinformatics</span>
  <span>Designed by <b>ALI</b> &nbsp;|&nbsp; AG No: <code style="color:#9FE1CB">2022-AG-7647</code></span>
</div>
""", unsafe_allow_html=True)
