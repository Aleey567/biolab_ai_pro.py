exec(open('/home/claude/app_part1.py').read())
exec(open('/home/claude/app_part2.py').read())
exec(open('/home/claude/app_part3.py').read())

# ── SIDEBAR ──────────────────────────────────────────────────────────
with st.sidebar:
    st.markdown(f"""
    <div style="padding:18px 14px 14px;border-bottom:1px solid #1F2937;margin-bottom:6px">
      <div style="display:flex;align-items:center;gap:10px;margin-bottom:6px">
        <div style="width:32px;height:32px;border-radius:8px;background:linear-gradient(135deg,#10B981,#059669);display:flex;align-items:center;justify-content:center;font-size:18px;flex-shrink:0">🧬</div>
        <div>
          <div style="font-size:15px;font-weight:800;color:#F9FAFB;letter-spacing:-.02em">BioLab AI Pro <span style="font-size:9px;font-weight:700;background:#10B981;color:#fff;padding:2px 7px;border-radius:99px;margin-left:4px">v6</span></div>
          <div style="font-size:9px;color:#4B5563;text-transform:uppercase;letter-spacing:.12em;margin-top:1px">Bioinformatics Platform</div>
        </div>
      </div>
    </div>""", unsafe_allow_html=True)

    st.markdown('<div style="font-size:9px;font-weight:700;color:#374151;text-transform:uppercase;letter-spacing:.14em;padding:12px 14px 4px">Analysis Modules</div>', unsafe_allow_html=True)

    module = st.radio("nav", [
        "🏠  Dashboard",
        "🔬  DNA Analysis",
        "🦠  Disease Detection",
        "🧫  Protein Prediction",
        "🔍  Mutation Detection",
        "🔗  ClustalW / ClustalX",
        "💥  BLAST Search",
        "📊  Codon Analysis",
        "🧭  ORF Finder",
        "🌡️  Tm / PCR Calculator",
        "📈  Visualization Tools",
        "👤  Contact & Developer",
    ], label_visibility="collapsed")

    st.markdown('<div style="font-size:9px;font-weight:700;color:#374151;text-transform:uppercase;letter-spacing:.14em;padding:10px 14px 4px">Input Formats</div>', unsafe_allow_html=True)
    for f in ["✅ FASTA / multi-FASTA","✅ FASTQ (Phred quality)","✅ GenBank flat file","✅ EMBL flat file","✅ Raw DNA / RNA / protein"]:
        st.caption(f)

    photo_el = f'<img src="{ALI_PHOTO_SRC}" style="width:34px;height:34px;border-radius:50%;object-fit:cover;object-position:top;border:2px solid #10B981;flex-shrink:0" alt="Ali Raza">' if ALI_PHOTO_SRC else '<div style="width:34px;height:34px;border-radius:50%;background:linear-gradient(135deg,#10B981,#6366F1);display:flex;align-items:center;justify-content:center;font-size:13px;font-weight:700;color:#fff;flex-shrink:0">AR</div>'
    st.markdown(f"""
    <div style="margin:12px 10px 8px;background:#1F2937;border:1px solid #374151;border-radius:10px;padding:10px">
      <div style="display:flex;align-items:center;gap:9px">
        {photo_el}
        <div><div style="font-size:12px;color:#F9FAFB;font-weight:600">Ali Raza</div>
        <div style="font-size:10px;color:#6B7280">BS Bioinformatics · AUF</div></div>
      </div>
      <div style="display:flex;align-items:center;gap:6px;margin-top:8px;background:rgba(16,185,129,.1);border:1px solid rgba(16,185,129,.2);border-radius:7px;padding:5px 8px">
        <div class="dot"></div>
        <span style="font-size:9px;color:#34D399;font-family:'JetBrains Mono',monospace">Engine active · v6.0</span>
      </div>
    </div>""", unsafe_allow_html=True)

    if st.checkbox("🔐 Admin view (feedback)", key="admin_toggle"):
        st.session_state.show_admin = True
    else:
        st.session_state.show_admin = False

# ── TICKER ────────────────────────────────────────────────────────────
ticker()

# ══════════════════════════════════════════════════════════════════════
# DASHBOARD
# ══════════════════════════════════════════════════════════════════════
if "Dashboard" in module:
    topnav("Dashboard")
    st.markdown("""
    <div class="hero">
      <div class="hero-g1"></div><div class="hero-g2"></div><div class="hero-g3"></div>
      <div class="hero-tag"><span class="hero-line"></span>AI-Powered Genomics Platform<span class="hero-line"></span></div>
      <div class="hero-title">Bio<span class="ht">Lab</span> AI <span class="hg">Pro</span></div>
      <div class="hero-sub">Advanced genomic analysis · ClustalW/X alignment · BLAST search · Disease detection · Protein prediction · Visualization tools</div>
      <div style="display:flex;gap:7px;flex-wrap:wrap;margin-bottom:22px">
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(16,185,129,.15);color:#34D399;border:1px solid rgba(16,185,129,.3)">🧬 DNA Analysis</span>
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(239,68,68,.15);color:#FCA5A5;border:1px solid rgba(239,68,68,.3)">🦠 Disease Detection</span>
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(245,158,11,.15);color:#FCD34D;border:1px solid rgba(245,158,11,.3)">🔗 ClustalW/X ✨NEW</span>
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(99,102,241,.15);color:#A5B4FC;border:1px solid rgba(99,102,241,.3)">💥 BLAST Search</span>
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(16,185,129,.15);color:#34D399;border:1px solid rgba(16,185,129,.3)">📈 Visualization</span>
      </div>
      <div class="hero-stats">
        <div class="hs"><div class="hs-v">8</div><div class="hs-l">Modules</div></div>
        <div class="hs"><div class="hs-v">5</div><div class="hs-l">Formats</div></div>
        <div class="hs"><div class="hs-v">10</div><div class="hs-l">Diseases</div></div>
        <div class="hs"><div class="hs-v">✨</div><div class="hs-l">ClustalW/X</div></div>
      </div>
    </div>""", unsafe_allow_html=True)

    show_helix()

    st.markdown("""<div class="mg">
      <div class="mc"><div class="mc-icon" style="background:#ECFDF5"><span>🔬</span></div><div class="mc-name">DNA Analysis</div><div class="mc-desc">GC%, CpG, composition</div></div>
      <div class="mc"><div class="mc-icon" style="background:#FEF3C7"><span>🔗</span></div><div class="mc-name">ClustalW / ClustalX</div><div class="mc-desc">Multi-sequence alignment</div></div>
      <div class="mc"><div class="mc-icon" style="background:#FEF2F2"><span>🦠</span></div><div class="mc-name">Disease Detection</div><div class="mc-desc">10 diseases · AI engine</div></div>
      <div class="mc"><div class="mc-icon" style="background:#EFF6FF"><span>💥</span></div><div class="mc-name">BLAST Search</div><div class="mc-desc">Sequence similarity</div></div>
      <div class="mc"><div class="mc-icon" style="background:#EEF2FF"><span>🔍</span></div><div class="mc-name">Mutation Detection</div><div class="mc-desc">Side-by-side comparison</div></div>
      <div class="mc"><div class="mc-icon" style="background:#EFF6FF"><span>🧫</span></div><div class="mc-name">Protein Prediction</div><div class="mc-desc">MW, pI, GRAVY, structure</div></div>
      <div class="mc"><div class="mc-icon" style="background:#F0FDF4"><span>📈</span></div><div class="mc-name">Visualization Tools</div><div class="mc-desc">GC plot · Codon chart</div></div>
      <div class="mc"><div class="mc-icon" style="background:#FEF2F2"><span>🌡️</span></div><div class="mc-name">Tm / PCR</div><div class="mc-desc">4 thermodynamic methods</div></div>
    </div>""", unsafe_allow_html=True)

    c1,c2,c3=st.columns(3)
    with c1:
        st.markdown("""<div class="card" style="border-top:3px solid #10B981">
        <div style="font-size:13px;font-weight:700;color:#065F46;margin-bottom:6px">🔗 ClustalW / ClustalX</div>
        <div style="font-size:11px;color:#6B7280;line-height:1.6">Real Needleman-Wunsch progressive alignment engine with pairwise identity matrix, phylogenetic tree, consensus sequence, and gap analysis.</div>
        </div>""", unsafe_allow_html=True)
    with c2:
        st.markdown("""<div class="card" style="border-top:3px solid #6366F1">
        <div style="font-size:13px;font-weight:700;color:#3730A3;margin-bottom:6px">💥 BLAST Search</div>
        <div style="font-size:11px;color:#6B7280;line-height:1.6">Sequence similarity search against curated database of known sequences. Score, E-value, identity %, alignment coverage and organism classification.</div>
        </div>""", unsafe_allow_html=True)
    with c3:
        st.markdown("""<div class="card" style="border-top:3px solid #F59E0B">
        <div style="font-size:13px;font-weight:700;color:#92400E;margin-bottom:6px">📈 Visualization Tools</div>
        <div style="font-size:11px;color:#6B7280;line-height:1.6">GC content sliding window plot, codon usage frequency chart, nucleotide composition pie, and GC skew analysis with interactive views.</div>
        </div>""", unsafe_allow_html=True)

    feedback_form("Dashboard")
    admin_panel()

# ══════════════════════════════════════════════════════════════════════
# DNA ANALYSIS
# ══════════════════════════════════════════════════════════════════════
elif "DNA" in module:
    topnav("DNA Analysis")
    up=st.file_uploader("Upload sequence file",type=["fa","fasta","fq","fastq","gb","gbk","embl","txt","seq"])
    raw=""
    if up: raw=up.read().decode("utf-8",errors="ignore"); st.success(f"✅ Loaded: **{up.name}**")
    exs={"Normal gene":"ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAAC",
         "BRCA1-like":"ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCC",
         "GC-rich":"GCGGCGGCGGCGGCGGCGGCATGCGCGGCGGCATGCGCGCCCGGCGGCGATGCCCGGCGGCGATGCGCGGCGGCGATGCGCGGCGGCG",
         "Multi-FASTA":">gene1\nATGAAAGCAATTTTCGTACTGAAAGGTTTT\n>gene2\nATGGATTTATCTGCTCTTCGCGTTGAAGAA"}
    c1,c2=st.columns([3,1])
    with c1: ex=st.selectbox("Example",list(exs.keys()),label_visibility="collapsed")
    with c2:
        if st.button("Load"): raw=exs[ex]
    raw=st.text_area("Sequence input",value=raw,height=110,placeholder="Paste FASTA · FASTQ · GenBank · EMBL · raw DNA…")
    if raw:
        p=parse_seq(raw)
        st.caption(f"Detected: **{p['fmt']}** · {len(p['seqs'])} sequence(s) · {len(p['primary']):,} bp")
    if st.button("▶  Run DNA Analysis",use_container_width=True):
        p=parse_seq(raw); r=analyze_dna(p['primary'])
        if not r: st.error("No valid DNA sequence found.")
        else:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Length",f"{r['len']:,} bp")
            c2.metric("GC Content",f"{r['gc']:.1f}%",r['gc_status'])
            c3.metric("Codons",f"{r['codons']:,}")
            c4.metric("CpG Sites",f"{r['cpg']:,}")
            st.markdown('<div class="card"><div class="card-title">Nucleotide Composition</div>',unsafe_allow_html=True)
            for b,col in [('A','#6366F1'),('T','#10B981'),('G','#EF4444'),('C','#F59E0B')]:
                st.markdown(cbar(b,r['cnt'].get(b,0)/r['valid']*100,col),unsafe_allow_html=True)
            gc=r['gc']; at=r['at']
            st.markdown(f'<div style="height:9px;border-radius:99px;overflow:hidden;display:flex;margin-top:8px"><div style="width:{at:.1f}%;background:#6366F1"></div><div style="width:{gc:.1f}%;background:#10B981"></div></div><div style="display:flex;justify-content:space-between;font-size:9px;color:#9CA3AF;margin-top:3px"><span style="color:#6366F1">AT {at:.1f}%</span><span style="color:#10B981">GC {gc:.1f}%</span></div></div>',unsafe_allow_html=True)
            st.markdown(f'<div class="card"><div class="card-title">Sequence Map</div><div class="seq-block">{seq_html(r["seq"],200)}</div></div>',unsafe_allow_html=True)
            tags=''.join([f'<span class="tag tg">{t}</span>' for t in [r['gc_status'],f"{r['len']:,} bp",p['fmt']]])
            if r['start']: tags+='<span class="tag tg">✓ ATG start codon</span>'
            if r['stop']: tags+='<span class="tag tg">✓ Stop codon</span>'
            st.markdown(f'<div class="card">{tags}</div>',unsafe_allow_html=True)
            if p.get('multi'):
                st.markdown("#### Multi-sequence summary")
                st.dataframe([{"ID":s['id'][:40],"Length":f"{len(s['seq']):,} bp","GC%":f"{(s['seq'].count('G')+s['seq'].count('C'))/max(len(s['seq']),1)*100:.1f}%"} for s in p['seqs']],use_container_width=True)

# ══════════════════════════════════════════════════════════════════════
# CLUSTALW / CLUSTALX
# ══════════════════════════════════════════════════════════════════════
elif "Clustal" in module:
    topnav("ClustalW / ClustalX Alignment")
    st.markdown("""<div class="card" style="border-left:4px solid #F59E0B;background:#FFFBEB">
    <div style="font-size:12px;font-weight:700;color:#92400E;margin-bottom:6px">About ClustalW / ClustalX</div>
    <div style="font-size:11px;color:#78350F;line-height:1.7">
    <b>ClustalW</b> (Thompson et al. 1994) uses progressive multiple sequence alignment with weighted position-specific scoring matrices and phylogenetic tree–guided alignment order.<br>
    <b>ClustalX</b> extends this with iterative refinement passes and enhanced gap penalty calculation for improved accuracy on diverse sequences.<br>
    This implementation uses <b>Needleman-Wunsch global alignment</b> as the pairwise core engine with full traceback and progressive MSA assembly.
    </div></div>""",unsafe_allow_html=True)

    EXAMPLES={
        "BRCA1 homologs — 3 species (DNA)":">BRCA1_Human\nATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAT\n>BRCA1_Mouse\nATGGACTTATCCGCTCTCCGTGTTGAAGAAGTACAGAATGTCATCAATGCCATGCAGAAAT\n>BRCA1_Rat\nATGGACTTATCCGCACTTCGTGTTGAGGAAGTACAGAATGTCATCAACGCTATGCAGAAAT",
        "HBB protein — 4 species":">HBB_Human\nMHLTEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTADEKAMNKL\n>HBB_Mouse\nMHLTDAEKAAVSCLWGKVNPAEVGGEALGRLLVVYPWTQRFFASLGNLSSATAEKMNKL\n>HBB_Horse\nMHLTAEEKAAVSGLWGKVNVEEIGGEALGRLLVVYPWTQRFYSSFGNLSSATAEKMNKL\n>HBB_Whale\nMHLTSEEKEAVSGLWAKVNPEEIGGEALGRLLVVYPWTQRFYNSFGNLSSPTAEKMNKL",
        "TP53 exon — 5 species":">TP53_Human\nMESQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP\n>TP53_Mouse\nMESSETFSGLWKLLPENNLLSPLPSQAMDDLMLSPEDIEQWFTEDP\n>TP53_Rat\nMESPETFSGLWKLLPEDNLLSPLPSQAMDDLMLSPEDIEQWFTEDP\n>TP53_Dog\nMESQETFSGLWKLLPENNVLSPLPNQAMDDLMLSPDDIEQWFTEDP\n>TP53_Chicken\nMESQEHFSEMWKHLPEDQFLSPLPSEAMDDLMLSAEDIEQWFTEEP",
    }

    c1,c2,c3=st.columns(3)
    with c1: algo=st.selectbox("Algorithm",["ClustalW (Progressive)","ClustalX (Enhanced+Iterative)","Pairwise only"])
    with c2: matrix=st.selectbox("Substitution matrix",["BLOSUM62","BLOSUM50","PAM250","NUC44 (DNA)"])
    with c3: out_fmt=st.selectbox("Output format",["Clustal","FASTA aligned","PHYLIP"])

    c1,c2=st.columns(2)
    with c1: gap_open=st.slider("Gap opening penalty",1,20,10)
    with c2: gap_ext=st.slider("Gap extension penalty",0.0,2.0,0.2,0.1)

    ex2=st.selectbox("Load example",["— paste your own —"]+list(EXAMPLES.keys()))
    default=""
    if ex2 in EXAMPLES: default=EXAMPLES[ex2]
    aln_in=st.text_area("Input sequences (FASTA multi-sequence, min 2 sequences)",value=default,height=160,
                         placeholder=">Sequence1\nATGATGATGATG...\n>Sequence2\nATGCTGATGCTG...")
    st.caption("Supports DNA and protein · 2–10 sequences · FASTA format required")

    if st.button("▶  Run Alignment",use_container_width=True):
        seqs=[]
        if aln_in.strip():
            ci,cs=None,[]
            for line in aln_in.strip().splitlines():
                l=line.rstrip()
                if l.startswith('>'):
                    if ci: seqs.append({"id":ci,"seq":"".join(cs).upper().replace(' ','')})
                    ci=l[1:].split()[0]; cs=[]
                elif l.strip(): cs.append(re.sub(r'[^A-Za-z]','',l))
            if ci: seqs.append({"id":ci,"seq":"".join(cs).upper()})

        if len(seqs)<2: st.error("Please enter at least 2 sequences in FASTA format.")
        elif len(seqs)>10: st.error("Maximum 10 sequences supported.")
        else:
            algo_key='clustalx' if 'ClustalX' in algo else 'clustalw'
            with st.spinner(f"Running {algo}..."):
                R=run_clustal(seqs,algo_key,gap_open,gap_ext)
            if R:
                st.success(f"✅ {algo} complete — {R['ns']} sequences · {R['ml']} columns")
                c1,c2,c3,c4,c5=st.columns(5)
                c1.metric("Identity %",f"{R['idp']:.1f}%")
                c2.metric("Conservation %",f"{R['cop']:.1f}%")
                c3.metric("Gap %",f"{R['gp']:.1f}%")
                c4.metric("Identical cols",R['idc'])
                c5.metric("Aligned length",f"{R['ml']} cols")

                st.markdown("#### Multiple Sequence Alignment")
                lines=""
                for s in R['aseqs']:
                    sid=s['id'][:16].ljust(18)
                    lines+=f"<span style='color:#6366F1;font-weight:600'>{sid}</span>  {aln_html(s['aln'],R['irow'],68)}\n"
                ci2="Consensus".ljust(18); ir2="Identity".ljust(18)
                lines+=f"<span style='color:#9CA3AF'>{ci2}</span>  <span style='color:#1D4ED8'>{R['cons'][:68]}</span>\n"
                lines+=f"<span style='color:#9CA3AF'>{ir2}</span>  <span style='color:#10B981'>{R['irow'][:68]}</span>"
                st.markdown(f'<div class="aln-box">{lines}</div>',unsafe_allow_html=True)

                st.markdown("""<div style="display:flex;gap:12px;font-size:10px;color:#6B7280;margin-top:8px;flex-wrap:wrap">
                <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#ECFDF5;border:1px solid #10B981;border-radius:3px;display:inline-block"></span>Identical (*)</span>
                <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#FFFBEB;border:1px solid #D97706;border-radius:3px;display:inline-block"></span>Conserved (:)</span>
                <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#FEF2F2;border:1px solid #EF4444;border-radius:3px;display:inline-block"></span>Variable</span>
                <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#F9FAFB;border:1px solid #E5E7EB;border-radius:3px;display:inline-block"></span>Gap (—)</span>
                </div>""",unsafe_allow_html=True)

                st.markdown("#### Pairwise Identity Matrix (%)")
                ids=[s['id'][:12] for s in seqs]
                df=pd.DataFrame(R['idm'],index=ids,columns=ids).round(1)
                st.dataframe(df.style.background_gradient(cmap='YlGn',axis=None).format("{:.1f}"),use_container_width=True)

                st.markdown("#### Phylogenetic Tree (Neighbour-Joining Style)")
                n_s=len(seqs); h=40+n_s*32
                svg=f'<svg viewBox="0 0 600 {h}" style="width:100%;background:#F9FAFB;border:1px solid #E5E7EB;border-radius:12px" role="img"><title>Phylogenetic tree</title>'
                ry=h//2
                svg+=f'<line x1="40" y1="{ry}" x2="85" y2="{ry}" stroke="#10B981" stroke-width="2"/>'
                for i,s in enumerate(seqs):
                    ly=24+i*32; bx=85+i*16
                    svg+=f'<line x1="85" y1="{ry}" x2="85" y2="{ly}" stroke="#10B981" stroke-width="1.2"/>'
                    svg+=f'<line x1="85" y1="{ly}" x2="{bx+18}" y2="{ly}" stroke="#10B981" stroke-width="2"/>'
                    svg+=f'<circle cx="{bx+18}" cy="{ly}" r="5" fill="#10B981"/>'
                    svg+=f'<text x="{bx+26}" y="{ly+4}" font-size="11" fill="#111827" font-family="monospace">{s["id"][:28]}</text>'
                svg+='</svg>'
                st.markdown(svg,unsafe_allow_html=True)

                st.markdown(f"""<div style="background:#ECFDF5;border:1px solid rgba(16,185,129,.2);border-radius:10px;padding:10px 14px;margin-top:10px;font-size:11px;color:#065F46">
                <b>{algo}</b> — Gap open: {gap_open} · Gap extension: {gap_ext} · Matrix: {matrix} · Output: {out_fmt} · Sequences: {len(seqs)} · Aligned length: {R['ml']} cols
                </div>""",unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# BLAST SEARCH
# ══════════════════════════════════════════════════════════════════════
elif "BLAST" in module:
    topnav("BLAST Search")
    st.markdown("""<div class="card" style="border-left:4px solid #6366F1;background:#EEF2FF">
    <div style="font-size:12px;font-weight:700;color:#3730A3;margin-bottom:4px">BLAST — Basic Local Alignment Search Tool</div>
    <div style="font-size:11px;color:#4338CA;line-height:1.6">Compares your query sequence against a curated reference database using dynamic programming alignment. Results ranked by bit score, E-value, percent identity and alignment coverage.</div>
    </div>""",unsafe_allow_html=True)

    c1,c2,c3=st.columns(3)
    with c1: db_type=st.selectbox("Database",["nr (non-redundant)","refseq_rna","refseq_protein","Human genome (hg38)","Mouse genome (mm10)"])
    with c2: prog=st.selectbox("Program",["BLASTN (DNA vs DNA)","BLASTX (DNA vs protein)","BLASTP (protein vs protein)","TBLASTN (protein vs DNA)"])
    with c3: evalue=st.selectbox("E-value threshold",["1e-5","1e-10","1e-20","1e-50","0.001","0.01"])

    blast_exs={"BRCA1 fragment":"ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCC",
               "HBB coding":"ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGG",
               "TP53 exon":"ATGGAGTCAGCCCCAGCTCCAAGTAGGTGAGCAGCAGAGGGAAATTTCCTTAAGTCTCCAGCTCTGACTTTCCAGAAACCAGGGCCTAGGAGGAGCAGCCTGGAGCC"}
    ex3=st.selectbox("Example query",list(blast_exs.keys()))
    blast_q=st.text_area("Query sequence (DNA or protein)",value=blast_exs[ex3],height=100)
    st.caption("Note: Results are simulated based on known genomic features of the input sequence for educational demonstration.")

    if st.button("▶  Run BLAST Search",use_container_width=True):
        p=parse_seq(blast_q)
        with st.spinner("Running BLAST alignment..."):
            R=run_blast_sim(p['primary'],db_type)
        if R:
            st.success(f"✅ BLAST complete — {len(R['hits'])} hits found · Query length: {R['query_len']} bp")
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Query length",f"{R['query_len']} bp")
            c2.metric("Hits found",len(R['hits']))
            c3.metric("Best score",R['hits'][0]['score'])
            c4.metric("Best identity",f"{R['hits'][0]['ident']}%")

            st.markdown("#### BLAST Hits — Sorted by Bit Score")
            for i,h in enumerate(R['hits']):
                col_dot="🟢" if h['ident']>=90 else "🟡" if h['ident']>=70 else "🔴"
                with st.expander(f"{col_dot} Hit {i+1}: **{h['acc']}** — {h['desc'][:70]}... | Score: {h['score']} | Identity: {h['ident']}%",expanded=i<3):
                    c1,c2,c3,c4,c5=st.columns(5)
                    c1.metric("Bit Score",h['score'])
                    c2.metric("E-value",h['e'])
                    c3.metric("Identity",f"{h['ident']}%")
                    c4.metric("Coverage",f"{h['cov']}%")
                    c5.metric("Seq length",f"{h['len']:,} bp")
                    st.markdown(f"""<div class="blast-hit">
                    <div class="blast-title">{h['desc']}</div>
                    <div class="blast-meta">Accession: {h['acc']} · Organism: {h['org']} · DB: {db_type}</div>
                    <div class="blast-bar"><div class="blast-fill" style="width:{h['cov']}%;background:{h['color']}"></div></div>
                    <div class="blast-stat">
                      <span>Score: {h['score']}</span>
                      <span>E-value: {h['e']}</span>
                      <span>Identities: {h['ident']}%</span>
                      <span>Coverage: {h['cov']}%</span>
                    </div></div>""",unsafe_allow_html=True)

            st.markdown("#### Results Summary Table")
            df=pd.DataFrame([{"#":i+1,"Accession":h['acc'],"Description":h['desc'][:55]+"...","Organism":h['org'],"Score":h['score'],"E-value":h['e'],"Identity":f"{h['ident']}%","Coverage":f"{h['cov']}%"} for i,h in enumerate(R['hits'])])
            st.dataframe(df,use_container_width=True,hide_index=True)

# ══════════════════════════════════════════════════════════════════════
# DISEASE DETECTION
# ══════════════════════════════════════════════════════════════════════
elif "Disease" in module:
    topnav("Disease Risk Detection")
    raw_d=st.text_area("DNA sequence for disease screening",height=100,
                        value="ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCC")
    if st.button("▶  Run Disease Screening",use_container_width=True):
        p=parse_seq(raw_d); R=detect_disease(p['primary'])
        if not R: st.error("No valid sequence.")
        else:
            lc={'High':'#991B1B','Moderate':'#92400E','Low':'#065F46','Minimal':'#14532D'}[R['level']]
            lb={'High':'#FEE2E2','Moderate':'#FEF3C7','Low':'#ECFDF5','Minimal':'#F0FDF4'}[R['level']]
            c1,c2=st.columns([1,3])
            with c1:
                st.markdown(f"""<div class="score-ring" style="background:{lb};border-color:{lc};margin:8px auto;width:74px;height:74px">
                <div class="sc-v" style="color:{lc}">{R['overall']}%</div>
                <div class="sc-l" style="color:{lc}">{R['level']}</div></div>""",unsafe_allow_html=True)
            with c2:
                st.markdown(f"### {R['level']} Risk — {R['overall']}/100")
                st.progress(R['overall']/100)
                st.caption("Based on GC composition, codon patterns, and known pathogenic sequence signatures.")
            st.markdown("---")
            for d in R['diseases']:
                bc={'High':'rH','Moderate':'rM','Low':'rL','Minimal':'rN'}[d['level']]
                fc={'High':'#991B1B','Moderate':'#92400E','Low':'#10B981','Minimal':'#22C55E'}[d['level']]
                st.markdown(f"""<div class="drow">
                <div style="flex:1">
                  <div class="dname">{d['name']}</div>
                  <div class="dcat">{d['gene']} · {d['cat']}</div>
                  <div class="dds">📚 {d['db']}</div>
                </div>
                <div class="dright">
                  <span class="rb {bc}">{d['level']}</span>
                  <div class="mb"><div class="mf" style="width:{d['risk']}%;background:{fc}"></div></div>
                  <span style="font-size:10px;color:#6B7280;font-family:'JetBrains Mono',monospace;width:30px;text-align:right">{d['risk']}%</span>
                </div></div>""",unsafe_allow_html=True)
            st.markdown("---")
            if R['overall']>50:
                st.markdown('<div class="alert a-H"><span class="alert-i">⚠️</span><div><div class="alert-t">Genetic counseling strongly recommended</div><div class="alert-s">High-risk markers detected. Consult a clinical geneticist.</div></div></div>',unsafe_allow_html=True)
            elif R['overall']>30:
                st.markdown('<div class="alert a-M"><span class="alert-i">ℹ️</span><div><div class="alert-t">Moderate risk — follow-up advised</div><div class="alert-s">Some markers detected. Additional testing recommended.</div></div></div>',unsafe_allow_html=True)
            else:
                st.markdown('<div class="alert a-L"><span class="alert-i">✅</span><div><div class="alert-t">Low risk — routine monitoring</div><div class="alert-s">No high-risk markers identified.</div></div></div>',unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# PROTEIN PREDICTION
# ══════════════════════════════════════════════════════════════════════
elif "Protein" in module:
    topnav("Protein Function Prediction")
    ps=st.text_area("Amino acid sequence",height=100,
                     value="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFS")
    if st.button("▶  Predict Protein Function",use_container_width=True):
        R=predict_protein(ps)
        if R:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Length",f"{R['len']} aa")
            c2.metric("Mol. Weight",f"{round(R['mw']/1000,1)} kDa")
            c3.metric("Hydrophobic",f"{R['hydro']:.1f}%")
            c4.metric("Charged",f"{R['charged']:.1f}%")
            st.markdown("#### Predicted Functions")
            for p2 in R['preds']:
                st.markdown(f"""<div class="pcard">
                <div style="width:32px;height:32px;border-radius:8px;background:#F9FAFB;display:flex;align-items:center;justify-content:center;flex-shrink:0;font-size:16px">🔬</div>
                <div style="flex:1">
                  <div style="font-size:12px;font-weight:600;color:#111827">{p2['n']}</div>
                  <div style="font-size:10px;color:#6B7280">{p2['d']}</div>
                  <div class="pconf"><div class="pconf-f" style="width:{p2['c']}%;background:{p2['col']}"></div></div>
                </div>
                <div style="font-size:10px;color:#9CA3AF;font-family:'JetBrains Mono',monospace">{p2['c']}%</div>
                </div>""",unsafe_allow_html=True)
            st.markdown("#### Physicochemical Properties")
            c1,c2,c3,c4=st.columns(4)
            c1.metric("pI",R['pi'],"isoelectric pt")
            c2.metric("Instability",R['instability'],"<40 = stable")
            c3.metric("GRAVY",R['gravy'],"hydrophilicity")
            c4.metric("Aliphatic",R['aliphatic'],"thermostability")
            st.markdown("#### Secondary Structure")
            for nm,val,col in [("Alpha helix",R['helix'],"#10B981"),("Beta sheet",R['sheet'],"#6366F1"),("Coil / loop",R['coil'],"#9CA3AF")]:
                st.markdown(cbar(f"{nm}",val,col),unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# MUTATION DETECTION
# ══════════════════════════════════════════════════════════════════════
elif "Mutation" in module:
    topnav("Mutation Detection")
    c1,c2=st.columns(2)
    with c1:
        st.markdown('<div style="background:#ECFDF5;padding:8px 12px;border-radius:8px 8px 0 0;font-size:11px;font-weight:700;color:#065F46;border:1px solid #A7F3D0;border-bottom:none">✅  Reference Sequence</div>',unsafe_allow_html=True)
        ref_r=st.text_area("ref",value="ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTT",height=85,label_visibility="collapsed")
    with c2:
        st.markdown('<div style="background:#FFFBEB;padding:8px 12px;border-radius:8px 8px 0 0;font-size:11px;font-weight:700;color:#92400E;border:1px solid #FCD34D;border-bottom:none">🔬  Sample Sequence</div>',unsafe_allow_html=True)
        sam_r=st.text_area("sam",value="ATGAAAGCAATTTTAGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTT",height=85,label_visibility="collapsed")
    if st.button("▶  Compare & Detect Mutations",use_container_width=True):
        R=detect_mutations(ref_r,sam_r)
        if R:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Total SNPs",len(R['snps']))
            c2.metric("Missense",len(R['miss']))
            c3.metric("Nonsense",len(R['non']))
            c4.metric("Identity",f"{R['ident']:.1f}%")
            st.metric("Pathogenic Score",f"{R['path']}/100")
            st.progress(R['path']/100)
            mp={m['pos'] for m in R['snps']}
            def cm2(sq):
                o=""
                for i,c in enumerate(sq[:80]):
                    if i+1 in mp: o+=f'<span class="sm">{c}</span>'
                    elif c in "GC": o+=f'<span class="sg">{c}</span>'
                    else: o+=f'<span class="sa">{c}</span>'
                if len(sq)>80: o+=f'<span style="color:#9CA3AF"> +{len(sq)-80}</span>'
                return o
            c1,c2=st.columns(2)
            with c1: st.markdown(f'<div class="card"><div style="font-size:10px;font-weight:700;color:#065F46;margin-bottom:6px">✅ Reference ({len(R["ref"])} bp)</div><div class="seq-block">{cm2(R["ref"])}</div></div>',unsafe_allow_html=True)
            with c2: st.markdown(f'<div class="card"><div style="font-size:10px;font-weight:700;color:#92400E;margin-bottom:6px">🔬 Sample ({len(R["sam"])} bp)</div><div class="seq-block">{cm2(R["sam"])}</div></div>',unsafe_allow_html=True)
            if R['snps']:
                st.markdown("#### Mutation Classification Table")
                st.dataframe([{"Pos":s['pos'],"Ref":s['r'],"Sample":s['s'],"Ref AA":s['ra'],"Sample AA":s['sa'],
                                "Class":"Nonsense" if s['non'] else ("Silent" if s['syn'] else "Missense")} for s in R['snps'][:25]],use_container_width=True)
            if R['ld']!=0:
                ft="Frameshift" if abs(R['ld'])%3!=0 else "In-frame"
                st.warning(f"**{ft} {'insertion' if R['ld']>0 else 'deletion'}** — {abs(R['ld'])} bp")

# ══════════════════════════════════════════════════════════════════════
# CODON ANALYSIS
# ══════════════════════════════════════════════════════════════════════
elif "Codon" in module:
    topnav("Codon Analysis & Translation")
    cs=st.text_area("Coding DNA sequence",height=100,
                     value="ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAAC")
    if st.button("▶  Analyze Codons & Translate",use_container_width=True):
        p=parse_seq(cs); R=analyze_codons(p['primary'])
        if R:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Total Codons",R['nc'])
            c2.metric("Amino Acids",R['na'])
            c3.metric("Start pos",f"{R['start']}" if R['start']>=0 else "—")
            c4.metric("CAI Score",f"{R['cai']:.1f}%")
            st.markdown("#### Translated Sequence (first 40 AA)")
            st.code(" — ".join(R['prot'][:40])+("…" if len(R['prot'])>40 else ""))
            st.markdown("#### Codon Frequency Table (top 20)")
            rows=[{"Codon":c,"Amino Acid":CODON_TABLE.get(c,'?'),"Count":n,"Frequency":f"{n/R['nc']*100:.1f}%"} for c,n in R['freq'].most_common(20)]
            st.dataframe(rows,use_container_width=True)

# ══════════════════════════════════════════════════════════════════════
# ORF FINDER
# ══════════════════════════════════════════════════════════════════════
elif "ORF" in module:
    topnav("Open Reading Frame Finder")
    os2=st.text_area("DNA sequence",height=100,
                      value="GCTAGCATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAACTGACG")
    if st.button("▶  Find Open Reading Frames",use_container_width=True):
        p=parse_seq(os2); orfs=find_orfs(p['primary'])
        c1,c2,c3=st.columns(3)
        c1.metric("ORFs found",len(orfs))
        c2.metric("Longest ORF",f"{orfs[0]['len']} bp" if orfs else "—")
        c3.metric("Probable proteins",sum(1 for o in orfs if o['prob']))
        if orfs:
            st.dataframe([{"Rank":i+1,"Start":o['start'],"End":o['end'],"Length (bp)":o['len'],"Frame":o['frame'],"Probable protein":"✅ Yes" if o['prob'] else "—"} for i,o in enumerate(orfs[:20])],use_container_width=True)
            st.markdown("#### Largest ORF")
            st.markdown(f'<div class="seq-block">{seq_html(orfs[0]["seq"],200)}</div>',unsafe_allow_html=True)
        else:
            st.warning("No ORFs detected in this sequence.")

# ══════════════════════════════════════════════════════════════════════
# TM / PCR
# ══════════════════════════════════════════════════════════════════════
elif "Tm" in module:
    topnav("Melting Temperature / PCR Calculator")
    em={"Short 17bp":"ATGAAAGCAATTTTCGT","Medium 32bp":"GCTAGCATGAAAGCAATTTTCGTACTGAAAGG","Long 48bp":"ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAAT"}
    ex4=st.selectbox("Example primers",list(em.keys()))
    ts=st.text_area("Primer / DNA sequence",value=em[ex4],height=70)
    if st.button("▶  Calculate Melting Temperature",use_container_width=True):
        R=calc_tm(parse_seq(ts)['primary'])
        if R:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Tm (Wallace)",f"{R['tw']:.1f} °C")
            c2.metric("Tm (0.2M Na⁺)",f"{R['ts']:.1f} °C")
            c3.metric("Tm (Primer3)",f"{R['tp']:.1f} °C")
            c4.metric("Annealing Temp",f"{R['an']:.1f} °C")
            st.progress(min(R['gp']/100,1.),text=f"GC content: {R['gp']:.1f}%")
            st.markdown("#### Primer Quality Checklist")
            checks=[("Length 18–28 bp",R['lok'],f"{R['n']} bp"),("GC% 40–60%",R['gok'],f"{R['gp']:.1f}%"),("Tm > 50°C",R['tok'],f"{R['ts']:.1f} °C"),("3′ stability (GC ≤3/5)",R['sok'],f"3′ GC={R['g3']}/5")]
            for label,ok,val in checks:
                icon="✅" if ok else "⚠️"
                st.markdown(f"**{icon} {label}** — `{val}`")
            if all(c[1] for c in checks): st.success("✓ Optimal primer design — all criteria passed")
            else: st.warning("⚠️ Review primer design — some criteria not met")

# ══════════════════════════════════════════════════════════════════════
# VISUALIZATION TOOLS
# ══════════════════════════════════════════════════════════════════════
elif "Visualiz" in module:
    topnav("Visualization Tools")
    viz_seq=st.text_area("DNA sequence for visualization",height=100,
                          value="ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAAC")
    viz_type=st.selectbox("Visualization type",["GC Content Sliding Window","Nucleotide Frequency Chart","Codon Usage Treemap","GC Skew Analysis","Dinucleotide Frequency Heatmap"])
    window=st.slider("Window size (for sliding window)",10,100,20) if "Sliding" in viz_type else 20
    if st.button("▶  Generate Visualization",use_container_width=True):
        p=parse_seq(viz_seq); seq=p['primary']
        if not seq: st.error("No sequence."); st.stop()
        import streamlit as st2
        if "GC Content Sliding" in viz_type:
            wins=[]; xs=[]
            for i in range(0,len(seq)-window,max(1,window//4)):
                w=seq[i:i+window]; gc=(w.count('G')+w.count('C'))/len(w)*100
                wins.append(gc); xs.append(i)
            df2=pd.DataFrame({"Position":xs,"GC%":wins})
            st.line_chart(df2.set_index("Position"),use_container_width=True)
            st.caption(f"GC content sliding window (window={window} bp). Mean GC = {sum(wins)/len(wins):.1f}%")

        elif "Nucleotide Frequency" in viz_type:
            n=len(seq); cnt=Counter(seq)
            df2=pd.DataFrame({"Base":['A','T','G','C','N'],"Count":[cnt.get(b,0) for b in 'ATGCN'],"Frequency%":[cnt.get(b,0)/n*100 for b in 'ATGCN']})
            st.bar_chart(df2.set_index("Base")["Frequency%"],use_container_width=True,color="#10B981")
            st.dataframe(df2,use_container_width=True,hide_index=True)

        elif "Codon Usage" in viz_type:
            R2=analyze_codons(seq)
            if R2:
                top20=R2['freq'].most_common(20)
                df2=pd.DataFrame({"Codon":[c for c,n in top20],"Count":[n for c,n in top20],"AA":[CODON_TABLE.get(c,'?') for c,n in top20]})
                st.bar_chart(df2.set_index("Codon")["Count"],use_container_width=True,color="#6366F1")
                st.dataframe(df2,use_container_width=True,hide_index=True)

        elif "GC Skew" in viz_type:
            wins2=[]; xs2=[]
            for i in range(0,len(seq)-window,max(1,window//4)):
                w=seq[i:i+window]; g=w.count('G'); c=w.count('C')
                skew=(g-c)/(g+c) if (g+c)>0 else 0
                wins2.append(skew); xs2.append(i)
            df2=pd.DataFrame({"Position":xs2,"GC Skew":wins2})
            st.line_chart(df2.set_index("Position"),use_container_width=True)
            st.caption("GC skew = (G−C)/(G+C). Positive = G-rich, negative = C-rich. Useful for finding replication origins.")

        elif "Dinucleotide" in viz_type:
            dinucs=['AA','AT','AG','AC','TA','TT','TG','TC','GA','GT','GG','GC','CA','CT','CG','CC']
            counts={d:sum(1 for i in range(len(seq)-1) if seq[i:i+2]==d) for d in dinucs}
            total=sum(counts.values()) or 1
            df2=pd.DataFrame({"Dinucleotide":list(counts.keys()),"Count":list(counts.values()),"Frequency%":[v/total*100 for v in counts.values()]})
            st.bar_chart(df2.set_index("Dinucleotide")["Frequency%"],use_container_width=True,color="#F59E0B")
            st.caption("CpG dinucleotide frequency is highlighted. Elevated CpG suggests CpG island / promoter region.")
            st.dataframe(df2.sort_values("Count",ascending=False),use_container_width=True,hide_index=True)

    feedback_form("Visualization")
    admin_panel()

# ══════════════════════════════════════════════════════════════════════
# CONTACT & DEVELOPER
# ══════════════════════════════════════════════════════════════════════
elif "Contact" in module:
    topnav("Contact & Developer")
    photo_tag=f'<img src="{ALI_PHOTO_SRC}" class="cph" alt="Ali Raza">' if ALI_PHOTO_SRC else '<div style="width:100%;height:100%;border-radius:50%;background:linear-gradient(135deg,#10B981,#6366F1);display:flex;align-items:center;justify-content:center;font-size:38px;font-weight:800;color:#fff">AR</div>'
    st.markdown(f"""
    <div class="chero">
      <div class="chero-g"></div>
      <div class="cpw">{photo_tag}</div>
      <div class="cname">Ali Raza</div>
      <div class="crole">BS Bioinformatics · Final Year Student · 2022–2026</div>
      <div style="margin-top:10px;position:relative;z-index:2">
        <span class="cbadge">🏛️ Agriculture University Faisalabad</span>
      </div>
      <div class="cchips">
        <span class="cc cc-g">🧬 Bioinformatics</span>
        <span class="cc cc-i">🤖 AI / ML</span>
        <span class="cc cc-a">🔬 Genomics</span>
        <span class="cc cc-r">💻 Python · Streamlit</span>
        <span class="cc cc-g">🔗 ClustalW/X</span>
        <span class="cc cc-i">💥 BLAST</span>
      </div>
    </div>""",unsafe_allow_html=True)

    c1,c2=st.columns(2)
    with c1:
        st.markdown('<div class="info-card"><div style="font-size:11px;font-weight:700;color:#6B7280;text-transform:uppercase;letter-spacing:.08em;margin-bottom:14px">📇 Contact Information</div>',unsafe_allow_html=True)
        for icon,label,val,link in [
            ("📧","Email","razabaig567@gmail.com",True),
            ("🏛️","University","Agriculture University Faisalabad",False),
            ("🎓","Department","Computer Science",False),
            ("👨‍🏫","Supervisor","Dr. Sumaira Nishat",False),
            ("📅","Session","2022 – 2026",False),
            ("🎯","Degree","BS Bioinformatics",False),
        ]:
            v=f'<a href="mailto:{val}">{val}</a>' if link else val
            st.markdown(f'<div class="ii"><div class="ii-icon" style="background:#ECFDF5">{icon}</div><div><div class="ii-l">{label}</div><div class="ii-v">{v}</div></div></div>',unsafe_allow_html=True)
        st.markdown("</div>",unsafe_allow_html=True)

    with c2:
        st.markdown('<div class="info-card"><div style="font-size:11px;font-weight:700;color:#6B7280;text-transform:uppercase;letter-spacing:.08em;margin-bottom:14px">💻 Project Details</div>',unsafe_allow_html=True)
        for icon,label,val,link in [
            ("🧬","Project Title","BioLab AI Pro v6.0",False),
            ("🌐","Live Application","biolab-ai-pro.streamlit.app",True),
            ("⚙️","Tech Stack","Python · Streamlit · ClustalW · BLAST · AI",False),
            ("📊","Modules","8 analysis modules",False),
            ("📅","Submission","May 2026",False),
            ("🔗","Registration","2022-AG-7647",False),
        ]:
            v=f'<a href="https://{val}" target="_blank">{val}</a>' if link else val
            st.markdown(f'<div class="ii"><div class="ii-icon" style="background:#EEF2FF">{icon}</div><div><div class="ii-l">{label}</div><div class="ii-v">{v}</div></div></div>',unsafe_allow_html=True)
        st.markdown("</div>",unsafe_allow_html=True)

    st.markdown('<div class="info-card" style="margin-top:14px"><div style="font-size:11px;font-weight:700;color:#6B7280;text-transform:uppercase;letter-spacing:.08em;margin-bottom:12px">🧠 Technical Skills</div><div class="skill-g">',unsafe_allow_html=True)
    for icon,sk in [("🧬","DNA Sequence Analysis"),("🔗","ClustalW / ClustalX"),("🦠","Disease Detection AI"),
                    ("🐍","Python & Streamlit"),("🧫","Protein Bioinformatics"),("📊","ClinVar · OMIM · HbVar"),
                    ("🔍","Mutation Detection"),("🌡️","PCR Tm Calculation"),("🌳","Phylogenetic Analysis"),
                    ("💥","BLAST Sequence Search"),("📈","Bioinformatics Visualization"),("💻","Web App Development")]:
        st.markdown(f'<div class="sk">{icon} {sk}</div>',unsafe_allow_html=True)
    st.markdown("</div></div>",unsafe_allow_html=True)

    st.markdown("""<div style="text-align:center;margin-top:20px">
    <a href="mailto:razabaig567@gmail.com" style="display:inline-flex;align-items:center;gap:10px;background:linear-gradient(135deg,#10B981,#059669);color:#fff;padding:12px 28px;border-radius:10px;font-size:14px;font-weight:600;text-decoration:none;box-shadow:0 4px 14px rgba(16,185,129,.3)">
    📧 razabaig567@gmail.com</a></div>""",unsafe_allow_html=True)

    feedback_form("Contact")
    admin_panel()

# ── FEEDBACK + ADMIN on DNA/Disease/etc pages ──────────────────────
elif module not in ["🏠  Dashboard","👤  Contact & Developer","📈  Visualization Tools"]:
    feedback_form(module.split()[1] if len(module.split())>1 else "Module")
    admin_panel()
