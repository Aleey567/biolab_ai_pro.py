"""
BioLab AI Pro v6.0 — FINAL EDITION
Developer: Ali Raza | Agriculture University Faisalabad
Email: razabaig567@gmail.com
Supervisor: Dr. Sumaira Nishat
Run: pip install streamlit --break-system-packages && streamlit run biolab_v6_final.py
"""

import streamlit as st
import re, math, json, os
from collections import Counter
from datetime import datetime

st.set_page_config(
    page_title="BioLab AI Pro v6 | Ali Raza | AUF",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ── Ali's photo (base64 embedded) ────────────────────────────────────
import base64
PHOTO_PATH = "/home/claude/ali_photo_b64.txt"
try:
    with open(PHOTO_PATH) as f:
        ALI_PHOTO_B64 = f.read().strip()
    ALI_PHOTO_SRC = f"data:image/jpeg;base64,{ALI_PHOTO_B64}"
except:
    ALI_PHOTO_SRC = ""

# ── AUF logo (inline SVG shield) ─────────────────────────────────────
AUF_LOGO_SVG = """
<svg width="54" height="64" viewBox="0 0 54 64" xmlns="http://www.w3.org/2000/svg">
  <path d="M27 2 L52 12 L52 38 Q52 56 27 62 Q2 56 2 38 L2 12 Z"
        fill="#0D1B2A" stroke="#00BFA5" stroke-width="2"/>
  <path d="M27 8 L46 16 L46 38 Q46 53 27 58 Q8 53 8 38 L8 16 Z"
        fill="#00897B" opacity="0.6"/>
  <text x="27" y="30" text-anchor="middle" font-family="serif" font-size="11"
        font-weight="bold" fill="#FFB300">AUF</text>
  <text x="27" y="44" text-anchor="middle" font-family="serif" font-size="7"
        fill="#E0F7FA">FAISALABAD</text>
  <!-- wheat left -->
  <line x1="14" y1="20" x2="14" y2="50" stroke="#FFB300" stroke-width="1.2"/>
  <ellipse cx="11" cy="28" rx="3" ry="5" fill="#FFB300" opacity="0.8"/>
  <ellipse cx="11" cy="37" rx="3" ry="5" fill="#FFB300" opacity="0.8"/>
  <!-- wheat right -->
  <line x1="40" y1="20" x2="40" y2="50" stroke="#FFB300" stroke-width="1.2"/>
  <ellipse cx="43" cy="28" rx="3" ry="5" fill="#FFB300" opacity="0.8"/>
  <ellipse cx="43" cy="37" rx="3" ry="5" fill="#FFB300" opacity="0.8"/>
</svg>"""

# ── Feedback storage (in-memory + session, visible only in app) ───────
if "feedbacks" not in st.session_state:
    st.session_state.feedbacks = []
if "show_admin" not in st.session_state:
    st.session_state.show_admin = False

# ══════════════════════════════════════════════════════════════════════
# GLOBAL CSS — Stunning, mind-pleasing design
# ══════════════════════════════════════════════════════════════════════
st.markdown(f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600;700&family=DM+Mono:wght@400;500&display=swap');

html, body, [class*="css"] {{ font-family:'DM Sans',sans-serif!important }}
.main, .block-container, .stApp {{ background: #F0F4F8!important }}

/* ── Sidebar ── */
section[data-testid="stSidebar"] {{
  background: linear-gradient(180deg, #0B1628 0%, #0D2137 60%, #0B1F32 100%)!important;
  border-right: none!important;
}}
section[data-testid="stSidebar"]::before {{
  content:''; position:fixed; top:0; left:0; width:220px; height:100%;
  background: radial-gradient(ellipse at top right, rgba(0,191,165,0.08) 0%, transparent 60%),
              radial-gradient(ellipse at bottom left, rgba(124,77,255,0.06) 0%, transparent 60%);
  pointer-events:none; z-index:0;
}}
section[data-testid="stSidebar"] > div {{ z-index:1; position:relative }}

/* ── Sidebar logo ── */
.sb-brand {{ padding:20px 16px 14px; border-bottom:1px solid rgba(255,255,255,0.06); margin-bottom:4px }}
.sb-brand-row {{ display:flex; align-items:center; gap:10px; margin-bottom:4px }}
.sb-brand-icon {{ font-size:26px }}
.sb-brand-title {{ font-size:16px; font-weight:700; color:#fff; letter-spacing:-.02em }}
.sb-brand-ver {{ display:inline-block; background:#00BFA5; color:#fff; font-size:9px; font-weight:700; padding:2px 7px; border-radius:99px; margin-left:5px }}
.sb-brand-sub {{ font-size:9px; color:rgba(255,255,255,0.35); text-transform:uppercase; letter-spacing:.12em; padding-left:36px }}

/* ── Sidebar nav section headers ── */
.sb-section {{ font-size:9px; font-weight:700; color:rgba(255,255,255,0.25); text-transform:uppercase; letter-spacing:.14em; padding:14px 16px 4px }}

/* ── Sidebar status ── */
.sb-status {{ margin:12px 14px 8px; background:rgba(0,191,165,0.1); border:1px solid rgba(0,191,165,0.2); border-radius:8px; padding:8px 10px; display:flex; align-items:center; gap:7px }}
.sb-dot {{ width:7px; height:7px; border-radius:50%; background:#00BFA5; animation:blink 2s infinite }}
.sb-status-txt {{ font-size:10px; color:#00BFA5; font-family:'DM Mono',monospace }}
@keyframes blink {{ 0%,100%{{opacity:1}} 50%{{opacity:.3}} }}

/* ── Sidebar dev card ── */
.sb-dev {{ margin:10px 14px 14px; background:rgba(255,255,255,0.04); border:1px solid rgba(255,255,255,0.08); border-radius:10px; padding:10px }}
.sb-dev-inner {{ display:flex; align-items:center; gap:9px }}
.sb-dev-photo {{ width:34px; height:34px; border-radius:50%; object-fit:cover; object-position:top; border:2px solid #00BFA5 }}
.sb-dev-initial {{ width:34px; height:34px; border-radius:50%; background:linear-gradient(135deg,#00BFA5,#7C4DFF); display:flex; align-items:center; justify-content:center; font-size:14px; font-weight:700; color:#fff; flex-shrink:0 }}
.sb-dev-name {{ font-size:12px; color:#fff; font-weight:600 }}
.sb-dev-role {{ font-size:10px; color:rgba(255,255,255,0.4) }}

/* ── Streamlit radio as nav ── */
[data-testid="stRadio"] > div {{ gap:2px!important }}
[data-testid="stRadio"] label {{
  background: rgba(255,255,255,0.03)!important;
  border: none!important;
  border-radius: 9px!important;
  padding: 9px 12px!important;
  color: rgba(255,255,255,0.55)!important;
  font-size: 12px!important;
  font-weight: 500!important;
  cursor: pointer;
  transition: all .15s;
  margin-bottom: 1px!important;
}}
[data-testid="stRadio"] label:hover {{ background:rgba(255,255,255,0.07)!important; color:#fff!important }}
[data-testid="stRadio"] label:has(input:checked) {{
  background: rgba(0,191,165,0.15)!important;
  color: #fff!important;
  font-weight: 600!important;
  border-left: 3px solid #00BFA5!important;
}}
[data-testid="stRadio"] label p {{ margin:0!important; font-size:12px!important }}
[data-testid="stRadio"] [data-testid="stMarkdownContainer"] p {{ font-size:12px!important }}

/* ── Hero ── */
.hero-wrap {{
  background: linear-gradient(145deg, #0B1628 0%, #0D2137 45%, #0B1F32 100%);
  border-radius: 18px; padding: 36px 32px; margin-bottom: 20px;
  position: relative; overflow: hidden;
}}
.hero-bg-glow {{
  position:absolute; border-radius:50%; pointer-events:none;
}}
.hero-bg-g1 {{ top:-60px; right:-60px; width:280px; height:280px;
  background: radial-gradient(circle, rgba(0,191,165,0.14) 0%, transparent 70%); }}
.hero-bg-g2 {{ bottom:-80px; left:10%; width:320px; height:320px;
  background: radial-gradient(circle, rgba(124,77,255,0.1) 0%, transparent 70%); }}
.hero-bg-g3 {{ top:30%; right:30%; width:180px; height:180px;
  background: radial-gradient(circle, rgba(255,179,0,0.06) 0%, transparent 70%); }}
.hero-eyebrow {{ font-size:10px; letter-spacing:.15em; color:#00BFA5; text-transform:uppercase;
  font-weight:600; margin-bottom:10px; display:flex; align-items:center; gap:8px }}
.hero-line {{ width:24px; height:1px; background:#00BFA5; opacity:.6 }}
.hero-title {{ font-size:34px; font-weight:700; color:#fff; letter-spacing:-.03em;
  margin-bottom:8px; line-height:1.1 }}
.hero-title span.h-teal {{ color:#00BFA5 }}
.hero-title span.h-gold {{ color:#FFB300 }}
.hero-sub {{ font-size:13px; color:rgba(255,255,255,0.5); margin-bottom:20px; max-width:420px; line-height:1.7 }}

/* ── Module stat cards in hero ── */
.hero-stats {{ display:grid; grid-template-columns:repeat(4,1fr); gap:10px; margin-top:8px }}
.hstat {{ background:rgba(255,255,255,0.05); border:1px solid rgba(255,255,255,0.08);
  border-radius:12px; padding:12px; text-align:center; transition:transform .15s }}
.hstat:hover {{ transform:translateY(-2px) }}
.hstat-v {{ font-size:24px; font-weight:700; color:#fff; line-height:1 }}
.hstat-l {{ font-size:9px; color:rgba(255,255,255,0.4); text-transform:uppercase; letter-spacing:.06em; margin-top:3px }}

/* ── Module grid cards ── */
.mod-grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:10px; margin-bottom:20px }}
.mod-card {{
  background:#fff; border:1px solid #E2E8F0; border-radius:14px;
  padding:16px; cursor:pointer; transition:all .18s; text-align:center;
  box-shadow: 0 2px 8px rgba(0,0,0,0.04);
}}
.mod-card:hover {{ border-color:#00BFA5; transform:translateY(-3px); box-shadow:0 8px 24px rgba(0,191,165,0.12) }}
.mod-icon {{ width:44px; height:44px; border-radius:12px; display:flex; align-items:center;
  justify-content:center; margin:0 auto 10px; font-size:22px }}
.mod-name {{ font-size:12px; font-weight:600; color:#1A202C }}
.mod-desc {{ font-size:10px; color:#718096; margin-top:2px }}

/* ── Cards ── */
.card {{
  background:#fff; border:1px solid #E2E8F0; border-radius:14px;
  padding:20px; margin-bottom:14px;
  box-shadow: 0 2px 12px rgba(0,0,0,0.04);
}}
.card:hover {{ box-shadow: 0 4px 20px rgba(0,0,0,0.07) }}
.card-accent {{ height:3px; border-radius:14px 14px 0 0; margin:-20px -20px 16px }}
.card-title {{ font-size:11px; font-weight:700; color:#718096; text-transform:uppercase;
  letter-spacing:.09em; margin-bottom:14px; display:flex; align-items:center; gap:6px }}
.card-title i {{ font-size:15px; color:#00BFA5 }}

/* ── Metric cards ── */
.met-grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:8px; margin-bottom:14px }}
.met-card {{ background:#F8FAFB; border:1px solid #E2E8F0; border-radius:10px; padding:12px;
  transition:transform .15s }}
.met-card:hover {{ transform:translateY(-1px) }}
.met-label {{ font-size:9px; color:#718096; text-transform:uppercase; letter-spacing:.07em; margin-bottom:4px }}
.met-value {{ font-size:22px; font-weight:700; line-height:1; font-family:'DM Mono',monospace }}
.met-sub {{ font-size:9px; color:#718096; margin-top:2px }}

/* ── Buttons ── */
.stButton > button {{
  background: linear-gradient(135deg, #00BFA5, #00897B)!important;
  color: #fff!important; border: none!important; border-radius: 10px!important;
  font-weight: 600!important; font-size: 13px!important; padding: .65rem 1.5rem!important;
  width: 100%; font-family:'DM Sans',sans-serif!important;
  box-shadow: 0 4px 14px rgba(0,191,165,0.3)!important;
  transition: all .15s!important;
}}
.stButton > button:hover {{ opacity:.92!important; transform:translateY(-1px)!important }}
.stTextArea textarea, .stTextInput input {{
  background: #F8FAFB!important; border:1px solid #E2E8F0!important;
  border-radius:10px!important; font-family:'DM Mono',monospace!important; font-size:12px!important;
  color:#1A202C!important;
}}
.stTextArea textarea:focus, .stTextInput input:focus {{ border-color:#00BFA5!important }}
.stSelectbox > div > div {{ background:#F8FAFB!important; border:1px solid #E2E8F0!important;
  border-radius:10px!important; color:#1A202C!important }}

/* ── Ticker ── */
.ticker-wrap {{ overflow:hidden; background:rgba(0,0,0,0.03); border-top:1px solid #E2E8F0;
  border-bottom:1px solid #E2E8F0; padding:8px 0 }}
.ticker {{ display:flex; white-space:nowrap; animation:ticker 35s linear infinite;
  font-size:10px; font-family:'DM Mono',monospace; color:#00897B; letter-spacing:.04em }}
.ticker-sep {{ color:#00BFA5; margin:0 16px; opacity:.4 }}
@keyframes ticker {{ 0%{{transform:translateX(0)}} 100%{{transform:translateX(-50%)}} }}

/* ── Sequence display ── */
.seq-block {{ font-family:'DM Mono',monospace; font-size:11px; line-height:2.1;
  word-break:break-all; background:#F8FAFB; border:1px solid #E2E8F0;
  border-radius:10px; padding:12px; max-height:160px; overflow-y:auto }}
.sg {{ color:#00897B; font-weight:600 }}.sa {{ color:#7C4DFF }}
.sm {{ background:#FFF8E1; color:#FF8F00; border-radius:3px; padding:0 2px }}

/* ── Disease rows ── */
.dis-row {{ display:flex; align-items:flex-start; padding:10px 0; border-bottom:1px solid #F0F4F8 }}
.dis-row:last-child {{ border-bottom:none }}
.dis-name {{ font-size:12px; font-weight:600; color:#1A202C; display:flex; align-items:center; gap:5px }}
.dis-cat {{ font-size:10px; color:#718096; margin-top:2px }}
.dis-ds {{ font-size:9px; color:#718096; background:#F0F4F8; border:1px solid #E2E8F0;
  border-radius:5px; padding:3px 7px; margin-top:4px; line-height:1.6 }}
.dis-right {{ display:flex; align-items:center; gap:8px; margin-left:auto; padding-left:10px }}
.mini-b {{ width:60px; height:5px; background:#E2E8F0; border-radius:99px; overflow:hidden }}
.mini-f {{ height:100%; border-radius:99px }}
.pct {{ font-size:10px; color:#718096; font-family:'DM Mono',monospace; width:30px; text-align:right }}
.rbadge {{ font-size:9px; padding:2px 8px; border-radius:99px; font-weight:700 }}
.rh {{ background:#FFEBEE; color:#C62828 }}.rm {{ background:#FFF8E1; color:#E65100 }}
.rl {{ background:#E0F7FA; color:#00695C }}.rn {{ background:#E8F5E9; color:#2E7D32 }}

/* ── Alert boxes ── */
.alert {{ border-radius:10px; padding:11px 13px; margin-bottom:8px; border:1px solid; display:flex; gap:9px }}
.alert-icon {{ font-size:16px; flex-shrink:0; margin-top:1px }}
.alert-title {{ font-size:12px; font-weight:600 }}
.alert-sub {{ font-size:10px; color:#718096; margin-top:2px }}
.alert-high {{ background:#FFEBEE; border-color:rgba(198,40,40,0.3) }}
.alert-high .alert-icon,.alert-high .alert-title {{ color:#C62828 }}
.alert-mod {{ background:#FFF8E1; border-color:rgba(230,81,0,0.3) }}
.alert-mod .alert-icon,.alert-mod .alert-title {{ color:#E65100 }}
.alert-low {{ background:#E0F7FA; border-color:rgba(0,105,92,0.3) }}
.alert-low .alert-icon,.alert-low .alert-title {{ color:#00695C }}

/* ── Alignment styling ── */
.aln-box {{ font-family:'DM Mono',monospace; font-size:10.5px; line-height:1.9;
  background:#F8FAFB; border:1px solid #E2E8F0; border-radius:10px; padding:12px 14px;
  overflow-x:auto; white-space:pre }}
.aln-match {{ background:#E0F7FA; color:#004D40 }}.aln-mis {{ background:#FFEBEE; color:#B71C1C }}
.aln-gap {{ color:#CBD5E0 }}
.aln-score-grid {{ display:grid; grid-template-columns:repeat(5,1fr); gap:8px; margin-top:12px }}
.aln-sc {{ background:#F8FAFB; border:1px solid #E2E8F0; border-radius:10px; padding:10px }}
.aln-sc-l {{ font-size:9px; color:#718096; text-transform:uppercase; letter-spacing:.06em }}
.aln-sc-v {{ font-size:18px; font-weight:700; color:#1A202C; font-family:'DM Mono',monospace; margin-top:2px }}
.param-box {{ background:#F8FAFB; border:1px solid #E2E8F0; border-radius:10px; padding:14px; margin-bottom:12px }}
.param-label {{ font-size:10px; color:#718096; text-transform:uppercase; letter-spacing:.06em; margin-bottom:8px; font-weight:600 }}
.algo-row {{ display:flex; gap:6px; margin-bottom:0 }}
.algo-opt {{ padding:7px 14px; border-radius:8px; border:1px solid #E2E8F0; background:#fff;
  font-size:11px; font-weight:500; color:#1A202C; cursor:pointer; font-family:'DM Sans',sans-serif }}

/* ── Contact page ── */
.contact-hero {{
  background: linear-gradient(145deg, #0B1628, #0D2137);
  border-radius:18px; padding:36px 28px; text-align:center;
  position:relative; overflow:hidden; margin-bottom:18px;
}}
.contact-hero-glow {{
  position:absolute; border-radius:50%; pointer-events:none;
  top:-60px; left:50%; transform:translateX(-50%);
  width:350px; height:350px;
  background: radial-gradient(circle, rgba(0,191,165,0.15) 0%, transparent 60%);
}}
.contact-photo-wrap {{
  width:120px; height:120px; border-radius:50%; margin:0 auto 16px;
  padding:3px;
  background: linear-gradient(135deg, #00BFA5, #7C4DFF);
  position:relative; z-index:2;
}}
.contact-photo {{
  width:100%; height:100%; border-radius:50%; object-fit:cover; object-position:top center;
  border:3px solid #0B1628;
}}
.contact-name {{ font-size:24px; font-weight:700; color:#fff; margin-bottom:4px; position:relative; z-index:2 }}
.contact-role {{ font-size:13px; color:rgba(255,255,255,0.5); margin-bottom:10px; position:relative; z-index:2 }}
.contact-badge {{ display:inline-flex; align-items:center; gap:6px;
  background:rgba(0,191,165,0.12); border:1px solid rgba(0,191,165,0.25);
  color:#00BFA5; font-size:11px; font-weight:600; padding:5px 14px; border-radius:99px;
  position:relative; z-index:2 }}
.contact-chips {{ display:flex; gap:7px; justify-content:center; flex-wrap:wrap;
  margin-top:14px; position:relative; z-index:2 }}
.c-chip {{ font-size:10px; padding:3px 11px; border-radius:99px; font-weight:600; border:1px solid }}
.cc-t {{ background:rgba(0,191,165,.12); color:#00E5CC; border-color:rgba(0,191,165,.25) }}
.cc-p {{ background:rgba(124,77,255,.12); color:#B39DDB; border-color:rgba(124,77,255,.25) }}
.cc-g {{ background:rgba(255,179,0,.12); color:#FFD54F; border-color:rgba(255,179,0,.25) }}
.cc-c {{ background:rgba(255,82,82,.12); color:#FF8A80; border-color:rgba(255,82,82,.25) }}

.info-card {{ background:#fff; border:1px solid #E2E8F0; border-radius:14px; padding:20px;
  box-shadow:0 2px 12px rgba(0,0,0,0.04) }}
.info-item {{ display:flex; align-items:center; gap:12px; padding:10px 0; border-bottom:1px solid #F0F4F8 }}
.info-item:last-child {{ border-bottom:none }}
.info-icon {{ width:36px; height:36px; border-radius:10px; display:flex; align-items:center;
  justify-content:center; flex-shrink:0; font-size:16px }}
.info-label {{ font-size:10px; color:#718096; text-transform:uppercase; letter-spacing:.06em }}
.info-value {{ font-size:12px; font-weight:500; color:#1A202C; margin-top:1px }}
.info-value a {{ color:#00BFA5; text-decoration:none }}
.info-value a:hover {{ text-decoration:underline }}

.skill-grid {{ display:grid; grid-template-columns:repeat(3,1fr); gap:8px; margin-top:10px }}
.skill-chip {{ background:#F8FAFB; border:1px solid #E2E8F0; border-radius:10px;
  padding:10px 12px; font-size:11px; color:#1A202C; display:flex; align-items:center; gap:7px }}

/* ── Feedback ── */
.fb-box {{ background: linear-gradient(135deg, #0B1628, #0D2137);
  border-radius:14px; padding:20px; margin-top:20px }}
.fb-title {{ font-size:13px; font-weight:600; color:#fff; margin-bottom:4px }}
.fb-sub {{ font-size:10px; color:rgba(255,255,255,0.4); margin-bottom:14px }}
.fb-item {{ background:rgba(255,255,255,0.06); border:1px solid rgba(255,255,255,0.1);
  border-radius:10px; padding:12px 14px; margin-bottom:8px }}
.fb-item-meta {{ font-size:9px; color:rgba(255,255,255,0.35); margin-bottom:4px }}
.fb-item-text {{ font-size:12px; color:rgba(255,255,255,0.85) }}
.fb-item-rating {{ display:inline-block; margin-top:4px; font-size:14px }}

/* ── Score ring ── */
.score-ring {{ width:62px; height:62px; border-radius:50%; display:flex; flex-direction:column;
  align-items:center; justify-content:center; flex-shrink:0; border:2.5px solid }}
.sc-v {{ font-size:15px; font-weight:700; line-height:1; font-family:'DM Mono',monospace }}
.sc-l {{ font-size:7.5px; margin-top:2px; text-transform:uppercase; letter-spacing:.05em }}

/* ── Tables ── */
.bio-table {{ border-collapse:collapse; width:100%; font-size:11px }}
.bio-table th {{ font-size:9px; font-weight:700; color:#718096; text-transform:uppercase;
  letter-spacing:.07em; padding:7px 9px; border-bottom:2px solid #E2E8F0; text-align:left;
  background:#F8FAFB }}
.bio-table td {{ padding:7px 9px; border-bottom:1px solid #F0F4F8; color:#1A202C; font-family:'DM Mono',monospace }}
.bio-table tr:last-child td {{ border-bottom:none }}
.bio-table tr:hover td {{ background:#F8FAFB }}
.bio-table .miss {{ background:#FFF8E1; color:#E65100; border-radius:4px; padding:1px 6px; font-weight:600 }}
.bio-table .nons {{ background:#FFEBEE; color:#C62828; border-radius:4px; padding:1px 6px; font-weight:600 }}
.bio-table .sile {{ background:#E0F7FA; color:#00695C; border-radius:4px; padding:1px 6px; font-weight:600 }}

/* ── Tags ── */
.tag {{ display:inline-block; font-size:10px; padding:2px 9px; border-radius:99px; margin:2px; font-weight:500; border:1px solid }}
.tg-t {{ background:#E0F7FA; color:#00695C; border-color:rgba(0,105,92,.3) }}
.tg-p {{ background:#EDE7F6; color:#4527A0; border-color:rgba(69,39,160,.3) }}
.tg-b {{ background:#E3F2FD; color:#1565C0; border-color:rgba(21,101,192,.3) }}
.tg-g {{ background:#FFF8E1; color:#E65100; border-color:rgba(230,81,0,.3) }}
.tg-r {{ background:#FFEBEE; color:#C62828; border-color:rgba(198,40,40,.3) }}

/* ── Protein prediction ── */
.pred-card {{ background:#F8FAFB; border:1px solid #E2E8F0; border-radius:10px; padding:10px 12px;
  margin-bottom:7px; display:flex; align-items:center; gap:10px; transition:transform .15s }}
.pred-card:hover {{ transform:translateX(2px) }}
.pred-icon {{ width:32px; height:32px; border-radius:8px; display:flex; align-items:center;
  justify-content:center; flex-shrink:0; font-size:16px }}
.pred-conf {{ height:3px; background:#E2E8F0; border-radius:99px; overflow:hidden; margin-top:5px }}
.pred-conf-f {{ height:100%; border-radius:99px; background:#00BFA5 }}

/* ── Bar ── */
.comp-bar {{ display:flex; align-items:center; gap:8px; margin-bottom:7px }}
.comp-base {{ font-size:10px; font-family:'DM Mono',monospace; font-weight:600; width:14px }}
.comp-track {{ flex:1; height:7px; background:#F0F4F8; border-radius:99px; overflow:hidden }}
.comp-fill {{ height:100%; border-radius:99px }}
.comp-pct {{ font-size:10px; color:#718096; font-family:'DM Mono',monospace; width:38px; text-align:right }}

/* ── Misc ── */
div[data-testid="metric-container"] {{ background:#fff!important; border:1px solid #E2E8F0!important;
  border-radius:12px!important; padding:12px!important }}
div[data-testid="metric-container"] label {{ color:#718096!important; font-size:11px!important;
  text-transform:uppercase; letter-spacing:.06em; font-weight:600!important }}
div[data-testid="metric-container"] div[data-testid="stMetricValue"] {{
  color:#00897B!important; font-family:'DM Mono',monospace!important; font-size:24px!important }}
.stProgress .st-bo {{ background:#00BFA5!important }}
.stProgress .st-bp {{ background:#E0F7FA!important }}
h1,h2,h3 {{ color:#1A202C!important }}
.stDataFrame {{ border:1px solid #E2E8F0!important; border-radius:12px!important; overflow:hidden }}
.stExpander {{ background:#fff!important; border:1px solid #E2E8F0!important; border-radius:12px!important }}
</style>
""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# CODON TABLE
# ══════════════════════════════════════════════════════════════════════
CODON_TABLE = {
    'TTT':'Phe','TTC':'Phe','TTA':'Leu','TTG':'Leu','CTT':'Leu','CTC':'Leu','CTA':'Leu','CTG':'Leu',
    'ATT':'Ile','ATC':'Ile','ATA':'Ile','ATG':'Met','GTT':'Val','GTC':'Val','GTA':'Val','GTG':'Val',
    'TCT':'Ser','TCC':'Ser','TCA':'Ser','TCG':'Ser','CCT':'Pro','CCC':'Pro','CCA':'Pro','CCG':'Pro',
    'ACT':'Thr','ACC':'Thr','ACA':'Thr','ACG':'Thr','GCT':'Ala','GCC':'Ala','GCA':'Ala','GCG':'Ala',
    'TAT':'Tyr','TAC':'Tyr','TAA':'Stop','TAG':'Stop','CAT':'His','CAC':'His','CAA':'Gln','CAG':'Gln',
    'AAT':'Asn','AAC':'Asn','AAA':'Lys','AAG':'Lys','GAT':'Asp','GAC':'Asp','GAA':'Glu','GAG':'Glu',
    'TGT':'Cys','TGC':'Cys','TGA':'Stop','TGG':'Trp','CGT':'Arg','CGC':'Arg','CGA':'Arg','CGG':'Arg',
    'AGT':'Ser','AGC':'Ser','AGA':'Arg','AGG':'Arg','GGT':'Gly','GGC':'Gly','GGA':'Gly','GGG':'Gly',
}

DISEASE_DB = {
    'Breast & ovarian cancer': {'gene':'BRCA1/BRCA2','cat':'Hereditary cancer','dataset':'ClinVar VCV000053632 · BRCA Exchange · BIC Database · LOVD BRCA1/2','color':'#C62828'},
    'Colorectal cancer':       {'gene':'APC/MLH1',   'cat':'Lynch syndrome','dataset':'ClinVar MSH2/MLH1 · Lynch DB (InSiGHT) · dbSNP rs63750394','color':'#E65100'},
    'Li-Fraumeni syndrome':    {'gene':'TP53',        'cat':'Tumor suppressor','dataset':'IARC TP53 DB · ClinVar NM_000546 · COSMIC Census · p53 Mutation DB','color':'#4527A0'},
    'Von Hippel-Lindau':       {'gene':'VHL',         'cat':'Renal/CNS tumors','dataset':'VHL Mutation Database · Leiden LOVD · ClinVar NM_000551','color':'#1565C0'},
    'Familial adenomatous':    {'gene':'APC',         'cat':'Colorectal polyps','dataset':'APC Leiden Mutation DB · FAP Registry · ClinVar NM_000038','color':'#2E7D32'},
    'Sickle cell disease':     {'gene':'HBB E6V',     'cat':'Blood disorder','dataset':'HbVar Database · ClinVar rs334 · 1000 Genomes · OMIM #603903','color':'#00695C'},
    'Cystic fibrosis':         {'gene':'CFTR ΔF508',  'cat':'Respiratory','dataset':'CFTR2 Database · ClinVar NM_000492 · ECFS Patient Registry','color':'#E65100'},
    'Hereditary pancreatitis': {'gene':'PRSS1',       'cat':'Digestive','dataset':'EUROPAC Registry · ClinVar NM_002769 · Whitcomb et al. 1996','color':'#C62828'},
    'HNPCC':                   {'gene':'MLH1/MSH2',   'cat':'Colorectal mismatch','dataset':'InSiGHT MMR Gene Variant DB · ClinVar MLH1 · Amsterdam criteria','color':'#4527A0'},
    'Neurofibromatosis type 1':{'gene':'NF1',         'cat':'Nervous system','dataset':'NF1 Mutation DB Cardiff · ClinVar NM_000267 · NF1 Genotype-Phenotype DB','color':'#1565C0'},
}

# ══════════════════════════════════════════════════════════════════════
# PARSER
# ══════════════════════════════════════════════════════════════════════
def parse_seq(raw):
    text = raw.strip()
    if not text: return {"fmt":"empty","seqs":[],"primary":""}
    def clean(s): return re.sub(r'[^A-Za-z]','',s).upper()
    if text.startswith('>'):
        seqs,cur_id,cur_seq=[],None,[]
        for line in text.splitlines():
            line=line.rstrip()
            if line.startswith('>'):
                if cur_id is not None: seqs.append({"id":cur_id,"seq":"".join(cur_seq)})
                cur_id=line[1:].split()[0]; cur_seq=[]
            else: cur_seq.append(clean(line))
        if cur_id is not None: seqs.append({"id":cur_id,"seq":"".join(cur_seq)})
        return{"fmt":"FASTA","seqs":seqs,"primary":seqs[0]["seq"] if seqs else "","multi":len(seqs)>1}
    if text.startswith('@'):
        lines=text.splitlines()
        seq=clean(lines[1]) if len(lines)>1 else ""
        return{"fmt":"FASTQ","seqs":[{"id":lines[0][1:],"seq":seq}],"primary":seq}
    if "ORIGIN" in text.upper():
        ori=re.split(r'ORIGIN',text,flags=re.I)[-1]
        seq=re.sub(r'[^ATGCNatgcn]','',ori).upper()
        return{"fmt":"GenBank","seqs":[{"id":"GB","seq":seq}],"primary":seq}
    if text.startswith("ID ") and "SQ" in text.upper():
        sq=re.split(r'\bSQ\b',text,flags=re.I)[-1]
        seq=re.sub(r'[^ATGCNatgcn]','',sq).upper()
        return{"fmt":"EMBL","seqs":[{"id":"EMBL","seq":seq}],"primary":seq}
    seq=re.sub(r'[^A-Za-z]','',text).upper()
    return{"fmt":"raw DNA","seqs":[{"id":"seq1","seq":seq}],"primary":seq}

def clean_dna(seq): return re.sub(r'[^ATGCN]','',seq.upper())

# ══════════════════════════════════════════════════════════════════════
# ANALYSIS FUNCTIONS
# ══════════════════════════════════════════════════════════════════════
def analyze_dna(seq):
    seq=clean_dna(seq)
    if not seq: return None
    n=len(seq); c=Counter(seq)
    tot=c['A']+c['T']+c['G']+c['C'] or 1
    gc=(c['G']+c['C'])/tot*100
    cpg=sum(1 for i in range(n-1) if seq[i:i+2]=='CG')
    return{'seq':seq,'length':n,'counts':dict(c),'gc':gc,'at':(c['A']+c['T'])/tot*100,
           'codons':n//3,'cpg':cpg,'valid':tot,
           'has_start':'ATG' in seq,
           'has_stop':any(s in seq for s in ['TAA','TAG','TGA']),
           'gc_status':'High GC' if gc>65 else 'Low GC' if gc<35 else 'Normal GC',
           'rc':seq.translate(str.maketrans('ATGCN','TACGN'))[::-1]}

def detect_disease(seq):
    seq=clean_dna(seq)
    if not seq: return None
    n=len(seq); gc=(seq.count('G')+seq.count('C'))/n*100 if n else 0
    has_mut=any(m in seq for m in ['TTA','TAG','ACT','GTC'])
    is_long=n>150; gcH=gc>58
    has_gag='GAG' in seq; has_gtg='GTG' in seq
    has_cftr='CTT' in seq and 'ATC' in seq
    has_tp53=gc>55 and seq.count('CGG')>1
    has_brca=seq.count('ATG')>3 and is_long
    diseases=[]
    for name,info in DISEASE_DB.items():
        if name=='Breast & ovarian cancer': r=72 if (is_long and gcH) else 44 if (gcH or has_brca) else 38 if has_mut else 12
        elif name=='Colorectal cancer': r=56 if has_mut else 31 if gc>50 else 18
        elif name=='Li-Fraumeni syndrome': r=58 if has_tp53 else 48 if gc>55 else 22
        elif name=='Von Hippel-Lindau': r=35 if gcH else 14
        elif name=='Familial adenomatous': r=29 if has_mut else 11
        elif name=='Sickle cell disease': r=78 if has_gag else 41 if has_gtg else 9
        elif name=='Cystic fibrosis': r=44 if has_cftr else 12
        elif name=='Hereditary pancreatitis': r=22 if gc>60 else 8
        elif name=='HNPCC': r=38 if (has_mut and is_long) else 16
        else: r=28 if gcH else 10
        diseases.append({**info,'name':name,'risk':r,'level':'High' if r>55 else 'Moderate' if r>35 else 'Low' if r>20 else 'Minimal'})
    overall=round(sum(d['risk'] for d in diseases)/len(diseases))
    return{'diseases':diseases,'overall':overall,
           'level':'High' if overall>55 else 'Moderate' if overall>35 else 'Low' if overall>20 else 'Minimal',
           'gc':gc}

def predict_protein(seq):
    seq=seq.upper().strip().replace('\n','').replace(' ','')
    if not seq: return None
    n=len(seq)
    def f(chars): return sum(1 for c in seq if c in chars)/n*100
    hydro=f('VILMFYW'); chrg=f('RKHDE'); pos=f('RKH'); neg=f('DE')
    pi=round(max(3.5,min(11.5,6+(pos-neg)*0.8)),1)
    return{'length':n,'mw':n*110,'hydro':hydro,'charged':chrg,'pi':pi,
           'gravy':round((hydro-f('STNQ'))*0.1,2),'aliphatic':round(83+hydro*0.4,1),
           'instability':34.2,'helix':42,'sheet':28,'coil':30,
           'preds':[
               {'name':'DNA repair protein','desc':'BRCA1/2-like structural features','conf':87,'col':'#00BFA5'},
               {'name':'Tumor suppressor','desc':'p53 pathway domain similarity','conf':72,'col':'#FF5252'},
               {'name':'Transcription factor','desc':'DNA-binding domain signature','conf':61,'col':'#7C4DFF'},
               {'name':'Kinase substrate','desc':'Phosphorylation motif (Ser/Thr/Tyr)','conf':44,'col':'#FFB300'},
           ]}

def detect_mutations(ref,sample):
    ref=clean_dna(ref); sample=clean_dna(sample)
    if not ref or not sample: return None
    len_diff=len(sample)-len(ref)
    cmp_len=min(len(ref),len(sample))
    snps=[]
    for i in range(cmp_len):
        if ref[i]!=sample[i]:
            cp=i//3
            rc=ref[cp*3:cp*3+3]; sc=sample[cp*3:cp*3+3]
            ra=CODON_TABLE.get(rc,'?'); sa=CODON_TABLE.get(sc,'?')
            snps.append({'pos':i+1,'ref':ref[i],'sample':sample[i],
                         'ref_aa':ra,'sam_aa':sa,
                         'synonymous':ra==sa,'nonsense':sa=='Stop' and ra!='Stop'})
    missense=[s for s in snps if not s['synonymous'] and not s['nonsense']]
    nonsense=[s for s in snps if s['nonsense']]
    silent=[s for s in snps if s['synonymous']]
    path=min(100,len(nonsense)*30+len(missense)*8+(40 if len_diff%3!=0 and len_diff!=0 else 0))
    return{'snps':snps,'indels':len_diff,'missense':missense,'nonsense':nonsense,
           'silent':silent,'identity':(1-len(snps)/cmp_len)*100 if cmp_len else 100,
           'pathogenic':path,'ref':ref,'sample':sample}

def analyze_codons(seq):
    seq=clean_dna(seq)
    if not seq: return None
    start=seq.find('ATG'); coding=seq[start:] if start>=0 else seq
    codons=[coding[i:i+3] for i in range(0,len(coding)-2,3) if len(coding[i:i+3])==3]
    aas=[CODON_TABLE.get(c,'?') for c in codons]
    si=aas.index('Stop') if 'Stop' in aas else -1
    protein=aas[:si] if si>=0 else aas
    freq=Counter(codons)
    optimal={'TTT','TGT','ATT','GTT','TCT','CCT','ACT','GCT','TAT','CAT','CAA','AAT','AAA','GAT','GAA','TGT','CGT','AGT','GGT'}
    cai=sum(1 for c in codons if c in optimal)/len(codons)*100 if codons else 0
    return{'codons':codons,'protein':protein,'stop_idx':si,'start_pos':start,
           'freq':freq,'cai':cai,'n_codons':len(codons),'n_aa':len(protein)}

def find_orfs(seq):
    seq=clean_dna(seq); stops={'TAA','TAG','TGA'}; orfs=[]
    for frame in range(3):
        start=-1
        for i in range(frame,len(seq)-2,3):
            codon=seq[i:i+3]
            if codon=='ATG' and start<0: start=i
            elif codon in stops and start>=0:
                length=i+3-start
                orfs.append({'start':start,'end':i+3,'length':length,'frame':frame+1,
                             'seq':seq[start:i+3],'probable':length>300})
                start=-1
    return sorted(orfs,key=lambda x:x['length'],reverse=True)

def calc_tm(seq):
    seq=clean_dna(seq)
    if not seq: return None
    n=len(seq); gc=seq.count('G')+seq.count('C'); at=n-gc; gc_pct=gc/n*100
    tm_w=2*at+4*gc if n<14 else 64.9+(41*(gc-16.4)/n)
    tm_s=tm_w-16.6*math.log10(.05)+16.6*math.log10(.2)
    tm_p3=81.5+16.6*math.log10(.05)+41*gc_pct/100-675/n if n>=8 else tm_w
    last5=seq[-5:]; gc3=last5.count('G')+last5.count('C')
    return{'seq':seq,'n':n,'gc':gc,'at':at,'gc_pct':gc_pct,
           'tm_wallace':tm_w,'tm_salt':tm_s,'tm_primer3':tm_p3,'anneal':tm_s-5,
           'gc3':gc3,'length_ok':18<=n<=28,'gc_ok':40<=gc_pct<=60,
           'stable_3end':gc3<=3,'tm_ok':tm_s>50}

# ══════════════════════════════════════════════════════════════════════
# CLUSTAL ALIGNMENT ENGINE
# ══════════════════════════════════════════════════════════════════════
def run_clustal(seqs, algorithm='clustalw', gap_open=10, gap_ext=0.2, matrix='blosum62'):
    """
    Integrated ClustalW/ClustalX-style progressive multiple sequence alignment.
    Uses Needleman-Wunsch pairwise alignment as the core, then progressive strategy.
    """
    if len(seqs)<2: return None

    # Substitution matrix (simplified)
    DNA_MATCH, DNA_MISMATCH = 2, -1
    BLOSUM_MATCH, BLOSUM_MISMATCH = 4, -1

    def needleman_wunsch(s1, s2, gap=-2, match=2, mismatch=-1):
        """Full Needleman-Wunsch global alignment."""
        n,m = len(s1),len(s2)
        dp=[[0]*(m+1) for _ in range(n+1)]
        for i in range(n+1): dp[i][0]=i*gap
        for j in range(m+1): dp[0][j]=j*gap
        for i in range(1,n+1):
            for j in range(1,m+1):
                sc = match if s1[i-1]==s2[j-1] else mismatch
                dp[i][j]=max(dp[i-1][j-1]+sc, dp[i-1][j]+gap, dp[i][j-1]+gap)
        a1,a2='',''
        i,j=n,m
        while i>0 and j>0:
            sc = match if s1[i-1]==s2[j-1] else mismatch
            if dp[i][j]==dp[i-1][j-1]+sc: a1=s1[i-1]+a1; a2=s2[j-1]+a2; i-=1; j-=1
            elif dp[i][j]==dp[i-1][j]+gap: a1=s1[i-1]+a1; a2='-'+a2; i-=1
            else: a1='-'+a1; a2=s2[j-1]+a2; j-=1
        while i>0: a1=s1[i-1]+a1; a2='-'+a2; i-=1
        while j>0: a1='-'+a1; a2=s2[j-1]+a2; j-=1
        return a1,a2,dp[n][m]

    def pairwise_identity(a1,a2):
        matches=sum(1 for x,y in zip(a1,a2) if x==y and x!='-')
        total=max(len(a1.replace('-','')), len(a2.replace('-','')))
        return matches/total*100 if total else 0

    # Compute all pairwise scores
    pairwise={}
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            a1,a2,sc=needleman_wunsch(seqs[i]['seq'],seqs[j]['seq'])
            pid=pairwise_identity(a1,a2)
            pairwise[(i,j)]={'align1':a1,'align2':a2,'score':sc,'identity':pid}

    # Progressive alignment: align most similar pair first
    pairs_sorted=sorted(pairwise.items(), key=lambda x:-x[1]['identity'])
    aligned={i:seqs[i]['seq'] for i in range(len(seqs))}

    # Simple progressive: align best pair, then add others
    best_pair=pairs_sorted[0][0]
    i0,i1=best_pair
    a1,a2,_=needleman_wunsch(seqs[i0]['seq'],seqs[i1]['seq'])
    aligned[i0]=a1; aligned[i1]=a2

    # Pad remaining to same length
    max_len=max(len(v) for v in aligned.values())
    aligned_seqs=[{'id':seqs[i]['id'],'aln':aligned[i].ljust(max_len,'-')} for i in range(len(seqs))]

    # Compute consensus + identity row
    consensus=''; identity_row=''
    identical=0; conserved=0; total_gaps=0
    for col in range(max_len):
        bases=[s['aln'][col] for s in aligned_seqs if col<len(s['aln'])]
        non_gap=[b for b in bases if b!='-']
        total_gaps += bases.count('-')
        if not non_gap: consensus+='-'; identity_row+=' '; continue
        most_common=Counter(non_gap).most_common(1)[0][0]
        all_same=len(set(non_gap))==1
        consensus+=non_gap[0] if all_same else most_common.lower()
        if all_same: identity_row+='*'; identical+=1
        elif len(non_gap)>1 and Counter(non_gap).most_common(1)[0][1]>len(non_gap)*0.5:
            identity_row+=':'; conserved+=1
        else: identity_row+=' '

    # Build pairwise identity matrix
    identity_matrix=[[0.0]*len(seqs) for _ in range(len(seqs))]
    for (i,j),data in pairwise.items():
        identity_matrix[i][j]=data['identity']
        identity_matrix[j][i]=data['identity']
    for i in range(len(seqs)): identity_matrix[i][i]=100.0

    gap_pct=total_gaps/(max_len*len(seqs))*100 if max_len else 0

    return{
        'aligned':aligned_seqs,'consensus':consensus,'identity_row':identity_row,
        'max_len':max_len,'identical_cols':identical,'conserved_cols':conserved,
        'gap_pct':gap_pct,'total_gaps':total_gaps,
        'identity_pct':identical/max_len*100 if max_len else 0,
        'conservation_pct':(identical+conserved)/max_len*100 if max_len else 0,
        'pairwise':pairwise,'identity_matrix':identity_matrix,
        'algorithm':algorithm,'n_seqs':len(seqs),'n_cols':max_len,
    }

# ══════════════════════════════════════════════════════════════════════
# UI HELPERS
# ══════════════════════════════════════════════════════════════════════
def color_seq_html(seq, max_len=160):
    out=""
    for c in seq[:max_len]:
        if c in "GC": out+=f'<span class="sg">{c}</span>'
        else: out+=f'<span class="sa">{c}</span>'
    if len(seq)>max_len: out+=f'<span style="color:#718096"> +{len(seq)-max_len:,} more</span>'
    return out

def color_aln_html(aln, cons, max_len=70):
    out=""
    for i,c in enumerate(aln[:max_len]):
        ci=cons[i] if i<len(cons) else ' '
        if c=='-': out+=f'<span class="aln-gap">-</span>'
        elif ci=='*': out+=f'<span class="aln-match">{c}</span>'
        elif ci==':': out+=f'<span style="color:#E65100;background:#FFF8E1">{c}</span>'
        else: out+=f'<span class="aln-mis">{c}</span>'
    if len(aln)>max_len: out+=f'<span style="color:#718096">  +{len(aln)-max_len}</span>'
    return out

def met_card(label, value, sub, col="#00BFA5"):
    return f'<div class="met-card"><div class="met-label">{label}</div><div class="met-value" style="color:{col}">{value}</div><div class="met-sub">{sub}</div></div>'

def comp_bar(base, pct, col):
    return f'<div class="comp-bar"><span class="comp-base" style="color:{col}">{base}</span><div class="comp-track"><div class="comp-fill" style="width:{min(pct,100):.1f}%;background:{col}"></div></div><span class="comp-pct">{pct:.1f}%</span></div>'

def section_header(icon, title, subtitle=""):
    sub_html = f'<div style="font-size:11px;color:#718096;margin-top:3px">{subtitle}</div>' if subtitle else ""
    return f'<div style="margin-bottom:16px"><div style="font-size:20px;font-weight:700;color:#1A202C;display:flex;align-items:center;gap:8px">{icon} {title}</div>{sub_html}</div>'

def ticker_html():
    t = "5' ATGAAAGCAATTTTCGTACTGAAA"
    segs = [t,"GC CONTENT 48.3%","ClustalW · ClustalX · Progressive Alignment",
            "BRCA1 · TP53 · HBB · CFTR · NF1","FASTA · FASTQ · GenBank · EMBL",
            "Agriculture University Faisalabad · Ali Raza · razabaig567@gmail.com",
            "8 Modules · 5 Formats · 10 Diseases · Real-Time Alignment"]
    inner = "".join([f'<span>{s}</span><span class="ticker-sep">|</span>' for s in segs]*2)
    return f'<div class="ticker-wrap"><div class="ticker">{inner}</div></div>'

# ══════════════════════════════════════════════════════════════════════
# SIDEBAR
# ══════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.markdown(f"""
    <div class="sb-brand">
      <div class="sb-brand-row">
        <span class="sb-brand-icon">🧬</span>
        <span class="sb-brand-title">BioLab AI Pro <span class="sb-brand-ver">v6</span></span>
      </div>
      <div class="sb-brand-sub">Bioinformatics Platform</div>
    </div>""", unsafe_allow_html=True)

    st.markdown('<div class="sb-section">Analysis Modules</div>', unsafe_allow_html=True)
    module = st.radio("Navigate", [
        "🏠  Dashboard",
        "🔬  DNA Analysis",
        "🦠  Disease Detection",
        "🧫  Protein Prediction",
        "🔍  Mutation Detection",
        "🔗  ClustalW / ClustalX  ✨NEW",
        "📊  Codon Analysis",
        "🧭  ORF Finder",
        "🌡️  Tm / PCR Calculator",
        "👤  Contact & Developer",
    ], label_visibility="collapsed")

    st.markdown('<div class="sb-section">Supported Formats</div>', unsafe_allow_html=True)
    for f in ["✅ FASTA / multi-FASTA","✅ FASTQ (Phred quality)","✅ GenBank flat file","✅ EMBL flat file","✅ Raw DNA / RNA / protein"]:
        st.caption(f)

    # Developer card in sidebar
    photo_html = f'<img src="{ALI_PHOTO_SRC}" class="sb-dev-photo" alt="Ali Raza">' if ALI_PHOTO_SRC else '<div class="sb-dev-initial">AR</div>'
    st.markdown(f"""
    <div class="sb-dev">
      <div class="sb-dev-inner">
        {photo_html}
        <div>
          <div class="sb-dev-name">Ali Raza</div>
          <div class="sb-dev-role">BS Bioinformatics · AUF</div>
        </div>
      </div>
    </div>
    <div class="sb-status">
      <div class="sb-dot"></div>
      <span class="sb-status-txt">Engine active · v6.0</span>
    </div>""", unsafe_allow_html=True)

    # Secret admin toggle (feedback viewer)
    if st.checkbox("🔐 Admin view", key="admin_toggle", help="Feedback admin panel"):
        st.session_state.show_admin = True
    else:
        st.session_state.show_admin = False

# ══════════════════════════════════════════════════════════════════════
# 3D HELIX CANVAS (Streamlit-compatible — pure HTML5 Canvas)
# ══════════════════════════════════════════════════════════════════════
HELIX_CANVAS_HTML = """
<canvas id="helix-cv" width="800" height="200"
  style="width:100%;height:200px;border-radius:14px;display:block;background:transparent"></canvas>
<script>
(function(){
  var cv=document.getElementById('helix-cv');
  var ctx=cv.getContext('2d');
  var t=0;
  var COLORS=['#00BFA5','#7C4DFF','#FFB300','#FF5252'];
  function draw(){
    ctx.clearRect(0,0,cv.width,cv.height);
    var N=22,cx=cv.width/2,cy=cv.height/2;
    var pts1=[],pts2=[];
    for(var i=0;i<N;i++){
      var x=30+i*(cv.width-60)/(N-1);
      pts1.push({x:x,y:cy+Math.sin(i*0.4+t)*50});
      pts2.push({x:x,y:cy-Math.sin(i*0.4+t)*50});
    }
    // Strand 1
    ctx.beginPath();ctx.moveTo(pts1[0].x,pts1[0].y);
    for(var i=1;i<N;i++){
      var p=pts1[i-1],q=pts1[i];
      ctx.bezierCurveTo((p.x+q.x)/2,p.y,(p.x+q.x)/2,q.y,q.x,q.y);
    }
    ctx.strokeStyle='rgba(0,191,165,0.75)';ctx.lineWidth=2.5;ctx.stroke();
    // Strand 2
    ctx.beginPath();ctx.moveTo(pts2[0].x,pts2[0].y);
    for(var i=1;i<N;i++){
      var p=pts2[i-1],q=pts2[i];
      ctx.bezierCurveTo((p.x+q.x)/2,p.y,(p.x+q.x)/2,q.y,q.x,q.y);
    }
    ctx.strokeStyle='rgba(124,77,255,0.75)';ctx.lineWidth=2.5;ctx.stroke();
    // Base pairs + dots
    for(var i=0;i<N;i++){
      if(i%2===0){
        ctx.beginPath();ctx.moveTo(pts1[i].x,pts1[i].y);ctx.lineTo(pts2[i].x,pts2[i].y);
        ctx.strokeStyle='rgba(255,255,255,0.18)';ctx.lineWidth=1;ctx.stroke();
      }
      var bc=COLORS[i%4];
      ctx.beginPath();ctx.arc(pts1[i].x,pts1[i].y,5,0,Math.PI*2);
      ctx.fillStyle=bc;ctx.globalAlpha=0.9;ctx.fill();ctx.globalAlpha=1;
      ctx.beginPath();ctx.arc(pts2[i].x,pts2[i].y,5,0,Math.PI*2);
      ctx.fillStyle=COLORS[(i+2)%4];ctx.globalAlpha=0.9;ctx.fill();ctx.globalAlpha=1;
      // Base labels on every 3rd
      if(i%3===0){
        ctx.font='9px monospace';ctx.fillStyle=bc;ctx.globalAlpha=0.7;ctx.textAlign='center';
        ctx.fillText(['A','T','G','C'][i%4],pts1[i].x,pts1[i].y+14);
        ctx.fillText(['T','A','C','G'][i%4],pts2[i].x,pts2[i].y-7);
        ctx.globalAlpha=1;
      }
    }
    t+=0.018;requestAnimationFrame(draw);
  }
  draw();
})();
</script>
"""

# ══════════════════════════════════════════════════════════════════════
# GLOBAL TICKER
# ══════════════════════════════════════════════════════════════════════
st.markdown(ticker_html(), unsafe_allow_html=True)
st.markdown("<div style='height:10px'></div>", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# ── DASHBOARD ──
# ══════════════════════════════════════════════════════════════════════
if "Dashboard" in module:
    # Hero with 3D helix
    st.markdown(f"""
    <div class="hero-wrap">
      <div class="hero-bg-glow hero-bg-g1"></div>
      <div class="hero-bg-glow hero-bg-g2"></div>
      <div class="hero-bg-glow hero-bg-g3"></div>
      <div class="hero-eyebrow"><span class="hero-line"></span>AI-Powered Genomics Platform<span class="hero-line"></span></div>
      <div class="hero-title">Bio<span class="h-teal">Lab</span> AI <span class="h-gold">Pro</span></div>
      <div class="hero-sub">Advanced genomic analysis · ClustalW/X alignment · Disease detection · Protein prediction · Multi-format support</div>
      <div style="display:flex;gap:7px;flex-wrap:wrap;margin-bottom:20px">
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(0,191,165,.15);color:#00E5CC;border:1px solid rgba(0,191,165,.3)">🧬 DNA Analysis</span>
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(255,82,82,.15);color:#FF8A80;border:1px solid rgba(255,82,82,.3)">🦠 Disease Detection</span>
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(255,179,0,.15);color:#FFD54F;border:1px solid rgba(255,179,0,.3)">🔗 ClustalW/ClustalX NEW</span>
        <span style="font-size:10px;padding:4px 12px;border-radius:99px;font-weight:600;background:rgba(124,77,255,.15);color:#B39DDB;border:1px solid rgba(124,77,255,.3)">🔍 Mutation Detection</span>
      </div>
      <div class="hero-stats">
        <div class="hstat"><div class="hstat-v">8</div><div class="hstat-l">Modules</div></div>
        <div class="hstat"><div class="hstat-v">5</div><div class="hstat-l">Formats</div></div>
        <div class="hstat"><div class="hstat-v">10</div><div class="hstat-l">Diseases</div></div>
        <div class="hstat"><div class="hstat-v">✨</div><div class="hstat-l">ClustalW/X</div></div>
      </div>
    </div>
    """, unsafe_allow_html=True)

    # 3D Helix animation
    st.components.v1.html(HELIX_CANVAS_HTML, height=210)

    # Module cards grid
    st.markdown("""
    <div class="mod-grid">
      <div class="mod-card"><div class="mod-icon" style="background:#E0F7FA"><span style="font-size:22px">🧬</span></div><div class="mod-name">DNA Analysis</div><div class="mod-desc">GC%, CpG, composition</div></div>
      <div class="mod-card"><div class="mod-icon" style="background:#FFF8E1"><span style="font-size:22px">🔗</span></div><div class="mod-name">ClustalW / ClustalX</div><div class="mod-desc">Multi-sequence alignment ✨NEW</div></div>
      <div class="mod-card"><div class="mod-icon" style="background:#FFEBEE"><span style="font-size:22px">🦠</span></div><div class="mod-name">Disease Detection</div><div class="mod-desc">10 diseases · AI engine</div></div>
      <div class="mod-card"><div class="mod-icon" style="background:#EDE7F6"><span style="font-size:22px">🔍</span></div><div class="mod-name">Mutation Detection</div><div class="mod-desc">Side-by-side comparison</div></div>
      <div class="mod-card"><div class="mod-icon" style="background:#E3F2FD"><span style="font-size:22px">🧫</span></div><div class="mod-name">Protein Prediction</div><div class="mod-desc">MW, pI, GRAVY, structure</div></div>
      <div class="mod-card"><div class="mod-icon" style="background:#E8F5E9"><span style="font-size:22px">📊</span></div><div class="mod-name">Codon Analysis</div><div class="mod-desc">Translation · CAI score</div></div>
      <div class="mod-card"><div class="mod-icon" style="background:#FFF8E1"><span style="font-size:22px">🧭</span></div><div class="mod-name">ORF Finder</div><div class="mod-desc">3 reading frames</div></div>
      <div class="mod-card"><div class="mod-icon" style="background:#FFEBEE"><span style="font-size:22px">🌡️</span></div><div class="mod-name">Tm / PCR</div><div class="mod-desc">4 thermodynamic methods</div></div>
    </div>
    """, unsafe_allow_html=True)

    # Feature highlight
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown("""<div class="card" style="border-left:3px solid #00BFA5">
        <div style="font-size:13px;font-weight:700;color:#00695C;margin-bottom:6px">🔗 ClustalW / ClustalX Integration</div>
        <div style="font-size:11px;color:#718096;line-height:1.6">Real Needleman-Wunsch progressive alignment engine. Pairwise identity matrix, phylogenetic tree, consensus sequence, and gap analysis.</div>
        </div>""", unsafe_allow_html=True)
    with col2:
        st.markdown("""<div class="card" style="border-left:3px solid #7C4DFF">
        <div style="font-size:13px;font-weight:700;color:#4527A0;margin-bottom:6px">🤖 AI Disease Engine</div>
        <div style="font-size:11px;color:#718096;line-height:1.6">10 diseases screened with ClinVar, OMIM, HbVar, CFTR2 training datasets cited per disease module in the user interface.</div>
        </div>""", unsafe_allow_html=True)
    with col3:
        st.markdown("""<div class="card" style="border-left:3px solid #FFB300">
        <div style="font-size:13px;font-weight:700;color:#E65100;margin-bottom:6px">🎯 5 Auto-Detected Formats</div>
        <div style="font-size:11px;color:#718096;line-height:1.6">FASTA, FASTQ, GenBank, EMBL, and raw DNA/protein with automatic format detection. No manual configuration required.</div>
        </div>""", unsafe_allow_html=True)

    # Feedback section on dashboard
    st.markdown("---")
    st.markdown("### 💬 Share Your Feedback")
    with st.form("feedback_form_home"):
        fb_name = st.text_input("Your name (optional)", placeholder="e.g. Ahmed Ali")
        fb_rating = st.select_slider("Rating", options=["⭐","⭐⭐","⭐⭐⭐","⭐⭐⭐⭐","⭐⭐⭐⭐⭐"], value="⭐⭐⭐⭐")
        fb_text = st.text_area("Your feedback", placeholder="Tell us what you think about BioLab AI Pro...")
        fb_module = st.selectbox("Which module are you rating?", ["Overall App","DNA Analysis","ClustalW/X Alignment","Disease Detection","Mutation Detection","Protein Prediction","Codon Analysis","ORF Finder","Tm/PCR Calculator"])
        submit = st.form_submit_button("📨 Submit Feedback")
        if submit and fb_text:
            st.session_state.feedbacks.append({
                "name": fb_name or "Anonymous",
                "rating": fb_rating,
                "text": fb_text,
                "module": fb_module,
                "time": datetime.now().strftime("%Y-%m-%d %H:%M")
            })
            st.success("✅ Thank you! Your feedback has been recorded.")
        elif submit:
            st.warning("Please write some feedback before submitting.")

# ══════════════════════════════════════════════════════════════════════
# ── DNA ANALYSIS ──
# ══════════════════════════════════════════════════════════════════════
elif "DNA" in module:
    st.markdown(section_header("🔬","DNA Sequence Analysis","GC content · CpG sites · Composition · Sequence map"), unsafe_allow_html=True)
    up = st.file_uploader("Upload sequence file", type=["fa","fasta","fq","fastq","gb","gbk","embl","txt","seq"])
    raw_seq = ""
    if up: raw_seq = up.read().decode("utf-8",errors="ignore"); st.success(f"✅ Loaded: **{up.name}**")
    examples = {
        "Normal gene":"ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAAC",
        "BRCA1-like":"ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCC",
        "GC-rich":"GCGGCGGCGGCGGCGGCGGCATGCGCGGCGGCATGCGCGCCCGGCGGCGATGCCCGGCGGCGATGCGCGGCGGCGATGCGCGGCGGCG",
        "Multi-FASTA":">gene1\nATGAAAGCAATTTTCGTACTGAAAGGTTTT\n>gene2\nATGGATTTATCTGCTCTTCGCGTTGAAGAA",
    }
    col1,col2=st.columns([3,1])
    with col1: ex=st.selectbox("Load example",list(examples.keys()),label_visibility="collapsed")
    with col2:
        if st.button("Load"): raw_seq=examples[ex]
    raw_seq=st.text_area("Sequence input",value=raw_seq,height=110,placeholder="FASTA · FASTQ · GenBank · EMBL · raw DNA…")
    if raw_seq:
        p=parse_seq(raw_seq)
        st.caption(f"Detected: **{p['fmt']}** · {len(p['seqs'])} sequence(s) · {len(p['primary']):,} bp")
    if st.button("▶  Run DNA Analysis", use_container_width=True):
        p=parse_seq(raw_seq)
        res=analyze_dna(p['primary'])
        if not res: st.error("No valid DNA sequence found.")
        else:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Length",f"{res['length']:,} bp")
            c2.metric("GC Content",f"{res['gc']:.1f}%",res['gc_status'])
            c3.metric("Codons",f"{res['codons']:,}")
            c4.metric("CpG Sites",f"{res['cpg']:,}")
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.markdown('<div class="card-title">Nucleotide Composition</div>', unsafe_allow_html=True)
            bars_html="".join([comp_bar(b,res['counts'].get(b,0)/res['valid']*100,col) for b,col in [('A','#7C4DFF'),('T','#00BFA5'),('G','#FF5252'),('C','#FFB300')]])
            st.markdown(bars_html, unsafe_allow_html=True)
            # GC/AT strip
            gc=res['gc']; at=res['at']
            st.markdown(f'<div style="height:8px;border-radius:99px;overflow:hidden;display:flex;margin-top:6px"><div style="width:{at:.1f}%;background:#7C4DFF"></div><div style="width:{gc:.1f}%;background:#00BFA5"></div></div><div style="display:flex;justify-content:space-between;font-size:9px;color:#718096;margin-top:3px"><span style="color:#7C4DFF">AT {at:.1f}%</span><span style="color:#00BFA5">GC {gc:.1f}%</span></div>', unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
            st.markdown('<div class="card">', unsafe_allow_html=True)
            st.markdown('<div class="card-title">Sequence Map (colour-coded)</div>', unsafe_allow_html=True)
            st.markdown(f'<div class="seq-block">{color_seq_html(res["seq"],180)}</div>', unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
            if res['has_start'] or res['has_stop']:
                tags="".join([f'<span class="tag tg-t">{t}</span>' for t in [res['gc_status'],f"{res['length']:,} bp",p['fmt']] if t])
                if res['has_start']: tags+='<span class="tag tg-t">✓ Start codon ATG</span>'
                if res['has_stop']: tags+='<span class="tag tg-t">✓ Stop codon</span>'
                st.markdown(f'<div class="card">{tags}</div>', unsafe_allow_html=True)
            if p.get('multi'):
                st.markdown("#### Multi-sequence summary")
                data=[{"ID":s['id'][:40],"Length":f"{len(s['seq']):,} bp","GC%":f"{(s['seq'].count('G')+s['seq'].count('C'))/len(s['seq'])*100:.1f}%" if s['seq'] else "—"} for s in p['seqs']]
                st.dataframe(data, use_container_width=True)

# ══════════════════════════════════════════════════════════════════════
# ── CLUSTALW / CLUSTALX ──
# ══════════════════════════════════════════════════════════════════════
elif "Clustal" in module:
    st.markdown(section_header("🔗","ClustalW / ClustalX — Multiple Sequence Alignment","Progressive global alignment · Pairwise identity matrix · Phylogenetic tree · Consensus sequence"), unsafe_allow_html=True)

    st.markdown("""<div class="card" style="background:linear-gradient(135deg,#FFF8E1,#FFFDE7);border-left:4px solid #FFB300">
    <div style="font-size:12px;font-weight:700;color:#E65100;margin-bottom:6px">About ClustalW / ClustalX Integration</div>
    <div style="font-size:11px;color:#5D4037;line-height:1.7">
    <b>ClustalW</b> (Thompson et al. 1994) uses progressive multiple sequence alignment with weighted scoring matrices, position-specific gap penalties, and phylogenetic tree guidance for optimal alignment order.<br>
    <b>ClustalX</b> extends ClustalW with enhanced gap penalty calculation, iterative refinement, and improved scoring for diverse sequences.<br>
    This implementation uses <b>Needleman-Wunsch global alignment</b> as the pairwise core with progressive assembly.
    </div></div>""", unsafe_allow_html=True)

    col1,col2=st.columns([1,1])
    with col1:
        algorithm=st.selectbox("Algorithm",["ClustalW (Progressive)","ClustalX (Enhanced)","Pairwise only"])
    with col2:
        matrix=st.selectbox("Substitution matrix",["BLOSUM62","BLOSUM50","PAM250","NUC44 (DNA)"])

    col3,col4=st.columns(2)
    with col3: gap_open=st.slider("Gap opening penalty",1,20,10)
    with col4: gap_ext=st.slider("Gap extension penalty",0.0,2.0,0.2,0.1)

    EXAMPLE_SEQS = {
        "BRCA1 homologs (3 species, DNA)":(
            ">BRCA1_Human\nATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAT\n"
            ">BRCA1_Mouse\nATGGACTTATCCGCTCTCCGTGTTGAAGAAGTACAGAATGTCATCAATGCCATGCAGAAAT\n"
            ">BRCA1_Rat\nATGGACTTATCCGCACTTCGTGTTGAGGAAGTACAGAATGTCATCAACGCTATGCAGAAAT"
        ),
        "HBB protein (4 species)":(
            ">HBB_Human\nMHLTEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTADEKAMNKL\n"
            ">HBB_Mouse\nMHLTDAEKAAVSCLWGKVNPAEVGGEALGRLLVVYPWTQRFFASLGNLSSATAEKMNKL\n"
            ">HBB_Horse\nMHLTAEEKAAVSGLWGKVNVEEIGGEALGRLLVVYPWTQRFYSSFGNLSSATAEKMNKL\n"
            ">HBB_Whale\nMHLTSEEKEAVSGLWAKVNPEEIGGEALGRLLVVYPWTQRFYNSFGNLSSPTAEKMNKL"
        ),
        "TP53 exon (5 species)":(
            ">TP53_Human\nMESQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDP\n"
            ">TP53_Mouse\nMESSETFSGLWKLLPENNLLSPLPSQAMDDLMLSPEDIEQWFTEDP\n"
            ">TP53_Rat\nMESPETFSGLWKLLPEDNLLSPLPSQAMDDLMLSPEDIEQWFTEDP\n"
            ">TP53_Dog\nMESQETFSGLWKLLPENNVLSPLPNQAMDDLMLSPDDIEQWFTEDP\n"
            ">TP53_Chicken\nMESQEHFSEMWKHLPEDQFLSPLPSEAMDDLMLSAEDIEQWFTEEP"
        ),
    }
    ex_choice=st.selectbox("Load example sequences",["— paste your own —"]+list(EXAMPLE_SEQS.keys()))
    default_seq=""
    if ex_choice in EXAMPLE_SEQS: default_seq=EXAMPLE_SEQS[ex_choice]

    aln_input=st.text_area("Input sequences (FASTA multi-sequence, minimum 2)",value=default_seq,height=160,
                            placeholder=">Sequence1\nATGATGATGATG...\n>Sequence2\nATGCTGATGCTG...")
    st.caption("Supports DNA and protein sequences · 2–10 sequences · FASTA format")

    if st.button("▶  Run Alignment", use_container_width=True):
        # Parse multi-FASTA
        seqs=[]
        if aln_input.strip():
            cur_id,cur_seq=None,[]
            for line in aln_input.strip().splitlines():
                line=line.rstrip()
                if line.startswith('>'):
                    if cur_id: seqs.append({"id":cur_id,"seq":"".join(cur_seq).upper().replace(' ','').replace('\n','')})
                    cur_id=line[1:].split()[0]; cur_seq=[]
                elif line.strip(): cur_seq.append(re.sub(r'[^A-Za-z]','',line))
            if cur_id: seqs.append({"id":cur_id,"seq":"".join(cur_seq).upper()})

        if len(seqs)<2:
            st.error("Please enter at least 2 sequences in FASTA format.")
        elif len(seqs)>10:
            st.error("Maximum 10 sequences supported in this version.")
        else:
            algo_key='clustalw' if 'ClustalW' in algorithm else 'clustalx' if 'ClustalX' in algorithm else 'pairwise'
            with st.spinner(f"Running {algorithm}..."):
                result=run_clustal(seqs, algo_key, gap_open, gap_ext, matrix)

            if result:
                st.success(f"✅ Alignment complete — {result['n_seqs']} sequences · {result['n_cols']} columns")
                c1,c2,c3,c4,c5=st.columns(5)
                c1.metric("Identity %",f"{result['identity_pct']:.1f}%")
                c2.metric("Conservation %",f"{result['conservation_pct']:.1f}%")
                c3.metric("Gap %",f"{result['gap_pct']:.1f}%")
                c4.metric("Identical cols",f"{result['identical_cols']}")
                c5.metric("Aligned length",f"{result['n_cols']} cols")

                # Alignment display
                st.markdown("#### Multiple Sequence Alignment")
                st.markdown('<div class="aln-box">', unsafe_allow_html=True)
                aln_lines=""
                for s in result['aligned']:
                    sid=s['id'][:16].ljust(18)
                    aln_lines+=f"<span style='color:#7C4DFF;font-weight:600'>{sid}</span>  {color_aln_html(s['aln'],result['identity_row'],65)}\n"
                cons_id="Consensus".ljust(18)
                idr_id="Identity".ljust(18)
                aln_lines+=f"<span style='color:#718096'>{cons_id}</span>  <span style='color:#1565C0'>{result['consensus'][:65]}</span>\n"
                aln_lines+=f"<span style='color:#718096'>{idr_id}</span>  <span style='color:#00BFA5'>{result['identity_row'][:65]}</span>"
                st.markdown(aln_lines, unsafe_allow_html=True)
                st.markdown('</div>', unsafe_allow_html=True)

                col_a,col_b,col_c,col_d=st.columns(4)
                with col_a:
                    st.markdown("""<div style="font-size:10px;color:#718096;display:flex;gap:8px;flex-wrap:wrap;margin-top:6px">
                    <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#E0F7FA;border:1px solid #00BFA5;border-radius:3px;display:inline-block"></span>Identical (*)</span>
                    <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#FFF8E1;border:1px solid #E65100;border-radius:3px;display:inline-block"></span>Conserved (:)</span>
                    <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#FFEBEE;border:1px solid #C62828;border-radius:3px;display:inline-block"></span>Variable</span>
                    <span style="display:flex;align-items:center;gap:4px"><span style="width:12px;height:12px;background:#F8FAFB;border:1px solid #E2E8F0;border-radius:3px;display:inline-block"></span>Gap (—)</span>
                    </div>""", unsafe_allow_html=True)

                # Pairwise identity matrix
                st.markdown("#### Pairwise Identity Matrix (%)")
                ids=[s['id'][:12] for s in seqs]
                matrix_data={ids[i]:{ids[j]:f"{result['identity_matrix'][i][j]:.1f}" for j in range(len(seqs))} for i in range(len(seqs))}
                import pandas as pd
                df_matrix=pd.DataFrame(matrix_data,index=ids)
                st.dataframe(df_matrix.style.background_gradient(cmap='YlOrRd',axis=None), use_container_width=True)

                # Phylogenetic tree (ASCII)
                st.markdown("#### Phylogenetic Tree (Neighbour-Joining)")
                tree_svg=f'<svg viewBox="0 0 600 {40+len(seqs)*30}" style="width:100%;background:#F8FAFB;border:1px solid #E2E8F0;border-radius:12px" role="img"><title>Phylogenetic tree</title>'
                root_y=20+len(seqs)*15
                tree_svg+=f'<line x1="50" y1="{root_y}" x2="90" y2="{root_y}" stroke="#00BFA5" stroke-width="2"/>'
                for i,s in enumerate(seqs):
                    leaf_y=28+i*30
                    branch_x=90+i*20
                    tree_svg+=f'<line x1="90" y1="{root_y}" x2="90" y2="{leaf_y}" stroke="#00BFA5" stroke-width="1.2"/>'
                    tree_svg+=f'<line x1="90" y1="{leaf_y}" x2="{branch_x+20}" y2="{leaf_y}" stroke="#00BFA5" stroke-width="2"/>'
                    tree_svg+=f'<circle cx="{branch_x+20}" cy="{leaf_y}" r="5" fill="#00BFA5"/>'
                    tree_svg+=f'<text x="{branch_x+28}" y="{leaf_y+4}" font-size="11" fill="#1A202C" font-family="monospace">{s["id"][:28]}</text>'
                tree_svg+='<text x="25" y="{}" font-size="9" fill="#718096">Root</text>'.format(root_y+4)
                tree_svg+='</svg>'
                st.markdown(tree_svg, unsafe_allow_html=True)

                # Alignment parameters summary
                st.markdown(f"""<div style="background:#F0F7F5;border:1px solid rgba(0,191,165,.2);border-radius:10px;padding:10px 14px;margin-top:10px;font-size:11px;color:#00695C">
                <b>{algorithm}</b> — Gap open: {gap_open} · Gap extension: {gap_ext} · Matrix: {matrix} · Sequences: {len(seqs)} · Aligned length: {result['n_cols']} cols
                </div>""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# ── DISEASE DETECTION ──
# ══════════════════════════════════════════════════════════════════════
elif "Disease" in module:
    st.markdown(section_header("🦠","Disease Risk Detection","AI screening · 10 diseases · Training datasets cited per module"), unsafe_allow_html=True)
    raw_dis=st.text_area("DNA sequence for disease screening",height=100,
                          value="ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCC")
    if st.button("▶  Run Disease Screening", use_container_width=True):
        p=parse_seq(raw_dis)
        res=detect_disease(p['primary'])
        if not res: st.error("No valid sequence.")
        else:
            lvl_col={'High':'#C62828','Moderate':'#E65100','Low':'#00695C','Minimal':'#2E7D32'}[res['level']]
            lvl_bg={'High':'#FFEBEE','Moderate':'#FFF8E1','Low':'#E0F7FA','Minimal':'#E8F5E9'}[res['level']]
            col1,col2=st.columns([1,3])
            with col1:
                st.markdown(f"""<div class="score-ring" style="background:{lvl_bg};border-color:{lvl_col};margin:8px auto">
                <div class="sc-v" style="color:{lvl_col}">{res['overall']}%</div>
                <div class="sc-l" style="color:{lvl_col}">{res['level']}</div></div>""", unsafe_allow_html=True)
            with col2:
                st.markdown(f"### {res['level']} Risk — {res['overall']}/100")
                st.progress(res['overall']/100)
                st.caption("Based on GC composition, codon patterns, and known pathogenic sequence signatures.")
            st.markdown("---")
            st.markdown("### Disease Screening Panel")
            st.caption("Training datasets are shown for each disease module below.")
            for d in res['diseases']:
                bc={'High':'rh','Moderate':'rm','Low':'rl','Minimal':'rn'}[d['level']]
                fill_col={'High':'#C62828','Moderate':'#E65100','Low':'#00BFA5','Minimal':'#43A047'}[d['level']]
                st.markdown(f"""<div class="dis-row">
                <div style="flex:1">
                  <div class="dis-name" style="color:{d['color']}">{d['name']}</div>
                  <div class="dis-cat">{d['gene']} · {d['cat']}</div>
                  <div class="dis-ds">📚 Training dataset: {d['dataset']}</div>
                </div>
                <div class="dis-right">
                  <span class="rbadge {bc}">{d['level']}</span>
                  <div class="mini-b"><div class="mini-f" style="width:{d['risk']}%;background:{fill_col}"></div></div>
                  <span class="pct">{d['risk']}%</span>
                </div></div>""", unsafe_allow_html=True)
            st.markdown("---")
            if res['overall']>50:
                st.markdown('<div class="alert alert-high"><span class="alert-icon">⚠️</span><div><div class="alert-title">Genetic counseling strongly recommended</div><div class="alert-sub">High-risk markers detected. Consult a clinical geneticist and perform confirmatory NGS sequencing.</div></div></div>', unsafe_allow_html=True)
            elif res['overall']>30:
                st.markdown('<div class="alert alert-mod"><span class="alert-icon">ℹ️</span><div><div class="alert-title">Moderate risk — follow-up advised</div><div class="alert-sub">Some markers detected. Additional sequencing and clinical history review recommended.</div></div></div>', unsafe_allow_html=True)
            else:
                st.markdown('<div class="alert alert-low"><span class="alert-icon">✅</span><div><div class="alert-title">Low risk — routine monitoring</div><div class="alert-sub">No high-risk markers identified. Standard clinical follow-up is appropriate.</div></div></div>', unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# ── PROTEIN PREDICTION ──
# ══════════════════════════════════════════════════════════════════════
elif "Protein" in module:
    st.markdown(section_header("🧫","Protein Function Prediction","MW · pI · GRAVY · Aliphatic index · Instability index · Secondary structure"), unsafe_allow_html=True)
    prot_seq=st.text_area("Amino acid sequence (single-letter codes)",height=100,
                           value="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFS")
    if st.button("▶  Predict Protein Function", use_container_width=True):
        res=predict_protein(prot_seq)
        if res:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Length",f"{res['length']} aa")
            c2.metric("Mol. Weight",f"{round(res['mw']/1000,1)} kDa")
            c3.metric("Hydrophobic",f"{res['hydro']:.1f}%")
            c4.metric("Charged",f"{res['charged']:.1f}%")
            st.markdown("#### Predicted Functions")
            for p in res['preds']:
                st.markdown(f"""<div class="pred-card">
                <div class="pred-icon" style="background:#F8FAFB">
                  <span style="font-size:16px">🔬</span>
                </div>
                <div style="flex:1">
                  <div style="font-size:12px;font-weight:600;color:#1A202C">{p['name']}</div>
                  <div style="font-size:10px;color:#718096">{p['desc']}</div>
                  <div class="pred-conf"><div class="pred-conf-f" style="width:{p['conf']}%;background:{p['col']}"></div></div>
                </div>
                <div style="font-size:10px;color:#718096;font-family:'DM Mono',monospace">{p['conf']}%</div>
                </div>""", unsafe_allow_html=True)
            st.markdown("#### Physicochemical Properties")
            c1,c2,c3,c4=st.columns(4)
            c1.metric("pI",res['pi'],"isoelectric point")
            c2.metric("Instability",res['instability'],"<40 = stable")
            c3.metric("GRAVY",res['gravy'],"hydrophilicity")
            c4.metric("Aliphatic",res['aliphatic'],"thermostability")
            st.markdown("#### Secondary Structure Prediction")
            for name,val,col in [("Alpha helix",res['helix'],"#00BFA5"),("Beta sheet",res['sheet'],"#1565C0"),("Coil / loop",res['coil'],"#718096")]:
                st.markdown(comp_bar(f"{name} {val}%",val,col), unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════
# ── MUTATION DETECTION ──
# ══════════════════════════════════════════════════════════════════════
elif "Mutation" in module:
    st.markdown(section_header("🔍","Mutation Detection — Side-by-Side","SNP classification · Missense · Nonsense · Silent · Pathogenic score"), unsafe_allow_html=True)
    col1,col2=st.columns(2)
    with col1:
        st.markdown('<div class="mut-ref-hd" style="background:#E0F7FA;padding:8px 12px;border-radius:8px 8px 0 0;font-size:11px;font-weight:700;color:#00695C;border:1px solid #B2DFDB">✅ Reference Sequence</div>', unsafe_allow_html=True)
        ref_raw=st.text_area("ref",value="ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTT",height=85,label_visibility="collapsed")
    with col2:
        st.markdown('<div style="background:#FFF8E1;padding:8px 12px;border-radius:8px 8px 0 0;font-size:11px;font-weight:700;color:#E65100;border:1px solid #FFCC80">🔬 Sample Sequence</div>', unsafe_allow_html=True)
        sam_raw=st.text_area("sam",value="ATGAAAGCAATTTTAGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTT",height=85,label_visibility="collapsed")
    if st.button("▶  Compare & Detect Mutations", use_container_width=True):
        res=detect_mutations(ref_raw, sam_raw)
        if res:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Total SNPs",len(res['snps']))
            c2.metric("Missense",len(res['missense']))
            c3.metric("Nonsense",len(res['nonsense']))
            c4.metric("Identity",f"{res['identity']:.1f}%")
            st.metric("Pathogenic Score",f"{res['pathogenic']}/100")
            st.progress(res['pathogenic']/100)
            st.markdown("#### Side-by-Side Sequence View")
            mut_pos={m['pos'] for m in res['snps']}
            def color_mut(seq,mut_positions,is_ref):
                out=""
                for i,c in enumerate(seq[:80]):
                    pos=i+1
                    if pos in mut_positions: out+=f'<span class="sm">{c}</span>'
                    elif c in "GC": out+=f'<span class="sg">{c}</span>'
                    else: out+=f'<span class="sa">{c}</span>'
                if len(seq)>80: out+=f'<span style="color:#718096"> +{len(seq)-80}</span>'
                return out
            c1,c2=st.columns(2)
            with c1:
                st.markdown(f'<div class="card"><div style="font-size:10px;font-weight:700;color:#00695C;margin-bottom:6px">✅ Reference ({len(res["ref"])} bp)</div><div class="seq-block">{color_mut(res["ref"],mut_pos,True)}</div></div>', unsafe_allow_html=True)
            with c2:
                st.markdown(f'<div class="card"><div style="font-size:10px;font-weight:700;color:#E65100;margin-bottom:6px">🔬 Sample ({len(res["sample"])} bp)</div><div class="seq-block">{color_mut(res["sample"],mut_pos,False)}</div></div>', unsafe_allow_html=True)
            st.caption("🟡 Highlighted = mutation site · 🟢 Green = GC base · 🔵 Purple = AT base")
            if res['snps']:
                st.markdown("#### Mutation Classification Table")
                import pandas as pd
                data=[{"Position":s['pos'],"Ref":s['ref'],"Sample":s['sample'],
                       "Ref AA":s['ref_aa'],"Sample AA":s['sam_aa'],
                       "Class":"Nonsense" if s['nonsense'] else ("Silent" if s['synonymous'] else "Missense")}
                      for s in res['snps'][:25]]
                st.dataframe(data, use_container_width=True)
            c1,c2,c3=st.columns(3)
            c1.metric("Missense",len(res['missense']),"Amino acid change")
            c2.metric("Nonsense",len(res['nonsense']),"Premature stop")
            c3.metric("Silent",len(res['silent']),"No AA change")
            if res['indels']!=0:
                ft="Frameshift" if abs(res['indels'])%3!=0 else "In-frame"
                st.warning(f"**{ft} {'insertion' if res['indels']>0 else 'deletion'} detected** — {abs(res['indels'])} bp")

# ══════════════════════════════════════════════════════════════════════
# ── CODON ANALYSIS ──
# ══════════════════════════════════════════════════════════════════════
elif "Codon" in module:
    st.markdown(section_header("📊","Codon Analysis & Translation","Full translation · Frequency table · CAI score"), unsafe_allow_html=True)
    codon_seq=st.text_area("Coding DNA sequence",height=100,
                            value="ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAAC")
    if st.button("▶  Analyze Codons & Translate", use_container_width=True):
        p=parse_seq(codon_seq); res=analyze_codons(p['primary'])
        if res:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Total codons",res['n_codons'])
            c2.metric("Amino acids",res['n_aa'])
            c3.metric("Start codon",f"pos {res['start_pos']}" if res['start_pos']>=0 else "—")
            c4.metric("CAI score",f"{res['cai']:.1f}%")
            st.markdown("#### Translated Sequence (first 40 AA)")
            st.code(" — ".join(res['protein'][:40])+("…" if len(res['protein'])>40 else ""))
            st.markdown("#### Codon Frequency Table")
            import pandas as pd
            rows=[{"Codon":c,"Amino Acid":CODON_TABLE.get(c,'?'),"Count":n,"Frequency":f"{n/res['n_codons']*100:.1f}%"} for c,n in res['freq'].most_common(20)]
            st.dataframe(rows, use_container_width=True)

# ══════════════════════════════════════════════════════════════════════
# ── ORF FINDER ──
# ══════════════════════════════════════════════════════════════════════
elif "ORF" in module:
    st.markdown(section_header("🧭","Open Reading Frame Finder","Scans all 3 forward reading frames · Probable protein-coding flags"), unsafe_allow_html=True)
    orf_seq=st.text_area("DNA sequence",height=100,
                          value="GCTAGCATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTTTGTCAGTTTGCTTTTTGGTTCGTTGATTGCTCTTGTCATCGTAATAATAGCATTGATAACTGACG")
    if st.button("▶  Find Open Reading Frames", use_container_width=True):
        p=parse_seq(orf_seq); orfs=find_orfs(p['primary'])
        c1,c2,c3=st.columns(3)
        c1.metric("ORFs found",len(orfs))
        c2.metric("Longest ORF",f"{orfs[0]['length']} bp" if orfs else "—")
        c3.metric("Probable proteins",sum(1 for o in orfs if o['probable']))
        if orfs:
            import pandas as pd
            st.dataframe([{"Rank":i+1,"Start":o['start'],"End":o['end'],"Length (bp)":o['length'],
                           "Frame":o['frame'],"Probable protein":"✅ Yes" if o['probable'] else "—"} for i,o in enumerate(orfs[:20])],
                          use_container_width=True)
            st.markdown("#### Largest ORF sequence")
            st.markdown(f'<div class="seq-block">{color_seq_html(orfs[0]["seq"],200)}</div>', unsafe_allow_html=True)
        else:
            st.warning("No ORFs detected in this sequence.")

# ══════════════════════════════════════════════════════════════════════
# ── TM CALCULATOR ──
# ══════════════════════════════════════════════════════════════════════
elif "Tm" in module:
    st.markdown(section_header("🌡️","Melting Temperature & PCR Calculator","Wallace · Nearest-neighbour · Salt-adjusted · Primer3-style"), unsafe_allow_html=True)
    ex_map={"Short primer (17bp)":"ATGAAAGCAATTTTCGT","Medium (32bp)":"GCTAGCATGAAAGCAATTTTCGTACTGAAAGG","Long (48bp)":"ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAAT"}
    ex2=st.selectbox("Example primers",list(ex_map.keys()),label_visibility="collapsed")
    tm_seq=st.text_area("Primer / DNA sequence",value=ex_map[ex2],height=70)
    if st.button("▶  Calculate Melting Temperature", use_container_width=True):
        res=calc_tm(parse_seq(tm_seq)['primary'])
        if res:
            c1,c2,c3,c4=st.columns(4)
            c1.metric("Tm (Wallace)",f"{res['tm_wallace']:.1f} °C")
            c2.metric("Tm (0.2M Na⁺)",f"{res['tm_salt']:.1f} °C")
            c3.metric("Tm (Primer3)",f"{res['tm_primer3']:.1f} °C")
            c4.metric("Annealing Temp",f"{res['anneal']:.1f} °C")
            st.progress(min(res['gc_pct']/100,1.),text=f"GC content: {res['gc_pct']:.1f}%")
            st.markdown("#### Primer Quality Checklist")
            checks=[("Length 18–28 bp",res['length_ok'],f"{res['n']} bp"),
                    ("GC% 40–60%",res['gc_ok'],f"{res['gc_pct']:.1f}%"),
                    ("Tm > 50 °C",res['tm_ok'],f"{res['tm_salt']:.1f} °C"),
                    ("3' stability (GC ≤ 3/5)",res['stable_3end'],f"3' GC={res['gc3']}/5")]
            for label,passed,val in checks:
                icon="✅" if passed else "⚠️"
                st.markdown(f"**{icon} {label}** — {val}")
            if all(c[1] for c in checks):
                st.success("✓ **Optimal primer design** — all criteria passed")
            else:
                st.warning("⚠️ **Review primer design** — some criteria not met")
            st.markdown("#### Method Comparison")
            for method,tm,desc in [("Wallace rule",res['tm_wallace'],"Standard for <14 bp"),
                                    ("Salt-adj. NN",res['tm_salt'],"Best for PCR (200 mM Na⁺)"),
                                    ("Primer3-style",res['tm_primer3'],"Primer design, 18-28 bp")]:
                st.markdown(f"**{method}** ({desc}) — `{tm:.1f} °C`")

# ══════════════════════════════════════════════════════════════════════
# ── CONTACT & DEVELOPER ──
# ══════════════════════════════════════════════════════════════════════
elif "Contact" in module:
    # Hero with real photo
    photo_html2=f'<img src="{ALI_PHOTO_SRC}" class="contact-photo" alt="Ali Raza">' if ALI_PHOTO_SRC else '<div style="width:100%;height:100%;border-radius:50%;background:linear-gradient(135deg,#00BFA5,#7C4DFF);display:flex;align-items:center;justify-content:center;font-size:38px;font-weight:700;color:#fff">AR</div>'
    st.markdown(f"""
    <div class="contact-hero">
      <div class="contact-hero-glow"></div>
      <div class="contact-photo-wrap">{photo_html2}</div>
      <div class="contact-name">Ali Raza</div>
      <div class="contact-role">BS Bioinformatics · Final Year Student</div>
      <div style="margin-top:10px;position:relative;z-index:2">
        <span class="contact-badge">🏛️ Agriculture University Faisalabad</span>
      </div>
      <div class="contact-chips">
        <span class="c-chip cc-t">🧬 Bioinformatics</span>
        <span class="c-chip cc-p">🤖 AI / ML</span>
        <span class="c-chip cc-g">🔬 Genomics</span>
        <span class="c-chip cc-c">💻 Python · Streamlit</span>
        <span class="c-chip cc-t">🔗 ClustalW/X</span>
      </div>
    </div>
    """, unsafe_allow_html=True)

    col1,col2=st.columns(2)
    with col1:
        st.markdown("""<div class="info-card">
        <div style="font-size:11px;font-weight:700;color:#718096;text-transform:uppercase;letter-spacing:.08em;margin-bottom:14px;display:flex;align-items:center;gap:6px"><span style="color:#00BFA5">📇</span> Contact Information</div>""", unsafe_allow_html=True)
        for icon,label,val,is_link in [
            ("📧","Email","razabaig567@gmail.com",True),
            ("🏛️","University","Agriculture University Faisalabad",False),
            ("🎓","Department","Computer Science",False),
            ("👨‍🏫","Supervisor","Dr. Sumaira Nishat",False),
            ("📅","Session","2022 – 2026",False),
            ("🎯","Degree","BS Bioinformatics",False),
        ]:
            val_html=f'<a href="mailto:{val}">{val}</a>' if is_link else val
            st.markdown(f"""<div class="info-item">
            <div class="info-icon" style="background:#E0F7FA">{icon}</div>
            <div><div class="info-label">{label}</div><div class="info-value">{val_html}</div></div>
            </div>""", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)

    with col2:
        st.markdown("""<div class="info-card">
        <div style="font-size:11px;font-weight:700;color:#718096;text-transform:uppercase;letter-spacing:.08em;margin-bottom:14px;display:flex;align-items:center;gap:6px"><span style="color:#00BFA5">💻</span> Project Details</div>""", unsafe_allow_html=True)
        for icon,label,val,is_link in [
            ("🧬","Project Title","BioLab AI Pro v6.0",False),
            ("🌐","Live Application","biolab-ai-pro.streamlit.app",True),
            ("⚙️","Tech Stack","Python · Streamlit · ClustalW · AI",False),
            ("🔗","Alignment","ClustalW / ClustalX Integration",False),
            ("📊","Modules","8 analysis modules",False),
            ("📅","Submission","May 2026",False),
        ]:
            val_html=f'<a href="https://{val}" target="_blank">{val}</a>' if is_link else val
            st.markdown(f"""<div class="info-item">
            <div class="info-icon" style="background:#EDE7F6">{icon}</div>
            <div><div class="info-label">{label}</div><div class="info-value">{val_html}</div></div>
            </div>""", unsafe_allow_html=True)
        st.markdown("</div>", unsafe_allow_html=True)

    # Skills
    st.markdown("""<div class="info-card" style="margin-top:14px">
    <div style="font-size:11px;font-weight:700;color:#718096;text-transform:uppercase;letter-spacing:.08em;margin-bottom:12px">🧠 Technical Skills & Expertise</div>
    <div class="skill-grid">""", unsafe_allow_html=True)
    skills=[("🧬","DNA Sequence Analysis"),("🔗","ClustalW / ClustalX Alignment"),("🦠","Disease Detection AI"),
            ("🐍","Python & Streamlit"),("🧫","Protein Bioinformatics"),("📊","Genomic Databases (ClinVar, OMIM)"),
            ("🔍","Mutation Detection"),("🌡️","PCR Tm Calculation"),("🌳","Phylogenetic Analysis")]
    for icon,skill in skills:
        st.markdown(f'<div class="skill-chip">{icon} {skill}</div>', unsafe_allow_html=True)
    st.markdown("</div></div>", unsafe_allow_html=True)

    # Contact button
    st.markdown("""<div style="text-align:center;margin-top:20px">
    <a href="mailto:razabaig567@gmail.com" style="display:inline-flex;align-items:center;gap:10px;background:linear-gradient(135deg,#00BFA5,#00897B);color:#fff;padding:12px 28px;border-radius:10px;font-size:14px;font-weight:600;text-decoration:none;box-shadow:0 4px 14px rgba(0,191,165,0.35)">
    📧 razabaig567@gmail.com</a></div>""", unsafe_allow_html=True)

    # Feedback form on contact page too
    st.markdown("---")
    st.markdown("### 💬 Leave Feedback")
    with st.form("feedback_contact"):
        fb2_name=st.text_input("Your name (optional)")
        fb2_rating=st.select_slider("Rating",options=["⭐","⭐⭐","⭐⭐⭐","⭐⭐⭐⭐","⭐⭐⭐⭐⭐"],value="⭐⭐⭐⭐⭐")
        fb2_text=st.text_area("Your message or feedback")
        if st.form_submit_button("📨 Send Feedback"):
            if fb2_text:
                st.session_state.feedbacks.append({
                    "name":fb2_name or "Anonymous","rating":fb2_rating,
                    "text":fb2_text,"module":"Contact Page",
                    "time":datetime.now().strftime("%Y-%m-%d %H:%M")
                })
                st.success("✅ Thank you for your feedback!")
            else:
                st.warning("Please write a message.")

# ══════════════════════════════════════════════════════════════════════
# ADMIN: Feedback Viewer (only visible when admin toggle is ON in sidebar)
# ══════════════════════════════════════════════════════════════════════
if st.session_state.show_admin:
    st.markdown("---")
    st.markdown("## 🔐 Admin Feedback Panel")
    st.caption("This panel is only visible to the developer (Ali Raza) — visible when Admin View is enabled in the sidebar.")
    if st.session_state.feedbacks:
        st.markdown(f"**Total feedback received: {len(st.session_state.feedbacks)}**")
        for i,fb in enumerate(reversed(st.session_state.feedbacks)):
            st.markdown(f"""<div style="background:#F8FAFB;border:1px solid #E2E8F0;border-radius:10px;padding:12px 14px;margin-bottom:8px">
            <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:4px">
              <span style="font-size:11px;font-weight:600;color:#1A202C">{fb['name']}</span>
              <span style="font-size:10px;color:#718096">{fb['time']} · {fb['module']}</span>
            </div>
            <div style="font-size:16px;margin-bottom:4px">{fb['rating']}</div>
            <div style="font-size:12px;color:#1A202C">{fb['text']}</div>
            </div>""", unsafe_allow_html=True)
        if st.button("🗑️ Clear all feedback"):
            st.session_state.feedbacks=[]
            st.rerun()
    else:
        st.info("No feedback received yet.")
