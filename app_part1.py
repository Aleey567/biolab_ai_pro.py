"""BioLab AI Pro v6 — Part 1: Config, CSS, Helpers"""
import streamlit as st, re, math, json
from collections import Counter
from datetime import datetime
import pandas as pd

st.set_page_config(page_title="BioLab AI Pro v6",page_icon="🧬",layout="wide",initial_sidebar_state="expanded")

ALI_PHOTO_B64 = open('/home/claude/ali_photo_b64.txt').read().strip()
ALI_PHOTO_SRC = f"data:image/jpeg;base64,{ALI_PHOTO_B64}"

# ── Session state init ─────────────────────────────────────────────
for k,v in [("feedbacks",[]),("show_admin",False),("module","🏠  Dashboard"),
             ("dna_result",None),("dis_result",None),("prot_result",None),
             ("mut_result",None),("aln_result",None),("tm_result",None),
             ("orf_result",None),("codon_result",None),("blast_result",None)]:
    if k not in st.session_state: st.session_state[k]=v

CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;500&display=swap');
html,body,[class*="css"]{font-family:'Inter',sans-serif!important}
.main,.block-container,.stApp{background:#F7F9FC!important}

/* ── SIDEBAR ── */
section[data-testid="stSidebar"]{background:#111827!important;border-right:1px solid #1F2937!important}
section[data-testid="stSidebar"] .stRadio label{
  display:flex;align-items:center;gap:8px;padding:8px 12px!important;
  border-radius:8px!important;color:#9CA3AF!important;font-size:12px!important;
  font-weight:500!important;cursor:pointer;transition:all .15s;margin-bottom:2px!important;
  background:transparent!important;border:none!important;
}
section[data-testid="stSidebar"] .stRadio label:hover{background:#1F2937!important;color:#F9FAFB!important}
section[data-testid="stSidebar"] .stRadio label:has(input:checked){
  background:#064E3B!important;color:#6EE7B7!important;font-weight:600!important;
  border-left:3px solid #10B981!important;
}
section[data-testid="stSidebar"] .stRadio label p{font-size:12px!important;margin:0!important}
[data-testid="stRadio"] div{gap:1px!important}

/* ── TOP NAV ── */
.topnav{background:#fff;border-bottom:1px solid #E5E7EB;padding:12px 0 10px;margin:-1rem 0 1.5rem;
  display:flex;align-items:center;justify-content:space-between}
.topnav-brand{display:flex;align-items:center;gap:10px}
.topnav-logo{width:32px;height:32px;border-radius:8px;background:linear-gradient(135deg,#10B981,#059669);
  display:flex;align-items:center;justify-content:center;font-size:18px;flex-shrink:0}
.topnav-title{font-size:17px;font-weight:800;color:#111827;letter-spacing:-.02em}
.topnav-ver{font-size:10px;font-weight:700;background:#ECFDF5;color:#059669;padding:2px 7px;border-radius:99px;margin-left:6px}
.topnav-crumb{font-size:12px;color:#6B7280}

/* ── CARDS ── */
.card{background:#fff;border:1px solid #E5E7EB;border-radius:14px;padding:20px;margin-bottom:14px;
  box-shadow:0 1px 6px rgba(0,0,0,.04)}
.card:hover{box-shadow:0 4px 18px rgba(0,0,0,.07);transition:box-shadow .2s}
.card-accent{height:4px;border-radius:14px 14px 0 0;margin:-20px -20px 18px}
.card-title{font-size:11px;font-weight:700;color:#6B7280;text-transform:uppercase;letter-spacing:.09em;
  margin-bottom:14px;display:flex;align-items:center;gap:7px}

/* ── HERO ── */
.hero{background:linear-gradient(135deg,#0F172A 0%,#1E293B 50%,#0F172A 100%);
  border-radius:16px;padding:32px;margin-bottom:20px;position:relative;overflow:hidden}
.hero-g1{position:absolute;top:-40px;right:-40px;width:250px;height:250px;border-radius:50%;
  background:radial-gradient(circle,rgba(16,185,129,.15) 0%,transparent 70%);pointer-events:none}
.hero-g2{position:absolute;bottom:-60px;left:15%;width:280px;height:280px;border-radius:50%;
  background:radial-gradient(circle,rgba(99,102,241,.1) 0%,transparent 70%);pointer-events:none}
.hero-g3{position:absolute;top:40%;right:25%;width:160px;height:160px;border-radius:50%;
  background:radial-gradient(circle,rgba(245,158,11,.06) 0%,transparent 70%);pointer-events:none}
.hero-tag{font-size:10px;letter-spacing:.14em;color:#10B981;text-transform:uppercase;font-weight:600;
  margin-bottom:10px;display:flex;align-items:center;gap:8px}
.hero-line{width:22px;height:1px;background:#10B981;opacity:.5}
.hero-title{font-size:32px;font-weight:800;color:#fff;letter-spacing:-.03em;line-height:1.1;margin-bottom:8px}
.hero-title .ht{color:#10B981}.hero-title .hg{color:#F59E0B}
.hero-sub{font-size:13px;color:rgba(255,255,255,.5);line-height:1.7;margin-bottom:20px;max-width:440px}
.hero-stats{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}
.hs{background:rgba(255,255,255,.06);border:1px solid rgba(255,255,255,.08);border-radius:11px;
  padding:12px;text-align:center}
.hs-v{font-size:22px;font-weight:800;color:#fff;line-height:1}
.hs-l{font-size:9px;color:rgba(255,255,255,.38);text-transform:uppercase;letter-spacing:.06em;margin-top:3px}

/* ── MODULE GRID ── */
.mg{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin-bottom:18px}
.mc{background:#fff;border:1px solid #E5E7EB;border-radius:13px;padding:16px;text-align:center;
  cursor:pointer;transition:all .18s;box-shadow:0 1px 4px rgba(0,0,0,.03)}
.mc:hover{border-color:#10B981;transform:translateY(-3px);box-shadow:0 8px 22px rgba(16,185,129,.12)}
.mc-icon{width:44px;height:44px;border-radius:11px;display:flex;align-items:center;
  justify-content:center;margin:0 auto 9px;font-size:22px}
.mc-name{font-size:12px;font-weight:600;color:#111827}
.mc-desc{font-size:10px;color:#6B7280;margin-top:2px}

/* ── METRIC CARDS ── */
.mets{display:grid;grid-template-columns:repeat(4,1fr);gap:8px;margin-bottom:14px}
.met{background:#F9FAFB;border:1px solid #E5E7EB;border-radius:10px;padding:12px}
.met-l{font-size:9px;color:#6B7280;text-transform:uppercase;letter-spacing:.07em;margin-bottom:4px}
.met-v{font-size:21px;font-weight:700;color:#111827;line-height:1;font-family:'JetBrains Mono',monospace}
.met-s{font-size:9px;color:#9CA3AF;margin-top:2px}

/* ── BUTTONS ── */
.stButton>button{
  background:linear-gradient(135deg,#10B981,#059669)!important;color:#fff!important;
  border:none!important;border-radius:10px!important;font-weight:600!important;
  font-size:13px!important;width:100%!important;font-family:'Inter',sans-serif!important;
  box-shadow:0 4px 12px rgba(16,185,129,.25)!important;padding:.65rem 1.5rem!important;
}
.stButton>button:hover{opacity:.92!important;transform:translateY(-1px)!important;transition:all .15s!important}

/* ── INPUTS ── */
.stTextArea textarea,.stTextInput input{
  background:#F9FAFB!important;border:1px solid #D1D5DB!important;border-radius:10px!important;
  font-family:'JetBrains Mono',monospace!important;font-size:12px!important;color:#111827!important;
}
.stTextArea textarea:focus,.stTextInput input:focus{border-color:#10B981!important;box-shadow:0 0 0 3px rgba(16,185,129,.1)!important}
.stSelectbox>div>div{background:#F9FAFB!important;border:1px solid #D1D5DB!important;border-radius:10px!important;color:#111827!important}
.stSlider>div>div>div{color:#10B981!important}

/* ── SEQUENCE DISPLAY ── */
.seq-block{font-family:'JetBrains Mono',monospace;font-size:11px;line-height:2.1;word-break:break-all;
  background:#F8FAFF;border:1px solid #E5E7EB;border-radius:10px;padding:12px;max-height:160px;overflow-y:auto}
.sg{color:#059669;font-weight:600}.sa{color:#6366F1}.sm{background:#FEF3C7;color:#92400E;border-radius:3px;padding:0 2px}

/* ── ALIGNMENT ── */
.aln-box{font-family:'JetBrains Mono',monospace;font-size:10.5px;line-height:1.9;background:#F8FAFF;
  border:1px solid #E5E7EB;border-radius:10px;padding:12px 14px;overflow-x:auto;white-space:pre}
.am{background:#ECFDF5;color:#065F46}.ami{background:#FEF2F2;color:#991B1B}.ag{color:#CBD5E1}

/* ── DISEASE ── */
.drow{display:flex;align-items:flex-start;padding:10px 0;border-bottom:1px solid #F3F4F6}
.drow:last-child{border-bottom:none}
.dname{font-size:12px;font-weight:600;color:#111827}
.dcat{font-size:10px;color:#6B7280;margin-top:1px}
.dds{font-size:9px;color:#6B7280;background:#F9FAFB;border:1px solid #E5E7EB;
  border-radius:5px;padding:3px 8px;margin-top:4px;line-height:1.6}
.dright{display:flex;align-items:center;gap:8px;margin-left:auto;padding-left:12px;flex-shrink:0}
.mb{width:58px;height:5px;background:#E5E7EB;border-radius:99px;overflow:hidden}
.mf{height:100%;border-radius:99px}
.rb{font-size:9px;padding:2px 8px;border-radius:99px;font-weight:700}
.rH{background:#FEE2E2;color:#991B1B}.rM{background:#FEF3C7;color:#92400E}
.rL{background:#ECFDF5;color:#065F46}.rN{background:#F0FDF4;color:#14532D}

/* ── ALERTS ── */
.alert{border-radius:10px;padding:11px 13px;margin-bottom:8px;border:1px solid;display:flex;gap:9px}
.alert-i{font-size:16px;flex-shrink:0}
.alert-t{font-size:12px;font-weight:600}
.alert-s{font-size:10px;color:#6B7280;margin-top:1px}
.a-H{background:#FEF2F2;border-color:rgba(153,27,27,.2)}.a-H .alert-i,.a-H .alert-t{color:#991B1B}
.a-M{background:#FFFBEB;border-color:rgba(146,64,14,.2)}.a-M .alert-i,.a-M .alert-t{color:#92400E}
.a-L{background:#ECFDF5;border-color:rgba(6,95,70,.2)}.a-L .alert-i,.a-L .alert-t{color:#065F46}
.a-I{background:#EFF6FF;border-color:rgba(29,78,216,.2)}.a-I .alert-i,.a-I .alert-t{color:#1D4ED8}

/* ── CONTACT ── */
.chero{background:linear-gradient(145deg,#0F172A,#1E293B);border-radius:16px;padding:36px 28px;
  text-align:center;position:relative;overflow:hidden;margin-bottom:18px}
.chero-g{position:absolute;border-radius:50%;pointer-events:none;top:-60px;left:50%;transform:translateX(-50%);
  width:380px;height:380px;background:radial-gradient(circle,rgba(16,185,129,.13) 0%,transparent 60%)}
.cpw{width:120px;height:120px;border-radius:50%;margin:0 auto 16px;padding:3px;
  background:linear-gradient(135deg,#10B981,#6366F1);position:relative;z-index:2}
.cph{width:100%;height:100%;border-radius:50%;object-fit:cover;object-position:top center;
  border:3px solid #0F172A}
.cname{font-size:24px;font-weight:800;color:#fff;margin-bottom:4px;position:relative;z-index:2}
.crole{font-size:13px;color:rgba(255,255,255,.5);margin-bottom:10px;position:relative;z-index:2}
.cbadge{display:inline-flex;align-items:center;gap:6px;background:rgba(16,185,129,.12);
  border:1px solid rgba(16,185,129,.25);color:#34D399;font-size:11px;font-weight:600;
  padding:5px 14px;border-radius:99px;position:relative;z-index:2}
.cchips{display:flex;gap:7px;justify-content:center;flex-wrap:wrap;margin-top:14px;position:relative;z-index:2}
.cc{font-size:10px;padding:3px 11px;border-radius:99px;font-weight:600;border:1px solid}
.cc-g{background:rgba(16,185,129,.12);color:#34D399;border-color:rgba(16,185,129,.25)}
.cc-i{background:rgba(99,102,241,.12);color:#A5B4FC;border-color:rgba(99,102,241,.25)}
.cc-a{background:rgba(245,158,11,.12);color:#FCD34D;border-color:rgba(245,158,11,.25)}
.cc-r{background:rgba(239,68,68,.12);color:#FCA5A5;border-color:rgba(239,68,68,.25)}

/* ── INFO TABLE ── */
.info-card{background:#fff;border:1px solid #E5E7EB;border-radius:14px;padding:20px;
  box-shadow:0 1px 4px rgba(0,0,0,.03)}
.ii{display:flex;align-items:center;gap:12px;padding:10px 0;border-bottom:1px solid #F3F4F6}
.ii:last-child{border-bottom:none}
.ii-icon{width:36px;height:36px;border-radius:9px;display:flex;align-items:center;
  justify-content:center;flex-shrink:0;font-size:16px}
.ii-l{font-size:10px;color:#9CA3AF;text-transform:uppercase;letter-spacing:.06em}
.ii-v{font-size:12px;font-weight:500;color:#111827;margin-top:1px}
.ii-v a{color:#10B981;text-decoration:none}.ii-v a:hover{text-decoration:underline}

/* ── SKILLS ── */
.skill-g{display:grid;grid-template-columns:repeat(3,1fr);gap:8px}
.sk{background:#F9FAFB;border:1px solid #E5E7EB;border-radius:9px;padding:9px 11px;
  font-size:11px;color:#111827;display:flex;align-items:center;gap:7px}

/* ── FEEDBACK ── */
.fb-panel{background:#111827;border-radius:14px;padding:20px;margin-top:16px}
.fb-title{font-size:14px;font-weight:700;color:#fff;margin-bottom:3px}
.fb-sub{font-size:10px;color:rgba(255,255,255,.4);margin-bottom:14px}
.fb-item{background:rgba(255,255,255,.06);border:1px solid rgba(255,255,255,.1);
  border-radius:10px;padding:12px 14px;margin-bottom:8px}
.fb-meta{font-size:9px;color:rgba(255,255,255,.3);margin-bottom:4px}
.fb-text{font-size:12px;color:rgba(255,255,255,.85)}
.fb-rating{font-size:14px;margin-top:4px}

/* ── COMP BARS ── */
.cbar{display:flex;align-items:center;gap:8px;margin-bottom:7px}
.cbar-b{font-size:10px;font-family:'JetBrains Mono',monospace;font-weight:700;width:14px}
.cbar-t{flex:1;height:7px;background:#F3F4F6;border-radius:99px;overflow:hidden}
.cbar-f{height:100%;border-radius:99px}
.cbar-p{font-size:10px;color:#6B7280;font-family:'JetBrains Mono',monospace;width:38px;text-align:right}

/* ── TAGS ── */
.tag{display:inline-block;font-size:10px;padding:2px 9px;border-radius:99px;margin:2px;font-weight:500;border:1px solid}
.tg{background:#ECFDF5;color:#065F46;border-color:rgba(6,95,70,.3)}
.tp{background:#EEF2FF;color:#3730A3;border-color:rgba(55,48,163,.3)}
.tb{background:#EFF6FF;color:#1D4ED8;border-color:rgba(29,78,216,.3)}
.ta{background:#FFFBEB;color:#92400E;border-color:rgba(146,64,14,.3)}
.tr{background:#FEF2F2;color:#991B1B;border-color:rgba(153,27,27,.3)}

/* ── PROTEIN PRED ── */
.pcard{background:#F9FAFB;border:1px solid #E5E7EB;border-radius:10px;padding:10px 12px;
  margin-bottom:7px;display:flex;align-items:center;gap:10px}
.pcard:hover{transform:translateX(2px);transition:transform .15s}
.pconf{height:3px;background:#E5E7EB;border-radius:99px;overflow:hidden;margin-top:5px}
.pconf-f{height:100%;border-radius:99px}

/* ── BLAST ── */
.blast-hit{background:#F9FAFB;border:1px solid #E5E7EB;border-radius:10px;padding:12px;margin-bottom:8px}
.blast-hit:hover{border-color:#10B981;box-shadow:0 2px 10px rgba(16,185,129,.1)}
.blast-title{font-size:12px;font-weight:600;color:#111827;margin-bottom:4px}
.blast-meta{font-size:10px;color:#6B7280;margin-bottom:6px}
.blast-bar{height:5px;background:#E5E7EB;border-radius:99px;overflow:hidden;margin-bottom:4px}
.blast-fill{height:100%;border-radius:99px;background:#10B981}
.blast-stat{display:flex;gap:10px;font-size:10px;color:#6B7280;font-family:'JetBrains Mono',monospace}

/* ── VIZ ── */
.viz-card{background:#fff;border:1px solid #E5E7EB;border-radius:14px;padding:16px;margin-bottom:12px}
.score-ring{width:64px;height:64px;border-radius:50%;display:flex;flex-direction:column;
  align-items:center;justify-content:center;flex-shrink:0;border:2.5px solid}
.sc-v{font-size:15px;font-weight:800;line-height:1;font-family:'JetBrains Mono',monospace}
.sc-l{font-size:7px;margin-top:2px;text-transform:uppercase;letter-spacing:.05em}

/* ── STREAMLIT OVERRIDES ── */
div[data-testid="metric-container"]{background:#fff!important;border:1px solid #E5E7EB!important;border-radius:12px!important;padding:12px!important}
div[data-testid="metric-container"] label{color:#6B7280!important;font-size:11px!important;text-transform:uppercase;letter-spacing:.06em;font-weight:600!important}
div[data-testid="metric-container"] div[data-testid="stMetricValue"]{color:#10B981!important;font-family:'JetBrains Mono',monospace!important;font-size:22px!important}
.stProgress .st-bo{background:#10B981!important}
.stProgress .st-bp{background:#ECFDF5!important}
.stDataFrame{border:1px solid #E5E7EB!important;border-radius:12px!important}
.stExpander{background:#fff!important;border:1px solid #E5E7EB!important;border-radius:12px!important}
h1,h2,h3{color:#111827!important}
.stSuccess{background:#ECFDF5!important;color:#065F46!important;border:1px solid rgba(6,95,70,.2)!important;border-radius:10px!important}
.stWarning{background:#FFFBEB!important;border:1px solid rgba(146,64,14,.2)!important;border-radius:10px!important}
.stError{background:#FEF2F2!important;border:1px solid rgba(153,27,27,.2)!important;border-radius:10px!important}
.stInfo{background:#EFF6FF!important;border:1px solid rgba(29,78,216,.2)!important;border-radius:10px!important}

/* ── TICKER ── */
.ticker-wrap{overflow:hidden;background:#ECFDF5;border-top:1px solid #D1FAE5;border-bottom:1px solid #D1FAE5;padding:7px 0;margin-bottom:16px}
.ticker{display:flex;white-space:nowrap;animation:tick 38s linear infinite;font-size:10px;font-family:'JetBrains Mono',monospace;color:#059669;letter-spacing:.04em}
.tsep{color:#10B981;margin:0 14px;opacity:.5}
@keyframes tick{0%{transform:translateX(0)}100%{transform:translateX(-50%)}}
@keyframes blink{0%,100%{opacity:1}50%{opacity:.2}}
.dot{display:inline-block;width:7px;height:7px;border-radius:50%;background:#10B981;animation:blink 2s infinite}
</style>"""

st.markdown(CSS, unsafe_allow_html=True)
print("Part 1 loaded")
