"""BioLab AI Pro v6 — Part 3: UI helpers + all module pages"""
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd, re, math
from collections import Counter

def ticker():
    segs=["5' ATGAAAGCAATTTTCGTACTGAAA","GC: 48.3% · AT: 51.7%",
          "ClustalW · ClustalX · BLAST · Progressive Alignment",
          "BRCA1 · TP53 · HBB · CFTR · NF1","FASTA · FASTQ · GenBank · EMBL · Raw DNA",
          "Agriculture University Faisalabad · Ali Raza · razabaig567@gmail.com",
          "8 Modules · 5 Formats · 10 Diseases · Real-time Alignment · BLAST Search"]
    inner="".join([f'<span>{s}</span><span class="tsep">|</span>' for s in segs]*2)
    st.markdown(f'<div class="ticker-wrap"><div class="ticker">{inner}</div></div>',unsafe_allow_html=True)

def topnav(page="Dashboard"):
    st.markdown(f"""<div class="topnav">
    <div class="topnav-brand">
      <div class="topnav-logo">🧬</div>
      <span class="topnav-title">BioLab AI Pro<span class="topnav-ver">v6</span></span>
    </div>
    <span class="topnav-crumb">→ {page}</span>
    </div>""",unsafe_allow_html=True)

def met_html(items):
    inner="".join([f'<div class="met"><div class="met-l">{l}</div><div class="met-v" style="color:{col}">{v}</div><div class="met-s">{s}</div></div>' for l,v,s,col in items])
    return f'<div class="mets">{inner}</div>'

def cbar(base,pct,col):
    return f'<div class="cbar"><span class="cbar-b" style="color:{col}">{base}</span><div class="cbar-t"><div class="cbar-f" style="width:{min(pct,100):.1f}%;background:{col}"></div></div><span class="cbar-p">{pct:.1f}%</span></div>'

def seq_html(seq,mx=160):
    out=""
    for c in seq[:mx]:
        if c in "GC": out+=f'<span class="sg">{c}</span>'
        else: out+=f'<span class="sa">{c}</span>'
    if len(seq)>mx: out+=f'<span style="color:#9CA3AF"> +{len(seq)-mx:,} more</span>'
    return out

def aln_html(aln,cons,mx=68):
    out=""
    for i,c in enumerate(aln[:mx]):
        ci=cons[i] if i<len(cons) else ' '
        if c=='-': out+=f'<span class="ag">-</span>'
        elif ci=='*': out+=f'<span class="am">{c}</span>'
        elif ci==':': out+=f'<span style="color:#D97706;background:#FFFBEB">{c}</span>'
        else: out+=f'<span class="ami">{c}</span>'
    if len(aln)>mx: out+=f'<span style="color:#9CA3AF">  +{len(aln)-mx}</span>'
    return out

# ── 3D HELIX (Streamlit-compatible, zero WebGL) ───────────────────
HELIX_HTML="""
<canvas id="hx" width="780" height="180"
  style="width:100%;height:180px;display:block;border-radius:14px;background:transparent"></canvas>
<script>
(function(){
var cv=document.getElementById('hx'),ctx=cv.getContext('2d'),t=0;
var C=['#10B981','#6366F1','#F59E0B','#EF4444'];
var B=['A','T','G','C'];
function draw(){
  ctx.clearRect(0,0,cv.width,cv.height);
  var N=24,cy=cv.height/2;
  var p1=[],p2=[];
  for(var i=0;i<N;i++){
    var x=24+i*(cv.width-48)/(N-1);
    p1.push({x:x,y:cy+Math.sin(i*.38+t)*52});
    p2.push({x:x,y:cy-Math.sin(i*.38+t)*52});
  }
  // Strand 1
  ctx.beginPath();ctx.moveTo(p1[0].x,p1[0].y);
  for(var i=1;i<N;i++){var a=p1[i-1],b=p1[i];ctx.bezierCurveTo((a.x+b.x)/2,a.y,(a.x+b.x)/2,b.y,b.x,b.y);}
  ctx.strokeStyle='rgba(16,185,129,.8)';ctx.lineWidth=2.5;ctx.stroke();
  // Strand 2
  ctx.beginPath();ctx.moveTo(p2[0].x,p2[0].y);
  for(var i=1;i<N;i++){var a=p2[i-1],b=p2[i];ctx.bezierCurveTo((a.x+b.x)/2,a.y,(a.x+b.x)/2,b.y,b.x,b.y);}
  ctx.strokeStyle='rgba(99,102,241,.8)';ctx.lineWidth=2.5;ctx.stroke();
  // Rungs + dots
  for(var i=0;i<N;i++){
    if(i%2===0){ctx.beginPath();ctx.moveTo(p1[i].x,p1[i].y);ctx.lineTo(p2[i].x,p2[i].y);ctx.strokeStyle='rgba(255,255,255,.14)';ctx.lineWidth=1;ctx.stroke();}
    var cl=C[i%4];
    ctx.beginPath();ctx.arc(p1[i].x,p1[i].y,5,0,6.28);ctx.fillStyle=cl;ctx.globalAlpha=.9;ctx.fill();ctx.globalAlpha=1;
    ctx.beginPath();ctx.arc(p2[i].x,p2[i].y,5,0,6.28);ctx.fillStyle=C[(i+2)%4];ctx.globalAlpha=.9;ctx.fill();ctx.globalAlpha=1;
    if(i%4===0){
      ctx.font='9px JetBrains Mono,monospace';ctx.textAlign='center';
      ctx.fillStyle=cl;ctx.globalAlpha=.65;
      ctx.fillText(B[i%4],p1[i].x,p1[i].y+14);
      ctx.fillText(B[(i+2)%4],p2[i].x,p2[i].y-7);
      ctx.globalAlpha=1;
    }
  }
  t+=0.016;requestAnimationFrame(draw);
}
draw();
})();
</script>"""

def show_helix():
    components.html(HELIX_HTML,height=190)

def add_feedback(name,rating,text,module):
    st.session_state.feedbacks.append({
        "name":name.strip() or "Anonymous",
        "rating":rating,"text":text.strip(),
        "module":module,
        "time":__import__('datetime').datetime.now().strftime("%Y-%m-%d %H:%M")
    })

def feedback_form(location="Dashboard"):
    st.markdown("---")
    st.markdown("### 💬 Leave Feedback")
    # Use unique key based on location to avoid conflict
    with st.form(key=f"fb_form_{location}", clear_on_submit=True):
        c1,c2=st.columns([2,1])
        with c1: fb_name=st.text_input("Your name (optional)",placeholder="e.g. Ahmed Ali",key=f"fb_name_{location}")
        with c2: fb_rating=st.select_slider("Rating ⭐",options=[1,2,3,4,5],value=5,key=f"fb_rat_{location}")
        fb_text=st.text_area("Your feedback / suggestions",placeholder="What do you think about BioLab AI Pro?",height=80,key=f"fb_txt_{location}")
        fb_module=st.selectbox("Module you're rating",["Overall App","Dashboard","DNA Analysis","ClustalW/X Alignment","Disease Detection","Mutation Detection","Protein Prediction","Codon Analysis","ORF Finder","Tm/PCR Calculator","BLAST Search","Visualization Tools"],key=f"fb_mod_{location}")
        submitted=st.form_submit_button("📨 Submit Feedback",use_container_width=True)
        if submitted:
            if fb_text.strip():
                add_feedback(fb_name,fb_rating,fb_text,fb_module)
                st.success(f"✅ Thank you for your feedback! ({len(st.session_state.feedbacks)} total received)")
                st.balloons()
            else:
                st.warning("Please write some feedback before submitting.")

def admin_panel():
    if not st.session_state.show_admin: return
    st.markdown("---")
    st.markdown("## 🔐 Admin Feedback Panel")
    st.caption(f"Visible only when Admin View is active in sidebar · {len(st.session_state.feedbacks)} total responses")
    if not st.session_state.feedbacks:
        st.info("No feedback received yet. Feedback from users will appear here.")
        return
    cols=st.columns(4)
    total=len(st.session_state.feedbacks)
    avg_rating=sum(f['rating'] for f in st.session_state.feedbacks)/total if total else 0
    cols[0].metric("Total Responses",total)
    cols[1].metric("Avg Rating",f"{avg_rating:.1f}/5 ⭐")
    cols[2].metric("Latest",st.session_state.feedbacks[-1]['time'] if st.session_state.feedbacks else "—")
    cols[3].metric("Modules",len(set(f['module'] for f in st.session_state.feedbacks)))
    # Sort: newest first
    for fb in reversed(st.session_state.feedbacks):
        stars="⭐"*fb['rating']
        st.markdown(f"""<div class="fb-item">
        <div class="fb-meta">{fb['name']} · {fb['module']} · {fb['time']}</div>
        <div class="fb-rating">{stars}</div>
        <div class="fb-text">{fb['text']}</div>
        </div>""",unsafe_allow_html=True)
    col1,col2=st.columns(2)
    with col1:
        if st.button("🗑️ Clear all feedback",key="clear_fb"):
            st.session_state.feedbacks=[]
            st.success("All feedback cleared.")
            st.rerun()
    with col2:
        if st.button("📊 Export feedback",key="exp_fb"):
            import json
            st.download_button("💾 Download JSON",data=json.dumps(st.session_state.feedbacks,indent=2),
                               file_name="biolab_feedback.json",mime="application/json")

print("Part 3 loaded")
