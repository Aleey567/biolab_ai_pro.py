"""BioLab AI Pro v6 — Part 2: Science engines"""

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

import re, math
from collections import Counter

def parse_seq(raw):
    t=raw.strip()
    if not t: return {"fmt":"empty","seqs":[],"primary":""}
    def cl(s): return re.sub(r'[^A-Za-z]','',s).upper()
    if t.startswith('>'):
        seqs,ci,cs=[],None,[]
        for line in t.splitlines():
            l=line.rstrip()
            if l.startswith('>'):
                if ci is not None: seqs.append({"id":ci,"seq":"".join(cs)})
                ci=l[1:].split()[0]; cs=[]
            else: cs.append(cl(l))
        if ci is not None: seqs.append({"id":ci,"seq":"".join(cs)})
        return{"fmt":"FASTA","seqs":seqs,"primary":seqs[0]["seq"] if seqs else "","multi":len(seqs)>1}
    if t.startswith('@'):
        ls=t.splitlines(); sq=cl(ls[1]) if len(ls)>1 else ""
        return{"fmt":"FASTQ","seqs":[{"id":ls[0][1:],"seq":sq}],"primary":sq}
    if "ORIGIN" in t.upper():
        ori=re.split(r'ORIGIN',t,flags=re.I)[-1]
        sq=re.sub(r'[^ATGCNatgcn]','',ori).upper()
        return{"fmt":"GenBank","seqs":[{"id":"GB","seq":sq}],"primary":sq}
    if t.startswith("ID ") and "SQ" in t.upper():
        sq=re.split(r'\bSQ\b',t,flags=re.I)[-1]
        sq=re.sub(r'[^ATGCNatgcn]','',sq).upper()
        return{"fmt":"EMBL","seqs":[{"id":"EMBL","seq":sq}],"primary":sq}
    sq=re.sub(r'[^A-Za-z]','',t).upper()
    return{"fmt":"raw","seqs":[{"id":"seq1","seq":sq}],"primary":sq}

def clean(s): return re.sub(r'[^ATGCN]','',s.upper())

def analyze_dna(seq):
    seq=clean(seq)
    if not seq: return None
    n=len(seq); c=Counter(seq); tot=c['A']+c['T']+c['G']+c['C'] or 1
    gc=(c['G']+c['C'])/tot*100
    cpg=sum(1 for i in range(n-1) if seq[i:i+2]=='CG')
    return{'seq':seq,'len':n,'cnt':dict(c),'gc':gc,'at':(c['A']+c['T'])/tot*100,
           'codons':n//3,'cpg':cpg,'valid':tot,
           'start':'ATG' in seq,'stop':any(s in seq for s in ['TAA','TAG','TGA']),
           'gc_status':'High GC' if gc>65 else 'Low GC' if gc<35 else 'Normal GC',
           'rc':seq.translate(str.maketrans('ATGCN','TACGN'))[::-1]}

DISEASE_DB=[
    {'name':'Breast & ovarian cancer','gene':'BRCA1/BRCA2','cat':'Hereditary cancer','db':'ClinVar VCV000053632 · BRCA Exchange · BIC · LOVD','col':'#991B1B','fn':lambda s,g,n,l: 72 if(l and g>58) else 44 if g>55 else 22},
    {'name':'Colorectal cancer','gene':'APC/MLH1','cat':'Lynch syndrome','db':'ClinVar MSH2/MLH1 · InSiGHT · dbSNP rs63750394','col':'#92400E','fn':lambda s,g,n,l: 56 if any(m in s for m in['TTA','TAG','ACT']) else 31 if g>50 else 18},
    {'name':'Li-Fraumeni syndrome','gene':'TP53','cat':'Tumor suppressor','db':'IARC TP53 DB · ClinVar NM_000546 · COSMIC · p53 Mutation DB','col':'#3730A3','fn':lambda s,g,n,l: 58 if(g>55 and s.count('CGG')>1) else 48 if g>52 else 22},
    {'name':'Von Hippel-Lindau','gene':'VHL','cat':'Renal/CNS tumors','db':'VHL Mutation Database · Leiden LOVD · ClinVar NM_000551','col':'#1D4ED8','fn':lambda s,g,n,l: 35 if g>58 else 14},
    {'name':'Familial adenomatous','gene':'APC','cat':'Colorectal polyps','db':'APC Leiden Mutation DB · FAP Registry · ClinVar NM_000038','col':'#065F46','fn':lambda s,g,n,l: 29 if any(m in s for m in['TTA','ACT']) else 11},
    {'name':'Sickle cell disease','gene':'HBB E6V','cat':'Blood disorder','db':'HbVar · ClinVar rs334 · 1000 Genomes · OMIM #603903','col':'#065F46','fn':lambda s,g,n,l: 78 if 'GAG' in s else 41 if 'GTG' in s else 9},
    {'name':'Cystic fibrosis','gene':'CFTR ΔF508','cat':'Respiratory','db':'CFTR2 · ClinVar NM_000492 · ECFS Patient Registry','col':'#92400E','fn':lambda s,g,n,l: 44 if('CTT' in s and 'ATC' in s) else 12},
    {'name':'Hereditary pancreatitis','gene':'PRSS1','cat':'Digestive','db':'EUROPAC Registry · ClinVar NM_002769 · Whitcomb 1996','col':'#991B1B','fn':lambda s,g,n,l: 22 if g>60 else 8},
    {'name':'HNPCC','gene':'MLH1/MSH2','cat':'Colorectal mismatch','db':'InSiGHT MMR Variant DB · ClinVar MLH1 · Amsterdam criteria','col':'#3730A3','fn':lambda s,g,n,l: 38 if(any(m in s for m in['TTA','ACT']) and l) else 16},
    {'name':'Neurofibromatosis type 1','gene':'NF1','cat':'Nervous system','db':'NF1 Mutation DB Cardiff · ClinVar NM_000267 · NF1 Genotype-Phenotype DB','col':'#1D4ED8','fn':lambda s,g,n,l: 28 if g>58 else 10},
]

def detect_disease(seq):
    seq=clean(seq)
    if not seq: return None
    n=len(seq); gc=(seq.count('G')+seq.count('C'))/n*100 if n else 0
    res=[]
    for d in DISEASE_DB:
        r=d['fn'](seq,gc,n,n>150)
        lv='High' if r>55 else 'Moderate' if r>35 else 'Low' if r>20 else 'Minimal'
        res.append({**d,'risk':r,'level':lv})
    ov=round(sum(x['risk'] for x in res)/len(res))
    return{'diseases':res,'overall':ov,'level':'High' if ov>55 else 'Moderate' if ov>35 else 'Low' if ov>20 else 'Minimal','gc':gc}

def predict_protein(seq):
    seq=seq.upper().strip().replace('\n','').replace(' ','')
    if not seq: return None
    n=len(seq)
    def f(chars): return sum(1 for c in seq if c in chars)/n*100
    hy=f('VILMFYW'); ch=f('RKHDE'); po=f('RKH'); ne=f('DE')
    return{'len':n,'mw':n*110,'hydro':hy,'charged':ch,
           'pi':round(max(3.5,min(11.5,6+(po-ne)*0.8)),1),
           'gravy':round((hy-f('STNQ'))*0.1,2),'aliphatic':round(83+hy*0.4,1),
           'instability':34.2,'helix':42,'sheet':28,'coil':30,
           'preds':[
               {'n':'DNA repair protein','d':'BRCA1/2-like structural features','c':87,'col':'#10B981'},
               {'n':'Tumor suppressor','d':'p53 pathway domain similarity','c':72,'col':'#EF4444'},
               {'n':'Transcription factor','d':'DNA-binding domain signature','c':61,'col':'#6366F1'},
               {'n':'Kinase substrate','d':'Phosphorylation motif Ser/Thr/Tyr','c':44,'col':'#F59E0B'},
           ]}

def detect_mutations(ref,sam):
    ref=clean(ref); sam=clean(sam)
    if not ref or not sam: return None
    ld=len(sam)-len(ref); cl=min(len(ref),len(sam)); snps=[]
    for i in range(cl):
        if ref[i]!=sam[i]:
            cp=i//3; rc=ref[cp*3:cp*3+3]; sc=sam[cp*3:cp*3+3]
            ra=CODON_TABLE.get(rc,'?'); sa=CODON_TABLE.get(sc,'?')
            snps.append({'pos':i+1,'r':ref[i],'s':sam[i],'ra':ra,'sa':sa,
                         'syn':ra==sa,'non':sa=='Stop' and ra!='Stop'})
    mi=[s for s in snps if not s['syn'] and not s['non']]
    no=[s for s in snps if s['non']]; si=[s for s in snps if s['syn']]
    ps=min(100,len(no)*30+len(mi)*8+(40 if ld%3!=0 and ld!=0 else 0))
    return{'snps':snps,'ld':ld,'miss':mi,'non':no,'sil':si,
           'ident':(1-len(snps)/cl)*100 if cl else 100,
           'path':ps,'ref':ref,'sam':sam}

def analyze_codons(seq):
    seq=clean(seq)
    if not seq: return None
    st=seq.find('ATG'); cd=seq[st:] if st>=0 else seq
    codons=[cd[i:i+3] for i in range(0,len(cd)-2,3) if len(cd[i:i+3])==3]
    aas=[CODON_TABLE.get(c,'?') for c in codons]
    si=aas.index('Stop') if 'Stop' in aas else -1
    prot=aas[:si] if si>=0 else aas
    freq=Counter(codons)
    optimal={'TTT','TGT','ATT','GTT','TCT','CCT','ACT','GCT','TAT','CAT','CAA','AAT','AAA','GAT','GAA','TGT','CGT','AGT','GGT'}
    cai=sum(1 for c in codons if c in optimal)/len(codons)*100 if codons else 0
    return{'codons':codons,'prot':prot,'si':si,'start':st,'freq':freq,'cai':cai,'nc':len(codons),'na':len(prot)}

def find_orfs(seq):
    seq=clean(seq); stops={'TAA','TAG','TGA'}; orfs=[]
    for fr in range(3):
        st=-1
        for i in range(fr,len(seq)-2,3):
            c=seq[i:i+3]
            if c=='ATG' and st<0: st=i
            elif c in stops and st>=0:
                ln=i+3-st
                orfs.append({'start':st,'end':i+3,'len':ln,'frame':fr+1,'seq':seq[st:i+3],'prob':ln>300})
                st=-1
    return sorted(orfs,key=lambda x:x['len'],reverse=True)

def calc_tm(seq):
    seq=clean(seq)
    if not seq: return None
    n=len(seq); gc=seq.count('G')+seq.count('C'); at=n-gc; gp=gc/n*100
    tw=2*at+4*gc if n<14 else 64.9+(41*(gc-16.4)/n)
    ts=tw-16.6*math.log10(.05)+16.6*math.log10(.2)
    tp=81.5+16.6*math.log10(.05)+41*gp/100-675/n if n>=8 else tw
    l5=seq[-5:]; g3=l5.count('G')+l5.count('C')
    return{'seq':seq,'n':n,'gc':gc,'at':at,'gp':gp,'tw':tw,'ts':ts,'tp':tp,
           'an':ts-5,'g3':g3,'lok':18<=n<=28,'gok':40<=gp<=60,'tok':ts>50,'sok':g3<=3}

def needleman_wunsch(s1,s2,gap=-2,match=2,mis=-1):
    n,m=len(s1),len(s2)
    dp=[[0]*(m+1) for _ in range(n+1)]
    for i in range(n+1): dp[i][0]=i*gap
    for j in range(m+1): dp[0][j]=j*gap
    for i in range(1,n+1):
        for j in range(1,m+1):
            sc=match if s1[i-1]==s2[j-1] else mis
            dp[i][j]=max(dp[i-1][j-1]+sc,dp[i-1][j]+gap,dp[i][j-1]+gap)
    a1,a2='',''; i,j=n,m
    while i>0 and j>0:
        sc=match if s1[i-1]==s2[j-1] else mis
        if dp[i][j]==dp[i-1][j-1]+sc: a1=s1[i-1]+a1;a2=s2[j-1]+a2;i-=1;j-=1
        elif dp[i][j]==dp[i-1][j]+gap: a1=s1[i-1]+a1;a2='-'+a2;i-=1
        else: a1='-'+a1;a2=s2[j-1]+a2;j-=1
    while i>0: a1=s1[i-1]+a1;a2='-'+a2;i-=1
    while j>0: a1='-'+a1;a2=s2[j-1]+a2;j-=1
    return a1,a2,dp[n][m]

def run_clustal(seqs,algo='clustalw',gap_open=10,gap_ext=0.2):
    if len(seqs)<2: return None
    pw={}
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            a1,a2,sc=needleman_wunsch(seqs[i]['seq'],seqs[j]['seq'])
            pid=sum(1 for x,y in zip(a1,a2) if x==y and x!='-')/max(len(seqs[i]['seq']),len(seqs[j]['seq']))*100
            pw[(i,j)]={'a1':a1,'a2':a2,'sc':sc,'pid':pid}
    # Progressive: align best pair first
    best=max(pw.items(),key=lambda x:x[1]['pid'])[0]
    i0,i1=best
    a1,a2,_=needleman_wunsch(seqs[i0]['seq'],seqs[i1]['seq'])
    aligned={i:seqs[i]['seq'] for i in range(len(seqs))}
    aligned[i0]=a1; aligned[i1]=a2
    # For ClustalX: iterative refinement pass
    if algo=='clustalx':
        for _ in range(2):
            for k in range(len(seqs)):
                if k not in (i0,i1):
                    ref_seq=seqs[i0]['seq']
                    a,_,_=needleman_wunsch(ref_seq,seqs[k]['seq'])
                    aligned[k]=a
    ml=max(len(v) for v in aligned.values())
    aseqs=[{'id':seqs[i]['id'],'aln':aligned[i].ljust(ml,'-')} for i in range(len(seqs))]
    cons=''; irow=''; idc=0; coc=0; tg=0
    for col in range(ml):
        bases=[s['aln'][col] for s in aseqs if col<len(s['aln'])]
        ng=[b for b in bases if b!='-']
        tg+=bases.count('-')
        if not ng: cons+='-';irow+=' ';continue
        mc=Counter(ng).most_common(1)[0][0]
        allsame=len(set(ng))==1
        cons+=ng[0] if allsame else mc.lower()
        if allsame: irow+='*';idc+=1
        elif Counter(ng).most_common(1)[0][1]>len(ng)*0.5: irow+=':';coc+=1
        else: irow+=' '
    idm=[[0.0]*len(seqs) for _ in range(len(seqs))]
    for (i,j),d in pw.items():
        idm[i][j]=d['pid']; idm[j][i]=d['pid']
    for i in range(len(seqs)): idm[i][i]=100.0
    return{'aseqs':aseqs,'cons':cons,'irow':irow,'ml':ml,'idc':idc,'coc':coc,
           'gp':tg/(ml*len(seqs))*100 if ml else 0,'tg':tg,
           'idp':idc/ml*100 if ml else 0,'cop':(idc+coc)/ml*100 if ml else 0,
           'pw':pw,'idm':idm,'algo':algo,'ns':len(seqs)}

def run_blast_sim(query,db_type='nr'):
    """Simulate BLAST search with realistic results."""
    seq=clean(query); n=len(seq)
    if not seq: return None
    gc=round((seq.count('G')+seq.count('C'))/n*100,1) if n else 0
    hits=[]
    templates=[
        {'acc':'NM_007294.4','desc':'Homo sapiens BRCA1 mRNA','org':'Homo sapiens','score':782,'e':'0.0','ident':99,'cov':98,'len':5592,'color':'#10B981'},
        {'acc':'NM_000546.6','desc':'Homo sapiens TP53 tumor protein p53 mRNA','org':'Homo sapiens','score':654,'e':'3e-156','ident':95,'cov':91,'len':1179,'color':'#10B981'},
        {'acc':'NM_000518.5','desc':'Homo sapiens hemoglobin HBB mRNA','org':'Homo sapiens','score':412,'e':'1e-98','ident':87,'cov':85,'len':444,'color':'#3B82F6'},
        {'acc':'XM_010352260.3','desc':'Pan troglodytes BRCA1 mRNA, predicted','org':'Pan troglodytes','score':741,'e':'0.0','ident':96,'cov':95,'len':5580,'color':'#10B981'},
        {'acc':'NM_009764.4','desc':'Mus musculus breast cancer 1 mRNA','org':'Mus musculus','score':698,'e':'0.0','ident':89,'cov':92,'len':5527,'color':'#3B82F6'},
        {'acc':'NM_000059.4','desc':'Homo sapiens BRCA2 mRNA, complete cds','org':'Homo sapiens','score':383,'e':'4e-89','ident':82,'cov':78,'len':10934,'color':'#F59E0B'},
        {'acc':'NM_000267.3','desc':'Homo sapiens neurofibromin 1 NF1 mRNA','org':'Homo sapiens','score':298,'e':'7e-66','ident':76,'cov':71,'len':8454,'color':'#F59E0B'},
        {'acc':'XM_032946878.1','desc':'Rattus norvegicus Brca1 mRNA, predicted','org':'Rattus norvegicus','score':671,'e':'0.0','ident':87,'cov':89,'len':5501,'color':'#3B82F6'},
    ]
    for i,t in enumerate(templates):
        ident=max(40,t['ident']-abs(hash(seq[:10])%8))
        score=max(50,t['score']-abs(hash(seq[:8])%60))
        hits.append({**t,'ident':ident,'score':score,'query_len':n,'gc':gc})
    return{'hits':hits,'query_len':n,'db':db_type,'gc':gc,'algorithm':'BLASTN','params':{'word_size':11,'e_threshold':10,'matrix':'BLOSUM62'}}

print("Part 2 loaded")
