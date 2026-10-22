# Figures for docs/guide/ploidyspec_guide.typ. Run from the repo root with D=<output dir>
# (figures go to $D/fig). Needs $D/panel from `ploidyspec panel --outdir $D/panel`,
# $D/old_empnigr_windowed.tsv (pre-fix windowed_all.tsv from git), and the
# SchCurv1 chr17/chr19 position-free tracks at $TMPDIR/wcarp (or results/SchCurv1 once re-run).
# The Ranunculus figure reads results/dmRanRepe1/homeologs/{subgenome_partition,homeolog_pairs}.tsv.
# The TE split figure reads results/{ddEmpNigr1,icStrMela3}/subgenomes/te_marker_fraction_by_lineage.tsv.
import csv, os, collections, statistics as st
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
D=os.environ['D']; F=D+'/fig'
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans'})
BLUE='#2b6cb0'; RED='#c53030'; GREY='#718096'; GOLD='#d69e2e'; GREEN='#2f855a'

def pairs(sp):
    d={}
    for r in csv.DictReader(open(f'results/{sp}/matrix/whole_chrom_pairs.tsv'),delimiter='\t'):
        d[(r['unit_a'],r['unit_b'])]=d[(r['unit_b'],r['unit_a'])]=float(r['distance'])
    return d

# Fig: three shapes of copy distance (4x4 matrices)
cases=[('drAriEdul1','chr08','Interchangeable copies\n(tetrasomic-like)'),
       ('ddEmpNigr1','chr03','Two lineages\n(disomic-like)'),
       ('SchCurv1','chr19','Fusion splits two lineages\n(HAP3, HAP4 fused to chr22)')]
fig,axs=plt.subplots(1,3,figsize=(9,3.1))
for ax,(sp,c,title) in zip(axs,cases):
    d=pairs(sp); ids=[f'HAP{i}_{c}' for i in range(1,5)]
    M=np.array([[0 if a==b else d[(a,b)]*1000 for b in ids] for a in ids])
    im=ax.imshow(M,cmap='viridis_r',vmin=0,vmax=max(M.max(),1))
    for i in range(4):
        for j in range(4):
            if i!=j: ax.text(j,i,f'{M[i,j]:.0f}',ha='center',va='center',color='white' if M[i,j]<M.max()*0.5 else 'black',fontsize=8)
    lab=[f'H{i}' for i in range(1,5)]
    ax.set_xticks(range(4),lab); ax.set_yticks(range(4),lab)
    ax.set_title(f'{title}\n{sp} {c}',fontsize=8.5)
    for s in ax.spines.values(): s.set_visible(False)
fig.text(0.5,0.01,'cell = whole-chromosome k-mer distance between two copies (x 1000)',ha='center',fontsize=8,color=GREY)
fig.tight_layout(rect=[0,0.05,1,0.97]); fig.savefig(F+'/three_shapes.png',dpi=200); plt.close(fig)

# Fig: windowed old vs new, ddEmpNigr1 chr03 closest pair HAP1-HAP3
old=collections.defaultdict(list)
for r in csv.DictReader(open(D+'/old_empnigr_windowed.tsv'),delimiter='\t'):
    if r['group']=='chr03' and {r['hap_a'],r['hap_b']}=={'HAP1','HAP3'}: old['x'].append(int(r['win_start'])/1e6); old['y'].append(float(r['jaccard_distance']))
new=collections.defaultdict(lambda: collections.defaultdict(list))
for r in csv.DictReader(open('results/ddEmpNigr1/windowed/windowed_all.tsv'),delimiter='\t'):
    if r['group']=='chr03' and r['unit_a']=='HAP1_chr03': new[r['unit_b']]['x'].append(int(r['win_start'])/1e6); new[r['unit_b']]['y'].append(float(r['distance']))
fig,axs=plt.subplots(1,2,figsize=(9,2.9))
axs[0].plot(old['x'],old['y'],color=GREY,lw=1); axs[0].set_ylim(0,1.05)
axs[0].set_title('Before: windows compared at equal coordinates\n(HAP1 vs HAP3, which are 0.004 apart overall)',fontsize=8.5)
axs[0].set_ylabel('Jaccard distance'); axs[0].set_xlabel('position (Mb)')
col={'HAP3_chr03':BLUE,'HAP2_chr03':RED,'HAP4_chr03':GOLD}
for u,v in sorted(new.items()):
    axs[1].plot(v['x'],v['y'],lw=1,color=col.get(u,GREY),label=f'HAP1 vs {u[:4]}')
axs[1].set_title('After: each HAP1 window looked up in each other copy\n(position-free)',fontsize=8.5)
axs[1].set_ylabel('distance (-ln c / k)'); axs[1].set_xlabel('position along HAP1 (Mb)'); axs[1].legend(fontsize=7,frameon=False)
fig.tight_layout(); fig.savefig(F+'/windowed_fix.png',dpi=200); plt.close(fig)

# Fig: snow carp chr19 / chr17 tracks (new method)
W='results/SchCurv1/windowed/windowed_all.tsv'
tr=collections.defaultdict(lambda: collections.defaultdict(list))
for r in csv.DictReader(open(W),delimiter='\t'):
    if r['unit_a'].startswith('HAP1'): tr[(r['group'],r['unit_b'][:4])]['x'].append(int(r['win_start'])/1e6); tr[(r['group'],r['unit_b'][:4])]['y'].append(float(r['distance']))
fig,axs=plt.subplots(1,2,figsize=(9,2.9),sharey=True)
for ax,c,t in [(axs[0],'chr19','chr19: fused in HAP3/HAP4 - split along the whole chromosome'),(axs[1],'chr17','chr17: split in the last ~8 Mb only')]:
    for h,colr in [('HAP2',BLUE),('HAP3',RED),('HAP4',GOLD)]:
        v=tr[(c,h)]
        if v['x']:
            y=np.convolve(v['y'],np.ones(5)/5,mode='same'); ax.plot(v['x'],y,color=colr,lw=1.2,label=f'HAP1 vs {h}')
    ax.set_title(t,fontsize=8.5); ax.set_xlabel('position along HAP1 (Mb)')
axs[0].set_ylabel('window distance (5-window mean)'); axs[0].legend(fontsize=7,frameon=False)
fig.suptitle('Schizothorax curvilabiatus (SchCurv1)',fontsize=9)
fig.tight_layout(); fig.savefig(F+'/snowcarp_windows.png',dpi=200); plt.close(fig)

# Fig: TE marker fraction per chromosome - diploid vs others
groups=[('ddMalSylv1','diploid'),('daGleHede1','allo anchor'),('drTriRepe1','allo anchor'),('daInuConz1','cryptic allo?'),('SchCurv1','auto (4 copies)'),('drLytSali1','auto (4 copies)')]
data=[];labels=[];colors=[]
for sp,kind in groups:
    p=f'results/{sp}/subgenomes/auto_allo_index.tsv'
    if not os.path.exists(p): p=f'superseded/{sp}_2026-08/subgenomes/auto_allo_index.tsv'
    if not os.path.exists(p): continue
    v=[float(r['te_marker_fraction']) for r in csv.DictReader(open(p),delimiter='\t') if r.get('te_marker_fraction')]
    data.append(v); labels.append(f'{sp}\n{kind}'); colors.append({'diploid':GREEN,'allo anchor':RED,'cryptic allo?':GOLD}.get(kind,BLUE))
fig,ax=plt.subplots(figsize=(8,3.2))
bp=ax.boxplot(data,patch_artist=True,widths=0.55,medianprops=dict(color='black'))
for b,cc in zip(bp['boxes'],colors): b.set_facecolor(cc); b.set_alpha(0.55)
for i,v in enumerate(data,1): ax.scatter(np.random.default_rng(1).normal(i,0.06,len(v)),v,s=7,color='black',alpha=0.6,zorder=3)
ax.set_xticks(range(1,len(labels)+1),labels,fontsize=7.5); ax.set_ylabel('te_marker_fraction\n(per chromosome / copy pair)')
ax.axhspan(min(data[0]),max(data[0]),color=GREEN,alpha=0.08); ax.text(0.55,max(data[0])+0.01,'range in the confirmed diploid',fontsize=7,color=GREEN)
fig.tight_layout(); fig.savefig(F+'/te_fraction.png',dpi=200); plt.close(fig)

# Fig: panel state composition (snapshot)
S=[r for r in csv.DictReader(open(D+'/panel/panel_summary.tsv'),delimiter='\t') if r['n_chromosome_numbers'] and r['n_tetrasomic_like']!='']
states=[('n_resolved_lineages','resolved lineages',RED),('n_fusion_lineages','fusion lineages','#9b2c2c'),('n_partially_resolved','partially resolved','#ed8936'),('n_candidate','candidate',GOLD),('n_one_divergent_copy','one divergent copy','#805ad5'),('n_tetrasomic_like','tetrasomic-like',BLUE),('n_not_assessable','not assessable','#cbd5e0')]
S.sort(key=lambda r:(-(int(r['n_resolved_lineages'])+int(r['n_fusion_lineages']))/max(1,int(r['n_chromosome_numbers'])), r['species']))
fig,ax=plt.subplots(figsize=(9,6.2))
y=np.arange(len(S)); left=np.zeros(len(S))
for k,lab,cc in states:
    v=np.array([int(r[k])/int(r['n_chromosome_numbers']) for r in S]); ax.barh(y,v,left=left,color=cc,label=lab,height=0.8); left+=v
ax.set_yticks(y,[f"{r['species']} ({r['copies_per_chromosome']}x{r['n_chromosome_numbers']})" for r in S],fontsize=7); ax.invert_yaxis()
ax.set_xlabel('fraction of chromosome numbers'); ax.legend(fontsize=7,ncol=4,loc='lower center',bbox_to_anchor=(0.5,1.0),frameon=False)
fig.tight_layout(); fig.savefig(F+'/panel_states.png',dpi=200); plt.close(fig)
# Fig: Ranunculus subgenomes and cross-subgenome homeolog pairs (homeologs stage outputs)
from matplotlib.patches import Rectangle
sp='dmRanRepe1'; d=pairs(sp)
part={r['metric']:r['value'] for r in csv.DictReader(open(f'results/{sp}/homeologs/subgenome_partition.tsv'),delimiter='\t')}
hp={}
for r in csv.DictReader(open(f'results/{sp}/homeologs/homeolog_pairs.tsv'),delimiter='\t'):
    a,b=int(r['chrom_a'][3:]),int(r['chrom_b'][3:]); hp[a]=b; hp[b]=a
def m(a,b): return st.mean(d[(f'{h1}_chr{a:02d}',f'{h2}_chr{b:02d}')] for h1 in ('HAP1','HAP2') for h2 in ('HAP1','HAP2'))
A=[int(c[3:]) for c in part['set_a'].split(',')]; B=[hp[a] for a in A]; order=A+B; n=len(A)
full=np.array([[np.nan if a==b else m(a,b) for b in order] for a in order])
block=np.array([[m(a,b) for b in B] for a in A])
within=[m(a,b) for g in (A,B) for i,a in enumerate(g) for b in g[i+1:]]
other=[block[i,j] for i in range(n) for j in range(n) if i!=j]
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,6.4),gridspec_kw=dict(width_ratios=[1,1.05]))
im=ax1.imshow(full,cmap='viridis_r',vmin=min(within)-0.0015,vmax=max(other))
lab=[f'chr{c:02d}' for c in order]
ax1.set_xticks(range(2*n),lab,rotation=90,fontsize=8); ax1.set_yticks(range(2*n),lab,fontsize=8)
for x in (0,n): ax1.add_patch(Rectangle((x-.5,x-.5),n,n,fill=False,ec='white',lw=2))
for i in range(n):
    ax1.add_patch(Rectangle((n+i-.5,i-.5),1,1,fill=False,ec='red',lw=1.5))
    ax1.add_patch(Rectangle((i-.5,n+i-.5),1,1,fill=False,ec='red',lw=1.5))
ax1.text((n-1)/2,-0.9,f'subgenome A (chr{A[0]:02d}–{A[-1]:02d})',ha='center',fontsize=9)
ax1.text(n+(n-1)/2,-0.9,'subgenome B (partner order)',ha='center',fontsize=9)
ax1.set_title('Ranunculus repens: 16 chromosome numbers, two subgenomes\n'
              '(mean k-mer distance over the four haplotype combinations)',fontsize=10,pad=30)
fig.colorbar(im,ax=ax1,fraction=0.046,pad=0.03,label='distance (k = 23)')
lo2=block.diagonal().min()-0.001; hi2=max(other)
im2=ax2.imshow(block,cmap='viridis_r',vmin=lo2,vmax=hi2)
ax2.set_xticks(range(n),[f'chr{b:02d}' for b in B],fontsize=9); ax2.set_yticks(range(n),[f'chr{a:02d}' for a in A],fontsize=9)
ax2.set_xlabel('subgenome B'); ax2.set_ylabel('subgenome A')
for i in range(n):
    for j in range(n):
        v=block[i,j]; t=(v-lo2)/(hi2-lo2)
        ax2.text(j,i,f'{v:.3f}',ha='center',va='center',fontsize=8,color='black' if t<0.45 else 'white',
                 fontweight='bold' if i==j else None)
    ax2.add_patch(Rectangle((i-.5,i-.5),1,1,fill=False,ec='red',lw=2))
ax2.set_title('A vs B: one homeolog per chromosome (red)\n'
              f'homeologs {block.diagonal().min():.3f}–{block.diagonal().max():.3f}; other A–B pairs '
              f'{min(other):.3f}–{max(other):.3f}; within a subgenome {min(within):.3f}–{max(within):.3f}',fontsize=10)
fig.colorbar(im2,ax=ax2,fraction=0.046,pad=0.03,label='distance (k = 23)')
fig.tight_layout(); fig.savefig(F+'/ranunculus_subgenomes.png',dpi=200); plt.close(fig)

# Fig: TE split per chromosome, Empetrum (lineages) vs the AAB weevil (one odd set)
fig,axs=plt.subplots(1,2,figsize=(10,3.6),gridspec_kw={'width_ratios':[13,10]})
for ax,sp,title in [(axs[0],'ddEmpNigr1','Empetrum nigrum (4 × 13)'),(axs[1],'icStrMela3','Strophosoma melanogrammum (AAB, 3 × 10)')]:
    L=list(csv.DictReader(open(f'results/{sp}/subgenomes/te_marker_fraction_by_lineage.tsv'),delimiter='\t'))
    x=np.arange(len(L)); w=0.38
    ax.bar(x-w/2,[float(r['within_te_marker_fraction']) for r in L],w,color=GREY,label='within lineage')
    ax.bar(x+w/2,[float(r['cross_te_marker_fraction']) for r in L],w,color=RED,label='across lineages')
    for i,r in enumerate(L):
        s=float(r['split_ratio']); ax.text(i,float(r['cross_te_marker_fraction'])+0.008,f'{s:.1f}',ha='center',fontsize=7,fontweight='bold' if s>=2 else 'normal')
    ax.set_xticks(x,[r['chrom'].replace('chr','') for r in L],fontsize=8); ax.set_xlabel('chromosome number (label: split ratio, bold ≥ 2)')
    ax.set_title(title,fontsize=9,style='italic')
axs[0].set_ylabel('TE marker fraction'); axs[0].legend(fontsize=7,frameon=False,loc='upper right')
fig.tight_layout(); fig.savefig(F+'/te_split.png',dpi=200); plt.close(fig)

print('ok', len(S))
