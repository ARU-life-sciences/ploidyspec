# Figures for docs/guide/ploidyspec_guide.typ. Run from the repo root with D=<output dir>
# (figures go to $D/fig). Needs $D/panel from `ploidyspec panel --outdir $D/panel`,
# $D/old_empnigr_windowed.tsv (pre-fix windowed_all.tsv from git), and the
# SchCurv1 chr17/chr19 position-free tracks at $TMPDIR/wcarp (or results/SchCurv1 once re-run).
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
W=os.environ['TMPDIR']+'/wcarp/windowed/windowed_all.tsv'
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
groups=[('ddMalSylv1','diploid'),('daGleHede1','allo anchor'),('drTriRepe1','allo anchor'),('daInuConz1','cryptic allo?'),('SchCurv1','auto (4 copies)'),('drLytSali1','auto (4 copies,\nprevious assembly)')]
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
print('ok', len(S))
