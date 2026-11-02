# E2 (LEAD-chi v1.2 delta): design numbers that change in v1.2. No MEGRE/chi values exist or are used.
# (1) S3 sample under v1.2 (existing pyAFQ SCP_R profiles from Paper 0), (2) S3a power with Paper 0 measured DWI reliability,
# (3) power of the new label- and age-free positive control PC2w, (4) outcome-B/H1 unchanged (see e2_power_v2_rerun.txt).
import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(20260925)
D = '/path/to/scratch/data/'
P0 = '/path/to/Penn_LEAD/paper0/data/'
W = '/path/to/Penn_LEAD/LEADchi/prereg/work/'
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)
inv = pd.read_csv(D + 'session_inventory_with_pennlinc_qc.csv')
p = pd.read_csv(D + 'participants.tsv', sep='\t')
inv = inv.merge(p[['participant_id', 'dx_adhd']], left_on='sub', right_on='participant_id', how='left')
inv['PS'] = (inv.study_group == 'PRO/CHR').astype(int); inv['AD'] = inv.dx_adhd.fillna(0).astype(int)
inv['cell'] = np.select([(inv.PS == 0) & (inv.AD == 0), (inv.PS == 0) & (inv.AD == 1), (inv.PS == 1) & (inv.AD == 0)], ['TD', 'ADHD-only', 'PS-only'], 'PS+ADHD')
inv['sesn'] = inv.ses.str[-1].astype(int)
t1 = pd.read_csv(D + 'qc/EF_desc-T1_qc.csv', encoding='utf-8-sig'); t1['t1qc'] = t1.qc_determination.str.strip()
inv = inv.merge(t1[['participant_id', 'session_id', 't1qc']], left_on=['sub', 'ses'], right_on=['participant_id', 'session_id'], how='left', suffixes=('', '_t1'))
dq = pd.read_csv(D + 'qc/EF_desc-diffusion_qc.csv', encoding='utf-8-sig')
dq['sub'] = 'sub-' + dq.participant_id.astype(str); dq['ses'] = 'ses-' + dq.session_id.astype(str)
inv = inv.merge(dq[['sub', 'ses', 'qc_determination_scalar_maps']], on=['sub', 'ses'], how='left')
m = inv[(inv.MEGRE_full == 1) & (inv.T1w == 1) & (inv.t1qc.fillna('missing') != 'fail')].copy()
first = m.sort_values('sesn').groupby('sub').head(1).copy()
order = ['TD', 'ADHD-only', 'PS-only', 'PS+ADHD']
def cells(df): c = df.cell.value_counts(); return ' / '.join(str(int(c.get(k, 0))) for k in order) + f'  (total {len(df)})'
P('H1 first-usable sample (pre-MEGRE-QC), TD / ADHD-only / PS-only / PS+ADHD:', cells(first))
# (1) S3 sample, v1.2: H1 session + same-session DWI with PennLINC scalar pass + SCP_R profile present in the Paper 0 files
geo = pd.read_csv(P0 + 'forensics/bundle_geometry.csv'); geo = geo[geo.bundle == 'SCPR']
prof = set(map(tuple, geo[['subject', 'session']].values))
f2 = first.copy(); f2['scal'] = f2.qc_determination_scalar_maps.fillna('').str.strip() == 'pass'
f2['prof'] = [(a, b) in prof for a, b in zip(f2['sub'], f2['ses'])]
s3 = f2[(f2.DWI == 1) & f2.scal & f2.prof]
P('S3 v1.2 (H1 session, same-session DWI scalar pass, SCP_R profile exists):', cells(s3))
P('  H1 subjects not in S3:', ', '.join(f"{r['sub']} {r.ses} (DWI={r.DWI}, scalar={r.qc_determination_scalar_maps}, profile={r.prof})" for _, r in f2[~f2.index.isin(s3.index)].iterrows()))
excl = [l.strip() for l in open(P0 + 'forensics/exclude_afq_qc.txt') if l.startswith('sub-')]
s3c = s3[~s3['sub'].isin(excl)]
P('S3 conservative (drop the 20-subject Paper 0 / round-1 AFQ-QC list):', cells(s3c))
# node-rule survival for the S3 sessions (label-free DWI values; ISOVF <= 0.30 in nodes 0-30; >= 10 nodes)
it = pd.read_csv(P0 + 'vta_noddi_profiles_long.csv', usecols=['subject', 'session', 'bundle', 'metric', 'node', 'value'], chunksize=2_000_000)
iso = pd.concat([c[(c.bundle == 'SCPR') & (c.metric == 'ISOVF') & (c.node <= 30)] for c in it])
g = geo.set_index(['subject', 'session']).mean_len_mm
iso['arc'] = iso.node.values * g.reindex(list(zip(iso.subject, iso.session))).values / 99.0
def survive(rule):
    ss = iso if rule == 'iso' else iso[iso.arc >= 3.0]
    n = ss[ss.value <= 0.30].groupby(['subject', 'session']).size()
    idx = list(zip(s3['sub'], s3.ses)); n = n.reindex(idx, fill_value=0)
    return n
for lab, rule in [('ISOVF<=0.30 only', 'iso'), ('ISOVF<=0.30 + arc>=3 mm from node 0 (fallback buffer)', 'arc')]:
    n = survive(rule); P(f'  S3 sessions keeping >=10 nodes in 0-30 [{lab}]: {int((n >= 10).sum())}/{len(n)}; >=8: {int((n >= 8).sum())}; median nodes {int(n.median())}')
P(f'  pooled SCP_R length (label-free): mean {geo.mean_len_mm.mean():.1f} mm (SD {geo.mean_len_mm.std():.1f}); node spacing ~{geo.mean_len_mm.mean()/99:.2f} mm')
# (2) S3a power: observed r = rho_true * sqrt(R_chi * R_dwi). R_dwi from Paper 0 profiles (e2_scpr_reliability.csv, label-free)
def pw_r(r, nn, kcov, a=.05):
    return stats.norm.sf(stats.norm.ppf(1 - a / 2) - np.arctanh(r) * np.sqrt(nn - kcov - 3))
R = pd.read_csv(W + 'e2_scpr_reliability.csv')
sel = R[(R.segment.str.startswith('0-30')) & (R.node_rule == 'ISOVF<=0.30')]
P('Paper 0 measured cross-session ICC(3,1), SCP_R nodes 0-30, ISOVF<=0.30 (label-free):')
for _, r in sel.iterrows(): P(f'  {r["sample"]:12s} {r.metric:9s} ICC {r.icc31:.2f} [{r.icc_lo95:.2f}, {r.icc_hi95:.2f}], n pairs {int(r.n_pairs)}')
kS3 = 20
P(f'S3a power (two-sided; k={kS3} covariates); rows: n; cols: true (disattenuated) partial rho; R_chi = split-echo dentate reliability; R_dwi from Paper 0')
for Rdwi_lab, Rdwi in [('twICVF point 0.55', .55), ('twICVF lower 95% bound 0.33', .33), ('FA_t_pc 0.26', .26)]:
    for Rchi in [.70, .85]:
        for nn in [108, 101, 94, 88]:
            row = []
            for rho in [.30, .40, .50]:
                ro = rho * np.sqrt(Rchi * Rdwi); row.append(f'rho {rho:.2f}->r_obs {ro:.3f}: {pw_r(ro, nn, kS3, .0167):.2f}/{pw_r(ro, nn, kS3):.2f}')
            P(f'  R_dwi={Rdwi_lab}, R_chi={Rchi}, n={nn}: ' + ' | '.join(row) + '   (Holm worst .0167 / .05)')
for nn in [108, 94]:
    P(f'  v1.1 reading (observed r = .30): n={nn}: {pw_r(.30, nn, kS3, .0167):.2f}/{pw_r(.30, nn, kS3):.2f}; implied true rho at R_chi .70, R_dwi .55: {.30/np.sqrt(.70*.55):.2f}')
# (3) PC2w: within-person increase of putamen and GP chi between the first two usable coded sessions (no age, no labels)
m2 = m.sort_values(['sub', 'sesn']).copy(); m2['k'] = m2.groupby('sub').cumcount()
pairs = m2[m2.k <= 1].groupby('sub').filter(lambda x: len(x) == 2)
dt = pairs.groupby('sub').age.agg(lambda a: a.iloc[1] - a.iloc[0])
P(f'PC2w design: {dt.size} returners with 2 T1-eligible MEGRE sessions; interval {dt.mean():.2f} +/- {dt.std():.2f} y (range {dt.min():.2f}-{dt.max():.2f})')
# expected annual change at the returners' mean first-session age from two published lifespan models (whole-brain referenced)
a0 = pairs[pairs.k == 0].age.mean()
def rate(alpha, beta, a): return alpha * beta * np.exp(-beta * a) * 1000   # ppb / y
mods = {'Li 2014 (HBM 35:2698)': {'PUT': (0.20, 0.01), 'GP': (0.13, 0.11)},
        'Zhang 2018 (NeuroImage 171:176)': {'PUT': (0.078, 0.023), 'GP': (0.086, 0.049)}}
for k_, v in mods.items():
    P(f'  {k_} at age {a0:.1f}: PUT {rate(*v["PUT"], a0):.2f} ppb/y, GP {rate(*v["GP"], a0):.2f} ppb/y -> expected change over {dt.mean():.2f} y: PUT {rate(*v["PUT"], a0)*dt.mean():.1f}, GP {rate(*v["GP"], a0)*dt.mean():.1f} ppb')
P('  PC2w rule: pooled change of (PUT+GP)/2 > 0, one-sided p < .05 (paired t), AND neither ROI mean change < 0.')
P('  Power by SD of the within-person 2-y change (sigma_D, ppb; includes measurement error of both sessions and slope heterogeneity);')
P('  columns n pairs = 53 / 48 / 44; ROI changes correlated 0.5; mean changes from the lower (Zhang) and higher (Li) model')
for lab, (mP, mG) in [('Zhang', (rate(.078, .023, a0)*dt.mean(), rate(.086, .049, a0)*dt.mean())), ('Li', (rate(.20, .01, a0)*dt.mean(), rate(.13, .11, a0)*dt.mean()))]:
    for sd in [4, 6, 8, 12, 16]:
        row = []
        for npair in [53, 48, 44]:
            C = np.array([[1, .5], [.5, 1]]) * sd ** 2; L = np.linalg.cholesky(C); ok = 0; nsim = 4000
            for _ in range(nsim):
                x = np.array([mP, mG]) + rng.standard_normal((npair, 2)) @ L.T
                pooled = x.mean(1); t = pooled.mean() / (pooled.std(ddof=1) / np.sqrt(npair))
                ok += (stats.t.sf(t, npair - 1) < .05) and (x[:, 0].mean() >= 0) and (x[:, 1].mean() >= 0)
            row.append(f'{ok/nsim:.2f}')
        P(f'  {lab:5s} (PUT {mP:.1f}, GP {mG:.1f} ppb) sigma_D={sd:2d}: ' + ' / '.join(row))
P('  False-pass rate if the true change is 0 (sigma_D=8, n=48): ' + f'{np.mean([(stats.ttest_1samp((z:=rng.standard_normal((48,2))@np.linalg.cholesky(np.array([[1,.5],[.5,1]])*64).T).mean(1), 0).statistic > stats.t.ppf(.95, 47)) and (z[:,0].mean()>=0) and (z[:,1].mean()>=0) for _ in range(4000)]):.3f}')
open(W + 'e2_design_v12.txt', 'w').write('\n'.join(out) + '\n')
