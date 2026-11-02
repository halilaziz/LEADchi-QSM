# E2 (LEAD-chi protocol v1.2 delta): LABEL-FREE reliability of right-SCP (pyAFQ SCP_R) proximal-segment DWI metrics
# from the existing Paper 0 profile files. No group labels are read (the 'group' column is never loaded).
# No MEGRE / chi data are touched. Outputs: e2_scpr_reliability.csv / .txt next to this folder.
import numpy as np, pandas as pd, os
P0 = '/path/to/Penn_LEAD/paper0/data/'
OUT = '/path/to/Penn_LEAD/LEADchi/prereg/work/'
cols = ['subject', 'session', 'bundle', 'metric', 'node', 'value']       # 'group' deliberately excluded
def load(f, metrics):
    it = pd.read_csv(P0 + f, usecols=cols, chunksize=2_000_000)
    return pd.concat([c[(c.bundle == 'SCPR') & (c.metric.isin(metrics))] for c in it])
d = pd.concat([load('vta_noddi_profiles_long.csv', ['ICVF_NDI', 'ISOVF']),
               load('noddi_mod_profiles_long.csv', ['ICVF_mod', 'TF']),
               load('fw_b1000_pc_profiles_long.csv', ['FA_t_pc', 'MD_t_pc', 'TW_pc'])])
w = d.pivot_table(index=['subject', 'session', 'node'], columns='metric', values='value').reset_index()
w['twICVF'] = w.ICVF_mod / w.TF                      # Parker 2021: tissue-weighted mean = sum(mod) / sum(TF)
lines = []
def P(s): print(s); lines.append(s)
P(f'SCPR profiles: {w.groupby(["subject","session"]).ngroups} sessions, {w.subject.nunique()} subjects')
# consistency checks of the metric definitions (label-free)
P(f'  node-level corr(TF, 1-ISOVF) = {np.corrcoef(w.TF, 1 - w.ISOVF)[0,1]:.4f}; median |TF-(1-ISOVF)| = {np.median(np.abs(w.TF-(1-w.ISOVF))):.4f}')
P(f'  node-level corr(ICVF_mod, ICVF_NDI*TF) = {np.corrcoef(w.ICVF_mod, w.ICVF_NDI*w.TF)[0,1]:.4f}; median rel. diff = {np.median(np.abs(w.ICVF_mod-w.ICVF_NDI*w.TF)/w.ICVF_mod):.4f}')
excl = [l.strip() for l in open(P0 + 'forensics/exclude_afq_qc.txt') if l.startswith('sub-')]
dq = pd.read_csv(P0 + 'pennlinc_qc/EF_desc-diffusion_qc.csv', encoding='utf-8-sig')
dq['subject'] = 'sub-' + dq.participant_id.astype(str); dq['session'] = 'ses-' + dq.session_id.astype(str)
scal_pass = set(map(tuple, dq.loc[dq.qc_determination_scalar_maps.str.strip() == 'pass', ['subject', 'session']].values))
def icc31(y1, y2):
    Y = np.column_stack([y1, y2]); n, k = Y.shape; gm = Y.mean()
    ss_r = k*((Y.mean(1)-gm)**2).sum(); ss_c = n*((Y.mean(0)-gm)**2).sum(); ss_e = ((Y-gm)**2).sum()-ss_r-ss_c
    ms_r = ss_r/(n-1); ms_e = ss_e/((n-1)*(k-1)); return (ms_r-ms_e)/(ms_r+(k-1)*ms_e)
def icc_ci(icc, n, k=2, a=.05):   # Shrout-Fleiss F-based CI for ICC(3,1)
    from scipy import stats
    F = (1+(k-1)*icc)/(1-icc); dfn, dfd = n-1, (n-1)*(k-1)
    FL = F/stats.f.ppf(1-a/2, dfn, dfd); FU = F*stats.f.ppf(1-a/2, dfd, dfn)
    return (FL-1)/(FL+k-1), (FU-1)/(FU+k-1)
rows = []
segs = {'0-30 (v1.2 candidate)': (0, 30), '10-30 (fallback buffer)': (10, 30), '0-16 (round-1 range)': (0, 16), '0-99 (whole)': (0, 99)}
for samp in ['P0_frozen', 'scalar_pass']:
    ww = w[~w.subject.isin(excl)] if samp == 'P0_frozen' else w[[ (s, e) in scal_pass for s, e in zip(w.subject, w.session)]]
    for seg, (a, b) in segs.items():
        s = ww[(ww.node >= a) & (ww.node <= b)].copy()
        for rule in ['plain', 'ISOVF<=0.30']:
            ss = s if rule == 'plain' else s[s.ISOVF <= 0.30]
            nn = ss.groupby(['subject', 'session']).size()
            keep = nn[nn >= 10].index if rule != 'plain' else nn.index
            m = ss.set_index(['subject', 'session']).loc[keep].reset_index()
            for met in ['twICVF', 'ICVF_mod', 'ICVF_NDI', 'FA_t_pc', 'MD_t_pc']:
                mm = m if not met.endswith('_pc') else m[m.TW_pc >= 0.2]
                v = mm.groupby(['subject', 'session'])[met].mean().unstack('session')
                if 'ses-1' not in v or 'ses-2' not in v: continue
                pr = v[['ses-1', 'ses-2']].dropna()
                if len(pr) < 10: continue
                icc = icc31(pr['ses-1'].values, pr['ses-2'].values); lo, hi = icc_ci(icc, len(pr))
                rows.append(dict(sample=samp, segment=seg, node_rule=rule, metric=met, n_pairs=len(pr),
                                 n_sessions_retained=len(keep), icc31=icc, icc_lo95=lo, icc_hi95=hi,
                                 r=np.corrcoef(pr['ses-1'], pr['ses-2'])[0, 1]))
R = pd.DataFrame(rows); R.to_csv(OUT + 'e2_scpr_reliability.csv', index=False)
P('\nICC(3,1) ses-1 vs ses-2 returners (raw values; includes ~1.7 y of development -> lower bound on single-scan reliability)')
P(R.round(3).to_string(index=False))
# node-rule survival (label-free): sessions with >= 10 nodes in 0-30 after ISOVF<=0.30
s = w[(w.node <= 30)]; nn = s[s.ISOVF <= 0.30].groupby(['subject', 'session']).size().reindex(s.groupby(['subject','session']).size().index, fill_value=0)
P(f'\nNode-rule survival, nodes 0-30, ISOVF<=0.30: {int((nn>=10).sum())}/{len(nn)} sessions keep >=10 nodes; median retained {int(nn.median())}')
s = w[(w.node >= 10) & (w.node <= 30)]; nn = s[s.ISOVF <= 0.30].groupby(['subject', 'session']).size().reindex(s.groupby(['subject','session']).size().index, fill_value=0)
P(f'Node-rule survival, nodes 10-30, ISOVF<=0.30: {int((nn>=10).sum())}/{len(nn)} sessions keep >=10 nodes; median retained {int(nn.median())}')
pool = w.groupby('node')[['ISOVF', 'TW_pc']].median()
P('\nPooled (label-free) median ISOVF / TW_pc at nodes 0,5,10,15,20,25,30,35,40,50: ' +
  '; '.join(f'{k}: {pool.loc[k,"ISOVF"]:.2f}/{pool.loc[k,"TW_pc"]:.2f}' for k in [0,5,10,15,20,25,30,35,40,50]))
open(OUT + 'e2_scpr_reliability.txt', 'w').write('\n'.join(lines) + '\n')
