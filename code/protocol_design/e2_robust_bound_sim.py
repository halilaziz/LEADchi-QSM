# E2: operating characteristics of the v1.2 outcome-B scale rule (Paper 0 lesson iv: dominant scans inflate the SD that
# converts an SD-unit bound into ppb). Registered v1.1: t_B = (b + dB*sigma_hat)/SE_HC3. v1.2: sigma_B = min(sigma_hat, s_rob),
# s_rob = Huber proposal-2 residual scale of the same model (statsmodels RLM). Real first-usable design (n=112), simulated errors.
import numpy as np, pandas as pd, patsy, warnings, statsmodels.api as sm
from scipy import stats
warnings.filterwarnings('ignore'); rng = np.random.default_rng(20260925)
D = '/path/to/scratch/data/'
inv = pd.read_csv(D + 'session_inventory_with_pennlinc_qc.csv'); p = pd.read_csv(D + 'participants.tsv', sep='\t')
inv = inv.merge(p[['participant_id', 'dx_adhd', 'race']], left_on='sub', right_on='participant_id', how='left')
inv['PS'] = (inv.study_group == 'PRO/CHR').astype(int); inv['AD'] = inv.dx_adhd.fillna(0).astype(int); inv['sesn'] = inv.ses.str[-1].astype(int)
t1 = pd.read_csv(D + 'qc/EF_desc-T1_qc.csv', encoding='utf-8-sig'); t1['t1qc'] = t1.qc_determination.str.strip()
inv = inv.merge(t1[['participant_id', 'session_id', 't1qc']], left_on=['sub', 'ses'], right_on=['participant_id', 'session_id'], how='left', suffixes=('', '_t1'))
m = inv[(inv.MEGRE_full == 1) & (inv.T1w == 1) & (inv.t1qc.fillna('missing') != 'fail')]
d = m.sort_values('sesn').groupby('sub').head(1); n = len(d)
ps = np.where(d.PS == 1, .5, -.5); ad = np.where(d.AD == 1, .5, -.5); sex = np.where(d.sex == 'M', .5, -.5)
rb = d.race.str.startswith('Black').astype(float).values; rw = (d.race == 'White').astype(float).values; ro = 1 - rb - rw
B = np.asarray(patsy.dmatrix('cr(x, df=4) - 1', {'x': d.age.values}))[:, 1:]
X = np.column_stack([np.ones(n), ps, ad, ps*ad, sex, rb-ro, rw-ro, B, rng.standard_normal(n), rng.standard_normal(n)])
k = X.shape[1]; dfr = n - k; XtXi = np.linalg.inv(X.T @ X); h = np.einsum('ij,jk,ik->i', X, XtXi, X); t95 = stats.t.ppf(.95, dfr)
dB = .6*np.sqrt(.7); out = []
def run(delta, contam, nsim=3000):
    reg = rob = 0
    for _ in range(nsim):
        e = rng.standard_normal(n)
        if contam:   # two sessions with implausible-but-retained values (+/-6 SD), as in Paper 0's MD_t composite
            j = rng.choice(n, 2, replace=False); e[j] = rng.choice([-6, 6], 2)
        y = delta*ps + e
        b = XtXi @ X.T @ y; r = y - X @ b; s = np.sqrt(r @ r / dfr)
        V = XtXi @ (X.T * (r/(1-h))**2) @ X @ XtXi; se = np.sqrt(V[1, 1])
        srob = sm.RLM(y, X, M=sm.robust.norms.HuberT()).fit(scale_est=sm.robust.scale.HuberScale()).scale
        reg += (b[1] + dB*s)/se > t95; rob += (b[1] + dB*min(s, srob))/se > t95
    return reg/nsim, rob/nsim
for contam in [False, True]:
    for delta in [0, -.2, -.5, -.6]:
        a, b_ = run(delta, contam)
        s = f'contamination={contam!s:5s} true delta={delta:+.1f}: P(outcome B) registered {a:.3f} | v1.2 min-scale rule {b_:.3f}'
        print(s); out.append(s)
open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_robust_bound_sim.txt', 'w').write('\n'.join(out) + '\n')
