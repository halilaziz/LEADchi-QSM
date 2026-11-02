# E2: Monte-Carlo power of outcomes A and B with the REGISTERED HC3 test on the real first-usable design (n=112),
# versus the analytic (model-based SE, sigma known) numbers of v1.1 power_v2.txt. Normal errors; 20,000 draws.
import numpy as np, pandas as pd, patsy, warnings
from scipy import stats
warnings.filterwarnings('ignore'); rng = np.random.default_rng(20260926)
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
out = []
def P(s): print(s); out.append(s)
res = {}
for rep in range(4):   # 4 draws of the simulated roi_vol / Q_phase columns
    X = np.column_stack([np.ones(n), ps, ad, ps*ad, sex, rb-ro, rw-ro, B, rng.standard_normal(n), rng.standard_normal(n)])
    k = X.shape[1]; dfr = n - k; XtXi = np.linalg.inv(X.T @ X); h = np.einsum('ij,jk,ik->i', X, XtXi, X); H = XtXi @ X.T
    t975 = stats.t.ppf(.975, dfr); t95 = stats.t.ppf(.95, dfr); nsim = 5000
    for delta in [0, -.1, -.2, -.3, -.4, -.5, -.6]:
        Y = delta*ps[:, None] + rng.standard_normal((n, nsim)); Bh = H @ Y; R = Y - X @ Bh; s = np.sqrt((R**2).sum(0)/dfr)
        W = (R/(1-h)[:, None])**2; se_hc3 = np.sqrt(np.einsum('j,jn,j->n', (XtXi @ X.T)[1], W, (XtXi @ X.T)[1]))
        se_ols = s*np.sqrt(XtXi[1, 1])
        for lab, se in [('HC3', se_hc3), ('OLS', se_ols)]:
            A = ((Bh[1]/se) < -t975).mean(); Bo = ((Bh[1] + .6*np.sqrt(.7)*s)/se > t95).mean(); B6 = ((Bh[1] + .6*s)/se > t95).mean()
            res.setdefault((lab, delta), []).append((A, Bo, B6))
for lab in ['HC3', 'OLS']:
    for delta in [0, -.1, -.2, -.3, -.4, -.5, -.6]:
        a = np.mean(res[(lab, delta)], 0)
        P(f'{lab}: true delta {delta:+.1f}: P(A) {a[0]:.3f}  P(B, bound -0.50) {a[1]:.3f}  P(B, bound -0.60) {a[2]:.3f}')
open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_hc3_power_check.txt', 'w').write('\n'.join(out) + '\n')
