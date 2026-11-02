# E2: claim-label probabilities (sec. 7.3) with the REGISTERED HC3 standard errors (v1.1 used model-based SEs).
import numpy as np, pandas as pd, patsy, warnings
from scipy import stats
warnings.filterwarnings('ignore'); rng = np.random.default_rng(20260927)
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
k = X.shape[1]; dfr = n - k; XtXi = np.linalg.inv(X.T @ X); h = np.einsum('ij,jk,ik->i', X, XtXi, X); A_ = XtXi @ X.T
tc = stats.t.ppf(.975, dfr); t1c = stats.t.ppf(.95, dfr)
c1 = np.zeros(k); c1[1] = 1; c1[3] = -.5; c2 = np.zeros(k); c2[1] = 1; c2[3] = .5
out = []
def P(s): print(s); out.append(s)
def sim(dPSo, dPSc, dAD, nsim=20000):
    mu = np.select([(ps < 0) & (ad < 0), (ps < 0) & (ad > 0), (ps > 0) & (ad < 0), (ps > 0) & (ad > 0)], [0, dAD, dPSo, dAD + dPSc])
    Y = mu[:, None] + rng.standard_normal((n, nsim)); Bh = A_ @ Y; R = Y - X @ Bh; s = np.sqrt((R**2).sum(0)/dfr)
    W = (R/(1-h)[:, None])**2
    def se(c): a = c @ A_; return np.sqrt(np.einsum('j,jn,j->n', a, W, a))
    e = np.eye(k); sps, sad, sint = se(e[1]), se(e[2]), se(e[3]); se1, se2 = se(c1), se(c2)
    A = Bh[1]/sps < -tc; e1 = Bh[1] - .5*Bh[3]; e2 = Bh[1] + .5*Bh[3]
    tost = ((Bh[2] + .5*s)/sad > t1c) & ((Bh[2] - .5*s)/sad < -t1c); ci0 = np.abs(Bh[3]/sint) < tc
    lab = dict(A=A.mean(), B=((Bh[1] + .6*np.sqrt(.7)*s)/sps > t1c).mean(), PS_effect=(A & (e1 < 0) & (e2 < 0)).mean(),
               PS_specific=(A & (e1 < 0) & (e2 < 0) & tost & ci0).mean(),
               comorbid=(A & (e2/se2 < -tc) & ~(e1/se1 < -tc) & (e1 > -.5*np.abs(e2))).mean(), TOST=tost.mean())
    return {kk: round(float(v), 2) for kk, v in lab.items()}
P('Claim-label probabilities with HC3 SEs (20,000 draws; normal errors; real n=112 design); B uses bound -0.50')
for sc in [(0, 0, 0), (-.4, -.4, 0), (-.6, -.6, 0), (0, -.8, 0)]:
    P('  PS-only %.2f, comorbid %.2f, ADHD %.2f -> ' % sc + str(sim(*sc)))
open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_labels_hc3.txt', 'w').write('\n'.join(out) + '\n')
