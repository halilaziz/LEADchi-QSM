# E2: HC3 Monte-Carlo power after random (group-independent) MEGRE-QC loss; A at delta -0.6, B(-0.50/-0.60) at delta 0.
import numpy as np
from scipy import stats
exec(open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_scripts/e2_labels_hc3.py').read().split("out = []")[0])
rng = np.random.default_rng(20260929); out = []
Xf = X.copy(); psf = ps.copy()
for keep in [104, 100, 96]:
    acc = []
    for rep in range(40):
        idx = rng.choice(n, keep, replace=False); Xk = Xf[idx]; pk = psf[idx]; kk = Xk.shape[1]; dfk = keep - kk
        XtXi_ = np.linalg.inv(Xk.T @ Xk); A2 = XtXi_ @ Xk.T; hk = np.einsum('ij,jk,ik->i', Xk, XtXi_, Xk)
        res = []
        for delta in [-.6, 0]:
            Y = delta*pk[:, None] + rng.standard_normal((keep, 2000)); Bh = A2 @ Y; R = Y - Xk @ Bh; s = np.sqrt((R**2).sum(0)/dfk)
            W = (R/(1-hk)[:, None])**2; se = np.sqrt(np.einsum('j,jn,j->n', A2[1], W, A2[1]))
            res.append(((Bh[1]/se) < -stats.t.ppf(.975, dfk)).mean() if delta < 0 else ((Bh[1] + .6*np.sqrt(.7)*s)/se > stats.t.ppf(.95, dfk)).mean())
            if delta == 0: res.append(((Bh[1] + .6*s)/se > stats.t.ppf(.95, dfk)).mean())
        acc.append(res)
    a = np.mean(acc, 0); s_ = f'n={keep}: P(A | delta -0.6) {a[0]:.2f}; P(B -0.50 | 0) {a[1]:.2f}; P(B -0.60 | 0) {a[2]:.2f}'; print(s_); out.append(s_)
open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_hc3_afterqc.txt', 'w').write('\n'.join(out) + '\n')
