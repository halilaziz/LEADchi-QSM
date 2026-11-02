# E2: S1 contrast L and "dissociation" label power with HC3 SEs (v1.1 used model-based SEs). Real n=112 design.
import numpy as np, pandas as pd, patsy, warnings
from scipy import stats
exec(open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_scripts/e2_labels_hc3.py').read().split("out = []")[0])
rng = np.random.default_rng(20260928); out = []
def P(s): print(s); out.append(s)
e = np.eye(k); cL = e[1] - e[2]
def se(c, W): a = c @ A_; return np.sqrt(np.einsum('j,jn,j->n', a, W, a))
P('S1 with HC3 (20,000 draws): L power at Holm worst .0167 / .05 ; P(dissociation label) at .0167 / .05')
for rho in [.3, .5, .7]:
    sdD = np.sqrt(2*(1-rho))
    for (psS, psT, adS, adT) in [(-.4, -.1, -.1, -.4), (-.5, -.1, 0, -.4)]:
        bps = psS-psT; bad = adS-adT; mu = bps*ps + bad*ad
        Y = mu[:, None] + sdD*rng.standard_normal((n, 20000)); Bh = A_ @ Y; R = Y - X @ Bh; W = (R/(1-h)[:, None])**2
        tL = (cL @ Bh)/se(cL, W); tps = Bh[1]/se(e[1], W); tad = Bh[2]/se(e[2], W); comp = (tps < -t1c) & (tad > t1c)
        r = [((tL < 0) & (np.abs(tL) > stats.t.ppf(1-a/2, dfr))) for a in (.0167, .05)]
        P(f'  rho={rho} L={bps-bad:.1f}: L {r[0].mean():.2f}/{r[1].mean():.2f}; label {(r[0]&comp).mean():.2f}/{(r[1]&comp).mean():.2f}')
open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_s1_hc3.txt', 'w').write('\n'.join(out) + '\n')
