# E2: type I (skewed, cell-heteroscedastic and normal errors) and power of HC2 vs HC3 for beta_PS on the real n=112 design.
import numpy as np
from scipy import stats
exec(open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_scripts/e2_labels_hc3.py').read().split("out = []")[0])
rng = np.random.default_rng(20260930); out = []
a = A_[1]; tc = stats.t.ppf(.975, dfr)
def tests(Y):
    Bh = A_ @ Y; R = Y - X @ Bh
    se3 = np.sqrt(np.einsum('j,jn,j->n', a, (R/(1-h)[:, None])**2, a)); se2 = np.sqrt(np.einsum('j,jn,j->n', a, R**2/(1-h)[:, None], a))
    return Bh[1]/se3, Bh[1]/se2
N = 40000
e_sk = ((rng.chisquare(4, (n, N)) - 4)/np.sqrt(8)) * np.where(X[:, 1] > 0, 1.3, 1.0)[:, None]
for lab, E in [('skewed+heteroscedastic', e_sk), ('normal', rng.standard_normal((n, N)))]:
    t3, t2 = tests(E); s = f'type I ({lab}): HC3 {(np.abs(t3) > tc).mean():.3f}  HC2 {(np.abs(t2) > tc).mean():.3f}'; print(s); out.append(s)
for d in [-.4, -.5, -.6]:
    t3, t2 = tests(d*ps[:, None] + rng.standard_normal((n, N))); s = f'power delta {d}: HC3 {(t3 < -tc).mean():.3f}  HC2 {(t2 < -tc).mean():.3f}'; print(s); out.append(s)
open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_hc2_vs_hc3.txt', 'w').write('\n'.join(out) + '\n')
