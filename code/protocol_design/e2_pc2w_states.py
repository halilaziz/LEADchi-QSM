# E2: three-state PC2w gate (PASS / INDETERMINATE / FAIL) operating characteristics; design-only simulation.
import numpy as np; from scipy import stats
rng = np.random.default_rng(20260925); out = []
def P(s): print(s); out.append(s)
C0 = np.array([[1, .5], [.5, 1]])
P('PC2w three-state rule: PASS = pooled (PUT+GP)/2 change > 0 with one-sided paired-t p < .05 and both ROI mean changes >= 0;')
P('FAIL = pooled mean change <= 0; INDETERMINATE otherwise. Mean changes: Zhang 2018 model (PUT 2.6, GP 4.4 ppb over 1.95 y); 0 = null.')
for mu_lab, mu in [('Zhang', (2.6, 4.4)), ('half-Zhang', (1.3, 2.2)), ('null', (0., 0.))]:
    for sd in [8, 12, 16]:
        for n in [48, 44]:
            L = np.linalg.cholesky(C0 * sd**2); res = np.zeros(3)
            for _ in range(10000):
                x = np.array(mu) + rng.standard_normal((n, 2)) @ L.T; pooled = x.mean(1); mb = pooled.mean()
                t = mb / (pooled.std(ddof=1) / np.sqrt(n)); p = stats.t.sf(t, n - 1)
                if mb <= 0: res[2] += 1
                elif p < .05 and x[:, 0].mean() >= 0 and x[:, 1].mean() >= 0: res[0] += 1
                else: res[1] += 1
            P(f'  {mu_lab:10s} sigma_D={sd:2d} n={n}: PASS {res[0]/1e4:.2f} / INDET {res[1]/1e4:.2f} / FAIL {res[2]/1e4:.2f}')
open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_pc2w_states.txt', 'w').write('\n'.join(out) + '\n')
