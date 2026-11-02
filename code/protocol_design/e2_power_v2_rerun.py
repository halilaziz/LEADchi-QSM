# v1.1 power recomputation: registered primary model (ns(age, df=3) = 3 columns, + roi_vol + Q_phase),
# S3 restricted to H1 session, T1 fixed-margin rule across ICC, S1 dissociation IUT label, B with attenuation-adjusted bound.
# No imaging values exist; roi_vol and Q_phase are simulated as standard-normal columns (their values do not exist yet).
import pandas as pd, numpy as np, patsy, warnings
from scipy import stats
warnings.filterwarnings('ignore')
rng = np.random.default_rng(20260924)
D = '/path/to/scratch/data/'
OUT = '/path/to/Penn_LEAD/LEADchi/prereg/work/e2_power_v2_rerun.txt'
inv = pd.read_csv(D + 'session_inventory_with_pennlinc_qc.csv')
p = pd.read_csv(D + 'participants.tsv', sep='\t')
inv = inv.merge(p[['participant_id', 'dx_adhd', 'race']], left_on='sub', right_on='participant_id', how='left')
inv['PS'] = (inv.study_group == 'PRO/CHR').astype(int); inv['AD'] = inv.dx_adhd.fillna(0).astype(int)
inv['sesn'] = inv.ses.str[-1].astype(int)
t1 = pd.read_csv(D + 'qc/EF_desc-T1_qc.csv', encoding='utf-8-sig'); t1['t1qc'] = t1.qc_determination.str.strip()
inv = inv.merge(t1[['participant_id', 'session_id', 't1qc']], left_on=['sub', 'ses'], right_on=['participant_id', 'session_id'], how='left', suffixes=('', '_t1'))
dq = pd.read_csv(D + 'qc/EF_desc-diffusion_qc.csv', encoding='utf-8-sig')
dq['sub'] = 'sub-' + dq.participant_id.astype(str); dq['ses'] = 'ses-' + dq.session_id.astype(str)
inv = inv.merge(dq[['sub', 'ses', 'qc_determination_scalar_maps']], on=['sub', 'ses'], how='left')
m = inv[(inv.MEGRE_full == 1) & (inv.T1w == 1) & (inv.t1qc.fillna('missing') != 'fail')].copy()
first = m.sort_values('sesn').groupby('sub').head(1).copy()
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.append(s)

def X_reg(d, extra):
    ps = np.where(d.PS == 1, .5, -.5); ad = np.where(d.AD == 1, .5, -.5); sex = np.where(d.sex == 'M', .5, -.5)
    rb = d.race.str.startswith('Black').astype(float).values; rw = (d.race == 'White').astype(float).values
    ro = 1 - rb - rw
    # deviation coding of race3 (Black, White, Other): two columns
    r1 = rb - ro; r2 = rw - ro
    B = np.asarray(patsy.dmatrix('cr(x, df=4) - 1', {'x': d.age.values}))[:, 1:]  # 3 df natural spline, not collinear with intercept
    cols = [np.ones(len(d)), ps, ad, ps * ad, sex, r1, r2, B[:, 0], B[:, 1], B[:, 2]] + list(extra)
    return np.column_stack(cols)

n = len(first)
ses_ = []; Xs = []
for i in range(400):
    X = X_reg(first, [rng.standard_normal(n), rng.standard_normal(n)])
    XtXi = np.linalg.inv(X.T @ X); ses_.append([np.sqrt(XtXi[1, 1]), np.sqrt(XtXi[2, 2]), np.sqrt(XtXi[3, 3])])
    Xs.append(X)
ses_ = np.array(ses_); s_ps, s_ad, s_int = ses_.mean(0)
k = Xs[0].shape[1]; dfree = n - k
X = Xs[0]
tc = stats.t.ppf(.975, dfree); t1c = stats.t.ppf(.95, dfree)
P(f'REGISTERED PRIMARY: n={n}, p={k}, df={dfree}; SE(beta_PS)={s_ps:.3f} (range over simulated roi_vol/Q_phase {ses_[:,0].min():.3f}-{ses_[:,0].max():.3f}); SE(ADHD)={s_ad:.3f}; SE(int)={s_int:.3f}')
for dlt in [.2, .3, .4, .5, .6]:
    P(f'  H1 power delta={dlt}: {1-stats.nct.cdf(tc, dfree, dlt/s_ps):.2f}')
P(f'  MDE80 = {(tc+stats.norm.ppf(.8))*s_ps:.3f}')
for ntar in [104, 100, 96]:
    s_q = s_ps*np.sqrt(n/ntar); dq_ = ntar-k
    P(f'  after QC n~{ntar}: power .6 = {1-stats.nct.cdf(stats.t.ppf(.975,dq_), dq_, .6/s_q):.2f}; B(-0.6) at 0 = {1-stats.nct.cdf(stats.t.ppf(.95,dq_), dq_, .6/s_q):.2f}; B(-0.5) at 0 = {1-stats.nct.cdf(stats.t.ppf(.95,dq_), dq_, .5/s_q):.2f}')
P('Outcome B power by bound (rows) and true delta (cols 0, -0.1, -0.15, -0.2, -0.3):')
for lab, Bd in [('-0.6', .6), ('-0.6*sqrt(.9)', .6*np.sqrt(.9)), ('-0.6*sqrt(.8)', .6*np.sqrt(.8)), ('-0.6*sqrt(.7)=-0.50', .6*np.sqrt(.7)), ('-0.5', .5)]:
    P(f'  bound {lab:22s} ' + ' / '.join(f'{1-stats.nct.cdf(t1c, dfree, (Bd-dt)/s_ps):.2f}' for dt in [0, .1, .15, .2, .3]))
c1 = np.zeros(k); c1[1] = 1; c1[3] = -.5; c2 = np.zeros(k); c2[1] = 1; c2[3] = .5
XtXi = np.linalg.inv(X.T @ X)
P(f'  SE simple effects PS-only-TD {np.sqrt(c1@XtXi@c1):.3f}; PS+ADHD-ADHD-only {np.sqrt(c2@XtXi@c2):.3f}')
for B_ in [.5, .6]:
    lo = -B_ + t1c*s_ad; hi = B_ - t1c*s_ad
    P(f'  ADHD TOST +/-{B_}: ' + ' / '.join(f'{max(0, stats.norm.cdf((hi-dt)/s_ad)-stats.norm.cdf((lo-dt)/s_ad)):.2f}' for dt in [0, -.2]))
P('Dilution (d_c=0.6): ' + '; '.join(f'p={pc}: delta={pc*.6:.2f} powerA={1-stats.nct.cdf(tc, dfree, pc*.6/s_ps):.2f} P(B)={1-stats.nct.cdf(t1c, dfree, (.6-pc*.6)/s_ps):.2f}' for pc in [.15, .2, .25, .3]))

def sim_labels(dPS_only, dPS_com, dAD, nsim=4000):
    ps = X[:, 1]; ad = X[:, 2]
    mu = np.select([(ps < 0) & (ad < 0), (ps < 0) & (ad > 0), (ps > 0) & (ad < 0), (ps > 0) & (ad > 0)], [0, dAD, dPS_only, dAD + dPS_com])
    H = XtXi @ X.T; lab = dict(A=0, PS_effect=0, PS_specific=0, comorbid=0, B=0)
    for _ in range(nsim):
        y = mu + rng.standard_normal(n); b = H @ y; r = y - X @ b; s2 = r @ r / dfree; sd = np.sqrt(s2)
        sps = np.sqrt(s2*XtXi[1, 1]); sad = np.sqrt(s2*XtXi[2, 2]); sint = np.sqrt(s2*XtXi[3, 3])
        A = b[1]/sps < -tc; e1 = b[1]-.5*b[3]; e2 = b[1]+.5*b[3]
        tost = ((b[2]+.5*sd)/sad > t1c) and ((b[2]-.5*sd)/sad < -t1c); ci0 = abs(b[3]/sint) < tc
        se1 = np.sqrt(s2*(c1@XtXi@c1)); se2 = np.sqrt(s2*(c2@XtXi@c2))
        lab['A'] += A; lab['B'] += (b[1]+.6*sd)/sps > t1c
        lab['PS_effect'] += A and e1 < 0 and e2 < 0
        lab['PS_specific'] += A and e1 < 0 and e2 < 0 and tost and ci0
        lab['comorbid'] += A and (e2/se2 < -tc) and not (e1/se1 < -tc) and e1 > -.5*abs(e2)
    return {kk: round(v/nsim, 2) for kk, v in lab.items()}
P('Claim labels:')
for sc in [(0, 0, 0), (-.4, -.4, 0), (-.6, -.6, 0), (0, -.8, 0)]:
    P('  PS-only %.2f, comorbid %.2f, ADHD %.2f ->' % sc, sim_labels(*sc))

def hc3_p(y):
    b = XtXi @ X.T @ y; r = y - X @ b; h = np.einsum('ij,jk,ik->i', X, XtXi, X)
    V = XtXi @ (X.T * (r/(1-h))**2) @ X @ XtXi; t = b[1]/np.sqrt(V[1, 1]); return 2*stats.t.sf(abs(t), dfree)
rej = [0, 0]
for i in range(4000):
    e = (rng.chisquare(4, n)-4)/np.sqrt(8)*np.where(X[:, 1] > 0, 1.3, 1.0); rej[0] += hc3_p(e) < .05
    rej[1] += hc3_p(rng.standard_normal(n)) < .05
P(f'HC3 type I: skewed+heteroscedastic={rej[0]/4000:.3f}; normal={rej[1]/4000:.3f}')

# S1: power of L and of the IUT "dissociation" label (L Holm-significant AND b_PS(D)<0 AND b_ADHD(D)>0, each one-sided .05)
P('S1 (per rho): L power at alpha .0167/.05 ; P(dissociation label) at alpha_L .0167/.05')
for rho in [.3, .5, .7]:
    sdD = np.sqrt(2*(1-rho))
    for (psS, psT, adS, adT) in [(-.4, -.1, -.1, -.4), (-.5, -.1, 0, -.4)]:
        bps = psS-psT; bad = adS-adT; L = bps-bad
        ps = X[:, 1]; ad = X[:, 2]; mu = bps*ps + bad*ad
        H = XtXi @ X.T; cL = np.zeros(k); cL[1] = 1; cL[2] = -1; vL = cL@XtXi@cL
        res = np.zeros(4); nsim = 4000
        for _ in range(nsim):
            y = mu + rng.standard_normal(n)*sdD; b = H @ y; r = y - X @ b; s2 = r @ r/dfree
            tL = (b[1]-b[2])/np.sqrt(s2*vL); tps = b[1]/np.sqrt(s2*XtXi[1, 1]); tad = b[2]/np.sqrt(s2*XtXi[2, 2])
            comp = (tps < -t1c) and (tad > t1c)
            for j, a in enumerate([.0167, .05]):
                sig = (tL < 0) and (abs(tL) > stats.t.ppf(1-a/2, dfree))
                res[j] += sig; res[2+j] += sig and comp
        P(f'  rho={rho} PS(SN,STR)=({psS},{psT}) ADHD=({adS},{adT}) L={L:.1f}: L power {res[0]/nsim:.2f}/{res[1]/nsim:.2f}; label {res[2]/nsim:.2f}/{res[3]/nsim:.2f}')

# S2
def pw_r(r, nn, kcov, a=.05):
    return stats.norm.sf(stats.norm.ppf(1-a/2) - np.arctanh(r)*np.sqrt(nn-kcov-3))
P(f'S2 (n=83, 12 covariates): power r=.25 {pw_r(.25,83,12):.2f}; r=.30 {pw_r(.30,83,12):.2f}; MDE r={np.tanh((1.96+.8416)/np.sqrt(83-12-3)):.3f}')
# S3 restricted to H1 session
dm = m[(m.DWI == 1) & (m.qc_determination_scalar_maps.fillna('').str.strip() == 'pass')]
s3 = first[first.set_index(['sub', 'ses']).index.isin(dm.set_index(['sub', 'ses']).index)]
P(f'S3 same-session sample n={len(s3)}')
kS3 = 20  # chi_DGM, b0, ns(age,3), sex, race3(2), PS, ADHD, PSxADHD, FD, NDC, log nSL, DN_R_vol, Q_phase, |pitch|, tract length (+1 slack)
for nn in [108, 101, 94, 92]:
    P(f'  S3a n={nn} ({kS3} covariates): power r=.30 at .0167 {pw_r(.3,nn,kS3,.0167):.2f}; at .05 {pw_r(.3,nn,kS3):.2f}; MDE r (alpha .05) {np.tanh((1.96+.8416)/np.sqrt(nn-kS3-3)):.3f}')

# T1: fixed-margin rule across ICC (critic design: linear age, no race/spline; GLS with known ICC)
m2 = m.sort_values(['sub', 'sesn']).copy(); m2['k'] = m2.groupby('sub').cumcount(); d = m2[m2.k <= 1].copy()
d['occ'] = d.k; d['PSe'] = np.where(d.PS == 1, .5, -.5); d['ADe'] = np.where(d.AD == 1, .5, -.5); d['sx'] = np.where(d.sex == 'M', .5, -.5)
d['a'] = (d.age-d.age.mean())/2
XT = np.column_stack([np.ones(len(d)), d.PSe, d.ADe, d.PSe*d.ADe, d.occ, d.PSe*d.occ, d.a, d.sx]); nT = len(d)
g = pd.factorize(d['sub'])[0]; zc = stats.norm.ppf(.975); z1 = stats.norm.ppf(.95)
P(f'T1 design: {nT} rows, {g.max()+1} subjects, {int(d.occ.sum())} returners. Fixed margin Delta=0.4 total SD; label = reject H0 b(PSxocc)>=Delta (one-sided .05) AND follow-up PS estimate <0, conditional on A')
for icc in [.5, .6, .7, .8]:
    V = np.eye(nT)*(1-icc) + icc*(g[:, None] == g[None, :]); Vi = np.linalg.inv(V); Ci = np.linalg.inv(XT.T@Vi@XT); H = Ci@XT.T@Vi
    se_ps = np.sqrt(Ci[1, 1]); se_int = np.sqrt(Ci[5, 5]); Lc = np.linalg.cholesky(V)
    row = []
    for d0, d1 in [(-.5, -.5), (-.4, -.4), (-.5, 0), (-.4, 0)]:
        mu = np.where(d.occ == 0, d0, d1)*(d.PSe.values > 0); A = ok = 0
        for _ in range(8000):
            y = mu + Lc@rng.standard_normal(nT); b = H@y
            if not (b[1]/se_ps < -zc): continue
            A += 1; ok += ((b[5]-.4)/se_int < -z1) and (b[1]+b[5] < 0)
        row.append(ok/max(A, 1))
    P(f'  ICC {icc}: SE(PSxocc)={se_int:.3f}; P(label|A): 0.5 persists {row[0]:.2f}, 0.4 persists {row[1]:.2f}, 0.5 vanishes {row[2]:.2f}, 0.4 vanishes {row[3]:.2f}')
open(OUT, 'w').write('\n'.join(out))
