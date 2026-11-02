# E2: youth susceptibility (ppb) implied by published lifespan models chi = a(1-exp(-b*age)) + g (ppm), ages 9/13/18.
import numpy as np
li = {'GP': (.13, .11, 0.), 'PUT': (.20, .01, -.01), 'CAU': (.06, .05, .01), 'SN': (.12, .08, -.02), 'RN': (.12, .06, -.03), 'DN': (.12, .04, -.01)}
zh = {'GP': (.086, .049, .017), 'PUT': (.078, .023, -.015), 'CAU': (.032, .039, .0054), 'SN': (.12, .028, -.0052), 'RN': (.10, .044, -.011)}
out = []
for lab, M in [('Li 2014 (LSQR, 1-83 y; HBM 35:2698, Table 2)', li), ('Zhang 2018 (STAR-QSM, whole-brain ref; NeuroImage 171:176, Table 2)', zh)]:
    out.append(lab)
    for r, (a, b, g) in M.items():
        v = [1000*(a*(1-np.exp(-b*t)) + g) for t in (9, 13, 18)]; asym = 1000*(a+g)
        out.append(f'  {r:4s} age 9/13/18: {v[0]:6.1f} / {v[1]:6.1f} / {v[2]:6.1f} ppb; asymptote {asym:6.1f} ppb')
print('\n'.join(out)); open('/path/to/Penn_LEAD/LEADchi/prereg/work/e2_lit_ranges.txt', 'w').write('\n'.join(out) + '\n')
