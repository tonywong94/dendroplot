from hist1d import hist1d_oneset as h1d1, hist1d_multiset as h1dm

lines = ['30dor_12', '30dor_13']
vals  = ['alpha', 'mlte', 'mlumco', 'mvir', 'siglum', 'siglte', 'sigvir', 'vrms_k']

for l in lines:
    for v in vals:
        h1d1(label = l, val = v)

h1dm(val = 'vrms_k')
