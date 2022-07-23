#!/usr/bin/env python

import os
from dendroplot.plotting import comp_props

## Plot the property correlations for all GMCs together
## Can do all dendros or restrict to SCIMES clusters only (clustonly=True)

# Input/output directories
magmacsv = 'islands.sageco.csv'
# Input
analdir  = 'dendro/'
# Output
compdir  = 'comp_props/'

# GMCs ordered by decreasing 8um brightness
doclouds  = ['30Dor', 'N59C', 'A439', 'GMC104', 'GMC1', 'PCC']
markers   = ['o', 'v', '^', 's', '<', 'D']
clustonly = False

# Which lines
dolines = ['12','13']

# Which quantities to color by in comp_props
colorcodes = ['sp8med', '8um_avg', 'siglte', 'sigvir']

# ----------------------------------------------------------------------------------

thisdir = os.getcwd() # returns absolute path
print(thisdir)
if not os.path.isdir(compdir):
    os.makedirs(compdir)

try:
    os.chdir(compdir)
except OSError:
    print('Could not change directory to {}'.format(compdir))
else:
    if clustonly:
        use = 'clusters'
        x = 'cl_'
    else:
        use = 'all'
        x = ''
    comp_props(dolines, dotypes=colorcodes, clouds=doclouds, binned=True, 
               markers=markers, include=use, analdir=os.path.join(thisdir, analdir), 
               magmacsv=os.path.join(thisdir, magmacsv), linefit=True,
               xplot=[  'rad_pc','siglum',  'siglte', '8um_avg', '8um_avg', '8um_avg'],
               yplot=[  'vrms_k','sigvir',  'sigvir',  'sigvir',  'siglum',  'siglte'],
               xlims=[[-0.5,1.5],  [-1,4],    [-1,4],  [-0.5,3],  [-0.5,3],  [-0.5,3]],
               ylims=[[-1.5,1.2],  [-1,4],    [-1,4],    [-1,4],    [-1,4],    [-1,4]],
               pltname=[ x+'rdv', x+'bnd',x+'bndlte',x+'svirsf',x+'slumsf',x+'sltesf'],
               slope=[       0.5,       1,         1,         0,         0,         0],
               inter=[     -0.14,       0,         0,      2.00,      2.00,      2.00],
               pad=[        0.04,    0.12,      0.12,      0.03,      0.03,      0.03])

# Exit statement
os.chdir(thisdir)
