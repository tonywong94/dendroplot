#!/usr/bin/env python

from pltprops import pltprops

for line in ['12', '13']:
    label = '30dor_'+line
    if line == '12':
        pltprops(label, dvkms=0.5, beam=2,
            xplot=['rad_pc', 'vrms_k', 'area_pc2', 'siglum',   'mlumco',   'mlte', 'siglte'],
            yplot=['vrms_k', 'mlumco',   'mlumco', 'sigvir',     'mvir',   'mvir', 'sigvir'],
            xlims=[[-1,1.5],   [-2,2],     [-1,3],   [-1,4], [-0.5,5.5], [-0.5,5.5], [-1,4]],
            ylims=[[-2,1.5],   [-1,5],     [-1,5],   [-1,4], [-0.5,5.5], [-0.5,5.5], [-1,4]],
            pltname=[ 'rdv', 'dvflux', 'areaflux',    'bnd',   'virial', 'ltevir', 'ltebnd'])
    else:
        pltprops(label, fghz=220.399, dvkms=0.5, beam=2,
            xplot=['rad_pc',   'vrms_k', 'area_pc2', 'siglum', 'mlte', 'siglte'],
            yplot=['vrms_k',   'mlumco',   'mlumco', 'sigvir', 'mvir', 'sigvir'],
            xlims=[[-1,1.5],     [-2,2],     [-1,3],   [-1,4], [-0.5,5.5], [-1,4]],
            ylims=[[-2,1.5], [-1.5,4.5],     [-1,4],   [-1,4], [-0.5,5.5], [-1,4]],
            pltname=[ 'rdv',   'dvflux', 'areaflux',    'bnd', 'ltevir', 'ltebnd'])

