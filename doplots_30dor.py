#!/usr/bin/env python

from pltprops import pltprops

for line in ['12', '13']:
    label = '30dor_'+line
    if line == '12':
        pltprops(label, dvkms=0.5, beam=2,
            xplot=['rad_pc', 'vrms_k', 'area_pc2', 'siglum',   'mlumco',   'mlte'],
            yplot=['vrms_k', 'mlumco',   'mlumco', 'sigvir',     'mvir',   'mvir'],
            xlims=[[-1,1.5],   [-2,2],     [-1,3],   [-1,4], [-0.5,5.5], [-0.5,5.5]],
            ylims=[[-2,1.5],   [-1,5],     [-1,5],   [-1,4], [-0.5,5.5], [-0.5,5.5]],
            pltname=[ 'rdv', 'dvflux', 'areaflux',    'bnd',   'virial', 'ltevir'])
    else:
        pltprops(label, fghz=220.399, dvkms=0.5, beam=2,
            xplot=['rad_pc',   'vrms_k', 'area_pc2',   'mlte'],
            yplot=['vrms_k',   'mlumco',   'mlumco',   'mvir'],
            xlims=[[-1,1.5],     [-2,2],     [-1,3], [-0.5,5.5]],
            ylims=[[-2,1.5], [-1.5,4.5],     [-1,4], [-0.5,5.5]],
            pltname=[ 'rdv',   'dvflux', 'areaflux', 'ltevir'])

