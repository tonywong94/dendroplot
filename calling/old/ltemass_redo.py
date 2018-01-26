from add_ltemass import add_ltemass

ncube = '../lte/30Dor_cube_n13cube.fits.gz'
ncube_uc = '../lte/30Dor_cube_n13cubeerr.fits.gz'
lines = ['12', '13']

for l in lines:
    label = '30dor_' + l
    add_ltemass(label, n13cub = ncube, n13cub_uc = ncube_uc)
