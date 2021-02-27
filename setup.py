from setuptools import setup

setup(name='dendroplot',
    version='1.0',
    description='Helper tasks for astrodendro',
    author='Tony Wong',
    author_email = 'tonywong94@gmail.com',
    install_requires = ['numpy',
                        'scipy',
                        'astropy',
                        'radio_beam',
                        'spectral_cube',
                        'kapteyn',
                        'astrodendro'],
    url='https://github.com/tonywong94/dendroplot',
    packages=['dendroplot', 'dendroplot.analysis', 
              'dendroplot.lte', 'dendroplot.plotting'],
    )
