from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='dendroplot',
    version='0.6.0',
    description='Helper tasks for astrodendro',
    long_description=long_description,
    long_description_content_type="text/markdown",
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
    download_url = 'https://github.com/tonywong94/dendroplot/archive/refs/tags/v0.5.0.tar.gz', 
    packages=['dendroplot', 'dendroplot.analysis', 'dendroplot.lte', 
              'dendroplot.plotting'],
    )
