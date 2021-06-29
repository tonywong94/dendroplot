from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='dendroplot',
    version='0.3.2',
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
    download_url = 'https://github.com/tonywong94/dendroplot/archive/refs/tags/v0.3.2.tar.gz', 
    packages=['dendroplot', 'dendroplot.analysis', 'dendroplot.lte', 
              'dendroplot.plotting'],
    )
