from setuptools import setup

setup(name='dendroplot',
    version='0.1.0',
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
    download_url = 'https://github.com/tonywong94/dendroplot/archive/refs/tags/v0.1.0.tar.gz', 
    packages=['dendroplot'],
    )
