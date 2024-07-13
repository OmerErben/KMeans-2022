#### `setup.py`
from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanssp',
    version='0.1.0',
    author="Omer&Miki",
    author_email="omer.erben@gmail.com",
    description="A sample C-API",
    install_requires=['numpy', 'pandas'],
    packages=find_packages(),  # find_packages(where='.', exclude=())
    license='GPL-2',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            'mykmeanssp',  # the qualified name of the extension module to build
            ['kmeans.c'],  # the C file to compile into our module
        ),
    ]
)