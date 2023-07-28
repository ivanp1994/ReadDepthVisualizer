# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 13:49:13 2023

@author: ivanp
"""

from setuptools import setup

with open("README.md", "r") as _f:
    long_description = _f.read()
    
setup(name="RDViz",
      version="0.0.2",
      description="GUI for visualization of RD signal from CNVPytor HDF files",
      package_dir = {"RDViz":"ReadDepthViz"},
      package_data = {"RDViz":["data/*","data/pytor/*"]},
      long_description=long_description,
      long_description_content_type='text/markdown',
      url="https://github.com/ivanp1994/ReadDepthVisualizer.git",
      author = "Ivan Pokrovac",
      author_email = "ivan.pokrovac.fbf@gmail.com",
      licence = "MIT",
      classifiers=["Development Status :: 3 - Alpha",
                   "License :: OSI Approved :: MIT License",
                   "Programming Language :: Python :: 3.8",
                   "Intended Audience :: Science/Research"],
      install_requires = ["requests >= 2.27.1",
                          "h5py >= 3.4.0",
                          "numpy >= 1.22.3",
                          "pandas >= 1.4.1",
                          "matplotlib >= 3.4.2",
                          ],
      python_requires=">=3.8.12"
      )
