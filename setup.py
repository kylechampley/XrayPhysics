################################################################################
# Copyright 2022-2023 Lawrence Livermore National Security, LLC and other 
# LEAP project developers. See the LICENSE file for details.
# SPDX-License-Identifier: MIT
#
# LivermorE AI Projector for Computed Tomography (LEAP)
# setup.py for pytorch module
################################################################################
import os
import pathlib

from setuptools import setup, find_packages
from setuptools.command.install import install

from sys import platform as _platform
if _platform == "linux" or _platform == "linux2":
    lib_fname = 'build/lib/libxrayphysics.so'
    os.system(r'sh ./etc/build.sh')
elif _platform == "win32":
    lib_fname = r'win_build\bin\Release\libxrayphysics.dll'
    os.system(r'.\etc\win_build.bat')
    
    import site
    copy_text = 'copy ' + str(lib_fname) + ' ' + str(os.path.join(site.getsitepackages()[1], 'libxrayphysics.dll'))
    os.system(copy_text)
    
elif _platform == "darwin":
    lib_fname = 'build/lib/libxrayphysics.dylib'
    os.system(r'sh ./etc/build.sh')

setup(
    name='xrayphysics',
    version='1.0', 
    author='Kyle Champley', 
    author_email='champley@gmail.com', 
    description='x-ray cross section tables, x-ray tube source simulation, beam hardening correction, and dual energy decompositon', 
    keywords='x-ray-physics, spectra, cross-sections, beamhardening, bremsstrahlung, x-ray-simulation', 
    python_requires='>=3.6', 
    packages=find_packages("src"), 
    package_dir={'': 'src'},
    install_requires=['numpy'], 
    py_modules=['xrayphysics'], 
    package_data={'': [lib_fname]},
)
