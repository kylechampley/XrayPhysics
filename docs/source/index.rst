.. XrayPhysics documentation master file, created by
   sphinx-quickstart on Tue Feb 20 18:59:53 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to XrayPhysics's documentation!
=======================================

This is a C/C++ library with Python bindings with the following capabilities:

1. x-ray cross sections (energies from 1 keV to 20 MeV and elements 1-100) specified by chemical formula or element/compound mass fractions
2. incoherent and coherent x-ray scattering angle distributions specified by chemical formula or element/compound mass fractions
3. x-ray source spectra modeling (for any source voltage and take-off angle, and the following anode types: Cu, Mo, W, Au)
4. one- and two-material beam hardening correction algorithms (theoretically exact and polynomial-based)
5. Calculation of the LLNL-defined effective atomic number; For more information see here
6. dual energy decomposition algorithms (and SIRZ if used in conjuction with LEAP)

There are very few dependencies and we don't use any specialized data structures, making it very easy to incorporate with other software packages. The cross section tables are hard-coded into C++ arrays, so queries of the database are instant.

The x-ray source models are at least as accurate as TASMICS and much more flexible. The cross section tables are based on `EPDL97 <https://www-nds.iaea.org/epdl97/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   xrayphysics


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
