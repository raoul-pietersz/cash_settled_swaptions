# Cash settled swaptions
A module for valuation of cash-settled swaptions

# Documentation
See https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3041228.
For an effective implementation in practice, speed is an important factor. The largest part of the computation time of the Python prototype is caused by the quadrature routine used for the numerical integration. If, instead, a C++-based component is used for the quadrature integration, we find this leads to a reduction of 98% of computational time.

# Dependencies

* mathplotlib
* numpy
* scipy

# Versions

Python 3.5
