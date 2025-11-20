# Why & how

SpectralFitting.jl is a package for fitting models to spectral data, similar to [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/), [Sherpa](https://sherpa.readthedocs.io/en/latest/ciao.html) or [ISIS](https://space.mit.edu/CXC/isis/).

The rationale for this package is to provide a unanimous interface for different model libraries, and to leverage advancements in the computional methods that are available in Julia, including the rich statistics ecosystem, with [automatic-differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) and [_speed_](https://julialang.org/benchmarks/).

Longer term ambitions include
- Multi-wavelength fits
- Radiative transfer embedded into the package
- Spectral and timing fits

SpectralFitting aims to provide highly optimised and flexible fitting algorithms, along with a library of spectral models, for use in any field of Astronomy that concerns itself with spectral data.
