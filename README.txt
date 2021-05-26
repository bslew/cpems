Cosmologicl Parameter Estimation and Map Simulation


Features:
- Healpix spehere pixelization system
- CMB map simulations
- MCMC parameter estimation
- smooth particle interpolations
- points density calculations
- spherical harmonics transformations (ccSHT)
- power spectrum calculations (fftw)
- 


TODO for ubuntu 20.04 CMakeLists should be updated as 
- CGAL version should be <5 or remove cgal libraries for linking (since v5 CGAL is a header only package)
- see some comments in install-requirements-ubuntu.sh
- update proj API or -DACCEPT_USE_OF_DEPRECATED_PROJ_API_H