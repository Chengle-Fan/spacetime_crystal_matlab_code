# Validation record

The toolbox was checked at three levels before release.

## Static syntax

Every `.m` file was parsed with the MATLAB Tree-sitter grammar. No `ERROR` or missing syntax nodes were reported.

## Independent numerical cross-checks

The governing equations were independently implemented with NumPy/SciPy rather than translated from the MATLAB functions.

| Check | Result |
| --- | --- |
| Sampled versus analytic Park-Min Fourier coefficients | maximum error \(7.81\times10^{-3}\) on \(m=-1:1,n=-3:3\); the residual is endpoint sampling of a discontinuity |
| Zero-duration intermediate temporal layer versus one direct interface | Frobenius error \(3.93\times10^{-17}\) |
| Periodic transparent temporal stack at \(f_0\) | \((|E^+|,|E^-|)=(1.0000,1.34\times10^{-15})\) |
| Amplifying quarter-wave temporal stack at \(f_0\) | \((|E^+|,|E^-|)=(48.0052,47.9948)\) |
| Binary time-crystal PWE versus exact TMM, \(M=19\) | median normalized error \(1.90\times10^{-3}\), maximum \(2.28\times10^{-2}\) across the test sweep |
| Temporal-interface FDTD versus exact Morgenthaler amplitudes | forward error \(1.12\times10^{-3}\), backward error \(1.07\times10^{-3}\) |
| FHS Rice-Mele pump | \(C=1.0000000000000002\) |
| Uniform periodic FDTD short-run energy drift | \(2.06\times10^{-6}\) |

The remaining PWE-TMM difference is concentrated near the discontinuous temporal waveform's band edge and decreases with the temporal Fourier cutoff.

## User-side MATLAB check

A MATLAB executable was not present in the build environment. Run:

```matlab
stm_init
stm_selftest
```

before starting a new research sweep. The self-test repeats the coefficient, matrix-composition, Chern-number, and uniform-FDTD checks in MATLAB.
