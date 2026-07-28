# Validation record — v2

The package was checked at three levels before packaging.

## Static syntax

All 47 `.m` files were parsed with the MATLAB Tree-sitter grammar. No `ERROR` or missing syntax nodes were reported.

## Independent numerical cross-checks

The governing equations were independently implemented with NumPy/SciPy rather than translated from the MATLAB functions.

| Check | Result |
| --- | --- |
| Sampled versus analytic Fig. 2 Fourier coefficients | maximum error $7.81\times10^{-3}$ on $m=-1:1,n=-3:3$; the residual is endpoint sampling of a discontinuity |
| Zero-duration intermediate temporal layer versus one direct interface | Frobenius error $3.93\times10^{-17}$ |
| Periodic transparent temporal stack at $f_0$ | $(|E^+|,|E^-|)=(1.0000,1.34\times10^{-15})$ |
| Amplifying quarter-wave temporal stack at $f_0$ | $(|E^+|,|E^-|)=(48.0052,47.9948)$ |
| Binary time-crystal PWE versus exact TMM, $M=19$ | median normalized error $1.90\times10^{-3}$, maximum $2.28\times10^{-2}$ across the test sweep |
| Temporal-interface FDTD versus exact Morgenthaler amplitudes | forward error $1.12\times10^{-3}$, backward error $1.07\times10^{-3}$ |
| FHS Rice-Mele pump | $C=1.0000000000000002$ |
| Uniform periodic FDTD short-run energy drift | $2.06\times10^{-6}$ |
| D/B jump-law preservation | $0$ in an independent complex-amplitude test |
| E/B jump-law preservation | $5.55\times10^{-17}$ |
| Two-input coherent cancellation | residual output $0$ |
| Finite-crystal helper versus explicit $U^5$ | state error $4.71\times10^{-16}$ |
| AB versus BA infinite-crystal spectrum | maximum trace difference $4.44\times10^{-16}$ |
| Reversed-cell temporal domain wall | $k/(2\pi/T)=0.595311$; refined eigenspace mismatch $4.25\times10^{-8}$; envelope maximum at the interface |

The remaining PWE-TMM difference is concentrated near the discontinuous temporal waveform's band edge and decreases with the temporal Fourier cutoff.

## User-side MATLAB check

A MATLAB executable was not present in the build environment. Run:

```matlab
startup_stm
test_smoke
```

before starting a new research sweep. The smoke test repeats the coefficient,
matrix-composition, jump-law, finite-crystal, domain-wall, Chern-number,
CFL/probe, and uniform-FDTD checks in MATLAB.

The new domain-wall test establishes a transfer-matrix localized state. It
does not substitute for the temporal Zak Wilson loop listed as the first
research extension in `RESEARCH_ROADMAP_CN.md`.
