# Time-Limited Balanced Truncation Model Order Reduction for Descriptor Systems

This repository contains demo files for time-limited balanced truncation for three structured descriptor system types. We consider semi-explicit systems of index 1, Stokes-like systems of index 2, and mechanical multibody systems of index 3.

## Abstract

Balanced truncation is a well-known model order reduction technique for large-scale systems. In recent years, time-limited balanced truncation, which restricts the system Gramians to finite time intervals, has been investigated for different system types. We extend the ideas to linear time-invariant continuous-time descriptor systems using the framework of projected generalized Lyapunov equations. The formulation of the resulting Lyapunov equations is challenging since the right-hand sides are unknown a priori. We propose Krylov subspace methods for the efficient computation of the right-hand sides for different system structures. Since the right-hand sides may become indefinite, we use an LDLT-factorization-based low-rank ADI iteration to solve the Lyapunov equations and obtain the system balancing transformation. Comparing the time-limited to the classical approach in numerical experiments, we observe a steeper decay of the Hankel singular values. This behavior renders useful, especially when employing low-rank approximation techniques. Further, we show that time-limited balanced truncation can deliver reduced models of smaller order with similar accuracy in the prescribed time domain.

## Dependencies

Requires M-M.E.S.S. toolbox version 2.0 - [available here](https://zenodo.org/record/3368844)

Includes modified B(FOM)^2 package - [available here](https://gitlab.com/katlund/bfomfom-main)

## Installation

Please copy the contents of this repository to the DEMO folder of the M-M.E.S.S toolbox. Then run `mess_path` to add the files to your Matlab path.

## Getting started

Navigate to either one of the subfolders (`TLBT_DAE1_BIPS`, `TLBT_DAE2_Stokes`, `TLBT_DAE3_Stykel`) and consult the respective `runme` files.

## License

The software uses a BSD 2-Clause license. See [LICENSE.md](LICENSE.md) for details.

## Citation

Please only cite this work if you are referring to time-limited balanced truncation for descriptor systems. Otherwise, consider the references given in the files and the references therein.

```
@THESIS{Base2021,
  author = {Bäse, Simon Michael},
  title = {Time-Limited Balanced Truncation Model Order Reduction for Descriptor Systems},
  year = {2021},
  institution = {Institut für Mathematik},
  location = {Technische Universität Berlin},
  type = {Master thesis}
}
```

## Feedback and Support

Please feel free to open an issue if you find a bug or seek support. Also, I encourage you to fork the project.

[Email](mailto:simonbaese@mailbox.tu-berlin.de)


## Release Notes

Version 1.0 - Initial Release
