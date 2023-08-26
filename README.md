# Synthetic-Turbulence-Generator

This repo provides code for our paper [An efficient and low-divergence method for generating inhomogeneous and anisotropic turbulence with arbitrary spectra](https://doi.org/10.1017/jfm.2023.548) ([arXiv version](https://doi.org/10.48550/arXiv.2301.05363)), implemented in the OpenFOAM framework.
* Authors: Hao Guo \[[Google Scholar](https://scholar.google.com/citations?hl=zh-CN&user=ZOhz0b8AAAAJ)\], Peixue Jiang, Lin Ye, Yinhai Zhu

<p align="center">
  <img align="center" width="800" src="/docs/effectiveness.png">
</p>
<p align="center" > Inverter-version method effectiveness in reducing divergence. </p>

<p align="center">
  <img align="center" width="800" src="/docs/performance.png">
</p>
<p align="center"> Real performance of inverter spectrum-based method against white noise maker. </p>

## Abstract:

In this article, we propose a divergence-free method for the generation of inhomogeneous and anisotropic turbulence. Based on the idea of correlation reconstruction, the method uses the Cholesky decomposition matrix to re-establish the turbulence correlation functions, which avoids the time-consuming procedure that solves eigenvalues and eigenvectors in every location needed by the coordinate transformation in the conventional method and thus reduces the computational complexity and improves the efficiency of generating synthetic turbulence. Through adjusting the generation strategy of specific random vectors, the proposed method, which is based on the classical spectrum-based method widely used to generate uniform isotropic turbulence, can obtain inhomogeneous and anisotropic turbulence with a relatively low divergence level in practice with almost no additional computational burden. There are two versions of this new method: the shifter version and the inverter version. Both versions of the method are highly efficient, easy to implement, and compatible with high-performance computing. Suitable for providing high-quality initial or boundary conditions for scale-resolving turbulence simulations with large grid numbers (such as direct numerical simulation or large eddy simulation), this method can be quickly implemented either based on various open-source CFD codes or common commercial CFD software.


## Code structure:

* `src/` contains the source code of different generators.
  * `whiteNoiseMaker` contains the code of a simple white-noise generator for comparison.
  * `spectralInverter` contains the code of inverter-version spectrum-based generator.
  * `auxiliaries` contains some auxiliary functions and classes used by spectralInverter.
    * `spectrumModel` includes a set of spectrum models.
    * `auxiliaryFunctions` includes random unit vector generation and Cholesky decomposition.
* `examples/` contains a set of examples of using the framework.
  * `boxTurbAnisotropy` includes a classical box turbulence of different anisotropy types.


## Compilation:

Run `Allwmake`.


## Usage:

Use `gitVerSpetralInverter -help`.


# Examples:

1. boxTurbAnisotropy: run the `buildInitialField` script.


## Requirements:

OpenFOAM v2106


## Citation

If you find this repo useful in your research, please consider citing our paper: [An efficient and low-divergence method for generating inhomogeneous and anisotropic turbulence with arbitrary spectra](https://doi.org/10.1017/jfm.2023.548).

```
@article{Guo2023,
    title={An efficient and low-divergence method for generating inhomogeneous and anisotropic turbulence with arbitrary spectra},
    author={Guo, Hao and Jiang, Peixue and Ye, Lin and Zhu, Yinhai},
    year={2023},
    journal={Journal of Fluid Mechanics},
    volume={970},
    pages={A2},
    DOI={10.1017/jfm.2023.548},
}
```


## License

Synthetic-Turbulence-Generator is published under the GNU GPL Version 3 license, which can be found in the LICENSE file.


## Problems

If you find any bugs in the code or have trouble in compiling/running syntheticTurbulenceGenerator in your machine, you are very welcome to [open an issue](https://github.com/Fracturist/syntheticTurbulenceGenerator/issues) in this repository.


## Future plan

A Fluent UDF version is in progress.

