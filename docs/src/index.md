# About the Package

[`SmoQyDEAC.jl`](https://github.com/SmoQySuite/SmoQyDEAC.jl.git) utilizes the Differential Evolution for Analytic Continuation algorithm developed by Nathan S. Nichols, Paul Sokol, and Adrian Del Maestro [arXiv:2201.04155](https://arxiv.org/abs/2201.04155).

This package takes imaginary time or Matsubara frequency correlation functions from condensed matter Monte Carlo simulations and provides the associated spectral function on the real axis. 

# Installation

**NOTE**: This package is in the experimental phase of development and is not yet published to the Julia [`General`](https://github.com/JuliaRegistries/General.git) registry.

Open a Julia REPL environment and run the following commands:
```
julia> ]
pkg> add SmoQyDEAC
```

# Running SmoQyDEAC

SmoQyDEAC has a simple API interface with a two callable functions `DEAC_Binned` and `DEAC_Std` for binned data and data with the standard error, respectively.

## API
- [`DEAC_Binned`](@ref)
- [`DEAC_Std`](@ref)

```@docs
DEAC_Binned
```
```@docs
DEAC_Std
```

To add additional mutations to the base DEAC algorithm you utilize the ```user_mutation!``` functionality. ```user_mutation!``` is a user defined function which takes three (3) parameters, ```(population_new, population_old, rng)```. The result in ```population_new``` is treated as a new trial population, and if the fitness improves for a population, the trial population is stored.

```population_new``` and ```population_old``` are of the shape ```[1:nω, 1:population_size]```. Note, there is no guarantee that ```population_new``` is a copy of ```population_old```, so the user should update the entire array. See `Example 3: User Mutation` for an example implementation.

## Output

### Default Keys
Both API functions return a ```Dict{String,Any}``` object as well as save that dictionary to the location specified in the parameter ```output_file```. The dictionary has the following keys
- `A`: 2D array of shape (n$\omega$,nFitness), where nFitness is the number of fitness targets in the run. 
- `fitness`: 1D array with all fitnesses associated with the run in descending order
- `σ`: The calculated standard error for the run. NOTE: DEAC is non-ergodic, so this does not correspond actual error bars!
- `ωs`: The $\omega$ values corresponding to the first dimenson of `A`.
- `avg_generations`: number of generations to reach lowest target fitness
- `runtime`: Total run time in seconds
- `zeroth_moment`: 1D array of the calculated zeroth moments for the reported range of `$\omega$s` for each fitness
- `zeroth_moment_σ`: Standard error for `zeroth_moment`.
- `full_eigenvalues`: Only pertinent for binned data. When finding the eigenbasis for the covariance matrix near-zero eigenvalues may arise due to linear correlations in your data. Those below the value eigenvalue_ratio_min*maximum(eigenvalues) will be ignored and their eigenvectors not used.

### Bin data
If `keep_bin_data==true` then the following keys are also in the output dictionary
- `bin_data`: 3D array of size (n$\omega$,nBins,nFitness) with data from each bin.
- `bin_zeroth_moment`: Same as `zeroth_moment` but on a bin by bin basis.

## Kernels
The following are the supported kernels
- `time_fermionic`$=\frac{e^{-\tau\omega}}{1+e^{-\beta\omega}}$
- `time_bosonic`$=\frac{e^{-\tau\omega}}{1-e^{-\beta\omega}}$
- `time_bosonic_symmetric`$=\frac{1}{2}\frac{e^{-\tau\omega}+e^{-(\beta-\tau)\omega}}{1-e^{-\beta\omega}}$
- `time_bosonic_symmetric_w`$=\frac{\omega}{2}\frac{e^{-\tau\omega}+e^{-(\beta-\tau)\omega}}{1-e^{-\beta\omega}}$
- `frequency_fermionic`$=\frac{1}{i\omega_n-\omega}$
- `frequency_bosonic`$=\frac{1}{i\omega_n-\omega}$
- `frequency_bosonic_symmetric`$=\frac{2 \omega}{\omega_n^2+\omega^2}$

# Multithreading
SmoQyDEAC utilizes Julia's `Threads.@threads` multithreading capability. To take advantage of this run your code using 
```$ julia --threads=auto yourscript.jl```. 
`auto` will automatically use any available cores or hyperthreads. You can set the value to a fixed number as you wish.

## Tips, tricks and caveats

- You will likely never need to adjust any of the Optional Algorithm Arguments from the default. Many are merely initial values and will be updated stochastically early and often in the code.
- The DEAC algorithm can have edge effects where it places spectral weight on the first or last ω point. This occurs when there is spectral weight just outside of your range of ωs. The solution is simply expanding the range of your output energies.
- For bosonic correlations SmoQyDEAC returns the spectral function, e.g. $\mathcal{A}(\omega)$ not $\frac{\mathcal{A}(\omega)}{\omega}$ as some MaxEnt codes do.
- Different simulation codes may report correlation functions slightly differently. E.g. for [`SmoQyDQMC`](https://github.com/SmoQySuite/SmoQyDQMC.jl) `phonon_greens` $=\langle X(\tau)X(0)\rangle$ not the actual phonon green's function of $-2\Omega_0\langle X(\tau)X(0)\rangle$. While the negative sign will cancel out by our choice of Kernel, you may need to postprocess the spectral function you recover. In this case $\mathcal{A}\rightarrow \dfrac{\mathcal{A}}{2\Omega_0}$. 

# Funding
This work was supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award Number DE-SC0022311. N.S.N. was supported by the Argonne Leadership Computing Facility, which is a U.S. Department of Energy Office of Science User Facility operated under contract DE-AC02-06CH11357. 

# Citation
If you found this library to be useful in the course of academic work, please consider citing us:

```bibtex
@Article{10.21468/SciPostPhysCodeb.39,
	title={{SmoQyDEAC.jl: A differential evolution package for the analytic continuation of imaginary time correlation functions}},
	author={James Neuhaus and Nathan S. Nichols and Debshikha Banerjee and Benjamin Cohen-Stead and Thomas A. Maier and Adrian Del Maestro and Steven Johnston},
	journal={SciPost Phys. Codebases},
	pages={39},
	year={2024},
	publisher={SciPost},
	doi={10.21468/SciPostPhysCodeb.39},
	url={https://scipost.org/10.21468/SciPostPhysCodeb.39},
}
```

## Contact Us

For question and comments regarding this package, please email either James Neuhaus at [jneuhau1@utk.edu](mailto:jneuhau1@utk.edu) or Professor Steven Johnston at [sjohn145@utk.edu](mailto:sjohn145@utk.edu).

# Publication List
```@bibliography
*
```
