# About the Package

SmoQyDEAC utilizes the Differential Evolution for Analytic Continuation algorithm developed by Nathan S. Nichols, Paul Sokol, and Adrian Del Maestro [arXiv:2201.04155](https://arxiv.org/abs/2201.04155).

This package takes imaginary time or Matsubara frequency correlation functions from condensed matter Monte Carlo simulations and provides the associated spectral function on the real axis. 

# Installation

**NOTE**: This package is in the experimental phase of development and is not yet published to the Julia [`General`](https://github.com/JuliaRegistries/General.git) registry.

Open a Julia REPL environment and run the following command:
```julia
] dev https://github.com/SmoQySuite/SmoQyDEAC.jl.git
```
This command clones the [`SmoQyDEAC.jl`](https://github.com/SmoQySuite/SmoQyDEAC.jl.git) repository to the hidden directory `.julia/dev` that exists in the same directory where Julia is installed.

# Running SmoQyDEAC

SmoQyDEAC has a simple API interface with a single callable function `DEAC`.

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

```population_new``` and ```population_old``` are of the shape ```[1:nω, 1:population_size]```. Note, there is no guarantee that ```population_new``` is a copy of ```population_old```, so the user should update the entire array. See `Example 3: User Mutation` for an implementation.

## Kernels
The following are the supported kernels
- `time_fermionic`$=\frac{e^{-\tau\omega}}{1+e^{-\beta\omega}}$
- `time_bosonic`$=\frac{e^{-\tau\omega}}{1-e^{-\beta\omega}}$
- `time_bosonic_symmetric`$=\frac{e^{-\tau\omega}+e^{-(\beta-\tau)\omega}}{1-e^{-\beta\omega}}$
- `frequency_fermionic`$=\frac{1}{i\omega_n-\omega}$
- `frequency_bosonic`$=\frac{1}{i\omega_n-\omega}$
- `frequency_bosonic_symmetric`$=\frac{2 \omega}{\omega_n^2+\omega^2}$

# Multithreading
SmoQyDEAC utilizes Julia's `Threads.@threads` multithreading capability. To take advantage of this run your code using 
```$ julia --threads=auto yourscript.jl```
`auto` will automatically use any available cores or hyperthreads. You can set the value to a fixed number as you wish.

## Tips, tricks and caveats

- You will likely never need to adjust any of the Optional Algorithm Arguments from the default. Many are merely initial values and will be updated stochastically early and often in the code.
- The DEAC algorithm can have edge effects where it places spectral weight on the first or last ω point. This occurs when there is spectral weight just outside of your range of ωs. The solution is simply expanding the range of your output energies.
- For bosonic correlations SmoQyDEAC returns the spectral function, e.g. $B(\omega)$ not $\frac{B(\omega)}{\omega}$ as some MaxEnt codes do.
- Different simulation codes may report correlation functions slightly differently. E.g. for [`SmoQyDQMC`](https://github.com/SmoQySuite/SmoQyDQMC.jl) `phonon_greens` $=\langle X(\tau)X(0)\rangle$ not the actual phonon green's function of $-2\Omega_0\langle X(\tau)X(0)\rangle$. While the negative sign will cancel out by our choice of Kernel, you may need to postprocess the spectral function you recover. In this case $B\rightarrow \dfrac{B}{2\Omega_0}$. 

# Funding
This work was supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award Number DE-SC0022311. N.S.N. was supported by the Argonne Leadership Computing Facility, which is a U.S. Department of Energy Office of Science User Facility operated under contract DE-AC02-06CH11357. 
