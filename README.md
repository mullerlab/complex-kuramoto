# Complex Kuramoto

This repository contains code and a demonstration implementing the complex-valued approach to the Kuramoto model introduced in:

> [Muller, Minac, Nguyen. Algebraic approach to the Kuramoto model. *Physical Review E*, 2021](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.104.L022201)

> [Budzinski\*, Nguyen\*, Đoàn, Mináč, Sejnowski, Muller. Geometry unites synchrony, chimeras, and waves in nonlinear oscillator networks. *Chaos*, 2022 (\*equal contribution)](https://aip.scitation.org/doi/full/10.1063/5.0078791)

In these works, we developed an analytical approach to the Kuramoto model. The fundamental idea is to analyze a complex-valued system closely related to the formulation of the original, nonlinear Kuramoto model defined in the real numbers. We find the complex-valued system can be exactly solved and, by iterating this explicit expression over short intervals, the trajectories of the two systems can precisely match for long times. The advantage of this approach is that it can describe the full nonlinear dynamics in the Kuramoto model in terms of a linear operator (the matrix exponential that fully specifies the evolution of the system on each iteration). Further, this approach allows us to directly connect the nonlinear dynamics in an individual simulation (either a single network taken from data, or a single realization of a random graph) to the spectrum of the adjacency matrix. This, in turn, provides an opportunity to describe sychronization dynamics in this system (phase synchronization, traveling waves, and chimeras) in terms of the eigenmodes of the system (Budzinski\*, Nguyen\* et al., *Chaos*, 2022).

## Demonstration

The file `budzinski2022_figure1.m` provides a demonstration for Figure 1 in (Budzinski\*, Nguyen\* et al., *Chaos*, 2022):

<p align="center">
	<img src="https://mullerlab.ca/assets/img/complex-kuramoto/fig1_chaos2022.jpg" width="600">
</p>

## Dependencies

The functions in this repository use [CircStat](https://github.com/circstat/circstat-matlab) by Philipp Berens and [Expokit](https://www.maths.uq.edu.au/expokit).

## Installation

First, download or clone the repository:

```
git clone https://github.com/mullerlab/complex-kuramoto
```

Then navigate to the installation. Once you have all dependencies ready on the MATLAB path (see above), you will be able to run the demonstrations in this repository.

## Developers

[Roberto Budzinski](https://scholar.google.com/citations?user=6rsul4YAAAAJ&hl=en), [Gabriel Benigno](https://scholar.google.com/citations?user=BsNdLCkAAAAJ&hl=en), and [Lyle Muller](http://mullerlab.ca) (Western Academy for Advanced Research)
