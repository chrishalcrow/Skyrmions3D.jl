---
title: 'Skyrmions3D: A Julia package to create and visualise 3D Skyrmions in the Skyrme model'
tags:
  - Julia
  - skyrmions
  - topological solitons
authors:
  - name: Chris Halcrow
    corresponding: true
    orcid: 0000-0002-0246-8075
    affiliation: 1
  - name: Linden Disney-Hogg
    orcid: 0000-0002-6597-2463
    affiliation: 2
  - name: Cas G. Chaudhuri
    orcid: 0009-0001-4281-2766
    affiliation: 2
affiliations:
 - name: University of Edinburgh, United Kingdom
   index: 1
 - name: University of Leeds, United Kingdom
   index: 2
date: 19 January 2026
bibliography: paper.bib
---

# Summary


Solitons are solutions to partial differential equations which behave like particles in certain ways: they have finite extent, a well defined centre of mass and momentum, and they retain their shape after scattering off of each other. Skyrmions are soliton solutions of the Skyrme model. They were introduced in the 1960s to model nuclei at low energies. Mathematically they are defined as smooth maps $U: \mathbb{R}^3 \to \mathrm{SU}(2)$ which minimise a certain integral (the Skyrme energy) subject to boundary conditions which topologically stabilise the field $U$, so as to prevent it decyaing to the vacuum solution $U=1$. Analytic solutions of the Skyrme model are not known, and even writing approximate solutions is difficult. Accordingly many approximation schemes, using for example rational maps or the ADHM data of instantons to construct Skyrme fields, have been proposed. 

`Skyrmions3D.jl` allows for the creation, transformation, and visualisation of skyrmions. The package allows the user to evolve the fields using different flows, such as a gradient flow used to find minimal energy solutions. Several methods of skyrmion generation have been implemented and these are written written flexibly, allowing users to investigate their own novel skyrmion approximations. The underlying algorithms are state of the art, and written to exploit efficiencies of Julia. 

# Statement of need

Computational tools have long been important to the study of topological solitons due to the many complicated high-dimensional integrals which must be executed in order to determine even fundamental properties of a soliton such as its energy. This difficulty is compounded in the study of skyrmions, where exact Skyrme fields are not available. Flowing these fields involves simulating a nonlinear partial differential equation, which is computationally expensive. Despite the necessity of these computation tools, aside from `Skyrmions3D.jl` there are no open-source packages available for the study of 3-dimensional skyrmions. This is a marked difference to the study of other solitons where a variety of tools exist, for example 2-dimensional magnetic skyrmions (e.g. [@Beg2022], [@CortesOrtuno2019], [@KanaszNagy2015]), monopoles (e.g. [@DisneyHogg2023], [@Garcia2025], [@Lang2020]), and vortices (e.g. [@GonzalezArroyo2004], [@GonzalezArroyo2007], [@Stoica2013]). 

As such `Skyrmions3D.jl` provides a tool which lowers barriers to skyrmion research presented by the computational skill required, prevents inefficies arising from duplication of software, supports rigorous science by making source code open, and aids reproducibility across publications.

# Background

The Skyrme model is a non-linear sigma model with field $U: \mathbb{R}^3 \to \mathrm{SU}(2)$, designed such that the topological solitons serve as models of baryons. An element in $SU(2)$ may be considered as a unit quaternion, and the coefficients of its imaginary parts are taken to be pion fields. Upon imposing a boundary condition that $\lim_{|x| \to \infty} U(x) = 1$, $U$ determines a map $S^3 \to S^3$ and the integer which determines the homotopy class of that map in $\pi_3(S^3) \cong \mathbb{Z}$ is interpreted as the baryon number. The skyrmions are the fields which minimise the Skyrme energy 
\begin{align*}
	E = \int_{\mathbb{R}^3}  & \left \lbrace -\frac{F_\pi^2}{16}\mathrm{Tr} \left ( R_i R_i \right ) - \frac{1}{32e^2} \mathrm{Tr}\left ( [R_i, R_j] [R_i, R_j] \right)   + m_\pi^2 \mathrm{Tr}(1-U) \right \rbrace d^3 x,
\end{align*}
in each homotopy class. Here $R_i = (\partial_i U) U^{-1}$, while constants $F_\pi$, $e$, and $m_\pi$ correspond to the pion decay constant, Skyrme parameter, and tree-level pion mass repsectivey. The first two of these may be removed using energy units $\frac{F_\pi}{4e}$ and length units $\frac{2}{e F_\pi}$ (so-called Skyrme units). It is known that for a skyrmion of baryon number $B \geq 0$, the energy in skyrme units with $m_\pi=0$ is bounded below by $12 \pi^2 B$, though this bound is not attained for $B \geq 1$ [@Manton1987]. 

The remarkable fact about the Skyrme model is that, despite its simplicity, the scattering and energy levels of skyrmions within the model qualitatively predict the observed properties of baryons when $F_\pi$, $e$ and $m_\pi$ are determined via experiment [@Adkins1983]. As a result, there is hope that studying the skyrme model (and simple modifications) can shed light into the theory of nuclei. For a comprehensive background on skyrmions and the surrounding literature see [@Manton2004], [@Manton2022]. 

`Skyrmions3D.jl` implements a structure in Julia to describe numerically a skyrmion. It has
 - the discrete grid $(x_i, y_j, z_k)$ of spatial values at which the Skyrme field will be given, 
 - the pion field value at the points in the spatial grid, 
 - the physical parameters $m_\pi$, $F_\pi$, $e$, and
 - boundary conditions determing how the Skyrme field is to be treated at the edge of the grid. 

Explicit formula for skyrmion fields which attain the minimal energy are not known. Moreover, though (local) minimisers of the energy functional are given by solutions to the Euler-Lagrange PDE, attempting to solve these using standard Julia PDE implementations (such as those from SciML) with topologically non-trivial boundary conditions is not feasible. Hence one generally starts with well-motivated approximations where the toplogical constraints are imposed from the outset, then either study these approximations or use them as initial data which can be flowed to a true minimum. Two common approaches to skyrmions arising from the study of other topological solitons are the rational map approximation (motivated by monopoles) and the Atiyah-Manton or ADHM approximation (motivated by instantons): both are implemented in `Skyrmions3D.jl`. In addition, `Skyrmions3D.jl` has been written in a flexible manner such that it is simple to implement new approximations within the existing framework, see for example [@Cork2025]. 

There are a variety of ways one can manipulate a skyrmion implemented in `Skyrmions3D.jl`.

- Translate, rotate and iso-rotate a skyrmion, for example by sending $U(x)$ to $U(x-x_0)$ when translating by fixed $x_0 \in \mathbb{R}^3$.
- Combine two skyrmions via the product ansatz.
- Evolve a skyrmion using a gradient flow, deforming it into a true energy minimiser.

Moreover there are many properties of skyrmions which one may naturally wish to compute, such as the total energy. `Skyrmions3D.jl` has the ability to compute a large number of such properties which are commonly used. It also supports plotting via `Makie`. Furthermore, the ability to export skyrmion fields is supported by saving the pion field in HDF5 format, meaning it is easy to share skyrmions between collaborators.

Comprehensive documentation for `Skyrmions3D.jl` is provided via a webpage, including an API, examples of how to use key features, and guidance on requesting features or raising bug reports. The package also has a comprehensive unit-test suite.

# Software design

Julia was chosen for the implementation due to its balance between speed of execution, required especially for the computation of high-dimensional numerical integrals, and accessibility, both in terms of simplicity of coding and reproducibility of scientific results. As an example of the latter, there exist many helpful Julia packages such as `Pkg.jl` and `Documenter.jl` which ease the process of creating high-quality packages with simple reproducible installation instructions and clear documentation; both of these examples were used in the creation of `Skyrmions3D.jl`. That this decision has been successful is evidenced by the fact that papers implementing code within the `Skyrmions3D.jl` framework are already appearing [@Cork2025]. 

By having the skyrmion structure implemented retain the information of the underlying spatial grid, pion field at all these points, as well as additional parameters, a single instance of a skyrmion often takes up a large amount of memory. This trade-off was accepted in order to allow the implementation of the saving and loading procedure of skymrions, a functionality which bolsters the reproducibilty of any scientific experiments performed with `Skyrmions3D.jl`. Additional steps were taken in the design of the package to boost its stability, such as the removal of interactive plotting functionality, which commonly interfered with other package depdencies.

The underlying code was deliberately modularised to support future development, as it is known that modularity makes code easier to read, test, and refactor in later instances. 

# State of the field

As previously written in the statement of need, there are no current alternatives to `Skyrmions3D.jl`, and this gives the software a unique relevance for researchers working on three dimensional skyrmions.

Moreover, while many of the underlying numerical methods in `Skyrmions3D.jl` are standard applications, there are notable exceptions which warrant special attention due to their wider significance for the topological solitions and differential geometry communities. The method used to approximate the holonomy of the instanton gauge field necessary for the Atiyah-Manton approximation is the only known implementation of [@Harland2023], developed specifically for the context of skyrmions, but which has wider applicability. Moreover, the ODE solved to flow an approximate skyrmion towards the (locally) minimal energy configuration is "arrested Newton flow", a modification of gradient flow common to the field of topological solitons but less widely known, see [@Battye2002], [@Gudnason2020]. The colouring used in plotting of skyrmion baryon density represents the dominant Skyrme field at that point in space, via the Runge colour sphere as introduced in [@Manton2012]. 

Future work on the package will focus on responding to the needs of the community. This may involve implementing new approximations for Skyrme fields, and by allowing for modifications of the standard Skyrme model. This further emphasizes the need to reproducible and refactorable code stated earlier. 

# Research Impact Statement

`Skyrmions3D.jl` is already having a tangible impact on the research of skyrmions, with a growing list of publications citing the package tracked in the software documentation. Moreover, the package has spurred future development with researchers independent of the core contributors creating forks of the code, and with new skyrmion software integrating `Skyrmions3D.jl` for its core functionality [@Cork2025]. 

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing of this manuscript, or the preparation of supporting materials. 

# Acknowledgements

CH is supported by the UKRI Biotechnology and Biological Sciences Research Council (BBSRC) grant number BB/X01861X/1.
The research of LDH is supported by a UK Engineering and Physical Sciences Research Council (EPSRC) doctoral prize fellowship.
The research of CGC is supported by the EPSRC Doctoral Training Partnership ["Geometry and Dynamics of Topological Solitons"](https://gtr.ukri.org/projects?ref=studentship-2650914).

# References
