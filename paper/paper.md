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
    orcird: 0009-0001-4281-2766
    affiliation: 2
affiliations:
 - name: University of Edinburgh, United Kingdom
   index: 1
 - name: University of Leeds, United Kingdom
   index: 2
date: 22 August 2025
bibliography: paper.bib
---

# Summary
*A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.*


Solitons are solutions to partial differential equations which behave like particles in certain ways: they have finite extent, a well defined centre of mass and momentum, and they retain their shape after scattering off each other. Skyrmions were introduced in the 1960s as topological solitons which serve as low-energy models of nuclei. Mathematically they are defined as smooth maps $U: \mathbb{R}^3 \to \mathrm{SU}(2)$ which minimise a certain integral (the Skyrme energy) subject to boundary conditions which topologically stabilise the field $U$ so as to prevent it only occupying the vacuum solution $U=1$. Analytically writing down skyrmions fields is difficult and accordingly many approximation of skyrmions fields, using for example rational maps or the ADHM data of instantons, have been proposed. 

`Skyrmions3D.jl` allows for the creation, transformation, and visualisation of skyrmions generated via various approximations already provided in the package. Moreover, the implementation is written flexibly so as to allow users to investigate their own novel skyrmion approximations. Underlying algorithms are state of the art, and written so as to exploit efficiencies of Julia. 

# Statement of need
*A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work.*


Computational tools have long been important to the study of topological solitons due to the many complicated high-dimensional integrals which must be executed in order to determine even fundamental properties of a soliton such as its energy. This difficulty is compounded in the study of skyrmions, where exact Skyrme fields are not available. Despite such a necessity, aside from `Skyrmions3D.jl` there are no open-source packages available for the study of 3-dimensional skyrmions. This is a marked difference to the study of other solitons where a variety of tools exist, for example 2-dimensional magnetic skyrmions (e.g. [@Beg2022, @CortesOrtuno019, @KanaszNagy2015]), monopoles (e.g. [@DisneyHogg2023, @Garcia2025, @Lang2020]), and vortices (e.g. [@GonzalezArroyo2004, @GonzalezArroyo2007, @Stoica2013]). 

As such `Skyrmions3D.jl` provides a tool which lowers barriers to skyrmion research presented by the computational skill required, prevents inefficies arising from duplication of software, supports rigorous science by making source code open, and aids reproducibility across publications. 

# Background
*A bit of background about skyrmions, their mathematics, and existing work.*
*A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.*
*Note that the paper begins by a metadata section (the enclosing — lines are mandatory) and ends with a References heading, and the references are built automatically from the content in the .bib file. You should enter in-text citations in the paper body following correct Markdown citation syntax. Also note that the references include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.*
*Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.*


The Skyrme model is a non-linear sigma model with field $U: \mathbb{R}^3 \to \mathrm{SU}(2)$, designed such that the topolocial solitons serve as models of baryons: an element in $SU(2)$ may be considered as a unit quaternion, and the coefficients of its imaginary parts are taken to be pion field. Upon imposing a boundary condition that $\lim_{|x| \to \infty} U(x) = 1$, $U$ determines a map $S^3 \to S^3$ and the integer which determines the homotopy class of that map in $\pi_3(S^3) \cong \mathbb{Z}$ is to be interpreted as the baryon number. The skyrmions are the fields which minimise the Skyrme energy 
$$
E = \int_{\mathbb{R}^3} \left ( \left \lvert U^{-1} dU \right \rvert^2 + \left \lvert U^{-1} dU \wedge U^{-1} dU \right\rvert^2 \right) d^3 x
$$
in each homotopy class. It is known that for a skyrmion of baryon number $B$, the energy is bounded below by $12 \pi^2 B$, though this bound is not attained [@Manton1987]. 

The remarkable fact about the Skyrme model is that, despite its simplicity, the scattering and energy levels of skyrmions within the model accurately predict the observed properties of baryons when scales $F_\pi$, $e$, and $m_\pi$ (corresponding to the pion decay constant, Skyrme parameter, and tree-level pion mass repsectivey, determined via experiment [@Adkins1983]) are introduced. As a result, there is hope that studying the skyrme model (and simple modifications thereof) can shed light into the theory of nuclei. For a comprehensive background on skyrmions and the surrounding literature see [@Manton2004, @Manton2022]. 

`Skyrmions3D.jl` implements a structure in Julia to describe numerically a skyrmion. It has
 - the discrete grid $(x_i, y_j, z_k)$ of spatial values at which the Skyrme field will be given, 
 - the pion field value at the points in the spatial grid, 
 - the physical parameters $m_\pi$, $F_\pi$, $e$, and
 - boundary conditions determing how the Skyrme field is to be treated at the edge of the grid. 

Explicit formula for skyrmion fields which attain the minimal energy are not known, and so one must work with well-motivated approximations. Two common approaches to skyrmions arising from the study of other topological solitons are the rational map approximation (motivated by monopoles) and the Atiyah-Manton or ADHM approximation (motivated by instantons): both are implemented in `Skyrmions3D.jl`. In addition, `Skyrmions3D.jl` has been written in a flexible manner such that it is simple to implement new approximation within the existing framework, see for example [@Cork2025]. 

Given skyrmions there are a variety of ways to manipulate them:
 - one can translate and (iso)rotate them, for example by sending $U(x)$ to $U(x-x_0)$ when translating by fixed $x_0 \in \mathbb{R}^3$, and
 - one can combine two skyrmions via the product ansatz.
Moreover, there are a variety of properties of skyrmions which one may naturally wish to compute, such as the total energy, or visualisations of the skyrmion baryon density. `Skyrmions3D.jl` has the ability to compute a large number of such properties which are commonly used, and supports interactive plotting via `Makie`. *Furthermore, the ability to export skyrmion fields for plotting with other software is supported*. 

Comprehensive documentation for `Skyrmions3D.jl` is provided via a webpage, including an API, examples of how to use key features, and guidance on requesting features or raising bug reports. In addition, this webpage provides a list of known publications which have used and cited `Skyrmions3D.jl`. 

While many of the underlying numerical methods are standard applications, there are notable exceptions which warrant special attention. The method used to approximate the holonomy of the instanton gauge field necessary for the Atiyah-Manton approximation is the only known implementation of [@Harland2023], developed specifically for the context of skyrmions, but having wider applicability. Moreover, the ODE solved to flow an approximate skyrmion towards the (locally) minimal energy configuration is "arrested Newton flow", a modification of gradient flow common to the field of topological solitons but less widely known, see [@Battye2002, @Gudnason2020]. Finally, the colouring used in plotting of skyrmion baryon density is not merely used to provide contrast, but colours the surface using the value of the Skyrme field at that point in space via the Runge colour sphere as introduced in [@Manton2012]. 

Future work on the package shall focus on allowing modification to the Skyrme action which make binding energies more realist, see for example [@Gudnason2020]. 

# Acknowledgements
*Acknowledgement of any financial support.*


The research of LDH is supported by a UK Engineering and Physical Sciences Research Council (EPSRC) doctoral prize fellowship.

# References
