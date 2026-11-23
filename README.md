# Skyrmions3D


<!--[![Build](url_build_status_badge)](url_build_status)-->
<!--[![Coverage](url_coverage_status_badge)](url_coverage_status)-->

A package to create and visualise 3D Skyrmions in the Skyrme model.

## Installation

You can install `Skyrmions3D.jl` using the Julia REPL.

```julia
julia> ]

pkg> add https://github.com/chrishalcrow/Skyrmions3D.jl.git

julia> using Skyrmions3D
```

A more comprehensive installation guide can be found on the [documentation webpage](https://chrishalcrow.github.io/Skyrmions3D.jl/dev/#Installation).
You can test whether the installation was successful by creating a Skyrmion and computing some invariants.

```julia
julia> skyrmion = Skyrmion(40, 0.2); p(z) = z^2; q(z) = 1;

julia> make_rational_map!(skyrmion, p, q; baryon=2);

julia> Baryon(skyrmion)
1.9910413040360426

julia> Energy(skyrmion)
2.496764502522567
```

The documentation also contains more detailed information on how to [create](https://chrishalcrow.github.io/Skyrmions3D.jl/dev/making_skyrmions/) and [visualise](https://chrishalcrow.github.io/Skyrmions3D.jl/dev/visualisation/) a Skyrmion.

## Tutorials and Documentation

A beginners tutorial written in Jupyter notebook can be downloaded from this repo.

The documentation can be found [here](https://chrishalcrow.github.io/Skyrmions3D.jl/dev/).

There are [one](https://youtu.be/HeSs7yVGXR4) or [two](https://youtu.be/TI5huk6Rqos) YouTube videos about the package.

## Features
- Create skyrmions using the Rational Map and ADHM approximations.
- Calculate skyrmion properties such as energy, baryon number, moments of inertia and more.
- Plot skyrmions.
- Transform and combine skyrmions.
- Flow skyrmions to find their energy minimisers or model scattering paths.
