# Make, save and load Skyrmions

Skyrmions in Skyrmions3D are represented by the [`Skyrmion`](@ref) struct.
This contains all the information needed to reproduce the skyrmion: it's pion field, grid, pion mass etc.
For example, we can create a vacuum skyrmion on a `30×30×30` grid with lattice spacing `0.2`.

```jldoctest basic_init
julia> using Skyrmions3D

julia> my_skyrmion = Skyrmion(30,0.2);

julia> overview(my_skyrmion)
This skyrmion is on a 30x30x30 grid, with lattice spacing [0.2, 0.2, 0.2].
The boundary conditions are dirichlet

            m = 0.0
Baryon number = 0.0
       Energy = 0.0
   Baryon rms = 0.0

With physical constants Fpi = 180.0 and e = 4.0,
the energy and length units are 11.25 MeV and 0.5472 fm.
So physical E = 0.0 MeV.
   Baryon rms = 0.0 fm.
```
The [`overview`](@ref) function computes and prints a brief description of a skyrmion's properties.
We'll look at these properties in more detail in [Properties](@ref tut/props).

The underlying grid for the skyrmion does not need to be uniform.
```jldoctest basic_init
julia> another_skyrmion = Skyrmion([60,60,120],[0.2,0.2,0.1]);

```

We can also set some of the skyrmions properties in the constructor.
```jldoctest basic_init
julia> massive_periodic_skyrmion = Skyrmion([30,30,30], [0.2,0.2,0.2], mpi=1.0, boundary_conditions="periodic", Fpi=184, ee=4.5,);

julia> overview(massive_periodic_skyrmion)
This skyrmion is on a 30x30x30 grid, with lattice spacing [0.2, 0.2, 0.2].
The boundary conditions are periodic

            m = 1.0
Baryon number = 0.0
       Energy = 0.0
   Baryon rms = 0.0

With physical constants Fpi = 184.0 and e = 4.5,
the energy and length units are 10.22 MeV and 0.4758 fm.
So physical E = 0.0 MeV.
   Baryon rms = 0.0 fm.
```
You can find a full overview of the properties exposed by the [`Skyrmion`](@ref) constructor on its [API page](@ref Skyrmion).

```@meta
DocTestSetup = nothing
DocTestTeardown = nothing
```

By default the pion field is set equal to `(0,0,0,1)`, which is that of a vacuum skyrmion.
There are various ways to add structure to it.

## Rational Maps

```@meta
DocTestSetup = quote
    using Skyrmions3D
end
DocTestFilters = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
```

A complex rational map is defined by two complex valued polynomials; we call these the numerator `p(z)` and the denominator `q(z)`.
Given these polynomials, we can create a skyrmion from the rational map approximation [houghtonRationalMapsMonopoles1998](@cite).
For example, we can make a new skyrmion and give it the structure of the baryon number 3 tetrahedral skyrmion as follows:
```jldoctest rat_map
julia> b_3_tet_skyrmion = Skyrmion(30,0.2);

julia> p3(z) = sqrt(3)*im*z^2 - 1; q3(z) = z*(z^2 - sqrt(3)*im);

julia> make_rational_map!(b_3_tet_skyrmion, p3, q3)
I think your baryon number is 3.0. If it is not, include '; baryon=B' in your argument.
```

Note: By convention the "bang" character `!` at the end of function means it is a _modifiying function_.
So we are modifiying `b_3_tet_skyrmion`, not creating a new skyrmion.
Indeed, every non-trivial skyrmion defined using Skyrmions3D is constructed this way: first creating a vacuum skyrmion on some underlying grid, then giving it some structure through a modifying function.

The [`make_rational_map!`](@ref) function tries to estimate the degree of the rational map and warns us that it is doing so.
Note that we can explicitly set the baryon number with the keyword argument `baryon = B` if it doesn't choose the baryon number you were expecting.
[`make_rational_map!`](@ref) then tries to find a reasonable profile function for the skyrmion.
If this all worked, you should get a sensible output when you compute the energy of the skyrmion.
```jldoctest rat_map
julia> Energy(b_3_tet_skyrmion)
3.684801341143756
```
[`make_rational_map!`](@ref) also accepts custom profile functions, though when doing this the baryon number is not automatically calculated.
For example:
```jldoctest rat_map_w_prof
julia> b_4_cube_skyrmion = Skyrmion(30,0.2);

julia> p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0; q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0;

julia> f4(r) = pi*exp(-(r .^ 3) ./ 12.0); # our custom profile function

julia> make_rational_map!(b_4_cube_skyrmion, p4, q4, f4, baryon = 4);

julia> Energy(b_4_cube_skyrmion)
5.238411177507044
```

```@meta
DocTestSetup = nothing
DocTestTeardown = nothing
DocTestFilters = nothing
```

## ADHM data

```@meta
DocTestSetup = quote
    using Skyrmions3D
end
DocTestFilters = r"(\d*)\.(\d{4})\d+" => s"\1.\2***"
```

ADHM skyrmions [corkADHMSkyrmions2022](@cite) are skyrmions generated from ADHM data.
The data consists of symmetric quaternionic matrices which satisfy a constraint.
Most highly symmetric skyrmions can be represented by ADHM data.
This package implements the very efficient parallel transport algorithm written by [Harland2023](@citet).
We can make the baryon number 2 toroidal skyrmion as follows:
```jldoctest adhm_data; output = false
using Quaternions

B=2
adhm_data = [ Quaternion(0.0,0.0,0.0,0.0) for a in 1:B+1, b in 1:B]

lam = 1.0
adhm_data[1,1] = Quaternion(0.0,0.0,0.0,1.0)*lam
adhm_data[1,2] = Quaternion(0.0,0.0,1.0,0.0)*lam

adhm_data[2,1] = Quaternion(1.0,0.0,0.0,0.0)*lam/sqrt(2)
adhm_data[2,2] = Quaternion(0.0,1.0,0.0,0.0)*lam/sqrt(2)

adhm_data[3,1] = Quaternion(0.0,1.0,0.0,0.0)*lam/sqrt(2)
adhm_data[3,2] = Quaternion(-1.0,0.0,0.0,0.0)*lam/sqrt(2)

b_2_tor_skyrmion = Skyrmion(30,0.2);
make_ADHM!(b_2_tor_skyrmion, adhm_data)

# output


```

Again, we can calculate its energy.
```jldoctest adhm_data
julia> Energy(b_2_tor_skyrmion)
2.3885880515156876
```

```@meta
ShareDefaultModule = false
DocTestSetup = nothing
DocTestTeardown = nothing
DocTestFilters = nothing
```

## Save your Skyrmion

Once you've created your skyrmion, you can save a copy of it in a folder.
This folder will contain the pion fields in the [`HDF5` format](https://en.wikipedia.org/wiki/Hierarchical_Data_Format).
These files can be read by all popular programming languages, and allows for excellent compression.
In addition we save a human-readable `metadata.toml` file, containing the metadata associated with the Skyrmion, such as it's grid properties.
Overall the folder structure is
``` bash
my_output_folder/
    pion_field.h5
    metadata.toml
```

To save your skyrmion, simply run
``` julia
save_skyrmion(my_skyrmion, path="path/to/my_output_folder")
```

By default, this function will not overwrite existing folders.
You can add additional metadata by passing a dictionary to the `additional_metadata` argument.
Read more in the API.
To load a save Skyrion, use the `load_skyrmion` function:
``` julia
loaded_skyrmion = load_skyrmion("path/to/my_output_folder")
```
