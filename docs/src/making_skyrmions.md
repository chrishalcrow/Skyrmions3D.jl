# Making Skyrmions

Skyrmions in Skyrmions3D are represented by the `Skyrmion` struct. This contains all the information needed to reproduce the skyrmion: it's pion field, grid, pion mass etc. You can create a Skyrmion by passing some of these parameters, for example

``` julia
 # A skyrmion on a 30^3 grid with lattice spacing 0.2
my_skyrmion = Skyrmion(30,0.2)

# A skyrmion on a 60x60x120 grid with lattice spacing 0.15
another_skyrmion = Skyrmion([60,60,120],[0.15,0.15,0.15]) 

# A skyrmion with pion mass 1 and periodic bounary conditions
massive_periodic_skyrmion = Skyrmion(
    [30,30,30],
    [0.2,0.2,0.2], 
    mpi=1.0, 
    boundary_conditions="periodic",
    Fpi=184, 
    ee=4.5,
)
```

Find out more about the Skyrmion constructor in the API. 


By default the pion field is set equal to `(0,0,0,1)`, which is that of a vacuum skyrmion. So we next want to add some field structure to it...

## Rational Maps

A complex rational map is defined by two complex valued polynomials; we call these the numerator `p(z)` and the denominator `q(z)`. Given these polynomials, we can create a skyrmion from the rational map approximation [houghtonRationalMapsMonopoles1998](@cite). For example, the baryon number 3 tetrahedral skyrmion can be constructed as follows:

``` julia
p3(z) = sqrt(3)*im*z^2 - 1
q3(z) = z*(z^2 - sqrt(3)*im)
make_rational_map!(my_skyrmion, p3, q3)
```

Note: By convention the "bang" character `!` at the end of function means it is a _modifiying function_. So we are modifiying `my_skyrmion`, not creating a new skyrmion.

The `make_rational_map!` function tries to estimate the degree of the rational map, and tries to find a reasonable profile function for it. If this all worked, you should get a sensible output when you compute the energy of the skyrmion.

``` julia
Energy(my_skyrmion)
>>> 4.106354608768089
```

Find out more about computing skyrmion properties in the (computing properties section)[].

The `make_rational_map!` function also accepts custom profile functions.

## ADHM data

ADHM skyrmions [corkADHMSkyrmions2022](@cite) are skyrmions generated from ADHM data. The data consists of symmetric quaternionic matrices which satisfy a constraint. Most highly symmetric skyrmions can be represented by ADHM data. This package implements the very efficient parallel transport algorithm developed by [Harland2023](@citet). We can make the baryon number 2 toroidal skyrmion as follows:

``` julia
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

make_ADHM!(my_skyrmion, adhm_data)
```







