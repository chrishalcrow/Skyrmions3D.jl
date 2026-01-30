# [Properties](@id tut/props)

A skyrmion has many properties.
In general, `Skyrmions3D` allows you to compute the integrand or integral of the property in question.
For all the following examples, if you pass `density=true` you'll get the integrand.
You can also get the `n`th moment by passing `moment=n`.

Suppose you have a ``B=4`` cubic skyrmion

```julia
using Skyrmions3D
skyrmion = Skyrmion(40,0.2)
p4(z) = z^4 + 2*sqrt(3)*im*z^2 + 1
q4(z) = z^4 - 2*sqrt(3)*im*z^2 + 0
make_rational_map!(skyrmion, p4, q4)
```

## Baryon number

Basic use:

```julia
Baryon(skyrmion)
>>> 3.968749058206481
```

Other examples:

```julia
baryon_density = Baryon(skyrmion, density=true)
rms_baryon = sqrt(Baryon(skyrmion, moment=2))
>>> 3.882435997469747
```

## Energy

Basic use:

```julia
Energy(skyrmion)
>>> 5.371145445834469
```

Note: you can turn on physical units to get the energy in those units:

```julia
skyrmion.Fpi = 174
skyrmion.ee = 3
skyrmion.physical = true
Energy(skyrmion)
>>> (6149.285364807477, "MeV")
```

Note that the energy of the skyrmion is dependent both on its underlying pion field (which determines the energy density at a point) but also the grid, as the energy is calculated as the sum of the energy density over the grid.
If the grid is not sufficiently large, the computed energy will be smaller than the 'true' value.

## Currents

For currents, we use the function `compute_current` and pass a label.
We list the currents and labels now:

- Rotational moment of intertia, `vMOI`
- Isorotational moment of intertia, `uMOI`
- Mixed moment of intertia, `wMOI`
- Iso-Axial, `uAxial`
- Rotational-Axial, `vAxial`
- Stress-tensor, `stress`
- Noether-iso, `NoetherIso`
- Noether-axial, `NoetherAxial`

Like the `Energy` and `Baryon` we can get the densities by passing `density=true` and the `n`th moment by passing `moment=n`.
For instance, the Isorotational moment of inertia of the cubic rational map skyrmion is

```julia
U = compute_current(skyrmion, label="uMOI")
```

The density of the stress-tensor is given by

```julia
stress_density = compute_current(skyrmion, label="stress", density=true)
```

The second moment of the Noether-iso current is

```julia
noether_moment = compute_current(skyrmion, label="NoetherIso", density=false, moment=2)
```
