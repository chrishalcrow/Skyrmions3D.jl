# Transforming & Combining Skyrmions

Once you have a skyrmion, there's a lot you can do with it. Suppose you've made a `B=2` torus like so

```julia
using Skyrmions3D
toroidal_skyrmion = Skyrmion(30,0.2)
p(z) = z^2
q(z) = 1
make_rational_map!(toroidal_skyrmion, p, q)
```

It's important here to remember that functions ending in `!` modify the skyrmion, while functions with no `!` will return a new skyrmion. With that in mind, you could now apply a...

## Translation

Translations take in a skyrmion and a 3-vector. You can make a new skyrmion

```julia
X = [0.0,0.0,1.0]
translated_skyrmion = translate_sk(toroidal_skyrmion, X=X)
```

or modify an existing skyrmion


``` julia
translate_sk!(toroidal_skyrmion, X=X)
```

You can also "center" your skyrmion. This computes the center of mass, then translates so that the skyrmion has zero center of mass:

``` julia
# compute CoM for the translated skyrmion
center_of_mass(toroidal_skyrmion)
>>>3-element Vector{Float64}:
 -7.075309546463488e-17
 -2.5354376531123876e-17
  0.976781807200457
```

``` julia
center_skyrmion!(toroidal_skyrmion)
center_of_mass(toroidal_skyrmion)
>>>3-element Vector{Float64}:
 -1.0516044516158566e-16
  8.101616160453117e-17
  1.0194091904102672e-8
```

## Rotation

Rotations take in a skyrmion, an angle and a rotation vector. You can make a new skyrmion

```julia
rotated_skyrmion = rotate_sk(toroidal_skyrmion, theta=pi/4, n=[1.0,0.0,0.0])
```

or modify an existing skyrmion

```julia
theta = pi/4
n = [1.0,0.0,0.0]
rotate_sk!(toroidal_skyrmion, theta=theta, n=n)
```

## Isorotation

Isorotations take in a skyrmion, an angle and a rotation vector. You can make a new skyrmion

```julia
isorotated_skyrmion = isorotate_sk(toroidal_skyrmion, theta=pi/4, n=[1.0,0.0,0.0])
```

or modify an existing skyrmion

```julia
theta = pi/4
n = [1.0,0.0,0.0]
isorotate_sk!(toroidal_skyrmion, theta=theta, n=n)
```

## Product Approximation

We can also combine skyrmions using the Product Approximation. This codebase always applies a symmetrised product approx. Let's first create two tori

```julia
p(z) = z^2
q(z) = 1
make_rational_map!(toroidal_skyrmion, p, q)

skyrmion_up = translate_sk(toroidal_skyrmion, X=[0.0,0.0,1.5])
skyrmion_down = translate_sk(toroidal_skyrmion, X=[0.0,0.0,-1.5])
```

Then apply a product approximation

```julia
product_skyrmion = product_approx(skyrmion_up, skyrmion_down)
Baryon(product_skyrmion)
>>> 3.938613862543325
```

Oop, the product skyrmion is losing some baryon density - the box might be a bit small!
