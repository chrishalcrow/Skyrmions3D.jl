# API

Here are all the exported functions of Skyrmions3D.

## Types

```@docs
Skyrmion
```

## Set up the system

```@docs
set_lattice!
set_dirichlet!
set_neumann!
set_periodic!
set_mpi!
set_physical!
set_Fpi!
set_ee!
```

## Getters and checkers
```@docs
get_grid
check_if_normalised
get_field
```

## Create

```@docs
make_rational_map!
make_RM_product!
make_ADHM!
product_approx!
product_approx
```

## Transform

```@docs
translate_sk!
translate_sk
rotate_sk!
rotate_sk
isorotate_sk!
isorotate_sk
center_skyrmion!
normer
normer!
```

## Probe

```@docs
overview
Energy
Baryon
center_of_mass
compute_current
rms_baryon
```

## Flow

```@docs
gradient_flow!
newton_flow!
arrested_newton_flow!
```

## Visualise

```@docs
plot_overview
plot_field
plot_baryon_density
interactive_flow
```