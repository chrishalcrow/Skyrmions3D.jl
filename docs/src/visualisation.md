# Visualisation

Suppose we've made a tetrahedral ``B=3`` rational map skyrmion:

``` julia
using Skyrmions3D

 # A skyrmion on a 30^3 grid with lattice spacing 0.2
my_skyrmion = Skyrmion(30,0.2)
p3(z) = sqrt(3)*im*z^2 - 1
q3(z) = z*(z^2 - sqrt(3)*im)
make_rational_map!(my_skyrmion, p3, q3)
```
It would be nice to visualise this. Depending on your set-up, visualisation will work in different ways. On this page, we'll focus on static plots from the terminal.

## Plotting fields

We can plot the field with component `3` and isosurface with `iso_value=0.3` as follows:

``` julia
plot_field(my_skyrmion, iso_value=0.3, component=3)
```

This might appear under your cell in Jupyter Notebook, pop-up in a seperate window, or something else. Hopefully whatever happens, you can see the following plot:

```@raw html
<img src="../../../src/images/visualisation/B3_field.png>
```

## Saving figures

The function `plot_field`, and all other plotting functions return a `Makie` `Figure`. We can save these

``` julia
using Makie
fig = plot_field(my_skyrmion, iso_value=0.3, component=3)
Make.save("my_figure.png", fig)
```

## Plotting the baryon density

We can also plot the baryon density 

``` julia
plot_baryon_density(my_skyrmion)
```

This function will automatically compute a reasonable iso_value, but you can insert your own if you'd like.

I also implemented a fairly hideous juggling colouring. This can be useful to see symmetries. If you'd like to make this more attractive, please submit a PR!!

``` julia
plot_baryon_density(my_skyrmion)
```

