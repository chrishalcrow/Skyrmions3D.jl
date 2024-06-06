# Welcome to the Skyrmions3D package. This is a package designed to help in the creation,
# manipulation, dynamics and plotting of Skyrmions in the Skyrme model of nuclear physics.

# # #      1. GETTING STARTED     # # #

# We first need to import the Skyrmion package and the GLMakie package, which is used for
# plotting. Additionally, we'll tell our editor to *use* GLMakie. This is a bit awkward,
# but the many commands below should mean that it is turned on on all OSs.

# Note: you need to type the next time into the REPL in Visual Studio code. That means:
# type ']',which will enter "package mode", then type
# "add https://github.com/chrishalcrow/Skyrmions3D.jl.git"

#]add https://github.com/Imzeep/Skyrmions3D.jl

using Skyrmions3D

using GLMakie
GLMakie.activate!()
Makie.inline!(false)

# Before getting started, let's try and make a "Makie" window. This is important for some
# more advanced features of Skyrmions3D. The following line should generate a _new window_.

display(plot( rand(3), rand(3) ) )

# If it does not, try this

display(GLMakie.Screen(), plot( rand(3), rand(3) ) )

# Later, when we plot, please use whichever line worked for you

# # #       2. Setting your grid        # # #

# In this code, your Skyrme field is contained in a structure which also contains the
# grid that your field lives on. The structure also contains the parameters of the model.
# Here's an example

my_skyrmion = Skyrmion( [40, 40, 40], [0.2, 0.2, 0.2] );

# This generates a Skyrme field on a 40^3 grid with 0.2 lattice spacing. Another:

my_skyrmion = Skyrmion( [40, 40, 40], [0.2, 0.2, 0.2], mpi = 1.0, Fpi=186, ee=4.0, boundary_conditions="periodic")

# This sets some parameters too and the periodicity of the fields
# You can check all the properties of the Skyrmion using overview

overview(my_skyrmion)

# You can change all the parameters using the set! functions

set_mpi!(my_skyrmion, 0.5)
set_Fpi!(my_skyrmion, 100)
set_ee!(my_skyrmion, 6.5)
set_lattice!(my_skyrmion, [60,60,60], [0.2,0.2,0.2])
set_periodic!(my_skyrmion)
# set_neumann!(my_skyrmion) # you can also make Neumann boundary conditions!
# set_dirichlet!(my_skyrmion)

overview(my_skyrmion)

# Note: whenever you are confused about a function you can use ?function_name

# # #       Initialising Rational Map skyrmions     # # #

# Let's now make a Rational Map skyrmion. To do this we can pass the numerator and
# denomenator of a Rational Map.

p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0;
q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0;
f4(r) = pi*exp( -(r.^3)./12.0 )

make_rational_map!(my_skyrmion, p4, q4, f4)

# We can calculate properties of the skyrmion using the following functions

Baryon(my_skyrmion)
Energy(my_skyrmion)

# Without a profile function, the make_rational_map! function will find an OK approximate

make_rational_map!(my_skyrmion, p4, q4; baryon=4)
Energy(my_skyrmion)

# ...and the energy of this was smaller than the other one. Excellent!

# We can take a look by plotting the baryon density. (Remember, this should pop up in a
# seperate window!)

plot_baryon_density(my_skyrmion)

# There is also a function for making a product of several RM skyrmions:

p1(z) = z; q1(z) = 1; f1(r) = 4*atan(exp(-r));
p2(z) = z^2; q2(z) = 1; f2(r) = 4*atan(exp(-0.7*r));
X_list = [ [ p1, q1, f1, [0.0,0.0,1.5], 0.0, [0.0,0.0,1.0], 0.0, [0.0,0.0,1.0] ], [ p2, q2, f2, [0.0,0.0,-1.5], pi, [1.0,0.0,0.0], 0.0, [0.0,0.0,1.0] ] ]

make_RM_product!(my_skyrmion, X_list)

plot_baryon_density(my_skyrmion, iso_value=1.0)

# # #       Calculating skyrmion properties     # # #

# As well as Energy and Baryon number, we can calculate many other things. 
# First, we can look at moments using the `moment` keyword

Energy(my_skyrmion, moment=2)

# It can be easier to compare to other work by using Physical units on. Let's do this

set_physical!(my_skyrmion,true)
U = compute_current(my_skyrmion, label="uMOI")

# ... and let's turn it off again

set_physical!(my_skyrmion,false)

# You can also get a density if you'd like

baryon_density = Baryon(my_skyrmion, density=true);

# Most of the properties needed for physics are in the `compute_current` function (thanks 
# to Alberto for the suggestion). You can find the moments of inertia

U = compute_current(my_skyrmion, label="uMOI")
W = compute_current(my_skyrmion, label="wMOI")
V = compute_current(my_skyrmion, label="vMOI")

# The other currents currently included are the u and w Axial currents

compute_current(my_skyrmion, label="uAxial")
compute_current(my_skyrmion, label="wAxial")

# as well as the stress tensor

compute_current(my_skyrmion, label="stress")

# and the Noether iso and Axial currents.

compute_current(my_skyrmion, label="NoetherIso")
compute_current(my_skyrmion, label="NoetherAxial")

# These are all defined in the technical notes (coming soon!) and will be updated to work
# with physical units soon...

# All these functions work with moment and density. Let's find the density of the (3,3) 
# component of the second moment of U

U33_density = compute_current(my_skyrmion, label="uMOI", indices=[3,3], density=true, moment=2);

# We can plot this...

volume(U33_density[1,1,:,:,:], algorithm=:iso, isorange=1, isovalue=10)

# Cool! Find out more about plotting by Googling `Makie'.

# # #       Initialising ADHM skyrmions     # # #

# Another way to create is to make ADHM skyrmion using Quaternions.

B=3

adhm_data = [ Quaternion(0.0,0.0,0.0,0.0) for a in 1:B+1, b in 1:B  ]

lam = 1.0

adhm_data[1,1] = Quaternion(lam, 0.0, 0.0, 0.0)
adhm_data[1,2] = Quaternion(0.0, lam, 0.0, 0.0)
adhm_data[1,3] = Quaternion(0.0, 0.0, lam, 0.0)

adhm_data[2,2] = Quaternion(0.0, 0.0, lam, 0.0)
adhm_data[2,3] = Quaternion(0.0, lam, 0.0, 0.0)

adhm_data[3,1] = Quaternion(0.0, 0.0, lam, 0.0)
adhm_data[3,3] = Quaternion(lam, 0.0, 0.0, 0.0)

adhm_data[4,1] = Quaternion(0.0, lam, 0.0, 0.0)
adhm_data[4,2] = Quaternion(lam, 0.0, 0.0, 0.0)

# The make_ADHM! function is similar to the make_rational_map! function:

make_ADHM!(my_skyrmion,adhm_data)

plot_baryon_density(my_skyrmion, iso_value=2.0)

# # #       Transforming data     # # #

# One of the most important things we can do is transform our skyrmions.
# That means translating, rotating and isorotating. Functions which end
# in ! modify the underlying field

translate_sk!(my_skyrmion, [1.0,0.0,0.0])
rotate_sk!(my_skyrmion, pi/8, [0.0,0.0,1.0])
isorotate_sk!(my_skyrmion, pi, [1.0,0.0,0.0])

plot_baryon_density(my_skyrmion, iso_value=2.0)

# Or if you don't use a !, you can make new skyrmions. This is good for
# making something using the product approximation.

skyrmion_1 = translate_sk(my_skyrmion, [1.0,0.0,0.0]);
skyrmion_2 = translate_sk(my_skyrmion, [-3.0,0.0,0.0]);

skyrmion_3 = product_approx(skyrmion_1, skyrmion_2);

plot_baryon_density(skyrmion_3, iso_value=2.0)

# Looks like the skyrmions are a bit squished. This is reflected in their baryon number:

Baryon(skyrmion_3)

# You'll need to try a bigger grid!

# Skyrmions are quite large objects. A 40^3 skyrmion takes up 4MB of RAM.
# You can clear them using `Nothing`:

skyrmion_1 = Nothing
skyrmion_2 = Nothing
skyrmion_3 = Nothing

# # #       Flowing data     # # #

# Let's start fresh with a new my_skyrmion:

my_skyrmion = Skyrmion( 40, 0.2 );

p4(z) = z^4 + 3.0*sqrt(3.0)*im*z^2 + 1.0;
q4(z) = z^4 - 3.0*sqrt(3.0)*im*z^2 + 1.0;

make_rational_map!(my_skyrmion, p4, q4)

plot_baryon_density(my_skyrmion)

# We can flow this using arrested Newton flow. You can flow for a set numbers of steps

arrested_newton_flow!(my_skyrmion, steps=100)

# The error is given by the maximum absolute value of the variation. You can also make the
# flow tell you the energy/error at some intervals using `checks`

arrested_newton_flow!(my_skyrmion, steps=100, checks=10)

# Another option is to stop after you reach some numerical tolerance

arrested_newton_flow!(my_skyrmion, tolerance=0.3, checks=10)

# The error in arrested Newton flow fluctuates a lot. This is in contrast to gradient flow

gradient_flow!(my_skyrmion, steps=500, checks=50)

# Again, you can set a tolerance instead

gradient_flow!(my_skyrmion, tolerance=0.0075, checks=50)

# Note that we've never set the time step. The algorithm picks a reasonable choice.
# But you might want to set it yourself, you can do this as follows:

gradient_flow!(my_skyrmion, tolerance=0.01, checks=10, dt=0.01)

# Oh! That time step was too big, and it has error'd. Maybe we can salvage this...

gradient_flow!(my_skyrmion, steps=100, checks=25)

# Not this time... We can start again.

make_rational_map!(my_skyrmion, p4, q4)

# We can have a bit more fun by using the interactive_flow feature

interactive_flow(my_skyrmion)

# Much more fun!

# That's all for this tutorial. There are a few features hidden in the ?help.
# If you have naming, features or usage feedback, please let me know!
# Once it's an official package (hopefully soon) lots of my decisions will be
# baked in for a long time!



