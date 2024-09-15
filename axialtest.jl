using Skyrmions3D
using GLMakie
using CSV
using DataFrames
using Serialization
using LinearAlgebra
GLMakie.activate!()
Makie.inline!(false)

nuc = Skyrmion([60, 60, 60], [0.2, 0.2, 0.2], mpi = 0, Fpi = 108, ee = 4.84, boundary_conditions = "dirichlet")

p1(z) = z
q1(z) = 1
f1(r) = 4 * atan(exp(-r))

make_rational_map!(nuc, p1, q1, f1)


u = [0, 0, 1]  
v = [1, 0, 0]  
a = 0
r = 10

axial_symmetry_plot(nuc, u, v, a, r, points = 50, field=4)
