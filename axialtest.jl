using GLMakie
GLMakie.activate!()
Makie.inline!(false)
using Skyrmions3D
using CSV
using DataFrames
using Serialization
using LinearAlgebra

nuc = Skyrmion([100, 100, 100], [0.1, 0.1, 0.1], mpi = 0, Fpi = 108, ee = 4.84, boundary_conditions = "dirichlet")

p1(z) = z
q1(z) = 1
f1(r) = 4 * atan(exp(-r))
make_rational_map!(nuc, p1, q1, f1)

arrested_newton_flow!(nuc, tolerance = 0.01, checks = 50, dt=0.01)


u = [0, 0, 1]  
v = [1, 1, 0]  
a = 0
r = 1

axial_symmetry_plot(nuc, u, v, a, r, points = 50, field=4)