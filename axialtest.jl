using Skyrmions3D
using GLMakie
using CSV
using DataFrames
using Serialization
using LinearAlgebra
GLMakie.activate!()
Makie.inline!(false)


nuc = deserialize("l120_B1M06.5")


u = [0, 0, 1]  
v = [1, 0, 0]  
a = 30
r = 10

axial_symmetry_plot(nuc, u, v, a, r, points = 50, field=4)
