
using Skyrmions3D

a_skyrmion = Skyrmion(5,0.2)

@test Skyrmions3D.getX(a_skyrmion, 1, 1, 1) == a_skyrmion.pion_field[1,1,1,:]

set_dirichlet!(a_skyrmion)

@test Skyrmions3D.getDP(a_skyrmion, 3, 3, 3) == zeros(3,4)

@test Skyrmions3D.getders_local_np(a_skyrmion,3,3,3) == (a_skyrmion.pion_field[3,3,3,:], zeros(3,4), zeros(3,4), zeros(3,4))

set_periodic!(a_skyrmion)

@test Skyrmions3D.getDP(a_skyrmion, 1, 1, 1) == zeros(3,4)

@test Skyrmions3D.getders_local_p(a_skyrmion,1,1,1) == (a_skyrmion.pion_field[1,1,1,:], zeros(3,4), zeros(3,4), zeros(3,4))


set_neumann!(a_skyrmion)

@test Skyrmions3D.getDP(a_skyrmion, 1, 1, 1) == zeros(3,4)

@test Skyrmions3D.getders_local_p(a_skyrmion,1,1,1) == (a_skyrmion.pion_field[1,1,1,:], zeros(3,4), zeros(3,4), zeros(3,4))





