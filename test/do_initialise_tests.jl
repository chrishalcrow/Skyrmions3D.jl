
using Skyrmions3D

a_skyrmion = Skyrmion(5,0.2)
the_same_skyrmion = Skyrmion([5,5,5],[0.2,0.2,0.2])

@test a_skyrmion.pion_field == the_same_skyrmion.pion_field
@test a_skyrmion.lp == the_same_skyrmion.lp
@test a_skyrmion.ls == the_same_skyrmion.ls
@test a_skyrmion.mpi == the_same_skyrmion.mpi
@test a_skyrmion.Fpi == the_same_skyrmion.Fpi
@test a_skyrmion.ee == the_same_skyrmion.ee
@test a_skyrmion.physical == the_same_skyrmion.physical
@test a_skyrmion.periodic == the_same_skyrmion.periodic
@test a_skyrmion.index_grid_x == the_same_skyrmion.index_grid_x
@test a_skyrmion.index_grid_y == the_same_skyrmion.index_grid_y
@test a_skyrmion.index_grid_z == the_same_skyrmion.index_grid_z
@test a_skyrmion.sum_grid == the_same_skyrmion.sum_grid

# setting stuff

set_mpi!(a_skyrmion, 1.0)
@test a_skyrmion.mpi == 1.0

set_Fpi!(a_skyrmion, 100)
@test a_skyrmion.Fpi == 100

set_ee!(a_skyrmion, 3.0)
@test a_skyrmion.ee == 3.0

set_lattice!(a_skyrmion, [6,4,7], [0.1, 0.15, 0.25])

@test a_skyrmion.lp[1] == 6
@test a_skyrmion.lp[2] == 4
@test a_skyrmion.lp[3] == 7

@test a_skyrmion.ls[1] == 0.1
@test a_skyrmion.ls[2] == 0.15
@test a_skyrmion.ls[3] == 0.25

set_periodic!(a_skyrmion, true)

@test a_skyrmion.periodic == true

@test a_skyrmion.sum_grid == [ 1:6, 1:4, 1:7 ] 

@test a_skyrmion.index_grid_x == [5, 6, 1, 2, 3, 4, 5, 6, 1, 2]
@test a_skyrmion.index_grid_y == [3, 4, 1, 2, 3, 4, 1, 2]
@test a_skyrmion.index_grid_z == [6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2]

set_periodic!(a_skyrmion, false)

@test a_skyrmion.periodic == false

@test a_skyrmion.sum_grid == [ 3:4, 3:2, 3:5 ] 

set_physical!(a_skyrmion, true)

@test a_skyrmion.physical == true

set_physical!(a_skyrmion, false)

@test a_skyrmion.physical == false

@test sum_grid([3,3,3],true) == sum_grid(3,true)

@test sum_grid([3,3,3],false) == sum_grid(3,false)

@test index_grid(6)[end]   == 2
@test index_grid(6)[end-1] == 1
@test index_grid(6)[1]     == 5
@test index_grid(6)[2]     == 6

a_skyrmion.pion_field[2,3,1,1] = 2.0
normer!(a_skyrmion)
check_if_normalised(a_skyrmion)

a_skyrmion.pion_field[2,3,1,1] = 2.0
check_if_normalised(normer(a_skyrmion))


