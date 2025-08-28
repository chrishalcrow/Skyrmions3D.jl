
using Skyrmions3D

a_skyrmion = Skyrmion(5, 0.2)
the_same_skyrmion = Skyrmion([5, 5, 5], [0.2, 0.2, 0.2])

@test a_skyrmion.pion_field == the_same_skyrmion.pion_field
@test a_skyrmion.grid.lp == the_same_skyrmion.grid.lp
@test a_skyrmion.grid.ls == the_same_skyrmion.grid.ls
@test a_skyrmion.mpi == the_same_skyrmion.mpi
@test a_skyrmion.Fpi == the_same_skyrmion.Fpi
@test a_skyrmion.ee == the_same_skyrmion.ee
@test a_skyrmion.physical == the_same_skyrmion.physical
@test a_skyrmion.grid.dirichlet == the_same_skyrmion.grid.dirichlet
@test a_skyrmion.grid.index_grid_x == the_same_skyrmion.grid.index_grid_x
@test a_skyrmion.grid.index_grid_y == the_same_skyrmion.grid.index_grid_y
@test a_skyrmion.grid.index_grid_z == the_same_skyrmion.grid.index_grid_z
@test a_skyrmion.grid.sum_grid == the_same_skyrmion.grid.sum_grid

# setting stuff

set_mpi!(a_skyrmion, 1.0)
@test a_skyrmion.mpi == 1.0

set_Fpi!(a_skyrmion, 100)
@test a_skyrmion.Fpi == 100

set_ee!(a_skyrmion, 3.0)
@test a_skyrmion.ee == 3.0

set_lattice!(a_skyrmion, [6, 4, 7], [0.1, 0.15, 0.25])

@test a_skyrmion.grid.lp[1] == 6
@test a_skyrmion.grid.lp[2] == 4
@test a_skyrmion.grid.lp[3] == 7

@test a_skyrmion.grid.ls[1] == 0.1
@test a_skyrmion.grid.ls[2] == 0.15
@test a_skyrmion.grid.ls[3] == 0.25

@test a_skyrmion.grid.dirichlet == true

set_neumann!(a_skyrmion)
@test a_skyrmion.grid.dirichlet == false
@test a_skyrmion.grid.boundary_conditions == "neumann"

@test a_skyrmion.grid.sum_grid == [1:6, 1:4, 1:7]

@test a_skyrmion.grid.index_grid_x == [2, 1, 1, 2, 3, 4, 5, 6, 6, 5]
@test a_skyrmion.grid.index_grid_y == [2, 1, 1, 2, 3, 4, 4, 3]
@test a_skyrmion.grid.index_grid_z == [2, 1, 1, 2, 3, 4, 5, 6, 7, 7, 6]

set_periodic!(a_skyrmion)
@test a_skyrmion.grid.dirichlet == false
@test a_skyrmion.grid.boundary_conditions == "periodic"

@test a_skyrmion.grid.sum_grid == [1:6, 1:4, 1:7]

@test a_skyrmion.grid.index_grid_x == [5, 6, 1, 2, 3, 4, 5, 6, 1, 2]
@test a_skyrmion.grid.index_grid_y == [3, 4, 1, 2, 3, 4, 1, 2]
@test a_skyrmion.grid.index_grid_z == [6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2]

set_dirichlet!(a_skyrmion)
@test a_skyrmion.grid.dirichlet == true
@test a_skyrmion.grid.boundary_conditions == "dirichlet"

@test a_skyrmion.grid.sum_grid == [3:4, 3:2, 3:5]

set_physical!(a_skyrmion, true)

@test a_skyrmion.physical == true

set_physical!(a_skyrmion, false)

@test a_skyrmion.physical == false

@test Skyrmions3D.sum_grid([3, 3, 3], "dirichlet") == Skyrmions3D.sum_grid(3, "dirichlet")

@test Skyrmions3D.sum_grid([3, 3, 3], "neumann") == Skyrmions3D.sum_grid(3, "neumann")

@test Skyrmions3D.sum_grid([3, 3, 3], "periodic") == Skyrmions3D.sum_grid(3, "periodic")

@test_logs (:warn, "Unrecognised boundary conditions: unexpected behaviour may occur") Skyrmions3D.sum_grid(
    [3, 3, 3],
    "nonsense_boundary_condition",
)

@test Skyrmions3D.index_grid(6, "periodic")[end] == 2
@test Skyrmions3D.index_grid(6, "periodic")[end-1] == 1
@test Skyrmions3D.index_grid(6, "periodic")[1] == 5
@test Skyrmions3D.index_grid(6, "periodic")[2] == 6

a_skyrmion.pion_field[2, 3, 1, 1] = 2.0

@test_throws AssertionError check_if_normalised(a_skyrmion)

normer!(a_skyrmion)
check_if_normalised(a_skyrmion)

a_skyrmion.pion_field[2, 3, 1, 1] = 2.0
check_if_normalised(normer(a_skyrmion))
