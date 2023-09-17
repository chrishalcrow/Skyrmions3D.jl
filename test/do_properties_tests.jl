
using Skyrmions3D

a_skyrmion = Skyrmion(5,0.2)
b_skyrmion = deepcopy(a_skyrmion)

overview(a_skyrmion)

@test Energy(a_skyrmion) ≈ 0.0

@test Energy(a_skyrmion, density=true) ≈ 0.0.*similar(a_skyrmion.pion_field[:,:,:,1])

@test Skyrmions3D.engpt(zeros(3,4),[0.0,0.0,0.0,1.0][4],rand()) ≈ 0.0


@test Baryon(a_skyrmion) ≈ 0.0
@test Baryon(a_skyrmion, density=true) ≈ 0.0.*similar(a_skyrmion.pion_field[:,:,:,1])


@test center_of_mass(a_skyrmion) == [0.0,0.0,0.0]

@test rms_baryon(a_skyrmion) == 0.0

@test compute_current(a_skyrmion, label="uMOI") == zeros(3,3)
@test compute_current(a_skyrmion, label="wMOI") == zeros(3,3)
@test compute_current(a_skyrmion, label="vMOI") == zeros(3,3)
@test compute_current(a_skyrmion, label="uAxial") == zeros(3,3)
@test compute_current(a_skyrmion, label="wAxial") == zeros(3,3)
@test compute_current(a_skyrmion, label="NoetherIso") == zeros(3,3)
@test compute_current(a_skyrmion, label="NoetherAxial") == zeros(3,3)
@test compute_current(a_skyrmion, label="stress") == zeros(3,3)

@test compute_current(a_skyrmion, label="uMOI", indices=[1,1])[1] == 0.0
@test compute_current(a_skyrmion, label="wMOI", indices=[1,2])[1] == 0.0
@test compute_current(a_skyrmion, label="vMOI", indices=[1,3])[1] == 0.0
@test compute_current(a_skyrmion, label="uAxial", indices=[2,1])[1] == 0.0
@test compute_current(a_skyrmion, label="wAxial", indices=[2,2])[1] == 0.0
@test compute_current(a_skyrmion, label="NoetherIso", indices=[2,3])[1] == 0.0
@test compute_current(a_skyrmion, label="NoetherAxial", indices=[3,1])[1] == 0.0
@test compute_current(a_skyrmion, label="stress", indices=[3,2])[1] == 0.0

@test compute_current(a_skyrmion, label="uMOI", indices=[3,3], density=true) == zeros(1,1,5,5,5)
@test compute_current(a_skyrmion, label="uMOI", density=true) == zeros(3,3,5,5,5)



