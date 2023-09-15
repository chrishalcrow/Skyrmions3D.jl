
using Skyrmions3D

a_skyrmion = Skyrmion(5,0.2)

@test R_from_axis_angle(0.0, [1.0,0.0,0.0]) == [ 1.0 0 0 ; 0 1.0 0 ; 0 0 1.0 ]
@test R_from_axis_angle(1.0, [0.0,0.0,0.0]) == [ 1.0 0 0 ; 0 1.0 0 ; 0 0 1.0 ]
@test  R_from_axis_angle(pi/4, [0.0,0.0,1.0]) ≈ 1/sqrt(2.0).*[ 1.0 -1.0  0.0 ; 1.0 1.0  0.0 ; 0.0 0.0 sqrt(2.0) ]