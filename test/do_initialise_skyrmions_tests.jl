
using Skyrmions3D

a_skyrmion = Skyrmion(5, 0.2)

@test Skyrmions3D.R_from_axis_angle(0.0, [1.0, 0.0, 0.0]) == [1.0 0 0; 0 1.0 0; 0 0 1.0]
@test Skyrmions3D.R_from_axis_angle(1.0, [0.0, 0.0, 0.0]) == [1.0 0 0; 0 1.0 0; 0 0 1.0]
@test Skyrmions3D.R_from_axis_angle(pi/4, [0.0, 0.0, 1.0]) ≈
      1/sqrt(2.0) .* [1.0 -1.0 0.0; 1.0 1.0 0.0; 0.0 0.0 sqrt(2.0)]




tsteps = 100
lstime = pi/tsteps

ctL = [cos(lstime*(tint-1)) for tint = 1:(tsteps+1)]
stL = [sin(lstime*(tint-1)) for tint = 1:(tsteps+1)]

Skpt = Skyrmions3D.ADHMpt(
    [0.0 1.0 0.0 0.0],
    zeros(1, 1, 4),
    [0.0, 0.0, 1.0],
    1,
    tsteps,
    ctL,
    stL,
)

@test abs(Skpt[4] + cos(pi/sqrt(2.0))) < 10^(-5)
@test abs(Skpt[3] + sin(pi/sqrt(2.0))) < 10^(-5)

B=1

L = [Quaternion(0.0, 0.0, 0.0, 0.0) for a = 1:B, b = 1:B]
M = [Quaternion(0.0, 0.0, 0.0, 0.0) for a = 1:B, b = 1:B]

lam = 1.0

L[1, 1] = Quaternion(lam, 0.0, 0.0, 0.0)
M[1, 1] = Quaternion(0.0, 0.0, 0.0, 0.0)

make_ADHM!(a_skyrmion, L, M)

@test a_skyrmion.pion_field[1, 1, 1, 1]^2 +
      a_skyrmion.pion_field[1, 1, 1, 2]^2 +
      a_skyrmion.pion_field[1, 1, 1, 3]^2 +
      a_skyrmion.pion_field[1, 1, 1, 4]^2 ≈ 1
