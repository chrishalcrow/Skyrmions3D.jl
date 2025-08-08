
using Skyrmions3D

a_skyrmion = Skyrmion(5, 0.2)
b_skyrmion = Skyrmion(6, 0.2)

@test_throws Exception product_approx(a_skyrmion, b_skyrmion)

b_skyrmion = Skyrmion(5, 0.2)

@test product_approx(a_skyrmion, b_skyrmion).pion_field == a_skyrmion.pion_field

product_approx!(a_skyrmion, b_skyrmion)
@test b_skyrmion.pion_field == a_skyrmion.pion_field

translate_sk!(a_skyrmion, X = [0.0, 0.0, -4.0])
@test a_skyrmion.pion_field == b_skyrmion.pion_field

rotate_sk!(a_skyrmion, theta = 0.5, n = [0.0, 0.0, 1.0])
@test a_skyrmion.pion_field == b_skyrmion.pion_field

isorotate_sk!(a_skyrmion, theta = 0.5, n = [0.0, 1.0, 0.0])
@test a_skyrmion.pion_field == b_skyrmion.pion_field


# Check if transformations work in the same way for rational maps and post-RM transforms.
# Need to actually make a B=1

a_skyrmion = Skyrmion(5, 0.2)
b_skyrmion = Skyrmion(5, 0.2)
set_neumann!(a_skyrmion)
set_neumann!(b_skyrmion)

p(z) = z;
q(z) = 1;
f(r) = pi*exp(-(r .^ 3) ./ 12.0)
make_rational_map!(a_skyrmion, p, q, f)

for X0 in [
    [0.2, 0, 0, 0, 0],
    [-0.2, 0.0, 0.0],
    [0.0, 0.2, 0.0],
    [0.0, -0.2, 0.0],
    [0.0, 0.0, 0.2],
    [0.0, 0.0, -0.2],
]

    make_rational_map!(a_skyrmion, p, q, f)
    make_rational_map!(b_skyrmion, p, q, f, X = X0)
    @test translate_sk(a_skyrmion, X = X0).pion_field[3, 3, 3, :] ≈
          b_skyrmion.pion_field[3, 3, 3, :]

end

theta = 2.0*pi*rand()
n = [rand(), rand(), rand()]

make_rational_map!(a_skyrmion, p, q, f)
make_rational_map!(b_skyrmion, p, q, f, iTH = theta, i_n = n)
@test isorotate_sk(a_skyrmion, theta = theta, n = n).pion_field[3, 3, 3, :] ≈
      b_skyrmion.pion_field[3, 3, 3, :]

make_rational_map!(a_skyrmion, p, q, f)
make_rational_map!(b_skyrmion, p, q, f, jTH = theta, j_n = n)
@test rotate_sk(a_skyrmion, theta = theta, n = n).pion_field[3, 3, 3, :] ≈
      b_skyrmion.pion_field[3, 3, 3, :]


b_skyrmion = Skyrmion(60, 0.2)
set_neumann!(b_skyrmion)
make_rational_map!(b_skyrmion, p, q, f, X = [0.2, 0.0, 0.0])
center_skyrmion!(b_skyrmion)
@test isapprox(center_of_mass(b_skyrmion), [0.0, 0.0, 0.0], atol=1e-10)

b_skyrmion = Skyrmion(6, 0.2)
set_neumann!(b_skyrmion)
make_rational_map!(b_skyrmion, p, q, f)
Skyrmions3D.set_dirichlet_boudary!(b_skyrmion, vac = [2.0, 0.2, -0.3, 0.5])
@test b_skyrmion.pion_field[1, 1, 1, :] == [2.0, 0.2, -0.3, 0.5]


p1(z) = z;
q1(z) = 1;
f1(r) = 4*atan(exp(-r));
p2(z) = z^2;
q2(z) = 1;
f2(r) = 4*atan(exp(-0.7*r));
X_list = [
    [p1, q1, f1, [0.0, 0.0, 1.5], 0.0, [0.0, 0.0, 1.0], 0.0, [0.0, 0.0, 1.0]],
    [p2, q2, f2, [0.0, 0.0, -1.5], pi, [1.0, 0.0, 0.0], 0.0, [0.0, 0.0, 1.0]],
]

make_RM_product!(a_skyrmion, X_list)

@test sum(a_skyrmion.pion_field[3, 3, 3, :] .^ 2) ≈ 1.0
