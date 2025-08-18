# Set of tests to check numerical results against known analytic results
# Always do this with a pion mass - a more stringent test than massless
# Good resources for these tests: https://arxiv.org/pdf/0707.0868

using Skyrmions3D

spherically_symmetric_skyrmion = Skyrmion(10, 0.2, mpi = 1.0)
p1(z) = z
q1(z) = 1
f1(r) = pi*exp(-(r .^ 3) ./ 12.0)
make_rational_map!(spherically_symmetric_skyrmion, p1, q1, f1)

U1 = compute_current(spherically_symmetric_skyrmion, label = "uMOI")
V1 = compute_current(spherically_symmetric_skyrmion, label = "vMOI")
W1 = compute_current(spherically_symmetric_skyrmion, label = "wMOI")

tol = 1e-10

@test isapprox(U1[1, 1], U1[2, 2], atol = tol)
@test isapprox(U1[2, 2], U1[3, 3], atol = tol)

@test isapprox(W1[1, 1], W1[2, 2], atol = tol)
@test isapprox(W1[2, 2], W1[3, 3], atol = tol)

@test isapprox(V1[1, 1], V1[2, 2], atol = tol)
@test isapprox(V1[2, 2], V1[3, 3], atol = tol)

for a = 1:3, b = 1:3
    if a != b
        @test isapprox(U1[a, b], 0, atol = tol)
        @test isapprox(W1[a, b], 0, atol = tol)
        @test isapprox(V1[a, b], 0, atol = tol)
    end
end

cubic_skyrmion = Skyrmion(10, 0.4, mpi = 1.0)
p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0;
q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0;
f4(r) = pi*exp(-(r .^ 3) ./ 12.0)

make_rational_map!(cubic_skyrmion, p4, q4, f4)

U4 = compute_current(cubic_skyrmion, label = "uMOI")
V4 = compute_current(cubic_skyrmion, label = "vMOI")
W4 = compute_current(cubic_skyrmion, label = "wMOI")

tol = 1e-7

@test isapprox(U4[1, 1], U4[2, 2], atol = tol)

@test isapprox(V4[1, 1], V4[2, 2], atol = tol)
@test isapprox(V4[2, 2], V4[3, 3], atol = tol)

for a = 1:3, b = 1:3
    if a != b
        @test isapprox(U4[a, b], 0, atol = tol)
        @test isapprox(V4[a, b], 0, atol = tol)
    end
    @test isapprox(W4[a, b], 0, atol = tol)
end

# tetrahedral ADHM skyrmion

using Quaternions

tetrahedral_skyrmion = Skyrmion(10, 0.4, mpi = 1.0)

B=3

adhm_data = [Quaternion(0.0, 0.0, 0.0, 0.0) for a = 1:(B+1), b = 1:B]

lam = 1.0

adhm_data[1, 1] = Quaternion(lam, 0.0, 0.0, 0.0)
adhm_data[1, 2] = Quaternion(0.0, lam, 0.0, 0.0)
adhm_data[1, 3] = Quaternion(0.0, 0.0, lam, 0.0)

adhm_data[2, 2] = Quaternion(0.0, 0.0, lam, 0.0)
adhm_data[2, 3] = Quaternion(0.0, lam, 0.0, 0.0)

adhm_data[3, 1] = Quaternion(0.0, 0.0, lam, 0.0)
adhm_data[3, 3] = Quaternion(lam, 0.0, 0.0, 0.0)

adhm_data[4, 1] = Quaternion(0.0, lam, 0.0, 0.0)
adhm_data[4, 2] = Quaternion(lam, 0.0, 0.0, 0.0)

# The make_ADHM! function is similar to the make_rational_map! function:

make_ADHM!(tetrahedral_skyrmion, adhm_data)

U3 = compute_current(tetrahedral_skyrmion, label = "uMOI")
V3 = compute_current(tetrahedral_skyrmion, label = "vMOI")
W3 = compute_current(tetrahedral_skyrmion, label = "wMOI")

@test isapprox(U3[1, 1], U3[2, 2], atol = tol)
@test isapprox(U3[2, 2], U3[3, 3], atol = tol)

@test isapprox(W3[1, 1], W3[2, 2], atol = tol)
@test isapprox(W3[2, 2], W3[3, 3], atol = tol)

@test isapprox(V3[1, 1], V3[2, 2], atol = tol)
@test isapprox(V3[2, 2], V3[3, 3], atol = tol)
