"""
    make_RM_product!(skyrmion, X_list) 

Makes a product approximation of many rational map skyrmions, determined through the  list `X_list`. The final field is written into `skyrmion`.

The formatting of the list is as follow:
`X_list = [ data_1, data_2, data_3, ... ]`
where
`data_1 = [ p(z), q(z), f(r), X, θiso, n_iso, θrot, n_rot ]`

See also [`product`]

# Example of list
```
p1(z) = z; q1(z) = 1; f1(r) = 4*atan(exp(-r));
p2(z) = z^2; q2(z) = 1; f2(r) = 4*atan(exp(-0.7*r));
X_list = [ [ p1, q1, f1, [0.0,0.0,1.5], 0.0, [0.0,0.0,1.0], 0.0, [0.0,0.0,1.0] ], [ p2, q2, f2, [0.0,0.0,-1.5], pi, [1.0,0.0,0.0], 0.0, [0.0,0.0,1.0] ] ]
```

# Technical details

The product is taken pairwise in order. E.g. for a list of 3 skyrmions, we first calculate the symmetrised product of the first and second skyrmions then calculate the symmtrised product with the third skyrmion. Hence the final solution is not symmetric under permutations.

"""
function make_RM_product!(sk, Xs)

    lp = sk.grid.lp
    ls = sk.grid.ls

    temp_sk = Skyrmion(lp, ls)

    a=1
    if size(Xs[a])[1] == 8
        make_rational_map!(
            sk,
            Xs[a][1],
            Xs[a][2],
            Xs[a][3],
            X = Xs[a][4],
            iTH = Xs[a][5],
            i_n = Xs[a][6],
            jTH = Xs[a][7],
            j_n = Xs[a][8],
        )
    else
        make_rational_map!(
            sk,
            Xs[a][1],
            Xs[a][2],
            X = Xs[a][3],
            iTH = Xs[a][4],
            i_n = Xs[a][5],
            jTH = Xs[a][6],
            j_n = Xs[a][7],
            print_things = false,
        )
    end

    for a = 2:size(Xs)[1]

        if size(Xs[a])[1] == 8
            make_rational_map!(
                temp_sk,
                Xs[a][1],
                Xs[a][2],
                Xs[a][3],
                X = Xs[a][4],
                iTH = Xs[a][5],
                i_n = Xs[a][6],
                jTH = Xs[a][7],
                j_n = Xs[a][8],
            )
        else
            make_rational_map!(
                temp_sk,
                Xs[a][1],
                Xs[a][2],
                X = Xs[a][3],
                iTH = Xs[a][4],
                i_n = Xs[a][5],
                jTH = Xs[a][6],
                j_n = Xs[a][7],
                print_things = false,
            )
        end

        product_approx!(sk, temp_sk)
    end

end

"""
    make_rational_map!(skyrmion, pfn, qfn, prof; kwargs... )
    
Writes a rational map skyrmion in to `skyrmion`. The rational map is given by the polynomials R(z) = p(z)/q(z) and the profile f(r).

If no `f` is given, the function will find an OK approximation for the profile.

# Optional arguments
-  `X=[0.0,0.0,0.0]`: translate the initial skyrmion by `X`
-  `iTH = 0.0`: isorotate by initial skyrmion by `iTH`
-  `i_n = 0.0`: isorotate initial skyrmion around `i_n`
-  `jTH = 0.0`: isorotate by initial skyrmion by `jTH`
-  `j_n = 0.0`: isorotate initial skyrmion around `j_n`

"""
function make_rational_map!(
    skyrmion,
    pfn,
    qfn,
    prof = nothing;
    baryon = nothing,
    X = [0.0, 0.0, 0.0],
    iTH = 0.0,
    i_n = [0.0, 0.0, 1.0],
    jTH = 0.0,
    j_n = [0.0, 0.0, 0.0],
    print_things = true,
)
    # If a profile function is not given choose a sensible choice
    if isnothing(prof)
        # If the baryon number has not been specified, deduce it from from the
        # polynomials provided
        if isnothing(baryon)
            baryon1 = abs((log(pfn(10000)) - log(pfn(1)))/log(10000))
            baryon2 = abs((log(qfn(10000)) - log(qfn(1)))/log(10000))
            baryon = max(round(baryon1), round(baryon2))
            if print_things == true
                println(
                    "I think your baryon number is ",
                    baryon,
                    ". If it is not, include '; baryon=B' in your argument.",
                )
            end
        end
        R(z) = pfn(z)/qfn(z)
        k1, k2 = getOKprofile(1.0, baryon, getI(R), skyrmion.mpi)
        # Need what seems like an extra line here to avoid "cannot add method
        # to function argument" error here. 
        nprof(r) = pi/(1 - tanh(-k2*k1))*(-tanh(k2*(r - k1)) + 1.0);
        prof = nprof
    end

    lp, x = skyrmion.grid.lp, skyrmion.grid.x

    RI = R_from_axis_angle(iTH, i_n)
    RJ = R_from_axis_angle(jTH, j_n)

    Threads.@threads for k = 1:lp[3]
        @inbounds for j = 1:lp[2], i = 1:lp[1]

            Xto = SVector{3,Float64}(x[1][i]-X[1], x[2][j]-X[2], x[3][k]-X[3])
            Xt = RJ*Xto;

            r = sqrt(Xt[1]^2 + Xt[2]^2 + Xt[3]^2)

            sine_of_prof_r = sin(prof(r))

            if r + Xt[3] == 0
                zRM = 0.0
            else
                zRM = (Xt[1] + 1.0im*Xt[2])/(r + Xt[3])
            end

            pRM = pfn(zRM)
            qRM = qfn(zRM)

            den = real(qRM*conj(qRM) + pRM*conj(pRM))

            skyrmion.pion_field[i, j, k, 1] =
                (sine_of_prof_r/den)*real(pRM*conj(qRM) + qRM*conj(pRM))
            skyrmion.pion_field[i, j, k, 2] =
                (sine_of_prof_r/den)*imag(pRM*conj(qRM) - qRM*conj(pRM))
            skyrmion.pion_field[i, j, k, 3] =
                (sine_of_prof_r/den)*real(qRM*conj(qRM) - pRM*conj(pRM))
            skyrmion.pion_field[i, j, k, 4] = cos(prof(r))

            for a = 1:4
                if isnan(skyrmion.pion_field[i, j, k, a])
                    println(r + Xt[3])
                end
            end

            if iTH != 0.0
                skyrmion.pion_field[i, j, k, 1:3] = RI*skyrmion.pion_field[i, j, k, 1:3]
            end

        end
    end

    if skyrmion.grid.boundary_conditions == "dirichlet"
        set_dirichlet_boudary!(skyrmion)
    end

end


function getI(R)

    I_tot = 0.0
    dz_r = 0.05
    dz_i = 0.05
    dx=0.0001;
    for z_real = (-10-dz_r/2):dz_r:(10+dz_r/2), z_imag = (-10-dz_i/2):dz_i:(10+dz_i/2)
        z = z_real + 1.0im*z_imag
        Rp = ((R(z+dx)-R(z-dx))/(2dx) + (R(z+dx*1.0im)-R(z-dx*1.0im))/(2im*dx))/2.0
        I_tot += real((Rp*conj(Rp))^2*(1 + z*conj(z))^2/(1.0 + R(z)*conj(R(z)))^4)
    end

    return I_tot*dz_r*dz_i/pi

end


function getOKprofile(k1, B, I, m)

    dk1=0.001;
    k2=1.0;

    for _ = 1:10

        dE =
            (
                energy_test(k1+dk1, k2, (B, I, m)) - energy_test(k1-dk1, k2, (B, I, m))
            )/(2*dk1)

        ddE =
            (
                energy_test(k1+dk1, k2, (B, I, m)) - 2.0*energy_test(k1, k2, (B, I, m)) +
                energy_test(k1-dk1, k2, (B, I, m))
            )/(dk1^2)

        change = dE/ddE

        k1 -= change[1]

    end

    return k1, 1.0



end


function energy_test(k1, k2, (B, I, m); lp = 500, ls = 0.05, test_prof = profile(lp, ls))

    test_prof.field .=
        pi/(1.0 - tanh(-k2*k1)) .* (-tanh.(k2*(test_prof.r_grid .- k1)) .+ 1.0);
    return energy(test_prof, (B, I, m))

end

function energy(p, (B, I, m))

    ED = zeros(p.lp)

    dp = getdpB(p, 2)
    ED[2] = Ept([p.field[2], dp, p.r_grid[2], B, I, m])
    for i = 3:(p.lp-2)

        dp = getdpD(p, i)
        ED[i] = Ept([p.field[i], dp, p.r_grid[i], B, I, m])

    end

    return sum(ED)*p.ls/(3pi)

end


# From symbolics code.
function Ept(ˍ₋arg1)
    begin
        (/)(
            (+)(
                (+)(
                    (+)(
                        (+)(
                            (+)(
                                (*)(ˍ₋arg1[5], (^)((sin)(ˍ₋arg1[1]), 4)),
                                (*)((^)(ˍ₋arg1[3], 4), (^)(ˍ₋arg1[2], 2)),
                            ),
                            (*)((*)(2, (^)(ˍ₋arg1[6], 2)), (^)(ˍ₋arg1[3], 4)),
                        ),
                        (*)(
                            (*)((*)(2, ˍ₋arg1[4]), (^)(ˍ₋arg1[3], 2)),
                            (^)((sin)(ˍ₋arg1[1]), 2),
                        ),
                    ),
                    (*)(
                        (*)((*)(-2, (^)(ˍ₋arg1[6], 2)), (^)(ˍ₋arg1[3], 4)),
                        (cos)(ˍ₋arg1[1]),
                    ),
                ),
                (*)(
                    (*)((*)((*)(2, ˍ₋arg1[4]), (^)(ˍ₋arg1[3], 2)), (^)(ˍ₋arg1[2], 2)),
                    (^)((sin)(ˍ₋arg1[1]), 2),
                ),
            ),
            (^)(ˍ₋arg1[3], 2),
        )
    end
end


function getdpD(p, i)
    return (-p.field[i+2] + 8.0*p.field[i+1] - 8.0*p.field[i-1] + p.field[i-2])/(12.0*p.ls)
end

function getdpB(p, i)
    return (p.field[i+1] - p.field[i-1])/(2.0*p.ls)
end

function R_from_axis_angle(th, n)

    if th == 0.0
        return [1.0 0 0; 0 1.0 0; 0 0 1.0]
    end

    normn=norm(n)

    if normn == 0.0
        println("ERROR: your vector is zero.")
        return [1.0 0 0; 0 1.0 0; 0 0 1.0]
    end

    n=n/normn
    normv=0.0

    while normv == 0.0

        v = rand(3)

        v=cross(n, v)

        normv=norm(v)

    end

    v=v/normv
    u=cross(v, n)
    u=u/norm(u)
    M=[u v n]
    R=[cos(th) sin(th) 0.0; -sin(th) cos(th) 0.0; 0.0 0.0 1.0]

    return M*R*transpose(M)

end


"""
    make_ADHM!(skyrmion, L, M )
    
Writes an ADHM skyrmion in to `skyrmion`. The ADHM data is given by L and M. L and M should be given by `B` and `BxB` arrays of Quaternions, from the `Quaternions` package.

# Example
```
using Quaternions

B=2

L = [ Quaternion(0.0,0.0,0.0,0.0) for a in 1:B ]
M = [ Quaternion(0.0,0.0,0.0,0.0) for a in 1:B, b in 1:B ]

L[1] = Quaternion(0.0, 0.0, 0.0, sqrt(2.0))
L[2] = Quaternion(0.0, 0.0, sqrt(2.0), 0.0)

M[1,1] = Quaternion(1.0, 0.0, 0.0, 0.0)
M[1,2] = Quaternion(0.0, 1.0, 0.0, 0.0)
M[2,1] = Quaternion(0.0, 1.0, 0.0, 0.0)
M[2,2] = Quaternion(-1.0, 0.0, 0.0, 0.0)

my_skyrmion = Skyrmion(30,0.2)
make_ADHM!(my_skyrmion, L, M)
```

"""
function make_ADHM!(an_ADHM_skyrmion, L, M = nothing; tsteps = 42)
    # If only one of L, M is provided, i.e. M is nothing, then we assume that
    # the user has provided the combined matrix LM, and we extract the two. 
    if isnothing(M)
        B = size(L)[2]
        M = L[2:(B+1), 1:B]
        L = L[1, 1:B]
    end 

    B = size(L)[1]

    L_final = zeros(B, 4)
    M_final = zeros(B, B, 4)

    #if typeof(L[end]) == Quaternion{Float64}

    #println("L[1]: ", typeof(L[1].v3))

    for a = 1:B
        L_final[a, 1] = L[a].v3
        L_final[a, 2] = L[a].s
        L_final[a, 3] = L[a].v1
        L_final[a, 4] = L[a].v2
    end

    for a = 1:B, b = 1:B
        M_final[a, b, 1] = M[a, b].v3
        M_final[a, b, 2] = M[a, b].s
        M_final[a, b, 3] = M[a, b].v1
        M_final[a, b, 4] = M[a, b].v2
    end

    lstime = pi/tsteps

    ctL = [cos(lstime*(tint-1)) for tint = 1:(tsteps+1)]
    stL = [sin(lstime*(tint-1)) for tint = 1:(tsteps+1)]

    x = an_ADHM_skyrmion.grid.x
    lp = an_ADHM_skyrmion.grid.lp

    Threads.@threads for k = 1:lp[3]
        for j = 1:lp[2], i = 1:lp[1]
            @inbounds an_ADHM_skyrmion.pion_field[i, j, k, :] =
                ADHMpt(L_final, M_final, [x[1][i], x[2][j], x[3][k]], B, tsteps, ctL, stL)
        end
    end

end

function ADHMpt(L, M, y, B, tsteps, ctL, stL)

    # We use pre-allocated matrices for speed so that matrix multiplication algorithms can use special methods for different sizes. Bizarrely, this was the only way I could figure out to make the compiler know which size of matrix it was getting. Rnm = zeros( MMatrix{B,B,Float64} ) does not work!
    if B == 2
        Rnm = zeros(MMatrix{2,2,Float64})
    elseif B == 3
        Rnm = zeros(MMatrix{3,3,Float64})
    elseif B == 4
        Rnm = zeros(MMatrix{4,4,Float64})
    elseif B == 5
        Rnm = zeros(MMatrix{5,5,Float64})
    elseif B == 6
        Rnm = zeros(MMatrix{6,6,Float64})
    elseif B == 7
        Rnm = zeros(MMatrix{7,7,Float64})
    elseif B == 8
        Rnm = zeros(MMatrix{8,8,Float64})
    else
        Rnm = zeros(B, B)
    end

    Ln = zeros(B, 4)
    Mn = zeros(B, B, 4)

    U = zeros(MVector{4,Float64})
    U[1] = 1.0

    p = zeros(B, 4)

    Ω1M = zeros(MMatrix{4,4,Float64})

    Ω1 = zeros(MVector{4,Float64})
    Ω2 = zeros(MVector{4,Float64})
    Ω3 = zeros(MVector{4,Float64})

    allN = make_intial_Nα(B, tsteps+1)

    for tint = 1:(tsteps+1)

        ct = ctL[tint]

        @inbounds x1t = y[1]*ct;
        x2t = y[2]*ct;
        x3t = y[3]*ct;
        Ln .= ct .* L
        Mn .= ct .* M

        Nfy!(allN, tint, Ln, Mn, [stL[tint], x1t, x2t, x3t], p, B, Rnm)

    end

    for tint = 1:2:(tsteps-1)

        getΩf!(Ω1, allN, tint+2, tint+1, B)
        getΩf!(Ω2, allN, tint+1, tint, B)
        getΩf!(Ω3, allN, tint+2, tint, B)

        third_order_update!(Ω1M, Ω1, Ω2, Ω3)

        U = Ω1M*U

    end

    normer = sqrt(U[1]^2 + U[2]^2 + U[3]^2 + U[4]^2)
    @simd for a = 1:4
        @inbounds U[a] /= normer
    end

    return [U[2], U[3], U[4], U[1]]

end


function Nfy!(Nα, tint, L, M, y, p, B, Rnm)

    for b = 1:B, a = 1:4
        M[b, b, a] -= y[a]
    end

    for a = 1:B, c = 1:4
        p[a, c] = 0.0
    end

    makeRnm!(Rnm, L, M, B)
    iRnm = inv(Rnm)

    @inbounds for a = 1:B, b = 1:B
        p[a, 1] += iRnm[a, b]*L[b, 1]
        for c = 2:4
            p[a, c] -= iRnm[a, b]*L[b, c]
        end
    end

    modify_N!(Nα, tint, B, p, L, M)

    normalise_N!(Nα, tint, B)

end

function modify_N!(Nα, tint, B, p, L, Mmdysp)

    @fastmath @inbounds for a = 1:B

        Nα[tint, 1, 1] -= q1(L, p, a)
        Nα[tint, 1, 2] -= q2(L, p, a)
        Nα[tint, 1, 3] -= q3(L, p, a)
        Nα[tint, 1, 4] -= q4(L, p, a)

        for b = 1:B

            Nα[tint, b+1, 1] -= q1(Mmdysp, p, b, a)
            Nα[tint, b+1, 2] -= q2(Mmdysp, p, b, a)
            Nα[tint, b+1, 3] -= q3(Mmdysp, p, b, a)
            Nα[tint, b+1, 4] -= q4(Mmdysp, p, b, a)

        end

    end

end

function normalise_N!(Nα, tint, B)

    @fastmath normer =
        Nα[tint, 1, 1]^2 + Nα[tint, 1, 2]^2 + Nα[tint, 1, 3]^2 + Nα[tint, 1, 4]^2
    @inbounds for a = 2:(B+1), b = 1:4
        @fastmath normer += Nα[tint, a, b]^2
    end
    @fastmath normer = sqrt(normer)

    @inbounds for a = 1:(B+1), b = 1:4
        Nα[tint, a, b] /= normer
    end

end

function make_intial_Nα(B, tsteps)

    allN = zeros(tsteps+1, B+1, 4)

    for i = 1:tsteps
        allN[i, 1, 1] = 1.0
    end

    return allN

end

function makeRnm!(Rnm, L, M, B)
    @inbounds @fastmath for a = 1:B, b = 1:B
        Rnm[a, b] = L[a, 1]*L[b, 1] + L[a, 2]*L[b, 2] + L[a, 3]*L[b, 3] + L[a, 4]*L[b, 4]
        for c = 1:B, d = 1:4
            Rnm[a, b] += M[a, c, d]*M[c, b, d]
        end
    end
end


function q1(l, p, a)
    @inbounds l[a, 1]*p[a, 1] - l[a, 2]*p[a, 2] - l[a, 3]*p[a, 3] - l[a, 4]*p[a, 4]
end
function q2(l, p, a)
    @inbounds l[a, 2]*p[a, 1] + l[a, 1]*p[a, 2] - l[a, 4]*p[a, 3] + l[a, 3]*p[a, 4]
end
function q3(l, p, a)
    @inbounds l[a, 3]*p[a, 1] + l[a, 4]*p[a, 2] + l[a, 1]*p[a, 3] - l[a, 2]*p[a, 4]
end
function q4(l, p, a)
    @inbounds l[a, 4]*p[a, 1] - l[a, 3]*p[a, 2] + l[a, 2]*p[a, 3] + l[a, 1]*p[a, 4]
end

function q1(l, p, a, b)
    @inbounds l[a, b, 1]*p[b, 1] - l[a, b, 2]*p[b, 2] - l[a, b, 3]*p[b, 3] -
              l[a, b, 4]*p[b, 4]
end
function q2(l, p, a, b)
    @inbounds l[a, b, 2]*p[b, 1] + l[a, b, 1]*p[b, 2] - l[a, b, 4]*p[b, 3] +
              l[a, b, 3]*p[b, 4]
end
function q3(l, p, a, b)
    @inbounds l[a, b, 3]*p[b, 1] + l[a, b, 4]*p[b, 2] + l[a, b, 1]*p[b, 3] -
              l[a, b, 2]*p[b, 4]
end
function q4(l, p, a, b)
    @inbounds l[a, b, 4]*p[b, 1] - l[a, b, 3]*p[b, 2] +
              l[a, b, 2]*p[b, 3] +
              l[a, b, 1]*p[b, 4]
end


function third_order_update!(Q, q1, q2, q3)

    # Q = 4/3 q1*q2 - 1/3 q3

    Q[1, 1] =
        (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0
    Q[1, 2] =
        (-4.0*q1[2]*q2[1] - 4.0*q1[1]*q2[2] + 4.0*q1[4]*q2[3] - 4.0*q1[3]*q2[4] + q3[2])/3.0
    Q[1, 3] =
        (-4.0*q1[3]*q2[1] - 4.0*q1[4]*q2[2] - 4.0*q1[1]*q2[3] + 4.0*q1[2]*q2[4] + q3[3])/3.0
    Q[1, 4] =
        (-4.0*q1[4]*q2[1] + 4.0*q1[3]*q2[2] - 4.0*q1[2]*q2[3] - 4.0*q1[1]*q2[4] + q3[4])/3.0
    Q[2, 1] =
        (4.0*q1[2]*q2[1] + 4.0*q1[1]*q2[2] - 4.0*q1[4]*q2[3] + 4.0*q1[3]*q2[4] - q3[2])/3.0
    Q[2, 2] =
        (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0
    Q[2, 3] =
        (-4.0*q1[4]*q2[1] + 4.0*q1[3]*q2[2] - 4.0*q1[2]*q2[3] - 4.0*q1[1]*q2[4] + q3[4])/3.0
    Q[2, 4] =
        (4.0*q1[3]*q2[1] + 4.0*q1[4]*q2[2] + 4.0*q1[1]*q2[3] - 4.0*q1[2]*q2[4] - q3[3])/3.0
    Q[3, 1] =
        (4.0*q1[3]*q2[1] + 4.0*q1[4]*q2[2] + 4.0*q1[1]*q2[3] - 4.0*q1[2]*q2[4] - q3[3])/3.0
    Q[3, 2] =
        (4.0*q1[4]*q2[1] - 4.0*q1[3]*q2[2] + 4.0*q1[2]*q2[3] + 4.0*q1[1]*q2[4] - q3[4])/3.0
    Q[3, 3] =
        (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0
    Q[3, 4] =
        (-4.0*q1[2]*q2[1] - 4.0*q1[1]*q2[2] + 4.0*q1[4]*q2[3] - 4.0*q1[3]*q2[4] + q3[2])/3.0
    Q[4, 1] =
        (4.0*q1[4]*q2[1] - 4.0*q1[3]*q2[2] + 4.0*q1[2]*q2[3] + 4.0*q1[1]*q2[4] - q3[4])/3.0
    Q[4, 2] =
        (-4.0*q1[3]*q2[1] - 4.0*q1[4]*q2[2] - 4.0*q1[1]*q2[3] + 4.0*q1[2]*q2[4] + q3[3])/3.0
    Q[4, 3] =
        (4.0*q1[2]*q2[1] + 4.0*q1[1]*q2[2] - 4.0*q1[4]*q2[3] + 4.0*q1[3]*q2[4] - q3[2])/3.0
    Q[4, 4] =
        (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0

end

function getΩf!(Ω1, allv, t1, t2, B)

    Ω1 .= 0.0;

    @simd for a = 1:(B+1)
        @inbounds Ω1[1] +=
            allv[t2, a, 1]*allv[t1, a, 1] +
            allv[t2, a, 2]*allv[t1, a, 2] +
            allv[t2, a, 3]*allv[t1, a, 3] +
            allv[t2, a, 4]*allv[t1, a, 4]
        @inbounds Ω1[2] +=
            allv[t2, a, 2]*allv[t1, a, 1] - allv[t2, a, 1]*allv[t1, a, 2] -
            allv[t2, a, 4]*allv[t1, a, 3] + allv[t2, a, 3]*allv[t1, a, 4]
        @inbounds Ω1[3] +=
            allv[t2, a, 3]*allv[t1, a, 1] + allv[t2, a, 4]*allv[t1, a, 2] -
            allv[t2, a, 1]*allv[t1, a, 3] - allv[t2, a, 2]*allv[t1, a, 4]
        @inbounds Ω1[4] +=
            allv[t2, a, 4]*allv[t1, a, 1] - allv[t2, a, 3]*allv[t1, a, 2] +
            allv[t2, a, 2]*allv[t1, a, 3] - allv[t2, a, 1]*allv[t1, a, 4]
    end

end
