"""
    gradient_flow!(skyrmion; steps = 1, tolerance = 0.0, dt = ls^2/100.0, checks = max(100, steps), verbose = true, dEdp = zero_array)
    
Applies a gradient flow to `skyrmion` with timestep `dt` (where `ls` is the average lattice spacing), either for `steps` steps or until the error falls below `tolerance`. The error is checked every `checks` steps.

`dEdp` is an array, initialised to be the correct shape and 0.0 everywhere, which stores the change in energy density. The flow is verbose, incrementally describing the error and new energy, if `verbose` is true.

See also [`arrested_newton_flow!`](@ref). 

"""
function gradient_flow!(
    ϕ;
    steps = 1,
    dt = ((ϕ.grid.ls[1]*ϕ.grid.ls[2]*ϕ.grid.ls[3])^(2/3))/100.0,
    tolerance = 0.0,
    checks = max(100, steps),
    verbose = true,
    dEdp = zeros(ϕ.grid.lp[1], ϕ.grid.lp[2], ϕ.grid.lp[3], 4)
)

    if tolerance == 0 && checks > steps
        checks = steps
    end

    if verbose
        println("initial: energy: ", Energy(ϕ))

    end

    counter = 0
    prev_error = 1.0e9

    while counter < steps

        gradient_flow_for_n_steps!(ϕ, dEdp, checks, dt)

        err = max_abs_err(dEdp)
        if err > 3*prev_error
            error(
                "Suspected numerical blowup. Please use a smaller dt. Currently, dt = ",
                dt,
            )
        end
        prev_error = err

        counter += checks

        if verbose
            println("after ", counter, " steps, error = ", round(err, sigdigits = 4))
        end

        if tolerance != 0.0    # => we are in tol mode    
            if err < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += checks    # => continue the while loop
            end
        end

    end

    if verbose
        println("final energy: ", Energy(ϕ))
    end

    return

end

function gradient_flow_for_n_steps!(ϕ, dEdp, n, dt)
    for _ = 1:n
        gradient_flow_1_step!(ϕ, dEdp, dt)
    end
end

function gradient_flow_1_step!(sk, dEdp, dt)

    getdEdp!(sk, dEdp)
    sk.pion_field .+= dt .* dEdp;
    normer!(sk)

end


# we split dEdp into Dirichlet and other options here, so that the if statement happens once, rather than inside a for loop.
# This creates some code duplication, for a ~10% performance boost.

function getdEdp!(sk, dEdp)
    if sk.grid.boundary_conditions == "dirichlet"
        getdEdp_np!(sk, dEdp)
    else
        getdEdp_p!(sk, dEdp)
    end
end


function getdEdp_np!(sk, dEdp)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_np(sk, i, j, k)
            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

        end
    end

end

function getdEdp_p!(sk, dEdp)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_p(sk, i, j, k)
            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

        end
    end

end

function getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, mpi, i, j, k)

    Aj = getAj(dp, ddp1, ddp2)
    Bj = getBj(dp)

    @inbounds for a = 1:4
        dEdp[i, j, k, a] =
            Aj[1]*dp[1, a] +
            Aj[2]*dp[2, a] +
            Aj[3]*dp[3, a] +
            Bj[1]*ddp1[1, a] +
            Bj[2]*ddp1[2, a] +
            Bj[3]*ddp1[3, a] +
            Bj[4]*ddp2[1, a] +
            Bj[5]*ddp2[2, a] +
            Bj[6]*ddp2[3, a]
    end
    dEdp[i, j, k, 4] += mpi^2

    @inbounds DEdotpion_field =
        dEdp[i, j, k, 1]*p[1] +
        dEdp[i, j, k, 2]*p[2] +
        dEdp[i, j, k, 3]*p[3] +
        dEdp[i, j, k, 4]*p[4]

    for a = 1:4
        dEdp[i, j, k, a] -= p[a]*DEdotpion_field
    end

end

function getAj(dp, ddp1, ddp2)

    return SVector{3,Float64}(
        -((ddp1[2, 1] + ddp1[3, 1])*dp[1, 1]) - (ddp1[2, 2] + ddp1[3, 2])*dp[1, 2] -
        (ddp1[2, 3] + ddp1[3, 3])*dp[1, 3] - (ddp1[2, 4] + ddp1[3, 4])*dp[1, 4] +
        ddp2[3, 1]*dp[2, 1] +
        ddp2[3, 2]*dp[2, 2] +
        ddp2[3, 3]*dp[2, 3] +
        ddp2[3, 4]*dp[2, 4] +
        ddp2[2, 1]*dp[3, 1] +
        ddp2[2, 2]*dp[3, 2] +
        ddp2[2, 3]*dp[3, 3] +
        ddp2[2, 4]*dp[3, 4],
        ddp2[3, 1]*dp[1, 1] +
        ddp2[3, 2]*dp[1, 2] +
        ddp2[3, 3]*dp[1, 3] +
        ddp2[3, 4]*dp[1, 4] - (ddp1[1, 1] + ddp1[3, 1])*dp[2, 1] -
        (ddp1[1, 2] + ddp1[3, 2])*dp[2, 2] - (ddp1[1, 3] + ddp1[3, 3])*dp[2, 3] -
        (ddp1[1, 4] + ddp1[3, 4])*dp[2, 4] +
        ddp2[1, 1]*dp[3, 1] +
        ddp2[1, 2]*dp[3, 2] +
        ddp2[1, 3]*dp[3, 3] +
        ddp2[1, 4]*dp[3, 4],
        ddp2[2, 1]*dp[1, 1] +
        ddp2[2, 2]*dp[1, 2] +
        ddp2[2, 3]*dp[1, 3] +
        ddp2[2, 4]*dp[1, 4] +
        ddp2[1, 1]*dp[2, 1] +
        ddp2[1, 2]*dp[2, 2] +
        ddp2[1, 3]*dp[2, 3] +
        ddp2[1, 4]*dp[2, 4] - (ddp1[1, 1] + ddp1[2, 1])*dp[3, 1] -
        (ddp1[1, 2] + ddp1[2, 2])*dp[3, 2] - (ddp1[1, 3] + ddp1[2, 3])*dp[3, 3] -
        (ddp1[1, 4] + ddp1[2, 4])*dp[3, 4],
    )

end

function getBj(dp)

    return SVector{6,Float64}(
        1 +
        dp[2, 1]^2 +
        dp[2, 2]^2 +
        dp[2, 3]^2 +
        dp[2, 4]^2 +
        dp[3, 1]^2 +
        dp[3, 2]^2 +
        dp[3, 3]^2 +
        dp[3, 4]^2,
        1 +
        dp[1, 1]^2 +
        dp[1, 2]^2 +
        dp[1, 3]^2 +
        dp[1, 4]^2 +
        dp[3, 1]^2 +
        dp[3, 2]^2 +
        dp[3, 3]^2 +
        dp[3, 4]^2,
        1 +
        dp[1, 1]^2 +
        dp[1, 2]^2 +
        dp[1, 3]^2 +
        dp[1, 4]^2 +
        dp[2, 1]^2 +
        dp[2, 2]^2 +
        dp[2, 3]^2 +
        dp[2, 4]^2,
        -2*dp[2, 1]*dp[3, 1] - 2*dp[2, 2]*dp[3, 2] - 2*dp[2, 3]*dp[3, 3] -
        2*dp[2, 4]*dp[3, 4],
        -2*dp[1, 1]*dp[3, 1] - 2*dp[1, 2]*dp[3, 2] - 2*dp[1, 3]*dp[3, 3] -
        2*dp[1, 4]*dp[3, 4],
        -2*dp[1, 1]*dp[2, 1] - 2*dp[1, 2]*dp[2, 2] - 2*dp[1, 3]*dp[2, 3] -
        2*dp[1, 4]*dp[2, 4],
    )

end

function EnergyANF(sk, ED)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            dp = getDP(sk, i, j, k)
            ED[i, j, k] = engpt(dp, sk.pion_field[i, j, k, 4], sk.mpi)

        end
    end

    return sum(ED)

end



"""
    arrested_newton_flow!(skyrmion; skyrmion_dot = zero_array, steps = 1, tolerance = 0.0, dt = ls/10.0, checks = max(100, steps), verbose = true, method = "RK4")
    
Applies an arrested Newton flow to `skyrmion` whose initial time derivative field is `skyrmion_dot` (an array of the correct shape initialised with value 0.0) with timestep `dt` (where `ls` is the lattice spacing in the first direction), either for `steps` steps or until the error falls below `tolerance`. The error is checked every `checks` steps.

The flow is verbose, incrementally describing the error and new energy, if `verbose` is true. `method` determines how each timestep is carried out: accepted values are "RK4" or "leapfrog". 

See also [`gradient_flow!`](@ref). 

"""
function arrested_newton_flow!(
    ϕ;
    ϕd = zeros(ϕ.grid.lp[1], ϕ.grid.lp[2], ϕ.grid.lp[3], 4),
    dt = ϕ.grid.ls[1]/10.0,
    steps = 1,
    tolerance = 0.0,
    checks = max(100, steps),
    verbose = true,
    method = "RK4"
)
    tolerance < 0 && @warn "Tolerance is negative"

    if tolerance == 0 && checks > steps
        checks = steps
    end

    energy_density = zeros(ϕ.grid.lp[1], ϕ.grid.lp[2], ϕ.grid.lp[3])
    old_pion_field = deepcopy(ϕ.pion_field);

    dEdp = zeros(ϕ.grid.lp[1], ϕ.grid.lp[2], ϕ.grid.lp[3], 4)
    dEdp2 = zeros(ϕ.grid.lp[1], ϕ.grid.lp[2], ϕ.grid.lp[3], 4)
    dEdp3 = zeros(ϕ.grid.lp[1], ϕ.grid.lp[2], ϕ.grid.lp[3], 4)
    dEdp4 = zeros(ϕ.grid.lp[1], ϕ.grid.lp[2], ϕ.grid.lp[3], 4)
    sk2 = deepcopy(ϕ)


    counter = 0
    while counter < steps

        arrested_newton_flow_for_n_steps!(
            ϕ,
            sk2,
            ϕd,
            old_pion_field,
            dEdp,
            dEdp2,
            dEdp3,
            dEdp4,
            dt,
            energy_density,
            checks,
            EnergyANF(ϕ, energy_density),
            method,
        )
        error = max_abs_err(dEdp)
        counter += checks

        if verbose
            println(
                "after ",
                counter,
                " steps, error = ",
                round(error, sigdigits = 4),
                " energy = ",
                round(
                    sum(energy_density)*ϕ.grid.ls[1]*ϕ.grid.ls[2]*ϕ.grid.ls[3]/(12.0*pi^2),
                    sigdigits = 8,
                ),
            )
            #println( round(error, sigdigits=8), "," )
        end

        if tolerance != 0.0    # => we are in tol mode
            if error < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += checks
            end
        end

    end

    return

end

function arrested_newton_flow_for_n_steps!(
    ϕ,
    sk2,
    ϕd,
    old_pion_field,
    dEdp1,
    dEdp2,
    dEdp3,
    dEdp4,
    dt,
    energy_density,
    n,
    initial_energy,
    method,
)

    new_energy = initial_energy

    for _ = 1:n

        old_energy = new_energy
        old_pion_field .= deepcopy(ϕ.pion_field)

        if method == "RK4"
            newton_flow_for_1_step!(ϕ, sk2, ϕd, dEdp1, dEdp2, dEdp3, dEdp4, dt)
        elseif method == "leapfrog"
            leapfrog_for_1_step!(ϕ, ϕd, dEdp1, dEdp2, dt)
        else
            throw(ArgumentError(method, "method must be either RK4 or leapfrog"))
        end

        new_energy = EnergyANF(ϕ, energy_density)

        if new_energy > old_energy

            fill!(ϕd, 0.0);
            ϕ.pion_field .= deepcopy(old_pion_field);

            if new_energy > 1.2*old_energy
                error(
                    "Suspected numerical blow-up. Please use smaller dt. Currently, dt = ",
                    dt,
                )
            end

        end

    end

end


function leapfrog_for_1_step!(sk, skd, dEdp1, dEdp2, dt)

    getdEdp!(sk, dEdp1)
    sk.pion_field .+= dt .* (skd .+ (0.5*dt) .* dEdp1);

    getdEdp!(sk, dEdp2)
    skd .+= (0.5*dt) .* (dEdp1 .+ dEdp2)

    orthog_skd_and_norm!(sk, skd)

end

function orthog_skd_and_norm!(sk, skd)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            sk_dot_skd = 0.0
            sk_dot_sk = 0.0

            for a = 1:4
                sk_dot_skd += sk.pion_field[i, j, k, a]*skd[i, j, k, a]
                sk_dot_sk += sk.pion_field[i, j, k, a]^2
            end

            sk_dot_sq = sqrt(sk_dot_sk)

            for a = 1:4
                skd[i, j, k, a] -= sk_dot_skd*sk.pion_field[i, j, k, a]
                sk.pion_field[i, j, k, a] /= sk_dot_sq
            end

        end
    end

end



# The newton flow code sacrifices beauty for optimization. Each RK4 step updates the 
# fields ready for the next step, meaning we only need two fields in memory.
# Additionally,  these updates happen within threaded loops, so that
# we only rethread 4 times for an RK4 method. This means we need four seperate
# update functions, and four more for different boundary conditions.

function newton_flow_for_1_step!(sk, sk2, skd, dEdp1, dEdp2, dEdp3, dEdp4, dt)

    if sk.grid.boundary_conditions == "dirichlet"
        getdEdp1!(sk, dEdp1, sk2, skd, dt)
        getdEdp2!(sk2, dEdp2, sk, dEdp1, dt)
        getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)
        getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)
    else
        getdEdp1_p!(sk, dEdp1, sk2, skd, dt)
        getdEdp2_p!(sk2, dEdp2, sk, dEdp1, dt)
        getdEdp3_p!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)
        getdEdp4_p!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)
    end


end

function getdEdp1!(sk, dEdp, sk2, skd, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_np(sk, i, j, k)

            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a = 1:4
                sk2.pion_field[i, j, k, a] = p[a] + (0.5*dt)*skd[i, j, k, a]
            end

        end
    end

end

function getdEdp2!(sk2, dEdp2, sk, dEdp1, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_np(sk2, i, j, k)

            getdEdp_pt!(dEdp2, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a = 1:4
                sk.pion_field[i, j, k, a] = p[a] + (0.5*dt)^2*dEdp1[i, j, k, a]
            end

        end
    end

end

function getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_np(sk, i, j, k)

            getdEdp_pt!(dEdp3, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a = 1:4
                sk2.pion_field[i, j, k, a] =
                    p[a] +
                    (0.5*dt) .* skd[i, j, k, a] +
                    (0.5*dt)^2*(4.0*dEdp2[i, j, k, a] - dEdp1[i, j, k, a])

            end

        end
    end

end

function getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_np(sk2, i, j, k)
            getdEdp_pt!(dEdp4, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            skd_dot_sk = 0.0
            sk_dot_sk = 0.0

            for a = 1:4
                sk.pion_field[i, j, k, a] =
                    p[a] -
                    dt*(
                        (5/6*dt)*dEdp2[i, j, k, a] -
                        (dt/6)*(dEdp1[i, j, k, a] .+ dEdp3[i, j, k, a])
                    )
                skd[i, j, k, a] +=
                    (
                        dt/6.0
                    )*(
                        dEdp1[i, j, k, a] +
                        2.0*dEdp2[i, j, k, a] +
                        2.0*dEdp3[i, j, k, a] +
                        dEdp4[i, j, k, a]
                    )

                skd_dot_sk += skd[i, j, k, a]*sk.pion_field[i, j, k, a]
                sk_dot_sk += sk.pion_field[i, j, k, a]^2
            end

            sk_dot_sk /= sqrt(sk_dot_sk)
            for a = 1:4
                skd[i, j, k, a] -= skd_dot_sk*sk.pion_field[i, j, k, a]
                sk.pion_field[i, j, k, a] /= sk_dot_sk
            end

        end
    end

end


function getdEdp1_p!(sk, dEdp, sk2, skd, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_p(sk, i, j, k)

            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a = 1:4
                sk2.pion_field[i, j, k, a] = p[a] + (0.5*dt)*skd[i, j, k, a]
            end

        end
    end

end

function getdEdp2_p!(sk2, dEdp2, sk, dEdp1, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_p(sk2, i, j, k)

            getdEdp_pt!(dEdp2, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a = 1:4
                sk.pion_field[i, j, k, a] = p[a] + (0.5*dt)^2*dEdp1[i, j, k, a]
            end

        end
    end

end

function getdEdp3_p!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_p(sk, i, j, k)

            getdEdp_pt!(dEdp3, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a = 1:4
                sk2.pion_field[i, j, k, a] =
                    p[a] +
                    (0.5*dt) .* skd[i, j, k, a] +
                    (0.5*dt)^2*(4.0*dEdp2[i, j, k, a] - dEdp1[i, j, k, a])

            end

        end
    end

end

function getdEdp4_p!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)

    Threads.@threads for k in sk.grid.sum_grid[3]
        @fastmath @inbounds for j in sk.grid.sum_grid[2], i in sk.grid.sum_grid[1]

            p, dp, ddp1, ddp2 = getders_local_p(sk2, i, j, k)
            getdEdp_pt!(dEdp4, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            skd_dot_sk = 0.0
            sk_dot_sk = 0.0

            for a = 1:4
                sk.pion_field[i, j, k, a] =
                    p[a] -
                    dt*(
                        (5/6*dt)*dEdp2[i, j, k, a] -
                        (dt/6)*(dEdp1[i, j, k, a] .+ dEdp3[i, j, k, a])
                    )
                skd[i, j, k, a] +=
                    (
                        dt/6.0
                    )*(
                        dEdp1[i, j, k, a] +
                        2.0*dEdp2[i, j, k, a] +
                        2.0*dEdp3[i, j, k, a] +
                        dEdp4[i, j, k, a]
                    )

                skd_dot_sk += skd[i, j, k, a]*sk.pion_field[i, j, k, a]
                sk_dot_sk += sk.pion_field[i, j, k, a]^2
            end

            sk_dot_sk /= sqrt(sk_dot_sk)
            for a = 1:4
                skd[i, j, k, a] -= skd_dot_sk*sk.pion_field[i, j, k, a]
                sk.pion_field[i, j, k, a] /= sk_dot_sk
            end

        end
    end

end

function max_abs_err(A)

    return maximum(abs, A)

end
