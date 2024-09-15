using StaticArrays
"""
    gradient_flow!(skyrmion; steps = n, tolerance = tol, dt=ls^2/80.0, checks = freq, print_stuff = true)
    
Applies a gradient flow to `skyrmion` with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `checks` steps.

See also [`newton_flow!`, `arrested_newton_flow!`]

"""
function gradient_flow!(ϕ; steps = 1, dt=((ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3])^(2/3))/100.0, tolerance = 0.0, checks = max(100,steps), print_stuff = true, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4), max_steps = Inf )

    if tolerance == 0 && checks > steps
        checks = steps
    end
    
    if print_stuff == true
        println("initial: energy: ", Energy(ϕ) )

    end

    counter = 0
    prev_error = 1.0e9
    
    while counter < steps && counter < max_steps
        
        gradient_flow_for_n_steps!(ϕ,dEdp,checks,dt)
        
        err = max_abs_err(dEdp)
        if err > 3*prev_error
            error("Suspected numerical blowup. Please use a smaller dt. Currently, dt = ", dt)
        end
        prev_error = err

        counter += checks
        
        if print_stuff == true
            println("after ", counter, " steps, error = ", round(err, sigdigits=4))
        end

        if tolerance != 0.0    # => we are in tol mode    
            if err < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += checks    # => continue the while loop
            end
        end

    end

    if print_stuff == true
        println("final energy: ", Energy(ϕ) )
    end

    return

end

function gradient_flow_for_n_steps!(ϕ,dEdp,n,dt)
    for _ in 1:n
        gradient_flow_1_step!(ϕ,dEdp,dt)
    end
end

function gradient_flow_1_step!(sk, dEdp, dt)

    getdEdp!(sk, dEdp)
    sk.pion_field .+= dt.*dEdp;
    normer!(sk)
   
end 


# we split dEdp into Dirichlet and other options here, so that the if statement happens once, rather than inside a for loop.
# This creates some code duplication, for a ~10% performance boost.

function getdEdp!(sk, dEdp)
    if sk.boundary_conditions == "dirichlet"
        getdEdp_np!(sk, dEdp)
    else
        getdEdp_p!(sk, dEdp)
    end
end


function getdEdp_np!(sk, dEdp)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
                    
            p, dp, ddp1, ddp2 = getders_local_np(sk,i,j,k)
            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k, sk.metric)

        end
    end

end

function getdEdp_p!(sk, dEdp)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
                    
            p, dp, ddp1, ddp2 = getders_local_p(sk,i,j,k)
            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k, sk.metric)

        end
    end

end

function getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, mpi, i, j, k, alpha)

    c = (alpha)^2 - 1
    Aj = getAj(dp,ddp1,ddp2)
    Bj = getBj(dp)
    b_t = get_berger_grad_e2_star(p,dp,ddp1)
    c_t = get_berger_grad_e4_star(dp, ddp1, ddp2)

    @inbounds for a in 1:4
        dEdp[i,j,k,a] = Aj[1]*dp[1,a] + Aj[2]*dp[2,a] + Aj[3]*dp[3,a] + Bj[1]*ddp1[1,a] + Bj[2]*ddp1[2,a] + Bj[3]*ddp1[3,a] + Bj[4]*ddp2[1,a] + Bj[5]*ddp2[2,a] + Bj[6]*ddp2[3,a] -0.5*(alpha^2 -1)*b_t[a] -0.5*(alpha^2 -1)*c_t[a]
    end
    dEdp[i,j,k,3] += - c*p[3]*mpi^2
    dEdp[i,j,k,4] += mpi^2 


    @inbounds DEdotpion_field = dEdp[i,j,k,1]*p[1] + dEdp[i,j,k,2]*p[2] + dEdp[i,j,k,3]*p[3] + dEdp[i,j,k,4]*p[4]

    for a in 1:4
        dEdp[i,j,k,a] -= p[a]*DEdotpion_field
    end

end

function get_berger_grad_e2_star(p, dp, ddp1)

    p1, p2, p3, p4 = p
    dp11, dp12, dp13, dp14 = dp[1,1], dp[1,2], dp[1,3], dp[1,4]
    dp21, dp22, dp23, dp24 = dp[2,1], dp[2,2], dp[2,3], dp[2,4]
    dp31, dp32, dp33, dp34 = dp[3,1], dp[3,2], dp[3,3], dp[3,4]

    L3_1 = (p4*dp13 - p3*dp14 + p1*dp12 - p2*dp11)
    L3_2 = (p4*dp23 - p3*dp24 + p1*dp22 - p2*dp21)
    L3_3 = (p4*dp33 - p3*dp34 + p1*dp32 - p2*dp31)

    sqd_term = (p4*(ddp1[1,3] + ddp1[2,3] + ddp1[3,3])) - (p3*(ddp1[1,4] + ddp1[2,4] + ddp1[3,4])) + (p1*(ddp1[1,2] + ddp1[2,2] + ddp1[3,2])) - (p2*(ddp1[1,1] + ddp1[2,1] + ddp1[3,1]))

    t1 = @SVector [dp12, -dp11, -dp14, dp13]
    t2 = @SVector [dp22, -dp21, -dp24, dp23]
    t3 = @SVector [dp32, -dp31, -dp34, dp33]

    wp =  @SVector [p2, -p1, -p4, p3]

    result = (4*(L3_1*t1 + L3_2*t2 + L3_3*t3) + 2*sqd_term*wp)
    
    return result
end

function get_berger_grad_e4_star(dp, ddp1, ddp2)

    dp11, dp12, dp13, dp14 = dp[1,1], dp[1,2], dp[1,3], dp[1,4]
    dp21, dp22, dp23, dp24 = dp[2,1], dp[2,2], dp[2,3], dp[2,4]
    dp31, dp32, dp33, dp34 = dp[3,1], dp[3,2], dp[3,3], dp[3,4]

    ddp11, ddp12, ddp13, ddp14 = ddp1[1,1], ddp1[1,2], ddp1[1,3], ddp1[1,4]
    ddp21, ddp22, ddp23, ddp24 = ddp1[2,1], ddp1[2,2], ddp1[2,3], ddp1[2,4]
    ddp31, ddp32, ddp33, ddp34 = ddp1[3,1], ddp1[3,2], ddp1[3,3], ddp1[3,4]

    dm11, dm12, dm13, dm14 = ddp2[1,1], ddp2[1,2], ddp2[1,3], ddp2[1,4]
    dm21, dm22, dm23, dm24 = ddp2[2,1], ddp2[2,2], ddp2[2,3], ddp2[2,4]
    dm31, dm32, dm33, dm34 = ddp2[3,1], ddp2[3,2], ddp2[3,3], ddp2[3,4]

    s1 = @SVector [dp12, -dp11, -dp14, dp13]
    s2 = @SVector [dp22, -dp21, -dp24, dp23]
    s3 = @SVector [dp32, -dp31, -dp34, dp33]

    t1_1 = (ddp14 + ddp24 + ddp34)*dp13 - (ddp13 + ddp23 + ddp33)*dp14 + (ddp11 + ddp21 + ddp31)*dp12 - (ddp12 + ddp22 + ddp32)* dp11
    t1_2 = (ddp14 + ddp24 + ddp34)*dp23 - (ddp13 + ddp23 + ddp33)*dp24 + (ddp11 + ddp21 + ddp31)*dp22 - (ddp12 + ddp22 + ddp32)* dp21
    t1_3 = (ddp14 + ddp24 + ddp34)*dp33 - (ddp13 + ddp23 + ddp33)*dp34 + (ddp11 + ddp21 + ddp31)*dp32 - (ddp12 + ddp22 + ddp32)* dp31

    v = t1_1*s1 + t1_2*s2 + t1_3*s3

    s11 = (dp14*ddp13 - dp13*ddp14 + dp11*ddp12 - dp12*ddp11) # j = 1 , i = 1
    s21 = (dp24*dm33 - dp23*dm34 + dp21*dm32 - dp22*dm31) # j = 2 , i = 1
    s31 = (dp34*dm23 - dp33*dm24 + dp31*dm22 - dp32*dm21) # j = 3 , i = 1

    s12 = (dp14*dm33 - dp13*dm34 + dp11*dm32 -dp12*dm31) # j = 1 , i = 2
    s22 = (dp24*ddp23 - dp23*ddp24 + dp21*ddp22 - dp22*ddp21) # j = 2 , i = 2
    s32 = (dp34*dm13 - dp33*dm14 + dp31*dm12 -dp32*dm11) # j = 3 , i = 2

    s13 = (dp14*dm23 - dp13*dm24 + dp11*dm22 - dp12*dm21) # j = 1 , i = 3
    s23 = (dp24*dm13 - dp23*dm14 + dp21*dm12 - dp22*dm11) # j = 2 , i = 3
    s33 = (dp34*ddp33 - dp33*ddp34 +dp31*ddp32 -dp32*ddp31) # j = 3, i = 3

    w = (s11 + s21 + s31)*s1+ (s12 + s22 + s32)*s2 + (s13 + s23 + s33)*s3

    return -2*(v+w)

end

function getAj(dp, ddp1, ddp2)

    return SVector{3,Float64}(
        -((ddp1[2,1] + ddp1[3,1])*dp[1,1]) - (ddp1[2,2] + ddp1[3,2])*dp[1,2] - (ddp1[2,3] + ddp1[3,3])*dp[1,3] - (ddp1[2,4] + ddp1[3,4])*dp[1,4] + ddp2[3,1]*dp[2,1] + ddp2[3,2]*dp[2,2] + ddp2[3,3]*dp[2,3] + ddp2[3,4]*dp[2,4] + ddp2[2,1]*dp[3,1] + ddp2[2,2]*dp[3,2] + ddp2[2,3]*dp[3,3] + ddp2[2,4]*dp[3,4],
         ddp2[3,1]*dp[1,1] + ddp2[3,2]*dp[1,2] + ddp2[3,3]*dp[1,3] + ddp2[3,4]*dp[1,4] - (ddp1[1,1] + ddp1[3,1])*dp[2,1] - (ddp1[1,2] + ddp1[3,2])*dp[2,2] - (ddp1[1,3] + ddp1[3,3])*dp[2,3] - (ddp1[1,4] + ddp1[3,4])*dp[2,4] + ddp2[1,1]*dp[3,1] + ddp2[1,2]*dp[3,2] + ddp2[1,3]*dp[3,3] + ddp2[1,4]*dp[3,4],
         ddp2[2,1]*dp[1,1] + ddp2[2,2]*dp[1,2] + ddp2[2,3]*dp[1,3] + ddp2[2,4]*dp[1,4] + ddp2[1,1]*dp[2,1] + ddp2[1,2]*dp[2,2] + ddp2[1,3]*dp[2,3] + ddp2[1,4]*dp[2,4] - (ddp1[1,1] + ddp1[2,1])*dp[3,1] - (ddp1[1,2] + ddp1[2,2])*dp[3,2] - (ddp1[1,3] + ddp1[2,3])*dp[3,3] - (ddp1[1,4] + ddp1[2,4])*dp[3,4]
    )
    
end

function getBj(dp)

    return SVector{6,Float64}(
     1 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2 + dp[3,4]^2,
        1 + dp[1,1]^2 + dp[1,2]^2 + dp[1,3]^2 + dp[1,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2 + dp[3,4]^2,
        1 + dp[1,1]^2 + dp[1,2]^2 + dp[1,3]^2 + dp[1,4]^2 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2 + dp[2,4]^2,
        -2*dp[2,1]*dp[3,1] - 2*dp[2,2]*dp[3,2] - 2*dp[2,3]*dp[3,3] - 2*dp[2,4]*dp[3,4],
        -2*dp[1,1]*dp[3,1] - 2*dp[1,2]*dp[3,2] - 2*dp[1,3]*dp[3,3] - 2*dp[1,4]*dp[3,4],
        -2*dp[1,1]*dp[2,1] - 2*dp[1,2]*dp[2,2] - 2*dp[1,3]*dp[2,3] - 2*dp[1,4]*dp[2,4]
    )

end

function EnergyANF(sk, ED)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]

            dp = getDP(sk ,i, j, k )
            ED[i,j,k] = engpt(dp,sk.pion_field[i,j,k,:],sk.mpi, sk.metric)

        end
    end    

    return sum(ED)

end 



"""
    arrested_newton_flow!(skyrmion; skyrmion_dot, steps = n, tolerance = tol, dt=ls^2/80.0, checks = freq, print_stuff = true)
    
Applies an arrested Newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `checks` steps.

See also [`gradient_flow!`, `newton_flow!`]
"""
function arrested_newton_flow!(ϕ; ϕd=zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4), dt=ϕ.ls[1]/10.0, steps=1, tolerance = 0.0, checks = max(100,steps), print_stuff = true, max_steps=Inf, method="RK4")

    if tolerance == 0 && checks > steps
        checks = steps
    end

    energy_density = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3])
    old_pion_field = deepcopy(ϕ.pion_field);

    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp2 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp3 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp4 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    sk2 = deepcopy(ϕ)


    counter = 0
    while counter < steps && counter < max_steps

        arrested_newton_flow_for_n_steps!(ϕ,sk2,ϕd,old_pion_field,dEdp,dEdp2,dEdp3,dEdp4,dt,energy_density,checks, EnergyANF(ϕ,energy_density), method)
        error = max_abs_err(dEdp)
        counter += checks

        if print_stuff == true 
            println("after ", counter, " steps, error = ", round(error, sigdigits=4), " energy = ", round(sum(energy_density)*ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3]/(12.0*pi^2), sigdigits=8) )
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

function arrested_newton_flow_for_n_steps!(ϕ,sk2,ϕd,old_pion_field,dEdp1,dEdp2,dEdp3,dEdp4,dt,energy_density,n, initial_energy, method)
    alpha = ϕ.metric
    new_energy = initial_energy
    
    for _ in 1:n

        old_energy = new_energy
        old_pion_field .= deepcopy(ϕ.pion_field)

        if method == "RK4"
            newton_flow_for_1_step!(ϕ,sk2,ϕd,dEdp1,dEdp2,dEdp3,dEdp4,dt,alpha)
        elseif method == "leapfrog"
            leapfrog_for_1_step!(ϕ,ϕd,dEdp1,dEdp2,dt)
        end

        new_energy = EnergyANF(ϕ,energy_density)

        if new_energy > old_energy

            fill!(ϕd, 0.0);
            ϕ.pion_field .= deepcopy(old_pion_field);

            if new_energy > 1.2*old_energy
                error("Suspected numerical blow-up. Please use smaller dt. Currently, dt = ", dt)
            end
    
        end

    end

end


function leapfrog_for_1_step!(sk,skd,dEdp1,dEdp2,dt)

    getdEdp!(sk, dEdp1)
    sk.pion_field .+= dt.*(skd .+ (0.5*dt).*dEdp1) ;

    getdEdp!(sk, dEdp2)
    skd .+= (0.5*dt).*(dEdp1 .+ dEdp2)

    orthog_skd_and_norm!(sk,skd)

end

function orthog_skd_and_norm!(sk, skd)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]

            sk_dot_skd = 0.0
            sk_dot_sk = 0.0

            for a in 1:4
                sk_dot_skd += sk.pion_field[i,j,k,a]*skd[i,j,k,a]
                sk_dot_sk += sk.pion_field[i,j,k,a]^2
            end

            sk_dot_sq = sqrt(sk_dot_sk)

            for a in 1:4
                skd[i,j,k,a] -= sk_dot_skd*sk.pion_field[i,j,k,a]
                sk.pion_field[i,j,k,a] /= sk_dot_sq
            end

        end
    end

end
                


# The newton flow code sacrifices beauty for optimization. Each RK4 step updates the 
# fields ready for the next step, meaning we only need two fields in memory.
# Additionally,  these updates happen within threaded loops, so that
# we only rethread 4 times for an RK4 method. This means we need four seperate
# update functions, and four more for different boundary conditions.

function newton_flow_for_1_step!(sk, sk2, skd ,dEdp1, dEdp2, dEdp3, dEdp4, dt, alpha)

    if sk.boundary_conditions == "dirichlet"
        getdEdp1!(sk, dEdp1, sk2, skd, dt, alpha)
        getdEdp2!(sk2, dEdp2, sk, dEdp1, dt, alpha)
        getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt, alpha)
        getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt, alpha)
    else
        getdEdp1_p!(sk, dEdp1, sk2, skd, dt, alpha)
        getdEdp2_p!(sk2, dEdp2, sk, dEdp1, dt, alpha)
        getdEdp3_p!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt, alpha)
        getdEdp4_p!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt, alpha)
    end

   
end

function getdEdp1!(sk, dEdp, sk2, skd, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_np(sk,i,j,k)

            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k, alpha)

            for a in 1:4
                sk2.pion_field[i,j,k,a] = p[a] + (0.5*dt)*skd[i,j,k,a]
            end

        end
    end

end

function getdEdp2!(sk2, dEdp2, sk, dEdp1, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_np(sk2,i,j,k)

            getdEdp_pt!(dEdp2, p, dp, ddp1, ddp2, sk.mpi, i, j, k, alpha)

            for a in 1:4
                sk.pion_field[i,j,k,a] = p[a] + (0.5*dt)^2*dEdp1[i,j,k,a]
            end

        end
    end

end

function getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_np(sk,i,j,k)

            getdEdp_pt!(dEdp3, p, dp, ddp1, ddp2, sk.mpi, i, j, k, alpha)

            for a in 1:4
                sk2.pion_field[i,j,k,a] = p[a] + (0.5*dt).*skd[i,j,k,a] + (0.5*dt)^2*(4.0*dEdp2[i,j,k,a] - dEdp1[i,j,k,a])
                
            end

        end
    end

end

function getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_np(sk2,i,j,k)
            getdEdp_pt!(dEdp4, p, dp, ddp1, ddp2, sk.mpi, i, j, k,alpha)

            skd_dot_sk = 0.0
            sk_dot_sk = 0.0

            for a in 1:4
                sk.pion_field[i,j,k,a] = p[a] - dt*( (5/6*dt)*dEdp2[i,j,k,a]  - (dt/6)*( dEdp1[i,j,k,a] .+ dEdp3[i,j,k,a] ) ) 
                skd[i,j,k,a] += (dt/6.0)*(dEdp1[i,j,k,a] + 2.0*dEdp2[i,j,k,a] + 2.0*dEdp3[i,j,k,a] + dEdp4[i,j,k,a])

                skd_dot_sk += skd[i,j,k,a]*sk.pion_field[i,j,k,a]
                sk_dot_sk += sk.pion_field[i,j,k,a]^2
            end

            sk_dot_sk /= sqrt( sk_dot_sk) 
            for a in 1:4
                skd[i,j,k,a] -=  skd_dot_sk*sk.pion_field[i,j,k,a]
                sk.pion_field[i,j,k,a] /=  sk_dot_sk 
            end

        end
    end

end


function getdEdp1_p!(sk, dEdp, sk2, skd, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_p(sk,i,j,k)

            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k, alpha)

            for a in 1:4
                sk2.pion_field[i,j,k,a] = p[a] + (0.5*dt)*skd[i,j,k,a]
            end

        end
    end

end

function getdEdp2_p!(sk2, dEdp2, sk, dEdp1, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_p(sk2,i,j,k)

            getdEdp_pt!(dEdp2, p, dp, ddp1, ddp2, sk.mpi, i, j, k, alpha)

            for a in 1:4
                sk.pion_field[i,j,k,a] = p[a] + (0.5*dt)^2*dEdp1[i,j,k,a]
            end

        end
    end

end

function getdEdp3_p!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_p(sk,i,j,k)

            getdEdp_pt!(dEdp3, p, dp, ddp1, ddp2, sk.mpi, i, j, k, alpha)

            for a in 1:4
                sk2.pion_field[i,j,k,a] = p[a] + (0.5*dt).*skd[i,j,k,a] + (0.5*dt)^2*(4.0*dEdp2[i,j,k,a] - dEdp1[i,j,k,a])
                
            end

        end
    end

end

function getdEdp4_p!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt, alpha)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local_p(sk2,i,j,k)
            getdEdp_pt!(dEdp4, p, dp, ddp1, ddp2, sk.mpi, i, j, k,alpha)

            skd_dot_sk = 0.0
            sk_dot_sk = 0.0

            for a in 1:4
                sk.pion_field[i,j,k,a] = p[a] - dt*( (5/6*dt)*dEdp2[i,j,k,a]  - (dt/6)*( dEdp1[i,j,k,a] .+ dEdp3[i,j,k,a] ) ) 
                skd[i,j,k,a] += (dt/6.0)*(dEdp1[i,j,k,a] + 2.0*dEdp2[i,j,k,a] + 2.0*dEdp3[i,j,k,a] + dEdp4[i,j,k,a])

                skd_dot_sk += skd[i,j,k,a]*sk.pion_field[i,j,k,a]
                sk_dot_sk += sk.pion_field[i,j,k,a]^2
            end

            sk_dot_sk /= sqrt( sk_dot_sk) 
            for a in 1:4
                skd[i,j,k,a] -=  skd_dot_sk*sk.pion_field[i,j,k,a]
                sk.pion_field[i,j,k,a] /=  sk_dot_sk 
            end

        end
    end

end





"""
    newton_flow!(skyrmion; skyrmion_dot, steps = n, dt=ls^2/80.0, frequency_of_printing = freq, print_stuff = true)
    
Applies a newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The energy is checked every `freq` steps.

See also [`gradient_flow!`, `arrested_newton_flow!`]
"""
function newton_flow!(ϕ; ϕd=zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4), dt=ϕ.ls[1]/20.0, steps=1, print_stuff = true, frequency_of_printing = steps)
    alpha = ϕ.metric
    if print_stuff == true
        println("intial energy: ", Energy(ϕ))
    end

    counter = 0
    while counter < steps

        newton_flow_for_n_steps!(ϕ,ϕd,dt,frequency_of_printing,alpha)
        counter += frequency_of_printing

        if print_stuff == true 
            println("after ", counter, " steps, energy = ", Energy(ϕ) )
        end

    end

    return

end

function newton_flow_for_n_steps!(ϕ,ϕd,dt,n,alpha)

    dEdp1 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp2 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp3 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp4 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    sk2 = deepcopy(ϕ)

    for _ in 1:n
        newton_flow_for_1_step!(ϕ,sk2,ϕd,dEdp1,dEdp2,dEdp3,dEdp4,dt,alpha)
    end

end

function max_abs_err(A)

    return maximum(abs, A)

end


function lin_interpolate(U, point)
    phi = U.pion_field

    x = point[1]
    y = point[2]
    z = point[3]

    x0 = floor(Int64,x)
    x1 = ceil(Int64,x)
    y0 = floor(Int64,y)
    y1 = ceil(Int64,y)
    z0 = floor(Int64,z)
    z1 = ceil(Int64,z)

    f_000 = phi[x0, y0, z0, :]
    f_100 = phi[x1, y0, z0, :]
    f_010 = phi[x0, y1, z0, :]
    f_001 = phi[x0, y0, z1, :]
    f_110 = phi[x1, y1, z0, :]
    f_101 = phi[x1, y0, z1, :]
    f_011 = phi[x0, y1, z1, :]
    f_111 = phi[x1, y1, z1, :]

    if x1 == x0
        x_d = 0.0
    else
        x_d = (x - x0) / (x1 - x0)
    end

    if y1 == y0
        y_d = 0.0
    else
        y_d = (y - y0) / (y1 - y0)
    end

    if z1 == z0
        z_d = 0.0
    else
        z_d = (z - z0) / (z1 - z0)
    end

    c00 = f_000 * (1 - x_d) + f_100 * x_d
    c01 = f_001 * (1 - x_d) + f_101 * x_d
    c10 = f_010 * (1 - x_d) + f_110 * x_d
    c11 = f_011 * (1 - x_d) + f_111 * x_d

    c0 = c00 * (1 - y_d) + c10 * y_d
    c1 = c01 * (1 - y_d) + c11 * y_d

    c = c0 * (1 - z_d) + c1 * z_d

    cdotphi = c[1]*phi[1] + c[2]*phi[2] + c[3]*phi[3] + c[4]*phi[4]

    for i in 1:4
        c[i] -= phi[i]*cdotphi
    end

    return c 

end

