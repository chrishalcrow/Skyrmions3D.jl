"""
    gradient_flow!(skyrmion; steps = n, tolerance = tol, dt=ls^2/80.0, frequency_of_checking_tolerance = freq, print_stuff = true)
    
Applies a gradient flow to `skyrmion` with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `freq` steps.

See also [`newton_flow!`, `arrested_newton_flow!`]

"""
function gradient_flow!(ϕ; steps = 1, dt=((ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3])^(2/3))/100.0, tolerance = 0.0, frequency_of_checking_tolerance = max(100,steps), print_stuff = true, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4) , error_function::Function=L2_err)

    if tolerance == 0 && frequency_of_checking_tolerance > steps
        frequency_of_checking_tolerance = steps
    end
    
    if print_stuff == true
        println(dt)
        println("initial: energy: ", Energy(ϕ) )

    end

    counter = 0
    prev_error = 1.0e9
    
    while counter < steps
        
        gradient_flow_for_n_steps!(ϕ,dEdp,frequency_of_checking_tolerance,dt)
        
        err = max_abs_err(dEdp)
        if err > 3*prev_error
            error("Suspected numerical blowup. Please use a smaller dt. Currently, dt = ", dt)
        end
        prev_error = err

        counter += frequency_of_checking_tolerance
        
        if print_stuff == true
            #println("after ", counter, " steps, error = ", round(error, sigdigits=4))
            println( round(err, sigdigits=8), "," )
        end

        if tolerance != 0.0    # => we are in tol mode    
            if err < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += frequency_of_checking_tolerance    # => continue the while loop
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

function getdEdp!(sk, dEdp)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            #p = getX(sk,i,j,k)
            #dp = getDP(sk ,i, j, k )
            #ddp1 = getDDX1(sk, i, j, k)
            #ddp2 = getDDX2(sk, i, j, k, ddp1)
            p, dp, ddp1, ddp2 = getders_local(sk,i,j,k)

            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

        end
    end

end

function getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, mpi, i, j, k)

    Aj = getAj(dp,ddp1,ddp2)
    Bj = getBj(dp)

    for a in 1:4
        dEdp[i,j,k,a] = Aj[1]*dp[1,a] + Aj[2]*dp[2,a] + Aj[3]*dp[3,a] + Bj[1]*ddp1[1,a] + Bj[2]*ddp1[2,a] + Bj[3]*ddp1[3,a] + Bj[4]*ddp2[1,a] + Bj[5]*ddp2[2,a] + Bj[6]*ddp2[3,a]
    end
    dEdp[i,j,k,4] += mpi^2

    @inbounds DEdotpion_field = dEdp[i,j,k,1]*p[1] + dEdp[i,j,k,2]*p[2] + dEdp[i,j,k,3]*p[3] + dEdp[i,j,k,4]*p[4]

    for a in 1:4
        dEdp[i,j,k,a] -= p[a]*DEdotpion_field
    end

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


"""
    arrested_newton_flow!(skyrmion, skyrmion_dot; steps = n, tolerance = tol, dt=ls^2/80.0, frequency_of_checking_tolerance = freq, print_stuff = true)
    
Applies an arrested Newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `freq` steps.

See also [`gradient_flow!`, `newton_flow!`]
"""
function arrested_newton_flow!(ϕ,ϕd; dt=ϕ.ls[1]/5.0, steps=1, tolerance = 0.0, frequency_of_checking_tolerance = max(100,steps), print_stuff = true, step_algorithm="RK4")

    if tolerance == 0 && frequency_of_checking_tolerance > steps
        frequency_of_checking_tolerance = steps
    end

    energy_density = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3])
    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    old_pion_field = deepcopy(ϕ.pion_field);

    counter = 0
    while counter < steps

        arrested_newton_flow_for_n_steps!(ϕ,ϕd,old_pion_field,dEdp,dt,energy_density,frequency_of_checking_tolerance,step_algorithm,initial_energy=EnergyANF(ϕ,energy_density))
        error = max_abs_err(dEdp)
        counter += frequency_of_checking_tolerance

        if print_stuff == true 
            println("after ", counter, " steps, error = ", round(error, sigdigits=4), " energy = ", round(sum(energy_density)*ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3]/(12.0*pi^2), sigdigits=8) )
        end

        if tolerance != 0.0    # => we are in tol mode
            if error < tolerance 
                counter = steps + 1    # => end the while loop
            else
                steps += frequency_of_checking_tolerance
            end
        end

    end

    return

end
 
function arrested_newton_flow_for_n_steps!(ϕ,ϕd,old_pion_field,dEdp1,dt,energy_density,n,step_algorithm::String; initial_energy)

    new_energy = initial_energy
    
    if step_algorithm == "RK4"
        dEdp2 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
        dEdp3 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
        dEdp4 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    end

    for _ in 1:n

        old_energy = new_energy
        old_pion_field .= ϕ.pion_field

        if step_algorithm == "Euler"
            newton_flow_for_1_step!(ϕ,ϕd,dEdp,dt)
        elseif step_algorithm == "RK4"
            newton_flow_for_1_step_RK4!(ϕ,ϕd,dEdp1,dEdp2,dEdp3,dEdp4,dt)
        end

        new_energy = EnergyANF(ϕ,energy_density)

        if new_energy > old_energy

            fill!(ϕd, 0.0);
            ϕ.pion_field .= old_pion_field;

            
    
        end

    end

end



function newton_flow_for_1_step_RK4!(sk, skd ,dEdp1, dEdp2, dEdp3, dEdp4, dt)

    getdEdp!(sk, dEdp1)
    # sk2 = sk + (0.5*dt).*skd
    sk.pion_field .+= (0.5*dt).*skd

    getdEdp!(sk, dEdp2)
    #getdEdp!(sk2, dEdp2), sk = sk2 + (0.5*dt)^2 .*dEdp1
    sk.pion_field .+= (0.5*dt)^2 .*dEdp1

    getdEdp!(sk, dEdp3)
    #getdEdp!(sk, dEdp2), sk2 = sk + (0.5*dt).*skd .+ (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)
    sk.pion_field .+= (0.5*dt).*skd .+ (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)


    getdEdp!(sk, dEdp4)
    # #getdEdp!(sk2, dEdp2), sk +=  dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )
    #RESET field to original value: sk.pion_field .-= (0.5*dt).*skd + (0.5*dt)^2 *(4.0 .* dEdp2 - dEdp1) + (0.5*dt)^2 .*dEdp1 + (0.5*dt).*skd
    #Then update: sk.pion_field .+= dt.*(skd + dt/6.0 .*( dEdp1 .+ dEdp2 .+ dEdp3 ) ), combined into:
    sk.pion_field .-=  dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )
    skd .+= (dt/6.0).*(dEdp1 .+ 2.0.*dEdp2 .+ 2.0.*dEdp3 .+ dEdp4)
   
    #orthog_skd_and_sk!(skd,sk)
    #normer!(sk)
    orthog_skd_and_sk_and_normer!(skd,sk)
   
end


function newton_flow_for_1_step!(sk, skd ,dEdp, dt)

    getdEdp!(sk, dEdp)
    sk.pion_field .+= dt.*skd
    skd .+= dt.*dEdp
    normer!(sk)
   
end

function orthog_skd_and_sk!(skd,sk)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]

            skd_dot_sk = 0.0
            for a in 1:4
                skd_dot_sk += skd[i,j,k,a]*sk.pion_field[i,j,k,a]
            end

            for a in 1:4
                skd[i,j,k,a] -=  skd_dot_sk*sk.pion_field[i,j,k,a]
            end

        end
    end

end

function orthog_skd_and_sk_and_normer!(skd,sk)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]

            skd_dot_sk = 0.0
            sk_dot_sk = 0.0
            for a in 1:4
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

function EnergyANF(sk, ED)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
    
            dp = getDP(sk ,i, j, k )
            
            ED[i,j,k] = engpt(dp,sk.pion_field[i,j,k,4],sk.mpi)

        end
    end    
        
    return sum(ED)

end 











function arrested_newton_flow_juggle!(ϕ,ϕd; dt=ϕ.ls[1]/10.0, steps=1, tolerance = 0.0, frequency_of_checking_tolerance = max(100,steps), print_stuff = true)

    println(dt)

    if tolerance == 0 && frequency_of_checking_tolerance > steps
        frequency_of_checking_tolerance = steps
    end

    energy_density = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3])
    old_pion_field = deepcopy(ϕ.pion_field);

    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp2 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp3 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dEdp4 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    sk2 = Skyrmion(ϕ.lp, ϕ.ls);


    counter = 0
    while counter < steps

        arrested_newton_flow_for_n_steps_juggle!(ϕ,sk2,ϕd,old_pion_field,dEdp,dEdp2,dEdp3,dEdp4,dt,energy_density,frequency_of_checking_tolerance, EnergyANF(ϕ,energy_density))
        error = max_abs_err(dEdp)
        counter += frequency_of_checking_tolerance

        if print_stuff == true 
            println("after ", counter, " steps, error = ", round(error, sigdigits=4), " energy = ", round(sum(energy_density)*ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3]/(12.0*pi^2), sigdigits=8) )
            #println( round(error, sigdigits=8), "," )
        end

        if tolerance != 0.0    # => we are in tol mode
            if error < tolerance 
                counter = steps + 1    # => end the while loop
            else
                steps += frequency_of_checking_tolerance
            end
        end

    end

    return

end
 
function arrested_newton_flow_for_n_steps_juggle!(ϕ,sk2,ϕd,old_pion_field,dEdp1,dEdp2,dEdp3,dEdp4,dt,energy_density,n, initial_energy)

    new_energy = initial_energy
    
    for _ in 1:n

        old_energy = new_energy
        old_pion_field .= ϕ.pion_field

        newton_flow_for_1_step_juggle!(ϕ,sk2,ϕd,dEdp1,dEdp2,dEdp3,dEdp4,dt)
        new_energy = EnergyANF(ϕ,energy_density)

        if new_energy > old_energy

            fill!(ϕd, 0.0);
            ϕ.pion_field .= old_pion_field;

            if new_energy > 1.2*old_energy
                error("Suspected numerical blow-up. Please use smaller dt. Currently, dt = ", dt)
            end
    
        end

    end

end

function newton_flow_for_1_step_juggle!(sk, sk2, skd ,dEdp1, dEdp2, dEdp3, dEdp4, dt)

    getdEdp1!(sk, dEdp1, sk2, skd, dt)
    # sk2 = sk + (0.5*dt).*skd

    #sk.pion_field .+= (0.5*dt).*skd
    getdEdp2!(sk2, dEdp2, sk, dEdp1, dt)
    #getdEdp!(sk2, dEdp2), sk = sk2 + (0.5*dt)^2 .*dEdp1

    ##sk.pion_field .+= (0.5*dt)^2 .*dEdp1
    getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)
    #getdEdp!(sk, dEdp3), sk2 = sk + (0.5*dt).*skd .+ (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)

    ##sk.pion_field .+= (0.5*dt).*skd .+ (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)
    getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)
    # #getdEdp!(sk2, dEdp4), sk += sk2 - dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )  

    #RESET field to original value: sk.pion_field .-= (0.5*dt).*skd + (0.5*dt)^2 *(4.0 .* dEdp2 - dEdp1) + (0.5*dt)^2 .*dEdp1 + (0.5*dt).*skd
    #Then update: sk.pion_field .+= dt.*(skd + dt/6.0 .*( dEdp1 .+ dEdp2 .+ dEdp3 ) ), combined into:
    ###sk.pion_field .-=  dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )
    
    #
   
    #orthog_skd_and_sk_and_normer!(skd,sk)
    #normer!(sk)

    ##orthog_skd_and_sk_and_normer!(skd,sk)
   
end


function getdEdp1!(sk, dEdp, sk2, skd, dt)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local(sk,i,j,k)

            getdEdp_pt!(dEdp, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a in 1:4
                sk2.pion_field[i,j,k,a] = p[a] + (0.5*dt)*skd[i,j,k,a]
            end

        end
    end

end

function getdEdp2!(sk2, dEdp2, sk, dEdp1, dt)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local(sk2,i,j,k)

            getdEdp_pt!(dEdp2, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a in 1:4
                sk.pion_field[i,j,k,a] = p[a] + (0.5*dt)^2*dEdp1[i,j,k,a]
            end

        end
    end

end

function getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local(sk,i,j,k)

            getdEdp_pt!(dEdp3, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

            for a in 1:4
                sk2.pion_field[i,j,k,a] = p[a] + (0.5*dt).*skd[i,j,k,a] + (0.5*dt)^2*(4.0*dEdp2[i,j,k,a] - dEdp1[i,j,k,a])
                
            end

        end
    end

end

function getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)

    Threads.@threads for k in sk.sum_grid[3]
        @fastmath @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            p, dp, ddp1, ddp2 = getders_local(sk2,i,j,k)
            getdEdp_pt!(dEdp4, p, dp, ddp1, ddp2, sk.mpi, i, j, k)

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
    newton_flow!(skyrmion, skyrmion_dot; steps = n, dt=ls^2/80.0, frequency_of_printing = freq, print_stuff = true)
    
Applies a newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The energy is checked every `freq` steps.

See also [`gradient_flow!`, `arrested_newton_flow!`]
"""
function newton_flow!(ϕ, ϕd; dt=ϕ.ls[1]/20.0, steps=1, print_stuff = true, frequency_of_printing = steps, step_algorithm::String)

    if print_stuff == true
        println("intial energy: ", Energy(ϕ))
    end

    counter = 0
    while counter < steps

        newton_flow_for_n_steps!(ϕ,ϕd,dt,frequency_of_printing, step_algorithm="RK4")
        counter += frequency_of_printing

        if print_stuff == true 
            println("after ", counter, " steps, energy = ", Energy(ϕ) )
        end

    end

    return

end


function newton_flow_for_n_steps!(ϕ,ϕd,dt,n;step_algorithm::String)

    
    if step_algorithm == "RK4"
        dEdp1 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
        dEdp2 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
        dEdp3 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
        dEdp4 = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    end

    for _ in 1:n

        if step_algorithm == "Euler"
            newton_flow_for_1_step!(ϕ,ϕd,dEdp1,dt)
        elseif step_algorithm == "RK4"
            newton_flow_for_1_step_RK4!(ϕ,ϕd,dEdp1,dEdp2,dEdp3,dEdp4,dt)
        end

    end

end

function L2_err(A)

    errtot = 0.0

    for i in eachindex(A)
        errtot += A[i]^2
    end

    return errtot

end

function abs_err(A)

    errtot = 0.0

    for i in eachindex(A)
        errtot += abs(A[i])
    end

    return errtot

end

function max_abs_err(A)

    return maximum(abs.(A)) 

end





#= RAK flow


"""
    flowDE!(skyrmion; steps = dt=0.0001)
    
Applies a gradient flow to `skyrmion` with timestep `dt` for `n` steps.

Within the code, an array which holds the variation is created. As such, it is significantly more efficient to use n=1000 than to loop the method 1000 times with `n=1`.

"""
function flowDE!(ϕ; time=1.0, print_stuff = false)




    #error = 2.0*tolerance

    tspan = (0.0,time)
    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    pars = (dEdp, ϕ.lp, ϕ.ls, ϕ.mpi)

    cb_normalise = DiscreteCallback(condition_always_true, normaliseDE!, save_positions=(false,false))
    prob = ODEProblem( dedp!, ϕ.pion_field, tspan, pars)
    sol = solve(prob,Tsit5(),save_everystep=false, callback=cb_normalise)

    ϕ.pion_field .= sol[end]
    
    #return sol[end]

end

function condition_always_true(u,t,integrator)
    true
end

function normaliseDE!(integrator)
    normer!(integrator.u)
end



#= Functions for calculating the functional gradients
function dedp!(dsk, sk, pars, t)

    (dEdp, lp, ls, mpi ) = pars

    Threads.@threads for i in 3:lp[1]-2
        @inbounds for j in 3:lp[2]-2, k in 3:lp[3]-2
        
            p = getX(sk,i,j,k)
            dp = getDX(sk ,i, j, k, ls )
            ddp = getDDX(sk, i, j, k, ls)

            dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            dEdp[i,j,k,4] = dedfpt4v(dp,ddp,mpi^2)

            DEdotpion_field = dEdp[i,j,k,1]*p[1] + dEdp[i,j,k,2]*p[2] + dEdp[i,j,k,3]*p[3] + dEdp[i,j,k,4]*p[4]

            for a in 1:4
                dEdp[i,j,k,a] -= p[a]*DEdotpion_field
                dsk[i,j,k,a] = dEdp[i,j,k,a]
            end


        end
    end

   
end =#


# Functions for calculating the functional gradients
function dedp!(dsk, sk, pars, t)

   # (dEdp, lp, ls, mpi ) = pars

    Threads.@threads for i in 3:pars[2][1]-2
        @inbounds for j in 3:pars[2][2]-2, k in 3:pars[2][3]-2
        
            p = getX(sk,i,j,k)
            dp = getDX(sk ,i, j, k, pars[3] )
            ddp = getDDX(sk, i, j, k, pars[3])

            pars[1][i,j,k,1] = dedfpt1v(dp,ddp)
            pars[1][i,j,k,2] = dedfpt2v(dp,ddp)
            pars[1][i,j,k,3] = dedfpt3v(dp,ddp)
            pars[1][i,j,k,4] = dedfpt4v(dp,ddp,pars[4]^2)

            DEdotpion_field = pars[1][i,j,k,1]*p[1] + pars[1][i,j,k,2]*p[2] + pars[1][i,j,k,3]*p[3] + pars[1][i,j,k,4]*p[4]

            for a in 1:4
                pars[1][i,j,k,a] -= p[a]*DEdotpion_field
                #dsk[i,j,k,a] = pars[1][i,j,k,a]
            end


        end
    end

    dsk .= pars[1]
    
    #dsk .= dt.*dEdp;
    #normer!(sk)
   
end 


GRAVEYARD


#=
function getdEdp!(sk, dEdp)

    Threads.@threads for i in sk.sum_grid[1]
        @inbounds for j in sk.sum_grid[2], k in sk.sum_grid[3]
        
            p = getX(sk,i,j,k)

            if sk.periodic == false
                dp = getDX(sk ,i, j, k )
                ddp = getDDX(sk, i, j, k)
            else
                dp = getDXp(sk ,i, j, k )
                ddp = getDDXp(sk, i, j, k)
            end

            getdEdp_pt!(dEdp, p, dp, ddp, sk.mpi, i, j, k)

        end
    end

end

function getdEdp_pt!(dEdp, p, dp, ddp, mpi, i, j, k)

    dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
    dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
    dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
    dEdp[i,j,k,4] = dedfpt4v(dp,ddp,mpi^2)

    DEdotpion_field = dEdp[i,j,k,1]*p[1] + dEdp[i,j,k,2]*p[2] + dEdp[i,j,k,3]*p[3] + dEdp[i,j,k,4]*p[4]

    for a in 1:4
        dEdp[i,j,k,a] -= p[a]*DEdotpion_field
    end
end


function dedfpt1v(dp,ddp)
    return ddp[3,1] + ddp[3,1]*dp[1,2]^2 + ddp[3,1]*dp[1,3]^2 + ddp[3,1]*dp[1,4]^2 - dp[1,1]*((ddp[2,2] + ddp[3,2])*dp[1,2] + (ddp[2,3] + ddp[3,3])*dp[1,3] + (ddp[2,4] + ddp[3,4])*dp[1,4]) + ddp[6,2]*dp[1,2]*dp[2,1] + ddp[6,3]*dp[1,3]*dp[2,1] + ddp[6,4]*dp[1,4]*dp[2,1] + ddp[6,2]*dp[1,1]*dp[2,2] - 2*ddp[6,1]*dp[1,2]*dp[2,2] - ddp[1,2]*dp[2,1]*dp[2,2] - ddp[3,2]*dp[2,1]*dp[2,2] + ddp[3,1]*dp[2,2]^2 + ddp[6,3]*dp[1,1]*dp[2,3] - 2*ddp[6,1]*dp[1,3]*dp[2,3] - ddp[1,3]*dp[2,1]*dp[2,3] - ddp[3,3]*dp[2,1]*dp[2,3] + ddp[3,1]*dp[2,3]^2 + ddp[6,4]*dp[1,1]*dp[2,4] - 2*ddp[6,1]*dp[1,4]*dp[2,4] - ddp[1,4]*dp[2,1]*dp[2,4] - ddp[3,4]*dp[2,1]*dp[2,4] + ddp[3,1]*dp[2,4]^2 + ddp[5,2]*dp[1,2]*dp[3,1] + ddp[5,3]*dp[1,3]*dp[3,1] + ddp[5,4]*dp[1,4]*dp[3,1] + ddp[4,2]*dp[2,2]*dp[3,1] + ddp[4,3]*dp[2,3]*dp[3,1] + ddp[4,4]*dp[2,4]*dp[3,1] + ddp[5,2]*dp[1,1]*dp[3,2] - 2*ddp[5,1]*dp[1,2]*dp[3,2] + ddp[4,2]*dp[2,1]*dp[3,2] - 2*ddp[4,1]*dp[2,2]*dp[3,2] - ddp[1,2]*dp[3,1]*dp[3,2] - ddp[2,2]*dp[3,1]*dp[3,2] + ddp[5,3]*dp[1,1]*dp[3,3] - 2*ddp[5,1]*dp[1,3]*dp[3,3] + ddp[4,3]*dp[2,1]*dp[3,3] - 2*ddp[4,1]*dp[2,3]*dp[3,3] - ddp[1,3]*dp[3,1]*dp[3,3] - ddp[2,3]*dp[3,1]*dp[3,3] + (ddp[5,4]*dp[1,1] - 2*ddp[5,1]*dp[1,4] + ddp[4,4]*dp[2,1] - 2*ddp[4,1]*dp[2,4] - (ddp[1,4] + ddp[2,4])*dp[3,1])*dp[3,4] + ddp[2,1]*(1 + dp[1,2]^2 + dp[1,3]^2 + dp[1,4]^2 + dp[3,2]^2 + dp[3,3]^2 + dp[3,4]^2) + ddp[1,1]*(1 + dp[2,2]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,2]^2 + dp[3,3]^2 + dp[3,4]^2)
end
function dedfpt2v(dp,ddp)
    return ddp[3,2] + ddp[3,2]*dp[1,1]^2 - ddp[2,1]*dp[1,1]*dp[1,2] - ddp[3,1]*dp[1,1]*dp[1,2] - ddp[2,3]*dp[1,2]*dp[1,3] - ddp[3,3]*dp[1,2]*dp[1,3] + ddp[3,2]*dp[1,3]^2 - ddp[2,4]*dp[1,2]*dp[1,4] - ddp[3,4]*dp[1,2]*dp[1,4] + ddp[3,2]*dp[1,4]^2 - 2*ddp[6,2]*dp[1,1]*dp[2,1] + ddp[6,1]*dp[1,2]*dp[2,1] + ddp[3,2]*dp[2,1]^2 + ddp[6,1]*dp[1,1]*dp[2,2] + ddp[6,3]*dp[1,3]*dp[2,2] + ddp[6,4]*dp[1,4]*dp[2,2] - ddp[1,1]*dp[2,1]*dp[2,2] - ddp[3,1]*dp[2,1]*dp[2,2] + ddp[6,3]*dp[1,2]*dp[2,3] - 2*ddp[6,2]*dp[1,3]*dp[2,3] - ddp[1,3]*dp[2,2]*dp[2,3] - ddp[3,3]*dp[2,2]*dp[2,3] + ddp[3,2]*dp[2,3]^2 + ddp[6,4]*dp[1,2]*dp[2,4] - 2*ddp[6,2]*dp[1,4]*dp[2,4] - ddp[1,4]*dp[2,2]*dp[2,4] - ddp[3,4]*dp[2,2]*dp[2,4] + ddp[3,2]*dp[2,4]^2 - 2*ddp[5,2]*dp[1,1]*dp[3,1] + ddp[5,1]*dp[1,2]*dp[3,1] - 2*ddp[4,2]*dp[2,1]*dp[3,1] + ddp[4,1]*dp[2,2]*dp[3,1] + ddp[5,1]*dp[1,1]*dp[3,2] + ddp[5,3]*dp[1,3]*dp[3,2] + ddp[5,4]*dp[1,4]*dp[3,2] + ddp[4,1]*dp[2,1]*dp[3,2] + ddp[4,3]*dp[2,3]*dp[3,2] + ddp[4,4]*dp[2,4]*dp[3,2] - ddp[1,1]*dp[3,1]*dp[3,2] - ddp[2,1]*dp[3,1]*dp[3,2] + ddp[5,3]*dp[1,2]*dp[3,3] - 2*ddp[5,2]*dp[1,3]*dp[3,3] + ddp[4,3]*dp[2,2]*dp[3,3] - 2*ddp[4,2]*dp[2,3]*dp[3,3] - ddp[1,3]*dp[3,2]*dp[3,3] - ddp[2,3]*dp[3,2]*dp[3,3] + (ddp[5,4]*dp[1,2] - 2*ddp[5,2]*dp[1,4] + ddp[4,4]*dp[2,2] - 2*ddp[4,2]*dp[2,4] - (ddp[1,4] + ddp[2,4])*dp[3,2])*dp[3,4] + ddp[2,2]*(1 + dp[1,1]^2 + dp[1,3]^2 + dp[1,4]^2 + dp[3,1]^2 + dp[3,3]^2 + dp[3,4]^2) + ddp[1,2]*(1 + dp[2,1]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,3]^2 + dp[3,4]^2)
end
function dedfpt3v(dp,ddp)
    return ddp[3,3] + ddp[3,3]*dp[1,1]^2 + ddp[3,3]*dp[1,2]^2 - ddp[2,1]*dp[1,1]*dp[1,3] - ddp[3,1]*dp[1,1]*dp[1,3] - ddp[2,2]*dp[1,2]*dp[1,3] - ddp[3,2]*dp[1,2]*dp[1,3] - ddp[2,4]*dp[1,3]*dp[1,4] - ddp[3,4]*dp[1,3]*dp[1,4] + ddp[3,3]*dp[1,4]^2 - 2*ddp[6,3]*dp[1,1]*dp[2,1] + ddp[6,1]*dp[1,3]*dp[2,1] + ddp[3,3]*dp[2,1]^2 - 2*ddp[6,3]*dp[1,2]*dp[2,2] + ddp[6,2]*dp[1,3]*dp[2,2] + ddp[3,3]*dp[2,2]^2 + ddp[6,1]*dp[1,1]*dp[2,3] + ddp[6,2]*dp[1,2]*dp[2,3] + ddp[6,4]*dp[1,4]*dp[2,3] - ddp[1,1]*dp[2,1]*dp[2,3] - ddp[3,1]*dp[2,1]*dp[2,3] - ddp[1,2]*dp[2,2]*dp[2,3] - ddp[3,2]*dp[2,2]*dp[2,3] + ddp[6,4]*dp[1,3]*dp[2,4] - 2*ddp[6,3]*dp[1,4]*dp[2,4] - ddp[1,4]*dp[2,3]*dp[2,4] - ddp[3,4]*dp[2,3]*dp[2,4] + ddp[3,3]*dp[2,4]^2 - 2*ddp[5,3]*dp[1,1]*dp[3,1] + ddp[5,1]*dp[1,3]*dp[3,1] - 2*ddp[4,3]*dp[2,1]*dp[3,1] + ddp[4,1]*dp[2,3]*dp[3,1] - 2*ddp[5,3]*dp[1,2]*dp[3,2] + ddp[5,2]*dp[1,3]*dp[3,2] - 2*ddp[4,3]*dp[2,2]*dp[3,2] + ddp[4,2]*dp[2,3]*dp[3,2] + ddp[5,1]*dp[1,1]*dp[3,3] + ddp[5,2]*dp[1,2]*dp[3,3] + ddp[5,4]*dp[1,4]*dp[3,3] + ddp[4,1]*dp[2,1]*dp[3,3] + ddp[4,2]*dp[2,2]*dp[3,3] + ddp[4,4]*dp[2,4]*dp[3,3] - ddp[1,1]*dp[3,1]*dp[3,3] - ddp[2,1]*dp[3,1]*dp[3,3] - ddp[1,2]*dp[3,2]*dp[3,3] - ddp[2,2]*dp[3,2]*dp[3,3] + (ddp[5,4]*dp[1,3] - 2*ddp[5,3]*dp[1,4] + ddp[4,4]*dp[2,3] - 2*ddp[4,3]*dp[2,4] - (ddp[1,4] + ddp[2,4])*dp[3,3])*dp[3,4] + ddp[2,3]*(1 + dp[1,1]^2 + dp[1,2]^2 + dp[1,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,4]^2) + ddp[1,3]*(1 + dp[2,1]^2 + dp[2,2]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,4]^2)
end
function dedfpt4v(dp,ddp,mpisq)
    return mpisq + ddp[3,4] + ddp[3,4]*dp[1,1]^2 + ddp[3,4]*dp[1,2]^2 + ddp[3,4]*dp[1,3]^2 - ddp[2,1]*dp[1,1]*dp[1,4] - ddp[3,1]*dp[1,1]*dp[1,4] - ddp[2,2]*dp[1,2]*dp[1,4] - ddp[3,2]*dp[1,2]*dp[1,4] - ddp[2,3]*dp[1,3]*dp[1,4] - ddp[3,3]*dp[1,3]*dp[1,4] - 2*ddp[6,4]*dp[1,1]*dp[2,1] + ddp[6,1]*dp[1,4]*dp[2,1] + ddp[3,4]*dp[2,1]^2 - 2*ddp[6,4]*dp[1,2]*dp[2,2] + ddp[6,2]*dp[1,4]*dp[2,2] + ddp[3,4]*dp[2,2]^2 - 2*ddp[6,4]*dp[1,3]*dp[2,3] + ddp[6,3]*dp[1,4]*dp[2,3] + ddp[3,4]*dp[2,3]^2 + ddp[6,1]*dp[1,1]*dp[2,4] + ddp[6,2]*dp[1,2]*dp[2,4] + ddp[6,3]*dp[1,3]*dp[2,4] - ddp[1,1]*dp[2,1]*dp[2,4] - ddp[3,1]*dp[2,1]*dp[2,4] - ddp[1,2]*dp[2,2]*dp[2,4] - ddp[3,2]*dp[2,2]*dp[2,4] - ddp[1,3]*dp[2,3]*dp[2,4] - ddp[3,3]*dp[2,3]*dp[2,4] - 2*ddp[5,4]*dp[1,1]*dp[3,1] + ddp[5,1]*dp[1,4]*dp[3,1] - 2*ddp[4,4]*dp[2,1]*dp[3,1] + ddp[4,1]*dp[2,4]*dp[3,1] - 2*ddp[5,4]*dp[1,2]*dp[3,2] + ddp[5,2]*dp[1,4]*dp[3,2] - 2*ddp[4,4]*dp[2,2]*dp[3,2] + ddp[4,2]*dp[2,4]*dp[3,2] - 2*ddp[5,4]*dp[1,3]*dp[3,3] + ddp[5,3]*dp[1,4]*dp[3,3] - 2*ddp[4,4]*dp[2,3]*dp[3,3] + ddp[4,3]*dp[2,4]*dp[3,3] + ddp[2,4]*(1 + dp[1,1]^2 + dp[1,2]^2 + dp[1,3]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2) + ddp[1,4]*(1 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2) + (ddp[5,1]*dp[1,1] + ddp[5,2]*dp[1,2] + ddp[5,3]*dp[1,3] + ddp[4,1]*dp[2,1] + ddp[4,2]*dp[2,2] + ddp[4,3]*dp[2,3] - (ddp[1,1] + ddp[2,1])*dp[3,1] - (ddp[1,2] + ddp[2,2])*dp[3,2] - (ddp[1,3] + ddp[2,3])*dp[3,3])*dp[3,4]
end




"""
    arrested_newton_flow!(skyrmion, skyrmion_dot; steps = n, tolerance = tol, dt=ls^2/80.0, frequency_of_checking_tolerance = freq, print_stuff = true)
    
Applies an arrested newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `freq` steps.

See also [`gradient_flow!`, `newton_flow!`]
"""
function arrested_newton_flow!(ϕ,ϕd; dt=ϕ.ls[1]/10.0, steps=1, tolerance = 0.0, frequency_of_checking_tolerance = max(100,steps), print_stuff = true, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4),error_function::Function=L2_err)

    println(dt)

    if tolerance == 0 && frequency_of_checking_tolerance > steps
        frequency_of_checking_tolerance = steps
    end

    energy_density = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3])
    old_pion_field = deepcopy(ϕ.pion_field);

    if print_stuff == true
        println("intial energy: ", Energy(ϕ))
    end

    counter = 0
    while counter < steps

        arrested_newton_flow_for_n_steps!(ϕ,ϕd,old_pion_field,dEdp,dt,energy_density,frequency_of_checking_tolerance)
        error = max_abs_err(dEdp)
        counter += frequency_of_checking_tolerance

        if print_stuff == true 
            println("after ", counter, " steps, error = ", round(error, sigdigits=4), " energy = ", round(sum(energy_density)*ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3]/(12.0*pi^2), sigdigits=8) )
        end

        if tolerance != 0.0    # => we are in tol mode
            if error < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += frequency_of_checking_tolerance
            end
        end

    end

    return

end

function arrested_newton_flow_for_n_steps!(ϕ,ϕd,old_pion_field,dEdp,dt,energy_density,n)

    new_energy = EnergyANF(ϕ,energy_density)

    for _ in 1:n

        old_energy = new_energy
        old_pion_field .= ϕ.pion_field

        newton_flow_for_1_step!(ϕ,ϕd,dEdp,dt)
        new_energy = EnergyANF(ϕ,energy_density)

        if new_energy > old_energy
            
            fill!(ϕd, 0.0);
            ϕ.pion_field .= old_pion_field;
            gradient_flow!(ϕ,steps=20,print_stuff=false,dEdp=dEdp)
            new_energy = EnergyANF(ϕ,energy_density)
    
        end

    end
    
end

function newton_flow_for_n_steps!(ϕ,ϕd,dEdp,dt,n)
        
    for _ in 1:n
        newton_flow_for_1_step!(ϕ,ϕd,dEdp,dt)
    end


end

=#



=#