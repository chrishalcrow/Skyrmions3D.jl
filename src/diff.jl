"""
    gradient_flow!(skyrmion; steps = n, tolerance = tol, dt=ls^2/80.0, frequency_of_checking_tolerance = freq, print_stuff = true)
    
Applies a gradient flow to `skyrmion` with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `freq` steps.

See also [`newton_flow!`, `arrested_newton_flow!`]

"""
function gradient_flow!(ϕ; steps = 1, dt=((ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3])^(2/3))/80.0, tolerance = 0.0, frequency_of_checking_tolerance = max(100,steps), print_stuff = true, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4) , error_function::Function=L2_err)

    if tolerance == 0 && frequency_of_checking_tolerance > steps
        frequency_of_checking_tolerance = steps
    end
    
    if print_stuff == true
        println("initial: energy: ", Energy(ϕ) )
    end

    counter = 0

    while counter < steps
        
        gradient_flow_for_n_steps!(ϕ,dEdp,frequency_of_checking_tolerance,dt)
        error = error_function(dEdp)
        counter += frequency_of_checking_tolerance
        
        if tolerance != 0.0    # => we are in tol mode
            steps += frequency_of_checking_tolerance
        end

        if print_stuff == true
            println("after ", counter, " steps, error = ", round(error; digits=6))
        end

        if error < tolerance
            counter = steps + 1    # => end the while loop
        end
            
    end

    if print_stuff == true
        println("final energy: ", Energy(ϕ) )
    end

    return

end

function gradient_flow_for_n_steps!(ϕ,dEdp,n,dt)
    
    if ϕ.periodic == true
        for _ in 1:n
            gradient_flow_1_step_periodic!(ϕ,dEdp,dt)
        end
    else 
        for _ in 1:n
            gradient_flow_1_step!(ϕ,dEdp,dt)
        end
    end

end


function gradient_flow_1_step!(sk, dEdp, dt)

    getdEdp!(sk, dEdp)
    sk.pion_field .+= dt.*dEdp;
    normer!(sk)
   
end 

function gradient_flow_1_step_periodic!(sk, dEdp, dt)

    getdEdp_periodic!(sk, dEdp)
    sk.pion_field .+= dt.*dEdp;
    normer!(sk)
   
end 

function getdEdp!(sk, dEdp)

    Threads.@threads for i in 3:sk.lp[1]-2
        @inbounds for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
            p = getX(sk,i,j,k)
            dp = getDX(sk ,i, j, k )
            ddp = getDDX(sk, i, j, k)

            getdEdp_pt!(dEdp, p, dp, ddp, sk.mpi, i, j, k)

        end
    end

end

function getdEdp_periodic!(sk, dEdp)

    Threads.@threads for i in 1:sk.lp[1]
        @inbounds for j in 1:sk.lp[2], k in 1:sk.lp[3]
        
            p = getX(sk,i,j,k)
            dp = getDXp(sk ,i, j, k )
            ddp = getDDXp(sk, i, j, k)

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

"""
    arrested_newton_flow!(skyrmion, skyrmion_dot; steps = n, tolerance = tol, dt=ls^2/80.0, frequency_of_checking_tolerance = freq, print_stuff = true)
    
Applies an arrested newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `freq` steps.

See also [`gradient_flow!`, `newton_flow!`]
"""
function arrested_newton_flow!(ϕ,ϕd; dt=ϕ.ls[1]/12.0, steps=1, tolerance = 0.0, frequency_of_checking_tolerance = max(100,steps), print_stuff = true, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4),error_function::Function=L2_err)

    if tolerance == 0 && frequency_of_checking_tolerance > steps
        frequency_of_checking_tolerance = steps
    end

    energy_density = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3])
    previous_pion_field = deepcopy(ϕ.pion_field);

    if print_stuff == true
        println("intial energy: ", Energy(ϕ))
    end

    counter = 0
    while counter < steps

        arrested_newton_flow_for_n_steps!(ϕ,ϕd,previous_pion_field,dEdp,dt,energy_density,frequency_of_checking_tolerance)
        error = error_function(dEdp)
        counter += frequency_of_checking_tolerance

        if tolerance != 0.0    # => we are in tol mode
            steps += frequency_of_checking_tolerance
        end

        if print_stuff == true 
            println("after ", counter, " steps, error = ", error, " energy = ", sum(energy_density)*ϕ.ls[1]*ϕ.ls[2]*ϕ.ls[3]/(12.0*pi^2) )
        end

        if error < tolerance
            counter = steps + 1    # => end the while loop
        end

    end

    return

end

function arrested_newton_flow_for_n_steps!(ϕ,ϕd,old_pion_field,dEdp,dt,energy_density,n)

    new_energy = EnergyANF(ϕ,energy_density)

    if ϕ.periodic == true
        for _ in 1:n

            old_energy = new_energy
            old_pion_field .= ϕ.pion_field

            newton_flow_for_1_step_periodic!(ϕ,ϕd,dEdp,dt)
            new_energy = EnergyANF_periodic(ϕ,energy_density)
    
            if new_energy > old_energy
                
                fill!(ϕd, 0.0);
                ϕ.pion_field .= old_pion_field;
      
            end
        end
    else 
        for _ in 1:n

            old_energy = new_energy
            old_pion_field .= ϕ.pion_field

            newton_flow_for_1_step!(ϕ,ϕd,dEdp,dt)
            new_energy = EnergyANF(ϕ,energy_density)
    
            if new_energy > old_energy
                
                fill!(ϕd, 0.0);
                ϕ.pion_field .= old_pion_field;
      
            end

        end
    end

end

function EnergyANF(sk, ED)

    Threads.@threads for i in 3:sk.lp[1]-2
        @inbounds for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
    
            dp = getDX(sk ,i, j, k )
            ED[i,j,k] = engpt(dp,sk.pion_field[i,j,k,4],sk.mpi)

        end
    end    
        
    return sum(ED)

end 


function EnergyANF_periodic(sk, ED)


    Threads.@threads for i in 1:sk.lp[1]
        @inbounds for j in 1:sk.lp[2], k in 1:sk.lp[3]
    
            dp = getDXp(sk ,i, j, k )
            ED[i,j,k] = engpt(dp,sk.pion_field[i,j,k,4],sk.mpi)

        end
    end
            
    return sum(ED)

end 

"""
    newton_flow!(skyrmion, skyrmion_dot; steps = n, dt=ls^2/80.0, frequency_of_printing = freq, print_stuff = true)
    
Applies a newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The energy is checked every `freq` steps.

See also [`gradient_flow!`, `arrested_newton_flow!`]
"""
function newton_flow!(ϕ, ϕd; dt=ϕ.ls[1]/20.0, steps=1, print_stuff = true, frequency_of_printing = steps, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4))

    if print_stuff == true
        println("intial energy: ", Energy(ϕ))
    end

    counter = 0
    while counter < steps

        newton_flow_for_n_steps!(ϕ,ϕd,dEdp,dt,frequency_of_printing)
        counter += frequency_of_printing

        if print_stuff == true 
            println("after ", counter, " steps, energy = ", Energy(ϕ) )
        end

    end

    return

end


function newton_flow_for_n_steps!(ϕ,ϕd,dEdp,dt,n)

    if ϕ.periodic == true
        for _ in 1:n

            newton_flow_for_1_step_periodic!(ϕ,ϕd,dEdp,dt)

        end
    else 
        for _ in 1:n

            newton_flow_for_1_step!(ϕ,ϕd,dEdp,dt)

        end
    end

end

function newton_flow_for_1_step!(sk, skd ,dEdp, dt)

    getdEdp!(sk, dEdp)
    sk.pion_field .-= dt.*skd
    skd .-= dt.*dEdp
    normer!(sk)
   
end

function newton_flow_for_1_step_periodic!(sk, skd ,dEdp, dt)

    getdEdp_periodic!(sk, dEdp)
    sk.pion_field .-= dt.*skd
    skd .-= dt.*dEdp
    normer!(sk)
   
end



# The gradient functions, in Voigt notation, generated externally in Mathematica.

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

function L2_err(A)

    errtot = 0.0

    for i in eachindex(A)
        errtot += A[i]^2
    end

    return sqrt(errtot)

end

function abs_err(A)

    errtot = 0.0

    for i in eachindex(A)
        errtot += abs(A[i])
    end

    return sqrt(errtot)

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

=#