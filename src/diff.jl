"""
    gradient_flow!(skyrmion, n; dt=0.0001)
    
Applies a gradient flow to `skyrmion` with timestep `dt` for `n` steps.

Within the code, an array which holds the variation is created. As such, it is significantly more efficient to use n=1000 than to loop the method 1000 times with `n=1`.

"""
function gradient_flow!(ϕ; steps = 0, dt=ϕ.ls[1]*ϕ.ls[2]/80.0, tolerance = 0.0, frequency_of_checking_tolerance = 100, print_stuff = false)

    if steps == 0 && tolerance == 0.0
        println("ERROR: no `steps` or `tolerance` given. Please either set the number of steps you want to flow, or the tolerance required")
        return
    end

    error = 2.0*tolerance

    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    
    if print_stuff == true
        println("initial: energy: ", Energy(ϕ) )
    end

    if steps != 0
        for i in 1:steps
            grad!(ϕ,dEdp,dt)
            if print_stuff == true 
                if i % frequency_of_checking_tolerance == 0
                    error = L2_err(dEdp)
                    println("error = ", error)
                end
            end
        end

    else
        while error > tolerance
            for _ in 1:frequency_of_checking_tolerance
                grad!(ϕ,dEdp,dt)
            end
            error = L2_err(dEdp)
            if print_stuff == true 
                println("error = ", error)
            end
        end
    end


    if print_stuff == true
        println("final energy: ", Energy(ϕ) )
        println("error: ", L2_err(dEdp) )
    end

end



# Chooses which gradient to apply, based on periodic BCs or not
function grad!(ϕ,dEdp,dt)
    if ϕ.periodic == true
        gradvDp!(ϕ,dEdp,dt)
    else
        gradvD!(ϕ,dEdp,dt)
    end
end



# Functions for calculating the functional gradients
function gradvD!(sk, dEdp, dt)

    Threads.@threads for i in 3:sk.lp[1]-2
        @inbounds for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
            p = getX(sk,i,j,k)
            dp = getDX(sk ,i, j, k )
            ddp = getDDX(sk, i, j, k)

            dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            dEdp[i,j,k,4] = dedfpt4v(dp,ddp,sk.mpi^2)

            DEdotpion_field = dEdp[i,j,k,1]*p[1] + dEdp[i,j,k,2]*p[2] + dEdp[i,j,k,3]*p[3] + dEdp[i,j,k,4]*p[4]

            for a in 1:4
                dEdp[i,j,k,a] -= p[a]*DEdotpion_field
            end

        end
    end
    
    sk.pion_field .+= dt.*dEdp;
    normer!(sk)
   
end 

# same as above, but periodic
function gradvDp!(sk, dEdp, dt)

    Threads.@threads for i in 1:sk.lp[1]
        @inbounds for j in 1:sk.lp[2], k in 1:sk.lp[3]
        
            p = getX(sk,i,j,k)
            dp = getDXp(sk ,i, j, k )
            ddp = getDDXp(sk, i, j, k)

            dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            dEdp[i,j,k,4] = dedfpt4v(dp,ddp,sk.mpi^2)

            DEdotpion_field = dEdp[i,j,k,1]*p[1] + dEdp[i,j,k,2]*p[2] + dEdp[i,j,k,3]*p[3] + dEdp[i,j,k,4]*p[4]

            for a in 1:4
                dEdp[i,j,k,a] -= p[a]*DEdotpion_field
            end

        end
    end
    
    sk.pion_field .+= dt.*dEdp;
    normer!(sk)
   
end 



function arrested_newton_flow!(ϕ,ϕd; dt=ϕ.ls[1]/12.0, steps=0, tolerance = 0.0,frequency_of_checking_tolerance = 100, print_stuff = false)

    if tolerance == 0.0 && steps == 0
        println("ERROR: no `steps` or `tolerance` given. Please either set the number of steps you want to flow, or the tolerance required")
        return
    end

    error = 2.0*tolerance

    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    energy_density = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3])

    previous_pion_field = zeros(Float64, ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    previous_pion_field .= ϕ.pion_field;

    old_energy = 10000000.0
    new_energy = 5000000.0

    arrestnumber = 0

    println("intial energy: ", Energy(ϕ))

    if steps != 0
        for _ in 1:steps
  
            old_energy = new_energy
            stepANF!(ϕ,ϕd,previous_pion_field,dEdp,dt)
            new_energy = EnergyANF(ϕ,energy_density)
    
            if new_energy > old_energy
            
                arrestnumber += 1
                fill!(ϕd, 0.0);
                ϕ.pion_field .= previous_pion_field;
      
            end        
    
        end

    else
        while error > tolerance

            for _ in 1:frequency_of_checking_tolerance

                old_energy = new_energy
                stepANF!(ϕ,ϕd,previous_pion_field,dEdp,dt)
                new_energy = EnergyANF(ϕ,energy_density)
    
                if new_energy > old_energy
            
                    arrestnumber += 1
                    fill!(ϕd, 0.0);
                    ϕ.pion_field .= previous_pion_field;
                    gradient_flow!(ϕ,steps=10,print_stuff=false)
                
                end
            
            end
            
            error = L2_err(dEdp)
            if print_stuff == true 
                println("error = ", error)
            end
        end
    end

    println("Final energy: ", Energy(ϕ))

    return arrestnumber

end

# Chooses which gradient to apply, based on periodic BCs or not
function stepANF!(ϕ,ϕd,previous_pion_field,dEdp,dt)
    if ϕ.periodic == true
        stepANFp!(ϕ,ϕd,previous_pion_field,dEdp,dt)
    else
        stepANFnp!(ϕ,ϕd,previous_pion_field,dEdp,dt)
    end
end


function stepANFnp!(sk, skd, previous_pion_field,dEdp, dt)

    Threads.@threads for i in 3:sk.lp[1]-2
        @inbounds for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2

            p = getX(sk, i, j, k)
            dp = getDX(sk ,i, j, k )
            ddp = getDDX(sk, i, j, k)
    
            dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            dEdp[i,j,k,4] = dedfpt4v(dp,ddp,sk.mpi^2)

            DEdotpion_field = Float64(0.0)
            for a in 1:4
                DEdotpion_field += dEdp[i,j,k,a]*p[a]
            end

            for a in 1:4
                dEdp[i,j,k,a] -= p[a]*DEdotpion_field
            end

        end
    end

    previous_pion_field .= sk.pion_field
    sk.pion_field .-= dt.*skd
    skd .-= dt.*dEdp
    normer!(sk)
   
end


function stepANFp!(sk, skd, previous_pion_field, dEdp, dt)

    Threads.@threads for i in 1:sk.lp[1]
        @inbounds for j in 1:sk.lp[2], k in 1:sk.lp[3]

            p = getX(sk, i, j, k)
            dp = getDXp(sk ,i, j, k )
            ddp = getDDXp(sk, i, j, k)
    
            dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            dEdp[i,j,k,4] = dedfpt4v(dp,ddp,sk.mpi^2)

            DEdotpion_field = Float64(0.0)
            for a in 1:4
                DEdotpion_field += dEdp[i,j,k,a]*p[a]
            end

            for a in 1:4
                dEdp[i,j,k,a] -= p[a]*DEdotpion_field
            end

        end
    end

    previous_pion_field .= sk.pion_field
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









### RAK flow


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