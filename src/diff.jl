
function dEdu!(du,u,p,t)
	
    println(t)

    dp = zeros(3,4)
    ddp = zeros(6,4)
    
    @simd for i in 3:p[2][1]-2
        @inbounds for j in 3:p[2][2]-2, k in 3:p[2][3]-2
        
            getDXf!(dp, u ,i, j, k, p[1])
            getDDXf!(ddp, u, i, j, k, p[1])

            p[4][i,j,k,1] = dedfpt1v(dp,ddp)
            p[4][i,j,k,2] = dedfpt2v(dp,ddp)
            p[4][i,j,k,3] = dedfpt3v(dp,ddp)
            p[4][i,j,k,4] = dedfpt4v(dp,ddp,p[3])

            for a in 1:4
                du[i,j,k,a] = p[4][i,j,k,a]
                for b in 1:4
                    du[i,j,k,a]  -= p[4][i,j,k,b]*u[i,j,k,b]*u[i,j,k,a]
                end
            end

                    
        end
    end
	
end




condition(u, t, integrator) = true#t ∈ integrator.p[4]

function affect!(integrator)
	
	normer!(integrator.u)
	
end

function flowRAK!(sk,mpi)
    
	#dosetimes = [ i for i=0.0:0.01:10.0 ]
	
    #(ls, lp, m) = p

    dEdp = zeros(sk.lp[1], sk.lp[2], sk.lp[3], 4)
    

	p=[sk.ls, sk.lp, mpi, dEdp ]

    tspan = (0.0,0.1)	
	
	cb = DiscreteCallback(condition,affect!, save_positions=(false,false))

    prob = ODEProblem(dEdu!,sk.phi,tspan,p)
	sol = solve(prob, Tsit5(), save_everystep=false, callback=cb)#,reltol=1e-8, abstol=1e-8,)

    for i in 3:p[2][1]-2, j in 3:p[2][2]-2, k in 3:p[2][3]-2, a in 1:4
        sk.phi[i,j,k,a] = sol[end][i,j,k,a]
    end

   # return sol[end]
    
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
function dedfpt4v(dp,ddp,mpi)
    return mpi + ddp[3,4] + ddp[3,4]*dp[1,1]^2 + ddp[3,4]*dp[1,2]^2 + ddp[3,4]*dp[1,3]^2 - ddp[2,1]*dp[1,1]*dp[1,4] - ddp[3,1]*dp[1,1]*dp[1,4] - ddp[2,2]*dp[1,2]*dp[1,4] - ddp[3,2]*dp[1,2]*dp[1,4] - ddp[2,3]*dp[1,3]*dp[1,4] - ddp[3,3]*dp[1,3]*dp[1,4] - 2*ddp[6,4]*dp[1,1]*dp[2,1] + ddp[6,1]*dp[1,4]*dp[2,1] + ddp[3,4]*dp[2,1]^2 - 2*ddp[6,4]*dp[1,2]*dp[2,2] + ddp[6,2]*dp[1,4]*dp[2,2] + ddp[3,4]*dp[2,2]^2 - 2*ddp[6,4]*dp[1,3]*dp[2,3] + ddp[6,3]*dp[1,4]*dp[2,3] + ddp[3,4]*dp[2,3]^2 + ddp[6,1]*dp[1,1]*dp[2,4] + ddp[6,2]*dp[1,2]*dp[2,4] + ddp[6,3]*dp[1,3]*dp[2,4] - ddp[1,1]*dp[2,1]*dp[2,4] - ddp[3,1]*dp[2,1]*dp[2,4] - ddp[1,2]*dp[2,2]*dp[2,4] - ddp[3,2]*dp[2,2]*dp[2,4] - ddp[1,3]*dp[2,3]*dp[2,4] - ddp[3,3]*dp[2,3]*dp[2,4] - 2*ddp[5,4]*dp[1,1]*dp[3,1] + ddp[5,1]*dp[1,4]*dp[3,1] - 2*ddp[4,4]*dp[2,1]*dp[3,1] + ddp[4,1]*dp[2,4]*dp[3,1] - 2*ddp[5,4]*dp[1,2]*dp[3,2] + ddp[5,2]*dp[1,4]*dp[3,2] - 2*ddp[4,4]*dp[2,2]*dp[3,2] + ddp[4,2]*dp[2,4]*dp[3,2] - 2*ddp[5,4]*dp[1,3]*dp[3,3] + ddp[5,3]*dp[1,4]*dp[3,3] - 2*ddp[4,4]*dp[2,3]*dp[3,3] + ddp[4,3]*dp[2,4]*dp[3,3] + ddp[2,4]*(1 + dp[1,1]^2 + dp[1,2]^2 + dp[1,3]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2) + ddp[1,4]*(1 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2) + (ddp[5,1]*dp[1,1] + ddp[5,2]*dp[1,2] + ddp[5,3]*dp[1,3] + ddp[4,1]*dp[2,1] + ddp[4,2]*dp[2,2] + ddp[4,3]*dp[2,3] - (ddp[1,1] + ddp[2,1])*dp[3,1] - (ddp[1,2] + ddp[2,2])*dp[3,2] - (ddp[1,3] + ddp[2,3])*dp[3,3])*dp[3,4]
end



function gradvD_SA!(sk, dEdp, dt)

    Threads.@threads for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
            dp = getDX(sk ,i, j, k )
            ddp = getDDX(sk, i, j, k)

            @inbounds dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            @inbounds dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            @inbounds dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            @inbounds dEdp[i,j,k,4] = dedfpt4v(dp,ddp,sk.mpi)
                    
        end
    end
    
    Threads.@threads for i in 3:sk.lp[1]-2
         for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
            
            DEdotPHI = Float64(0.0)
            for a in 1:4
                @inbounds DEdotPHI += dEdp[i,j,k,a]*sk.phi[i,j,k,a]
            end

            for a in 1:4
                @inbounds sk.phi[i,j,k,a] += dt*(dEdp[i,j,k,a] - sk.phi[i,j,k,a]*DEdotPHI)
            end

            # reusing DEdotPHI as the normalising constant, to reduce allocations
            @inbounds DEdotPHI = 1.0/sqrt( sk.phi[i,j,k,1]^2 + sk.phi[i,j,k,2]^2 + sk.phi[i,j,k,3]^2 + sk.phi[i,j,k,4]^2 )
			for a in 1:4
				@inbounds sk.phi[i,j,k,a] *= DEdotPHI
			end
        end
    end
    

   
end 



function gradvD!(sk, dEdp, dt, dp, ddp)

    mpi = sk.mpi

    @inbounds for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
            getDX!(dp, sk ,i, j, k )
            getDDX!(ddp, sk, i, j, k)

            dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            dEdp[i,j,k,4] = dedfpt4v(dp,ddp,mpi)
                    
        end
    end
    
    @inbounds for i in 3:sk.lp[1]-2
         for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
            
            DEdotPHI = 0.0
            for a in 1:4
                DEdotPHI += dEdp[i,j,k,a]*sk.phi[i,j,k,a]
            end

            for a in 1:4
                sk.phi[i,j,k,a] += dt*(dEdp[i,j,k,a] - sk.phi[i,j,k,a]*DEdotPHI)
            end
        end
    end
    
    normer!(sk)
   
end 


"""
    flow!(skyrmion, n; dt=0.0001)
    
Applies a gradient flow to `skyrmion` with timestep `dt` for `n` steps.

Within the code, an array which holds the variation is created. As such, it is significantly more efficient to use n=1000 than to loop the method 1000 times with `n=1`.

"""
function flow!(ϕ, n; dt=0.0001)

    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    

    println("initial: ", Energy(ϕ) )

    if Threads.nthreads() == 1
        
        dp = zeros(3,4)
        ddp = zeros(6,4)
        
        for _ in 1:n
            gradvD!(ϕ,dEdp,dt,dp,ddp)
        end

    else
        for _ in 1:n
            gradvD_SA!(ϕ,dEdp,dt)
        end
    end

    println("  final: ", Energy(ϕ) )
end




function ANFflow!(ϕ,ϕd,dt,n)

    println("intial energy: ", Energy(ϕ))

    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    ED = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3])

    previous_phi = zeros(Float64, ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    previous_phi .= ϕ.phi;

    old_energy = 10000000.0
    new_energy = 5000000.0

    arrestnumber = Int64(0)

    for _ in 1:n
  
        old_energy = new_energy
        stepANF!(ϕ,ϕd,previous_phi,dEdp,dt)
        new_energy = EnergyANF(ϕ,ED)
        #println("old: ", old_energy, ", new: ", new_energy)

        if new_energy > old_energy

            arrestnumber += 1

            #println("ARREST!")

            fill!(ϕd, 0.0);
            ϕ.phi .= previous_phi;
  
        end        

    end

    println("Final energy: ", Energy(ϕ))

    return arrestnumber

end



function stepANF!(sk, skd, previous_phi,dEdp, dt)


    Threads.@threads for i in 3:sk.lp[1]-2
    #for i in 3:sk.lp[1]-2
        @inbounds for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2

            
        
            dp = getDX(sk ,i, j, k )
            ddp = getDDX(sk, i, j, k)
    
            dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            dEdp[i,j,k,4] = dedfpt4v(dp,ddp,sk.mpi)

            DEdotPHI = Float64(0.0)
            for a in 1:4
                previous_phi[i,j,k,a] = sk.phi[i,j,k,a]
                @inbounds DEdotPHI += dEdp[i,j,k,a]*sk.phi[i,j,k,a]
            end

            for a in 1:4
                dEdp[i,j,k,a] -= sk.phi[i,j,k,a]*DEdotPHI
            end


        end
    end

    Threads.@threads for i in 3:sk.lp[1]-2
    #for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
           
           for a in 1:4
               @inbounds sk.phi[i,j,k,a] -= dt*skd[i,j,k,a]
           end
           for a in 1:4

               @inbounds skd[i,j,k,a] -= dt*dEdp[i,j,k,a]    
           end
   
           @inbounds normer = 1.0/sqrt( sk.phi[i,j,k,1]^2 + sk.phi[i,j,k,2]^2 + sk.phi[i,j,k,3]^2 + sk.phi[i,j,k,4]^2 )
           for a in 1:4
               @inbounds sk.phi[i,j,k,a] *= normer
           end
       end
   end
   
end

























