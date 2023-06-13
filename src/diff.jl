
function dEdu!(du,u,p,t)
	
    println(t)

    dp = zeros(3,4)
    ddp = zeros(6,4)
    
    for i in 3:p[2][1]-2
        for j in 3:p[2][2]-2, k in 3:p[2][3]-2
        
            getDXf!(dp, u ,i, j, k, p[1])
            getDDXf!(ddp, u, i, j, k, p[1])


            #@inbounds du[i,j,k,1] = dedfpt1v(dp,ddp)
            #@inbounds du[i,j,k,2] = dedfpt2v(dp,ddp)
            #@inbounds du[i,j,k,3] = dedfpt3v(dp,ddp)
            #@inbounds du[i,j,k,4] = dedfpt4v(dp,ddp,p[3])

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
	
	normer!(integrator.u, integrator.p[2])
	
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


function gradvD!(sk, dEdp,mpi, dt, dp, ddp)

    @simd for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
            getDX!(dp, sk ,i, j, k )
            getDDX!(ddp, sk, i, j, k)

            @inbounds dEdp[i,j,k,1] = dedfpt1v(dp,ddp)
            @inbounds dEdp[i,j,k,2] = dedfpt2v(dp,ddp)
            @inbounds dEdp[i,j,k,3] = dedfpt3v(dp,ddp)
            @inbounds dEdp[i,j,k,4] = dedfpt4v(dp,ddp,mpi)
                    
        end
    end
    
        
    @simd for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2, a in 1:4
            @inbounds sk.phi[i,j,k,a] += dt*dEdp[i,j,k,a]
            for b in 1:4
                sk.phi[i,j,k,a] -= dt*( dEdp[i,j,k,b]*sk.phi[i,j,k,b]*sk.phi[i,j,k,a] ) 
            end
        end
    end
    
    normer(sk)
   
end 

function flow!(ϕ,mpi,dt,n)

    dEdp = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dp = zeros(3,4)
    ddp = zeros(6,4)

    println("initial: ", Energy(ϕ, mpi) )

    for _ in 1:n
        gradvD!(ϕ,dEdp,mpi,dt,dp,ddp)
    end

    println("  final: ", Energy(ϕ, mpi) )
end



function ANFflow!(ϕ,ϕd,mpi,dt,n)

    dEdp = zeros(2,ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dp = zeros(3,4)
    ddp = zeros(6,4)

    previous_phi = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    previous_phi .= ϕ.phi;

    old_energy = 1000.0
    new_energy = 500.0

    for _ in 1:n
        previous_phi .= ϕ.phi;
        old_energy = new_energy
        stepANF!(ϕ,ϕd,dEdp,mpi,dt,dp,ddp)
        new_energy = Energy(ϕ,mpi)

        if new_energy > old_energy

            println("ARREST!")

            ϕd = zeros(ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
            ϕ.phi .= previous_phi

            #new_energy = old_energy


        end

        println(Energy(ϕ, mpi))

    end
end


function stepANF!(sk, skd, dEdp, mpi, dt, dp, ddp)

    dE_dot_mom = 0.0

    for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
            getDX!(dp, sk ,i, j, k )
            getDDX!(ddp, sk, i, j, k)

            @inbounds dEdp[2,i,j,k,1] = dedfpt1v(dp,ddp)
            @inbounds dEdp[2,i,j,k,2] = dedfpt2v(dp,ddp)
            @inbounds dEdp[2,i,j,k,3] = dedfpt3v(dp,ddp)
            @inbounds dEdp[2,i,j,k,4] = dedfpt4v(dp,ddp,mpi)

            for a in 1:4
                dEdp[1,i,j,k,a] = skd[i,j,k,a]
            end

            for a in 1:4
                dE_dot_mom += dEdp[2,i,j,k,a]*skd[i,j,k,a]
            end

        end
    end

    #if dE_dot_mom < 0 
    #    for i in 3:sk.lp[1]-2, j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2, a in 1:4
    #        dEdp[1,i,j,k,a] = 0.0
    #    end
    #end


    #println("dot: ", dE_dot_mom)
    
        
    @simd for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2, a in 1:4

            @inbounds sk.phi[i,j,k,a] += dt*dEdp[1,i,j,k,a]


            @inbounds skd[i,j,k,a] += dt*dEdp[2,i,j,k,a]


            for b in 1:4
                skd[i,j,k,a] -= dt*( dEdp[2,i,j,k,b]*sk.phi[i,j,k,b]*sk.phi[i,j,k,a] ) 
            end



        end
    end
    
    normer(sk)
   
end







function momflow!(ϕ,ϕd,mpi,dt,n;α=1.0, β=1.0)

    dEdp = zeros(2,ϕ.lp[1], ϕ.lp[2], ϕ.lp[3], 4)
    dp = zeros(3,4)
    ddp = zeros(6,4)


    for _ in 1:n

        stepMOM!(ϕ,ϕd,dEdp,mpi,dt,dp,ddp,α, β)
        println(Energy(ϕ, mpi))

    end
end


function stepMOM!(sk, skd, dEdp, mpi, dt, dp, ddp, α, β)


    for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
            getDX!(dp, sk ,i, j, k )
            getDDX!(ddp, sk, i, j, k)

            @inbounds dEdp[2,i,j,k,1] = dedfpt1v(dp,ddp)
            @inbounds dEdp[2,i,j,k,2] = dedfpt2v(dp,ddp)
            @inbounds dEdp[2,i,j,k,3] = dedfpt3v(dp,ddp)
            @inbounds dEdp[2,i,j,k,4] = dedfpt4v(dp,ddp,mpi)

            for a in 1:4
                dEdp[1,i,j,k,a] = skd[i,j,k,a]
            end

            #for a in 1:4
            #    dE_dot_mom += dEdp[2,i,j,k,a]*skd[i,j,k,a]
            #end

        end
    end


    #println("dot: ", dE_dot_mom)
    
        
    @simd for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2, a in 1:4

            @inbounds sk.phi[i,j,k,a] = sk.phi[i,j,k,a] + dt*α*skd[i,j,k,a]


            @inbounds skd[i,j,k,a] = β*skd[i,j,k,a] + dt*dEdp[2,i,j,k,a] 


            for b in 1:4
                skd[i,j,k,a] -= dt*( dEdp[2,i,j,k,b]*sk.phi[i,j,k,b]*sk.phi[i,j,k,a] ) 
            end

            #dt*dEdp[1,i,j,k,a]



        end
    end
    
    normer(sk)
   
end 





















