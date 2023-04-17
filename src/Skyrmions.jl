module Skyrmions

	using Makie
	using GLMakie
	using DifferentialEquations, DiffEqCallbacks

	export setgrid
	export makeRM, makeVac, SkyrIso, MakeProduct, SkyrShift, multicubes!
	export Skyr
	export Energy, Baryon
	export plot_skyrmion
	export flow!
	export flowRAK
	export checkunit, getMOI
	
	include("initialise.jl")
	include("plotting.jl")
	include("derivatives.jl")
	include("properties.jl")
	include("diff.jl")
	
	mutable struct Skyr
	    phi::Array{Float64, 4}
	    dEdp::Array{Float64, 4}
	    ED::Array{Float64, 3}
	    BD::Array{Float64, 3}
	    lp::Vector{Int64}
	    ls::Vector{Float64}
	end
	
	function checkunit(phi,lp)
	    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
	        @assert  phi[1,i,j,k]^2 + phi[2,i,j,k]^2 + phi[3,i,j,k]^2 + phi[4,i,j,k]^2 ≈ 1.0 "nooo"
	    end
	end
	
	
	
	function setgrid(lp,ls)
	
		x = [zeros(lp[1]), zeros(lp[2]), zeros(lp[3])]
		for a in 1:3
		    x[a] = [-0.5*ls[a]*(lp[a]-1) + n*ls[a] for n=0:lp[a]-1]
		end

		return x
	
	end
	
	function normer(sk)
    
	    @simd for i in 3:sk.lp[1]-2
	        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
                
	        @inbounds normer = 1.0/sqrt( sk.phi[1,i,j,k]^2 + sk.phi[2,i,j,k]^2 + sk.phi[3,i,j,k]^2 + sk.phi[4,i,j,k]^2 )
	        for a in 1:4; @inbounds sk.phi[a,i,j,k] *= normer; end;
       
	        end
	    end
	end
	
	function normer!(phi,lp)
    
	    @simd for i in 3:lp[1]-2
	        for j in 3:lp[2]-2, k in 3:lp[3]-2
                
	        @inbounds normer = 1.0/sqrt( phi[1,i,j,k]^2 + phi[2,i,j,k]^2 + phi[3,i,j,k]^2 + phi[4,i,j,k]^2 )
	        for a in 1:4; @inbounds phi[a,i,j,k] *= normer; end;
       
	        end
	    end
	end

end
