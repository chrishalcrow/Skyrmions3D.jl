module Skyrmions

using Makie
using GLMakie
using DifferentialEquations, DiffEqCallbacks, Quaternions

using Meshing, GeometryBasics, Interpolations, Colors

	export Skyrmion, getmesh,  Nfy!, ADHMpt2, checkunitpt, makeADHM!

	mutable struct Skyrmion
	    phi::Array{Float64, 4}
	    lp::Vector{Int64}
	    ls::Vector{Float64}
	    x::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
	end


	Skyrmion(lp::Int64, ls::Float64; vac = [0.0,0.0,0.0,1.0] ) = Skyrmion(zeros(lp,lp,lp,4) ,[lp,lp,lp],[ls,ls,ls], [ -ls*(lp - 1)/2.0 : ls : ls*(lp - 1)/2.0 for a in 1:3 ] )

	Skyrmion(lp::Vector{Int64}, ls::Vector{Float64}; vac = [0.0,0.0,0.0,1.0] ) = Skyrmion(zeros(lp[1],lp[2],lp[3],4) ,lp, ls, [ -ls[a]*(lp[a] - 1)/2.0 : ls[a] : ls[a]*(lp[a] - 1)./2.0 for a in 1:3 ] )


	include("initialise.jl")
	export makeRM!

	export makeRM, SkyrIso, MakeProduct, SkyrShift, multicubes!
	export Skyr
	export Energy, Baryon

	export flow!
	export flowRAK!
	export checkunit, getMOI

	export imusingnotebook, imusingterminal
	

	include("plotting.jl")
		export plot_skyrmion, plot_field, plot_baryon_density
	include("derivatives.jl")
	include("properties.jl")
	export EnergyD, BaryonD

	include("diff.jl")
	

	function checkunit(phi,lp)
	    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
	        @assert  phi[i,j,k,1]^2 + phi[i,j,k,2]^2 + phi[i,j,k,3]^2 + phi[i,j,k,4]^2 ≈ 1.0 "nooo"
	    end
	end
	
	function checkunitpt(phi,i,j,k)
		return phi[i,j,k,1]^2 + phi[i,j,k,2]^2 + phi[i,j,k,3]^2 + phi[i,j,k,4]^2
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
                
			@inbounds normer = 1.0/sqrt( sk.phi[i,j,k,1]^2 + sk.phi[i,j,k,2]^2 + sk.phi[i,j,k,3]^2 + sk.phi[i,j,k,4]^2 )
	        for a in 1:4
				@inbounds sk.phi[i,j,k,a] *= normer
			end
       
	        end
	    end
	end
	
	function normer!(phi,lp)
    
	    @simd for i in 3:lp[1]-2
	        for j in 3:lp[2]-2, k in 3:lp[3]-2
                
	        @inbounds normer = 1.0/sqrt( phi[i,j,k,1]^2 + phi[i,j,k,2]^2 + phi[i,j,k,3]^2 + phi[i,j,k,4]^2 )
	        for a in 1:4
				phi[i,j,k,a] *= normer
			end
       
	        end
	    end
	end

end
