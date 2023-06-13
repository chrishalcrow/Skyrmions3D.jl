module Skyrmions

using Makie
using GLMakie, WGLMakie
using DifferentialEquations, DiffEqCallbacks

using Meshing, GeometryBasics, Interpolations, Colors, StaticArrays, Quaternionic

export Skyrmion, check_if_normalised, makeADHM!, normer!, R_from_axis_angle

include("transform.jl")
export translate, translate!, isorotate, isorotate!, rotate!, rotate, product, product!, make_RM_product!

include("properties.jl")
export EnergyD, BaryonD, Energy, Baryon, getMOI, center_of_mass


mutable struct Skyrmion
	phi::Array{Float64, 4}
	lp::Vector{Int64}
	ls::Vector{Float64}
	x::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
	mpi::Float64
end


Skyrmion(lp::Int64, ls::Float64; vac = [0.0,0.0,0.0,1.0], mpi = 0.0 ) = Skyrmion(zeros(lp,lp,lp,4) ,[lp,lp,lp],[ls,ls,ls], [ -ls*(lp - 1)/2.0 : ls : ls*(lp - 1)/2.0 for a in 1:3 ] , mpi)
Skyrmion(lp::Vector{Int64}, ls::Vector{Float64}; vac = [0.0,0.0,0.0,1.0], mpi = 0.0 ) = Skyrmion(zeros(lp[1],lp[2],lp[3],4) ,lp, ls, [ -ls[a]*(lp[a] - 1)/2.0 : ls[a] : ls[a]*(lp[a] - 1)./2.0 for a in 1:3 ], mpi )


include("initialise.jl")
export makeRM!, R_from_axis_angle

export makeRM, SkyrIso, MakeProduct, SkyrShift, multicubes!
export Skyr, ANFflow!,  momflow!, array2, B3_tet_data, B4_cube_data


export flow!
export flowRAK!

export imusingnotebook, imusingterminal, imusingJupyter


include("plotting.jl")
export plot_skyrmion, plot_field, plot_baryon_density, getmesh

include("derivatives.jl")



include("diff.jl")


function imusingJupyter()
	Makie.inline!(true) 
	GLMakie.activate!()
end


"""
    check_if_normalised(skyrmion)

Check if skyrmion is normalised.

Throws an error if *any* point is not normalised

"""
function check_if_normalised(skyrmion)
	for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
		@assert  skyrmion.phi[i,j,k,1]^2 + skyrmion.phi[i,j,k,2]^2 + skyrmion.phi[i,j,k,3]^2 + skyrmion.phi[i,j,k,4]^2 ≈ 1.0 "nooo"
	end
end



"""
    setgrid(lp, ls)

Compute a Cartesian lattice with `lp` lattice points and `ls` lattice spacing.

`lp` and `ls` should be 3-vectors

"""
function setgrid(lp,ls)

	x = [zeros(lp[1]), zeros(lp[2]), zeros(lp[3])]
	for a in 1:3
		x[a] = [-0.5*ls[a]*(lp[a]-1) + n*ls[a] for n=0:lp[a]-1]
	end

	return x

end

"""
    normer!(skyrmion)

Normalise a skyrmion.

See also [`normer`]
"""
function normer!(sk)

	for i in 3:sk.lp[1]-2
		for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
			
			@inbounds normer = 1.0/sqrt( sk.phi[i,j,k,1]^2 + sk.phi[i,j,k,2]^2 + sk.phi[i,j,k,3]^2 + sk.phi[i,j,k,4]^2 )
			for a in 1:4
				@inbounds sk.phi[i,j,k,a] *= normer
			end
	
		end
	end
end


"""
    normer(skyrmion)

Returns the normalised skyrmion

See also [`normer!`]

"""
function normer(sk)

    lp = sk.lp
    ls = sk.ls

    sk_new = Skyrmion(lp,ls)

	@simd for i in 3:lp[1]-2
		for j in 3:lp[2]-2, k in 3:lp[3]-2
			
			@inbounds normer = 1.0/sqrt( sk.phi[i,j,k,1]^2 + sk.phi[i,j,k,2]^2 + sk.phi[i,j,k,3]^2 + sk.phi[i,j,k,4]^2 )
			for a in 1:4
				@inbounds sk_new.phi[i,j,k,a] = sk.phi[i,j,k,a]*normer
			end
	
		end
	end

	return sk_new

end

end