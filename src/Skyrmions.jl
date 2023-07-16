module Skyrmions

using Makie
using GLMakie, WGLMakie, CairoMakie

using DifferentialEquations, DiffEqCallbacks

using Meshing, GeometryBasics, Interpolations, Colors, StaticArrays, LinearAlgebra

export Skyrmion, check_if_normalised, makeADHM!, normer!, normer_SA!, R_from_axis_angle, turn_on_physical!,  turn_off_physical!, stepANF!, compute_current, center_skyrmion!, resize_lattice!, flowDE!, an_flow!

export dxD, dyD, dzD, d2xD, d2yD, d2zD, dxdyD, dxdzD, dydzD



include("transform.jl")
export translate_sk, translate_sk!, isorotate_sk, isorotate_sk!, rotate_sk!, rotate_sk, product_approx, product_approx!, make_RM_product!, set_dirichlet!

include("properties.jl")
export EnergyD, BaryonD, Energy, Baryon, getMOI, center_of_mass, rms_baryon, Baryon_SA


include("initialise.jl")
export makeRationalMap!, R_from_axis_angle

export makeRM, SkyrIso, MakeProduct, SkyrShift, multicubes!
export Skyr, ANFflow!,  momflow!, array2, B3_tet_data, B4_cube_data


export gradient_flow!, arrested_newton_flow!




include("plotting.jl")
export plot_field, plot_baryon_density, interactive_flow
 
include("derivatives.jl")


include("diff.jl")
export newton_flow!


"""
    Skyrmion(lp::Int64, ls::Float64)
	Skyrmion([lpx,lpy,lpx], [lsx,lsy,lsz])
    
Create a skyrme field with `lp` lattice points and `ls` lattice spacing. 

# Optional arguments
- `mpi = 0.0`: sets the pion mass for this Skyrme field

"""
mutable struct Skyrmion
	pion_field::Array{Float64, 4}
	lp::Vector{Int64}
	ls::Vector{Float64}
	x::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
	mpi::Float64
	Fpi::Float64
	ee::Float64
	physical::Bool
	periodic::Bool
	index_grid_x::Vector{Int64}
	index_grid_y::Vector{Int64}
	index_grid_z::Vector{Int64}
end


Skyrmion(lp::Int64, ls::Float64; vac = [0.0,0.0,0.0,1.0], mpi = 0.0, periodic=false ) = Skyrmion(vacuum_skyrmion(lp,lp,lp,vac) ,[lp,lp,lp],[ls,ls,ls], [ -ls*(lp - 1)/2.0 : ls : ls*(lp - 1)/2.0 for a in 1:3 ] , mpi, 180.0, 4.0, false, periodic,index_grid(lp), index_grid(lp), index_grid(lp) )

Skyrmion(lp::Vector{Int64}, ls::Vector{Float64}; vac = [0.0,0.0,0.0,1.0], mpi = 0.0 , periodic=false) = Skyrmion(vacuum_skyrmion(lp[1],lp[2],lp[3],vac) ,lp, ls, [ -ls[a]*(lp[a] - 1)/2.0 : ls[a] : ls[a]*(lp[a] - 1)./2.0 for a in 1:3 ], mpi ,180.0, 4.0, false, periodic,index_grid(lp[1]), index_grid(lp[2]), index_grid(lp[3]) )

function vacuum_skyrmion(lpx,lpy,lpz,vac)

	vac_sk = zeros(lpx,lpy,lpz,4)

	for i in 1:lpx, j in 1:lpy, k in 1:lpz
		vac_sk[i,j,k,:] = vac
	end

	return vac_sk

end

function index_grid(lp)

	index_grid_array = zeros(lp+4)

	for i in 1:lp+4
		index_grid_array[i] = mod1(i-2,lp)
	end

	return index_grid_array

end


"""
	turn_on_physical!(skyrmion)

Turns on physical units. Output of energy etc will now be displayed in units of MeV and fm.

"""
function turn_on_physical!(skyrmion)
	
	skyrmion.physical = true

	println("Fpi = ", skyrmion.Fpi, ",  e = ", skyrmion.ee, " and m = ", skyrmion.mpi)
	println("Hence, mpi = ", skyrmion.Fpi*skyrmion.ee*skyrmion.mpi/2.0, ", length unit = ", 197.327*2.0/(skyrmion.ee*skyrmion.Fpi), "and energy unit = ", skyrmion.Fpi/(4.0*skyrmion.ee))

end

function turn_off_physical!(skyrmion)
	
	skyrmion.physical = false

end



"""
    check_if_normalised(skyrmion)

Check if skyrmion is normalised.

Throws an error if any point is not normalised

"""
function check_if_normalised(skyrmion)
	for i in 1:skyrmion.lp[1], j in 1:skyrmion.lp[2], k in 1:skyrmion.lp[3]
		@assert  skyrmion.pion_field[i,j,k,1]^2 + skyrmion.pion_field[i,j,k,2]^2 + skyrmion.pion_field[i,j,k,3]^2 + skyrmion.pion_field[i,j,k,4]^2 ≈ 1.0 "nooo"
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

	Threads.@threads for i in 1:sk.lp[1]
		for j in 1:sk.lp[2], k in 1:sk.lp[3]
			
			@inbounds normer = 1.0/sqrt( sk.pion_field[i,j,k,1]^2 + sk.pion_field[i,j,k,2]^2 + sk.pion_field[i,j,k,3]^2 + sk.pion_field[i,j,k,4]^2 )
			for a in 1:4
				@inbounds sk.pion_field[i,j,k,a] *= normer
			end
	
		end
	end
end


function normer!(pion_field::Array{Float64, 4})

	lp = size(pion_field)[1:3]

	Threads.@threads for i in 1:lp[1]
		for j in 1:lp[2], k in 1:lp[3]
			
			@inbounds normer = 1.0/sqrt( pion_field[i,j,k,1]^2 + pion_field[i,j,k,2]^2 + pion_field[i,j,k,3]^2 + pion_field[i,j,k,4]^2 )
			for a in 1:4
				@inbounds pion_field[i,j,k,a] *= normer
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

	Threads.@threads for i in 1:sk.lp[1]
		for j in 1:sk.lp[2], k in 1:sk.lp[3]
			
			@inbounds normer = 1.0/sqrt( sk.pion_field[i,j,k,1]^2 + sk.pion_field[i,j,k,2]^2 + sk.pion_field[i,j,k,3]^2 + sk.pion_field[i,j,k,4]^2 )
			for a in 1:4
				@inbounds sk_new.pion_field[i,j,k,a] = sk.pion_field[i,j,k,a]*normer
			end
	
		end
	end

	return sk_new

end

end

