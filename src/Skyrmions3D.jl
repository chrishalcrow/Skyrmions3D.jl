module Skyrmions3D

#= TO DO LIST 

- Scanning plot
- B=4 Alberto stuff

=#


using Makie
using GLMakie, WGLMakie, CairoMakie
using Optimization, OptimizationOptimJL, ModelingToolkit, Symbolics

using Meshing, GeometryBasics, Interpolations, Colors, StaticArrays, LinearAlgebra

export Skyrmion, set_mpi!, set_periodic!, set_lattice!, set_Fpi!, set_ee!, set_physical!, set_lattice!
export check_if_normalised, normer!

include("transform.jl")
export translate_sk, translate_sk!, isorotate_sk, isorotate_sk!, rotate_sk!, rotate_sk
export product_approx, product_approx!, center_skyrmion!

include("properties.jl")
export Energy, Baryon, center_of_mass, rms_baryon, compute_current, overview

include("initialise.jl")
export make_rational_map!, make_RM_product!, make_ADHM!, get_close_ADHM_data

include("plotting.jl")
export plot_field, plot_baryon_density, interactive_flow, plot_overview, plot_scan
 
include("derivatives.jl")


include("diff.jl")
export gradient_flow!, arrested_newton_flow!


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
	sum_grid::Vector{UnitRange{Int64}}
	derivative_functions::Vector{Function}
end


Skyrmion(lp::Int64, ls::Float64; Fpi = 180, ee = 4.0, vac = [0.0,0.0,0.0,1.0], mpi = 0.0, periodic=false ) = Skyrmion(vacuum_skyrmion(lp,lp,lp,vac) ,[lp,lp,lp],[ls,ls,ls], [ -ls*(lp - 1)/2.0 : ls : ls*(lp - 1)/2.0 for a in 1:3 ] , mpi, Fpi, ee, false, periodic,index_grid(lp), index_grid(lp), index_grid(lp), sum_grid(lp, periodic), [getDX, getDDX] )

Skyrmion(lp::Vector{Int64}, ls::Vector{Float64}; Fpi = 180, ee = 4.0, vac = [0.0,0.0,0.0,1.0], mpi = 0.0 , periodic=false) = Skyrmion(vacuum_skyrmion(lp[1],lp[2],lp[3],vac) ,lp, ls, [ -ls[a]*(lp[a] - 1)/2.0 : ls[a] : ls[a]*(lp[a] - 1)./2.0 for a in 1:3 ], mpi ,Fpi, ee, false, periodic,index_grid(lp[1]), index_grid(lp[2]), index_grid(lp[3]), sum_grid(lp,periodic), [getDX, getDDX] )

mutable struct profile
    field::Vector{Float64}
    lp::Int64
    ls::Float64
    r_grid::Vector{Float64}
end


profile(lp,ls) = profile( zeros(Float64,lp), lp, ls,  [ ls*i for i in 0:(lp-1) ] )



"""
    set_mpi!(skyrmion::Skyrmion, mpi)

Set the pion mass of `skyrmion` to `mpi`.
"""
function set_mpi!(sk::Skyrmion, mpi)
	sk.mpi = mpi
end

"""
	set_periodic!(skyrmion::Skyrmion, is_periodic)

Sets the `skyrmion` to have periodic boundary conditions if `is_periodic` is `true`, and Dirichlet boundary conditions if `is_periodic` is `false.
"""
function set_periodic!(sk::Skyrmion, periodic::Bool)
	
	sk.periodic = periodic
	sk.sum_grid = sum_grid(sk.lp, sk.periodic)

	if periodic == true
		println("Periodic boundary conditions activated")
	else
		set_dirichlet!(sk)
		println("Dirichlet boundary conditions activated")
	end

end

"""
	set_Fpi!(skyrmion::Skyrmion, Fpi)

Sets the pion decay constant of `skyrmion` to `Fpi`. 
"""
function set_Fpi!(sk::Skyrmion, Fpi)
	
	sk.Fpi = Fpi

end


"""
	set_ee!(skyrmion::Skyrmion, ee)

Sets the Skyrme coupling constant of `skyrmion` to `ee`. 
"""
function set_ee!(sk::Skyrmion, ee)
	
	sk.ee = ee

end


"""
    set_physical!(skyrmion::Skyrmion, is_physical; Fpi=Fpi, ee=ee)

Sets `skyrmion` to use physical units with `Fpi` MeV and skyrme coupling `ee`, when `is_physical` is `true`.
Also used to turn off physical units by setting is_physical=false
"""
function set_physical!(skyrmion::Skyrmion, physical::Bool; Fpi=skyrmion.Fpi, ee=skyrmion.ee)
	
	skyrmion.physical = physical

	if skyrmion.physical == true
		println("Fpi = ", skyrmion.Fpi, ",  e = ", skyrmion.ee, " and m = ", skyrmion.mpi)
		println("Hence, mpi = ", skyrmion.Fpi*skyrmion.ee*skyrmion.mpi/2.0, ", length unit = ", 197.327*2.0/(skyrmion.ee*skyrmion.Fpi), "and energy unit = ", skyrmion.Fpi/(4.0*skyrmion.ee))
	end

end

"""
    set_lattice!(skyrmion, lp = [lpx, lpy, lpz], ls = [lsx, lsy, lsz])

Sets the underlying lattice to one with `lpx`x`lpy`x`lpz` points and `lsx`x`lsy`x`lsz` spacings, and reinterpolates `skyrmion` on the new grid.

"""
function set_lattice!(skyrmion, lp, ls)

    old_x = skyrmion.x
    x = setgrid(lp,ls)

    sky_temp = Skyrmion(lp, ls,  mpi = skyrmion.mpi , periodic = skyrmion.periodic)
    vac = [0.0,0.0,0.0,1.0]

    ϕinterp = [ extrapolate(scale(interpolate( skyrmion.pion_field[:,:,:,a] , BSpline(Quadratic()) ), (old_x[1],old_x[2],old_x[3]) ), Throw()) for a in 1:4 ]

    for k in 1:lp[3], j in 1:lp[2], i in 1:lp[1]

        if old_x[1][1] < x[1][i] < old_x[1][end] && old_x[2][1] < x[2][j] < old_x[2][end] && old_x[3][1] < x[3][k] < old_x[3][end]
            for a in 1:4
                sky_temp.pion_field[i,j,k,a] = ϕinterp[a](x[1][i], x[2][j], x[3][k])
            end
        else
            sky_temp.pion_field[i,j,k,:] .= vac
        end
        
    end
    
    skyrmion.lp = sky_temp.lp
    skyrmion.ls = sky_temp.ls
    skyrmion.x = sky_temp.x
    skyrmion.pion_field = zeros(lp[1],lp[2],lp[3],4)
    skyrmion.pion_field .= sky_temp.pion_field

    skyrmion.index_grid_x = index_grid(lp[1])
    skyrmion.index_grid_y = index_grid(lp[2])
    skyrmion.index_grid_z = index_grid(lp[3])
    
    skyrmion.sum_grid = sum_grid(lp,skyrmion.periodic)
    
    normer!(skyrmion)
    if skyrmion.periodic == false
        set_dirichlet!(skyrmion)
    end

	println("Your new lattice has ", lp[1],"*",lp[2],"*",lp[3]," points with lattice spacing [",ls[1],", ",ls[2],", ",ls[3],"].")

end


function vacuum_skyrmion(lpx,lpy,lpz,vac)

	vac_sk = zeros(lpx,lpy,lpz,4)

	for i in 1:lpx, j in 1:lpy, k in 1:lpz
		vac_sk[i,j,k,:] = vac
	end

	return vac_sk

end

function sum_grid(lp::Integer,periodic::Bool)

	if periodic == true
		return [ 1:lp, 1:lp, 1:lp ]
	else
		return [ 3:lp-2, 3:lp-2, 3:lp-2]
	end

end

function sum_grid(lp::Vector{Int64},periodic::Bool)

	if periodic == true
		return [ 1:lp[1], 1:lp[2], 1:lp[3] ]
	else
		return [ 3:lp[1]-2, 3:lp[2]-2, 3:lp[3]-2]
	end
	
end

function index_grid(lp)

	index_grid_array = zeros(lp+4)

	for i in 1:lp+4
		index_grid_array[i] = mod1(i-2,lp)
	end

	return index_grid_array

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

Normalises `skyrmion`.

See also [`normer`]
"""
function normer!(sk)

	Threads.@threads for k in 1:sk.lp[3]
		for j in 1:sk.lp[2], i in 1:sk.lp[1]
			
			@inbounds normer = 1.0/sqrt( sk.pion_field[i,j,k,1]^2 + sk.pion_field[i,j,k,2]^2 + sk.pion_field[i,j,k,3]^2 + sk.pion_field[i,j,k,4]^2 )
			for a in 1:4
				@inbounds sk.pion_field[i,j,k,a] *= normer
			end
	
		end
	end
end


function normer!(pion_field::Array{Float64, 4})

	lp = size(pion_field)[1:3]

	Threads.@threads for k in 1:lp[3]
		for j in 1:lp[2], i in 1:lp[1]
			
			@inbounds normer = 1.0/sqrt( pion_field[i,j,k,1]^2 + pion_field[i,j,k,2]^2 + pion_field[i,j,k,3]^2 + pion_field[i,j,k,4]^2 )
			for a in 1:4
				@inbounds pion_field[i,j,k,a] *= normer
			end
	
		end
	end
end



"""
    normer(skyrmion)

Returns normalised `skyrmion`.

See also [`normer!`]

"""
function normer(sk)

    lp = sk.lp
    ls = sk.ls

    sk_new = Skyrmion(lp,ls)

	Threads.@threads for k in 1:sk.lp[3]
		for j in 1:sk.lp[2], i in 1:sk.lp[1]
			
			@inbounds normer = 1.0/sqrt( sk.pion_field[i,j,k,1]^2 + sk.pion_field[i,j,k,2]^2 + sk.pion_field[i,j,k,3]^2 + sk.pion_field[i,j,k,4]^2 )
			for a in 1:4
				@inbounds sk_new.pion_field[i,j,k,a] = sk.pion_field[i,j,k,a]*normer
			end
	
		end
	end

	return sk_new

end

end
#=
function Base.:+(q1::Quaternion{Float64}, q2::Quaternion{Float64})
    return Quaternion(q1[1] + q2[1], q1[2] + q2[2], q1[3] + q2[3], q1[4] + q2[4])
end 

function Base.:-(q1::Quaternion{Float64}, q2::Quaternion{Float64})
    return Quaternion(q1[1] - q2[1], q1[2] - q2[2], q1[3] - q2[3], q1[4] - q2[4])
end 
=#