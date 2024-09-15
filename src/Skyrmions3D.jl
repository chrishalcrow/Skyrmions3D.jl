module Skyrmions3D

# Plotting
using Makie, CairoMakie, Requires, Meshing, GeometryBasics, Colors

# Functionality
using StaticArrays, LinearAlgebra, Interpolations

export Skyrmion, get_grid, get_field, set_metric!, set_mpi!,  set_lattice!, set_Fpi!, set_ee!, set_physical!
export set_periodic!, set_dirichlet!, set_neumann!
export check_if_normalised, normer!, normer

include("transform.jl")
export translate_sk, translate_sk!, isorotate_sk, isorotate_sk!, rotate_sk!, rotate_sk
export product_approx, product_approx!, center_skyrmion!

include("properties.jl")
export Energy, Baryon, center_of_mass, rms_baryon, compute_current, overview, sphericity, Berger_Isospin

include("initialise.jl")
export make_rational_map!, make_RM_product!, make_ADHM!

include("plotting.jl")
export activate_CairoMakie, plot_field, plot_baryon_density, plot_overview, plot_scan, axial_symmetry_plot
 
include("derivatives.jl")

include("diff.jl")
export gradient_flow!, arrested_newton_flow!, newton_flow!, lin_interpolate

function __init__()
	CairoMakie.activate!()
	@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("plottingGPU.jl")
	@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" export interactive_flow, activate_GLMakie
	@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" println("You have GLMakie installed. Interactive plotting is supported.")
	@require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" GLMakie.activate!(); 
end


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
    metric::Float64
    mpi::Float64
    Fpi::Float64
    ee::Float64
    physical::Bool
    dirichlet::Bool
    index_grid_x::Vector{Int64}
    index_grid_y::Vector{Int64}
    index_grid_z::Vector{Int64}
    sum_grid::Vector{UnitRange{Int64}}
    boundary_conditions::String
end

function Skyrmion(lp::Int64, ls::Float64; metric=1.0, Fpi=180, ee=4.0, vac=[0.0,0.0,0.0,1.0], mpi=0.0, boundary_conditions="dirichlet")
    return Skyrmion(vacuum_skyrmion(lp, lp, lp, vac),
                    [lp, lp, lp],
                    [ls, ls, ls],
                    [-ls*(lp - 1)/2.0 : ls : ls*(lp - 1)/2.0 for a in 1:3],
                    metric,
                    mpi,
                    Fpi,
                    ee,
                    false,
                    is_dirichlet(boundary_conditions),
                    index_grid(lp, boundary_conditions),
                    index_grid(lp, boundary_conditions),
                    index_grid(lp, boundary_conditions),
                    sum_grid(lp, boundary_conditions),
                    boundary_conditions)
end

function Skyrmion(lp::Vector{Int64}, ls::Vector{Float64}; metric=1.0, Fpi=180, ee=4.0, vac=[0.0, 0.0, 0.0, 1.0], mpi=0.0, boundary_conditions="dirichlet")
    return Skyrmion(vacuum_skyrmion(lp[1], lp[2], lp[3], vac),
                    lp,
                    ls,
                    [-ls[a]*(lp[a] - 1)/2.0 : ls[a] : ls[a]*(lp[a] - 1)/2.0 for a in 1:3],
                    metric,
                    mpi,
                    Fpi,
                    ee,
                    false,
                    is_dirichlet(boundary_conditions),
                    index_grid(lp[1], boundary_conditions),
                    index_grid(lp[2], boundary_conditions),
                    index_grid(lp[3], boundary_conditions),
                    sum_grid(lp, boundary_conditions),
                    boundary_conditions)
end


"""
    get_field(skyrmion::Skyrmion)

Returns an array of pion fields `[π1, π2, π3, π0]`, which can be used in integrals.
"""
function get_field(skyrmion::Skyrmion)

	return [ skyrmion.pion_field[:,:,:,a] for a in 1:4 ]

end

"""
    get_grid(skyrmion::Skyrmion)

Returns an array of 3D arrays `[x, y, z]`, which can be used in integrals.
"""
function get_grid(skyrmion::Skyrmion)

	x_grid = [ skyrmion.x[1][i] for i in 1:skyrmion.lp[1], j in 1:skyrmion.lp[2], k in 1:skyrmion.lp[3] ]
	y_grid = [ skyrmion.x[2][j] for i in 1:skyrmion.lp[1], j in 1:skyrmion.lp[2], k in 1:skyrmion.lp[3] ]
	z_grid = [ skyrmion.x[3][k] for i in 1:skyrmion.lp[1], j in 1:skyrmion.lp[2], k in 1:skyrmion.lp[3] ]

	return (x_grid, y_grid, z_grid)

end

mutable struct profile
    field::Vector{Float64}
    lp::Int64
    ls::Float64
    r_grid::Vector{Float64}
end

profile(lp,ls) = profile( zeros(Float64,lp), lp, ls,  [ ls*i for i in 0:(lp-1) ] )

function is_dirichlet(boundary_conditions)

	if boundary_conditions == "dirichlet"
		return true
	else
		return false
	end
end

"""
    set_metric_var!(sk::Skyrmion, metric_var)

Sets the metric variation to `metric_var`.
"""

function set_metric!(sk::Skyrmion, metric)
    sk.metric = metric
end

"""
    set_mpi!(skyrmion::Skyrmion, mpi)

Set the pion mass of `skyrmion` to `mpi`.
"""
function set_mpi!(sk::Skyrmion, mpi)
	sk.mpi = mpi
end

function set_periodic!(sk::Skyrmion, periodic::Bool)
	
	sk.periodic = periodic
	sk.sum_grid = sum_grid(sk.lp, sk.periodic)

	if dirichlet == false
		println("Periodic boundary conditions activated")
	else
		set_dirichlet!(sk)
		println("Dirichlet boundary conditions activated")
	end

end


"""
	set_periodic!(skyrmion::Skyrmion)

Sets the `skyrmion` to have periodic boundary conditions.
"""
function set_periodic!(sk::Skyrmion)
	
	sk.dirichlet = false

	sk.boundary_conditions = "periodic"
	sk.sum_grid = sum_grid(sk.lp, sk.boundary_conditions)

	sk.index_grid_x = index_grid(sk.lp[1], sk.boundary_conditions)
	sk.index_grid_y = index_grid(sk.lp[2], sk.boundary_conditions)
	sk.index_grid_z = index_grid(sk.lp[3], sk.boundary_conditions)

	println("Periodic boundary conditions activated")

end

"""
	set_dirichlet!(skyrmion::Skyrmion)

Sets the `skyrmion` to have Dirichlet boundary conditions.
"""
function set_neumann!(sk::Skyrmion)

	sk.dirichlet = false
	
	sk.boundary_conditions = "neumann"

	sk.sum_grid = sum_grid(sk.lp, sk.boundary_conditions)
	sk.index_grid_x = index_grid(sk.lp[1], sk.boundary_conditions)
	sk.index_grid_y = index_grid(sk.lp[2], sk.boundary_conditions)
	sk.index_grid_z = index_grid(sk.lp[3], sk.boundary_conditions)

	println("Neumann boundary conditions activated")

end

"""
	set_dirichlet!(skyrmion::Skyrmion)

Sets the `skyrmion` to have periodic boundary conditions.
"""
function set_dirichlet!(sk::Skyrmion)

	sk.dirichlet = true
	
	sk.boundary_conditions = "dirichlet"
	sk.sum_grid = sum_grid(sk.lp, sk.boundary_conditions)

	set_dirichlet_boudary!(sk)
	println("Dirichlet boundary conditions activated")

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

Sets `skyrmion` to use physical units, when `is_physical` is `true`.
Also used to turn off physical units by setting is_physical=false
"""
function set_physical!(skyrmion::Skyrmion, physical::Bool; Fpi=skyrmion.Fpi, ee=skyrmion.ee)
	
	skyrmion.physical = physical

	if skyrmion.physical == true
		println("Fpi = ", skyrmion.Fpi, ", e = ", skyrmion.ee, " and m = ", skyrmion.mpi)
		println("Hence, mpi = ", skyrmion.Fpi*skyrmion.ee*skyrmion.mpi/2.0, ", length unit = ", 197.327*2.0/(skyrmion.ee*skyrmion.Fpi), " and energy unit = ", skyrmion.Fpi/(4.0*skyrmion.ee))
	end

end

"""
    set_lattice!(skyrmion, lp = [lpx, lpy, lpz], ls = [lsx, lsy, lsz])

Sets the underlying lattice to one with `lpx`x`lpy`x`lpz` points and `lsx`x`lsy`x`lsz` spacings, and reinterpolates `skyrmion` on the new grid.

"""
function set_lattice!(skyrmion, lp, ls)

    old_x = skyrmion.x
    x = setgrid(lp,ls)

    sky_temp = Skyrmion(lp, ls,  mpi = skyrmion.mpi , boundary_conditions = skyrmion.boundary_conditions)
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

    skyrmion.index_grid_x = index_grid(lp[1], sky_temp.boundary_conditions)
    skyrmion.index_grid_y = index_grid(lp[2], sky_temp.boundary_conditions)
    skyrmion.index_grid_z = index_grid(lp[3], sky_temp.boundary_conditions)
    
    skyrmion.sum_grid = sum_grid(lp, sky_temp.boundary_conditions)
    
    normer!(skyrmion)
    if skyrmion.boundary_conditions == "dirichlet"
        set_dirichlet!(skyrmion)
    end

	println("Your new lattice has ", lp[1],"*",lp[2],"*",lp[3]," points with lattice spacing [",ls[1],", ",ls[2],", ",ls[3],"].")

	sky_temp = nothing

end


function vacuum_skyrmion(lpx,lpy,lpz,vac)

	vac_sk = zeros(lpx,lpy,lpz,4)

	for i in 1:lpx, j in 1:lpy, k in 1:lpz
		vac_sk[i,j,k,:] = vac
	end

	return vac_sk

end

function sum_grid(lp::Integer, boundary_conditions::String)

	if boundary_conditions == "dirichlet"
		return [ 3:lp-2, 3:lp-2, 3:lp-2]
	else
		return [ 1:lp, 1:lp, 1:lp ]
	end

end

function sum_grid(lp::Vector{Int64}, boundary_conditions::String)

	if boundary_conditions == "dirichlet"
		return [ 3:lp[1]-2, 3:lp[2]-2, 3:lp[3]-2]
	else
		return [ 1:lp[1], 1:lp[2], 1:lp[3] ]
	end
	
end

function index_grid(lp, boundary_conditions::String)

	index_grid_array = zeros(Int64, lp+4)

	if boundary_conditions == "periodic"

		for i in 1:lp+4
			index_grid_array[i] = mod1(i-2,lp)
		end

	else

		for i in 3:lp+2
			index_grid_array[i] = i-2
		end
		
		index_grid_array[1] = 2
		index_grid_array[2] = 1
		
		index_grid_array[lp+3] = lp
		index_grid_array[lp+4] = lp-1

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

    normer!(sk.pion_field, sk.lp)

end

"""
    normer(skyrmion)

Returns normalised `skyrmion`.

See also [`normer!`]

"""
function normer(sk)

    sk_new = deepcopy(sk)
    normer!(sk_new.pion_field, sk_new.lp)

	return sk_new

end


function normer!(pion_field::Array{Float64, 4}, lp)

	Threads.@threads for k in 1:lp[3]
		for j in 1:lp[2], i in 1:lp[1]
			
			@inbounds normer = 1.0/sqrt( pion_field[i,j,k,1]^2 + pion_field[i,j,k,2]^2 + pion_field[i,j,k,3]^2 + pion_field[i,j,k,4]^2 )
			for a in 1:4
				@inbounds pion_field[i,j,k,a] *= normer
			end
	
		end
	end
end

end