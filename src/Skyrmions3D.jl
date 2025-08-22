module Skyrmions3D

# Plotting
using Makie, CairoMakie, Requires, Meshing, GeometryBasics, Colors

# Functionality
using StaticArrays, LinearAlgebra, Interpolations

export Skyrmion,
    get_grid, get_field, set_mpi!, set_lattice!, set_Fpi!, set_ee!, set_physical!
export set_periodic!, set_dirichlet!, set_neumann!, set_bounary_conditions!
export check_if_normalised, normer!, normer

include("transform.jl")
export translate_sk, translate_sk!, isorotate_sk, isorotate_sk!, rotate_sk!, rotate_sk
export product_approx, product_approx!, center_skyrmion!, evaluate_sk, quadratic_spline_interpolation

include("properties.jl")
export Energy, get_energy_density!, Baryon, get_baryon_density!, center_of_mass, rms_baryon, compute_current, overview, sphericity

include("initialise.jl")
export make_rational_map!, make_RM_product!, make_ADHM!
export R_from_axis_angle

include("plotting.jl")
export activate_CairoMakie, plot_field, plot_baryon_density, plot_overview, plot_scan

include("derivatives.jl")

include("grid.jl")
export Grid

include("diff.jl")
export gradient_flow!, arrested_newton_flow!, newton_flow!

"""
    Skyrmion(lp::Int64, ls::Float64; kwargs...)
	Skyrmion([lpx, lpy, lpx], [lsx, lsy, lsz]; kwargs...)
    
Create a skyrme field with `lp` lattice points and `ls` lattice spacing. 

# Optional arguments
- `mpi = 0.0`: sets the pion mass for this Skyrme field
- `Fpi = 180`: sets the pion decay constant for this Skyrme field
- `ee = 4.0`: sets the Skyrme constant for this Skyrme field
- `physical = false`: whether the Skyrmion is using physical units

"""
mutable struct Skyrmion
    pion_field::Array{Float64,4}
    grid::Grid
    mpi::Float64
    Fpi::Float64
    ee::Float64
    physical::Bool
end


Skyrmion(
    lp::Int64,
    ls::Float64;
    Fpi = 180,
    ee = 4.0,
    vac = [0.0, 0.0, 0.0, 1.0],
    mpi = 0.0,
    boundary_conditions = "dirichlet",
) = Skyrmion(
    vacuum_skyrmion(lp, lp, lp, vac),
    Grid([lp, lp, lp], [ls, ls, ls], boundary_conditions),
    mpi,
    Fpi,
    ee,
    false,
)

Skyrmion(
    lp::Vector{Int64},
    ls::Vector{Float64};
    Fpi = 180,
    ee = 4.0,
    vac = [0.0, 0.0, 0.0, 1.0],
    mpi = 0.0,
    boundary_conditions = "dirichlet",
) = Skyrmion(
    vacuum_skyrmion(lp[1], lp[2], lp[3], vac),
    Grid([lp[1], lp[2], lp[3]], [ls[1], ls[2], ls[3]], boundary_conditions),
    mpi,
    Fpi,
    ee,
    false,
)


"""
    get_field(skyrmion)

Returns the array of pion fields `[π1, π2, π3, π0]` of `skyrmion`, which can be used in integrals.

"""
function get_field(skyrmion)

    return [skyrmion.pion_field[:, :, :, a] for a = 1:4]

end

""" 
    get_grid(skyrmion)

Returns an array of 3D arrays `[x, y, z]`, which can be used in integrals.

"""
function get_grid(skyrmion)

    x_grid = [
        skyrmion.x[1][i] for
        i = 1:skyrmion.grid.lp[1], j = 1:skyrmion.grid.lp[2], k = 1:skyrmion.grid.lp[3]
    ]
    y_grid = [
        skyrmion.x[2][j] for
        i = 1:skyrmion.grid.lp[1], j = 1:skyrmion.grid.lp[2], k = 1:skyrmion.grid.lp[3]
    ]
    z_grid = [
        skyrmion.x[3][k] for
        i = 1:skyrmion.grid.lp[1], j = 1:skyrmion.grid.lp[2], k = 1:skyrmion.grid.lp[3]
    ]

    return (x_grid, y_grid, z_grid)

end

mutable struct profile
    field::Vector{Float64}
    lp::Int64
    ls::Float64
    r_grid::Vector{Float64}
end

profile(lp, ls) = profile(zeros(Float64, lp), lp, ls, [ls*i for i = 0:(lp-1)])

function is_dirichlet(boundary_conditions)

    if boundary_conditions == "dirichlet"
        return true
    else
        return false
    end
end


"""
    set_mpi!(skyrmion, mpi)

Set the pion mass of `skyrmion` to `mpi`.

"""
function set_mpi!(sk, mpi)
    sk.mpi = mpi
end



function set_bounary_conditions!(sk, boundary_conditions)

    sk.grid.boundary_conditions = boundary_conditions
    sk.grid.sum_grid = sum_grid(sk.grid.lp, boundary_conditions)
    sk.grid.index_grid_x = index_grid(sk.grid.lp[1], boundary_conditions)
    sk.grid.index_grid_y = index_grid(sk.grid.lp[2], boundary_conditions)
    sk.grid.index_grid_z = index_grid(sk.grid.lp[3], boundary_conditions)

end

"""
    set_periodic!(skyrmion)

Sets the `skyrmion` to have periodic boundary conditions.

"""
function set_periodic!(sk)

    sk.grid.dirichlet = false

    set_bounary_conditions!(sk, "periodic")

    println("Periodic boundary conditions activated")

end

"""
    set_dirichlet!(skyrmion)

Sets the `skyrmion` to have Dirichlet boundary conditions.

"""
function set_neumann!(sk)

    sk.grid.dirichlet = false

    set_bounary_conditions!(sk, "neumann")

    println("Neumann boundary conditions activated")

end

"""
    set_dirichlet!(skyrmion)

Sets the `skyrmion` to have periodic boundary conditions.

"""
function set_dirichlet!(sk)

    sk.grid.dirichlet = true

    sk.grid.boundary_conditions = "dirichlet"
    sk.grid.sum_grid = sum_grid(sk.grid.lp, sk.grid.boundary_conditions)

    set_dirichlet_boudary!(sk)
    println("Dirichlet boundary conditions activated")

end


"""
    set_Fpi!(skyrmion, Fpi)

Sets the pion decay constant of `skyrmion` to `Fpi`. 

"""
function set_Fpi!(sk, Fpi)

    sk.Fpi = Fpi

end


"""
    set_ee!(skyrmion, ee)

Sets the Skyrme coupling constant of `skyrmion` to `ee`. 

"""
function set_ee!(sk, ee)

    sk.ee = ee

end


"""
    set_physical!(skyrmion, is_physical; Fpi = Fpi, ee = ee)

Sets `skyrmion` to use physical units, when `is_physical` is `true`.

Also used to turn off physical units by setting `is_physical=false`.

"""
function set_physical!(
    skyrmion,
    physical;
    Fpi = skyrmion.Fpi,
    ee = skyrmion.ee,
)

    skyrmion.physical = physical

    if skyrmion.physical == true
        println("Fpi = ", skyrmion.Fpi, ", e = ", skyrmion.ee, " and m = ", skyrmion.mpi)
        println(
            "Hence, mpi = ",
            skyrmion.Fpi*skyrmion.ee*skyrmion.mpi/2.0,
            ", length unit = ",
            197.327*2.0/(skyrmion.ee*skyrmion.Fpi),
            " and energy unit = ",
            skyrmion.Fpi/(4.0*skyrmion.ee),
        )
    end

end

"""
    set_lattice!(skyrmion, lp = [lpx, lpy, lpz], ls = [lsx, lsy, lsz])

Sets the underlying lattice to one with `lpx`x`lpy`x`lpz` points and `lsx`x`lsy`x`lsz` spacings, and reinterpolates `skyrmion` on the new grid.

"""
function set_lattice!(skyrmion, lp, ls)

    old_x = skyrmion.grid.x
    x = setgrid(lp, ls)

    sky_temp = Skyrmion(
        lp,
        ls,
        mpi = skyrmion.mpi,
        boundary_conditions = skyrmion.grid.boundary_conditions,
    )
    vac = [0.0, 0.0, 0.0, 1.0]

    ϕinterp = quadratic_spline_interpolation(skyrmion.pion_field, old_x)

    for k = 1:lp[3], j = 1:lp[2], i = 1:lp[1]

        if old_x[1][1] < x[1][i] < old_x[1][end] &&
           old_x[2][1] < x[2][j] < old_x[2][end] &&
           old_x[3][1] < x[3][k] < old_x[3][end]
            for a = 1:4
                sky_temp.pion_field[i, j, k, a] = ϕinterp[a](x[1][i], x[2][j], x[3][k])
            end
        else
            sky_temp.pion_field[i, j, k, :] .= vac
        end

    end

    skyrmion.grid.lp = sky_temp.grid.lp
    skyrmion.grid.ls = sky_temp.grid.ls
    skyrmion.grid.x = sky_temp.grid.x
    skyrmion.pion_field = zeros(lp[1], lp[2], lp[3], 4)
    skyrmion.pion_field .= sky_temp.pion_field

    skyrmion.grid.index_grid_x = index_grid(lp[1], sky_temp.grid.boundary_conditions)
    skyrmion.grid.index_grid_y = index_grid(lp[2], sky_temp.grid.boundary_conditions)
    skyrmion.grid.index_grid_z = index_grid(lp[3], sky_temp.grid.boundary_conditions)

    skyrmion.grid.sum_grid = sum_grid(lp, sky_temp.grid.boundary_conditions)

    normer!(skyrmion)
    if skyrmion.grid.boundary_conditions == "dirichlet"
        set_dirichlet!(skyrmion)
    end

    println(
        "Your new lattice has ",
        lp[1],
        "*",
        lp[2],
        "*",
        lp[3],
        " points with lattice spacing [",
        ls[1],
        ", ",
        ls[2],
        ", ",
        ls[3],
        "].",
    )

    sky_temp = nothing

end


function vacuum_skyrmion(lpx, lpy, lpz, vac)

    vac_sk = zeros(lpx, lpy, lpz, 4)

    for i = 1:lpx, j = 1:lpy, k = 1:lpz
        vac_sk[i, j, k, :] = vac
    end

    return vac_sk

end



function sum_grid(lp, boundary_conditions)
    # We allow for lp to be given as a single integer, in which case we set
    # the number of lattice points in each direction to be lp. 
    if isa(lp, Integer)
        lp = [lp, lp, lp]
    end

    if boundary_conditions == "dirichlet"
        return [3:(lp[1]-2), 3:(lp[2]-2), 3:(lp[3]-2)]
    else
        return [1:lp[1], 1:lp[2], 1:lp[3]]
    end

end

function index_grid(lp, boundary_conditions)

    index_grid_array = zeros(Int64, lp+4)

    if boundary_conditions == "periodic"

        for i = 1:(lp+4)
            index_grid_array[i] = mod1(i-2, lp)
        end

    else

        for i = 3:(lp+2)
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

Check if `skyrmion` is normalised.

Throws an error if any point is not normalised, i.e. the pion field does not have norm 1.

"""
function check_if_normalised(skyrmion)
    for i = 1:skyrmion.grid.lp[1], j = 1:skyrmion.grid.lp[2], k = 1:skyrmion.grid.lp[3]
        @assert skyrmion.pion_field[i, j, k, 1]^2 +
                skyrmion.pion_field[i, j, k, 2]^2 +
                skyrmion.pion_field[i, j, k, 3]^2 +
                skyrmion.pion_field[i, j, k, 4]^2 ≈ 1.0 "nooo"
    end
end



function setgrid(lp, ls)

    x = [zeros(lp[1]), zeros(lp[2]), zeros(lp[3])]
    for a = 1:3
        x[a] = [-0.5*ls[a]*(lp[a]-1) + n*ls[a] for n = 0:(lp[a]-1)]
    end

    return x

end


"""
    normer!(skyrmion)

Normalises `skyrmion`.

See also [`normer`](@ref). 

"""
function normer!(sk)

    normer!(sk.pion_field, sk.grid.lp)

end

"""
    normer(skyrmion)

Returns normalised `skyrmion`.

See also [`normer!`](@ref). 

"""
function normer(sk)

    sk_new = deepcopy(sk)
    normer!(sk_new.pion_field, sk_new.grid.lp)

    return sk_new

end


function normer!(pion_field::Array{Float64,4}, lp)

    Threads.@threads for k = 1:lp[3]
        for j = 1:lp[2], i = 1:lp[1]

            @inbounds normer =
                1.0/sqrt(
                    pion_field[i, j, k, 1]^2 +
                    pion_field[i, j, k, 2]^2 +
                    pion_field[i, j, k, 3]^2 +
                    pion_field[i, j, k, 4]^2,
                )
            for a = 1:4
                @inbounds pion_field[i, j, k, a] *= normer
            end

        end
    end
end

end
