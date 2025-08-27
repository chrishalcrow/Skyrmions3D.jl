
"""
	Grid([lpx, lpy, lpx], [lsx, lsy, lsz], boundary_conditions::String)
    
Create a Grid field with `lp` lattice points and `ls` lattice spacing. 

The current accepted values for `boundary_conditions` are `"dirichlet"`, `"neumann"`, and `"periodic"`. If a different value is provided, a warning is given. 

"""
mutable struct Grid
    lp::Vector{Int64}
    ls::Vector{Float64}
    x::Vector{
        StepRangeLen{
            Float64,
            Base.TwicePrecision{Float64},
            Base.TwicePrecision{Float64},
            Int64,
        },
    }
    dirichlet::Bool
    index_grid_x::Vector{Int64}
    index_grid_y::Vector{Int64}
    index_grid_z::Vector{Int64}
    sum_grid::Vector{UnitRange{Int64}}
    boundary_conditions::String
end

Grid(lp::Vector{Int64}, ls::Vector{Float64}, boundary_conditions::String) = Grid(
    lp,
    ls,
    [(-ls[a]*(lp[a]-1)/2.0):ls[a]:(ls[a]*(lp[a]-1) ./ 2.0) for a = 1:3],
    is_dirichlet(boundary_conditions),
    index_grid(lp[1], boundary_conditions),
    index_grid(lp[2], boundary_conditions),
    index_grid(lp[3], boundary_conditions),
    sum_grid(lp, boundary_conditions),
    boundary_conditions,
)
