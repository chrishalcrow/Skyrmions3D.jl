function activate_CairoMakie()
    CairoMakie.activate!()
end


"""
    plot_field(skyrmion; component = 3, iso_value = 0.5, kwargs...)
    
Plots an isosurface of constant value, `skyrme_field[component] = iso_value`. 

# Optional argument

Can accept any arguments used in `Axis3` from the `Makie` package.

"""
function plot_field(skyrmion; component = 3, iso_value = 0.5, kwargs...)


    pion_field_to_be_plotted = skyrmion.pion_field[:, :, :, component]

    println(
        minimum(pion_field_to_be_plotted),
        ", ",
        iso_value,
        ", ",
        maximum(pion_field_to_be_plotted),
    )

    if minimum(pion_field_to_be_plotted) < iso_value < maximum(pion_field_to_be_plotted)

        fig = Figure()

        field_mesh = getmesh(pion_field_to_be_plotted, iso_value, skyrmion.grid.x)


        ax = Axis3(
            fig[1, 1],
            title = "Isosurface ϕ" * string(component) * " = ±" * string(iso_value),
            aspect = (skyrmion.grid.lp[1], skyrmion.grid.lp[2], skyrmion.grid.lp[3]),
            kwargs...,
        )

        Makie.mesh!(ax, field_mesh, shading = NoShading)

        return fig
    else
        error("Your iso value is out of range")
    end

end


"""
    plot_overview(skyrmion; iso_value = 0.5)
    
Plots the pion fields and a baryon density of `skyrmion`.

Pion field isosurfaces are taken at the value `iso_value`, while the baryon density isosurface is taken at halfway between the minimum and maximum density value. 

"""
function plot_overview(skyrmion; iso_value = 0.5)

    x = skyrmion.grid.x

    BD = Baryon(skyrmion, density = true)
    (bdmin, bdmax) = extrema(BD)
    Biso_value = (bdmax - bdmin)/2.0
    if Biso_value == 0
        error("Cannot plot. You are likely plotting the vacuum state.")
    end

    BDmesh = getmesh(BD, Biso_value, x)

    skcolormap = make_color_map(skyrmion, BDmesh)

    the_aspect = get_aspect_based_on_density(BDmesh)


    fig = Figure()

    ax = Axis3(fig[1, 1:2], aspect = the_aspect, title = "Baryon density")

    Makie.mesh!(ax, BDmesh, color = skcolormap, shading = NoShading)

    ax11 = Axis3(
        fig[2, 1],
        title = "Isosurface ϕ" * string(1) * " = ±" * string(iso_value),
        aspect = (skyrmion.grid.lp[1], skyrmion.grid.lp[2], skyrmion.grid.lp[3]),
    )
    ax12 = Axis3(
        fig[2, 2],
        title = "Isosurface ϕ" * string(2) * " = " * string(iso_value),
        aspect = (skyrmion.grid.lp[1], skyrmion.grid.lp[2], skyrmion.grid.lp[3]),
    )
    ax21 = Axis3(
        fig[3, 1],
        title = "Isosurface ϕ" * string(3) * " = " * string(iso_value),
        aspect = (skyrmion.grid.lp[1], skyrmion.grid.lp[2], skyrmion.grid.lp[3]),
    )
    ax22 = Axis3(
        fig[3, 2],
        title = "Isosurface ϕ" * string(4) * " = " * string(iso_value),
        aspect = (skyrmion.grid.lp[1], skyrmion.grid.lp[2], skyrmion.grid.lp[3]),
    )

    if Makie.current_backend() == CairoMakie

        field_mesh = [
            getmesh(skyrmion.pion_field[:, :, :, a], iso_value, skyrmion.grid.x) for a = 1:4
        ]
        Makie.mesh!(ax11, field_mesh[1], shading = NoShading)
        Makie.mesh!(ax12, field_mesh[2], shading = NoShading)
        Makie.mesh!(ax21, field_mesh[3], shading = NoShading)
        Makie.mesh!(ax22, field_mesh[4], shading = NoShading)

    else

        volume!(
            ax11,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 1],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = iso_value,
        )
        volume!(
            ax11,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 1],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = -iso_value,
        )

        volume!(
            ax12,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 2],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = iso_value,
        )
        volume!(
            ax12,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 2],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = -iso_value,
        )

        volume!(
            ax21,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 3],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = iso_value,
        )
        volume!(
            ax21,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 3],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = -iso_value,
        )

        volume!(
            ax22,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 4],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = iso_value,
        )
        volume!(
            ax22,
            x[1],
            x[2],
            x[3],
            skyrmion.pion_field[:, :, :, 4],
            algorithm = :iso,
            isorange = 0.2,
            isovalue = -iso_value,
        )

    end

    return fig

end



"""
    plot_baryon_density(skyrmion; iso_value = 0.25*(max(BD) + min(BD)), juggling = false, kwargs...)
    
Plots an isosurface of constant baryon density, with value `iso_value` (where `BD` which determines the default value is the baryon density of `skyrmion` over its grid). The isosurface is coloured to reveal the pion field structure as originally described in [Manton2012](@cite).

Can use a _juggling ball_ colouring scheme by setting `juggling = true`.

# Optional argument

Can accept any arguments used in `Axis3` from the `Makie` package.

# References

* [Manton2012](@cite) Manton. _Classical Skyrmions—static solutions and dynamics_ (2012).

"""
function plot_baryon_density(skyrmion; juggling = false, iso_value = 0.0, kwargs...)

    x = skyrmion.grid.x

    BD = Baryon(skyrmion, density = true)
    (bdmin, bdmax) = extrema(BD)

    if bdmax == bdmin == 0
        error("There is nothing to plot. You are likely trying to plot the vacuum.")
    end

    if iso_value == 0.0
        iso_value = (bdmax + bdmin)/4.0
    end

    if iso_value > bdmax || iso_value < bdmin
        error(
            "Your iso_value, ",
            iso_value,
            " is out of range. The baryon density of your skymion has a minimum ",
            bdmin,
            " and maximum ",
            bdmax,
        )
        return
    end
    println("iso_value is ", iso_value)

    BDmesh = getmesh(BD, iso_value, x)

    if juggling
        skcolormap = make_color_map_juggle(skyrmion, BDmesh)
    else
        skcolormap = make_color_map(skyrmion, BDmesh)
    end

    the_aspect = get_aspect_based_on_density(BDmesh)

    fig = Figure()

    ax = Axis3(fig[1, 1], aspect = the_aspect; kwargs...)

    Makie.mesh!(ax, BDmesh, color = skcolormap, shading = NoShading)


    return fig

end

function get_aspect_based_on_density(mesh)

    the_extrema = [extrema([mesh[a][1][b] for a = 1:size(mesh.faces)[1]]) for b = 1:3]
    return ntuple(a -> the_extrema[a][2] - the_extrema[a][1], Val(3))

end

function getmesh(a_density, iso_value, x)

    points, faces = isosurface(a_density, MarchingCubes(iso = iso_value), -1:1, -1:1, -1:1)

    Npts = size(points)[1]
    Nfaces = size(faces)[1]

    pointsaspoints = fill(Point3f(0.0, 0.0, 0.0), Npts);
    facesagain = TriangleFace{Int}[faces[a] for a = 1:Nfaces]

    for i = 1:Npts
        pointsaspoints[i] = Point3f(points[i][1], points[i][2], points[i][3])
    end

    return GeometryBasics.Mesh(pointsaspoints, facesagain)

end

function make_color_map_juggle(skyrmion, BDmesh)

    x = skyrmion.grid.x

    ϕinterp = [
        linear_interpolation((x[1], x[2], x[3]), skyrmion.pion_field[:, :, :, a]) for
        a = 1:3
    ]

    Npts = size(coordinates(BDmesh))[1]

    the_color_function = zeros(RGB{Float64}, Npts)

    cube_faces = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]
    face_color_map = [
        Colors.RGB(1.0, 0.0, 0.0),
        Colors.RGB(0.0, 1.0, 0.0),
        Colors.RGB(0.0, 0.0, 1.0),
        Colors.RGB(1.0, 1.0, 0.0),
        Colors.RGB(0.0, 1.0, 1.0),
        Colors.RGB(1.0, 0.0, 1.0),
    ]

    for (a, vertex) in enumerate(coordinates(BDmesh))

        pion_field_on_mesh = [ϕinterp[b](vertex[1], vertex[2], vertex[3]) for b = 1:3]

        distances = zeros(6)
        for c = 1:6, b = 1:3
            distances[c] += (pion_field_on_mesh[b] - cube_faces[c][b])^2
        end

        the_color_function[a] = face_color_map[findmin(distances)[2]]

    end

    return the_color_function
end

function make_color_map(skyrmion, BDmesh)

    x = skyrmion.grid.x

    ϕinterp = [
        linear_interpolation((x[1], x[2], x[3]), skyrmion.pion_field[:, :, :, a]) for
        a = 1:3
    ]

    Npts = size(coordinates(BDmesh))[1]

    pion_field_on_mesh = zeros(3, Npts)
    for (a, vertex) in enumerate(coordinates(BDmesh)), b = 1:3
        pion_field_on_mesh[b, a] = ϕinterp[b](vertex[1], vertex[2], vertex[3])
    end

    (minp, maxp) = extrema(pion_field_on_mesh[3, :])
    p3color = (pion_field_on_mesh[3, :] .- minp) ./ (maxp - minp)

    return [
        Colors.HSL(
            360*(atan.(-pion_field_on_mesh[2, a], -pion_field_on_mesh[1, a]) .+ pi) ./
            (2pi),
            1,
            p3color[a],
        ) for a = 1:Npts
    ]

end
