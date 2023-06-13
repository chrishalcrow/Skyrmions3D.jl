function imusingnotebook()
	CairoMakie.activate!()
end



function imusingterminal()
	Makie.inline!(false) 
	GLMakie.activate!()
end

function plot_skyrmion(phi)
    
	Makie.volume(phi.ED,algorithm = :iso, isorange = 0.5, isovalue = 7.0)

end

function plot_field(skyrmion; component=3, iso_value = 0.5)
    
    fig = Figure()



    ax = Axis3(fig[1,1], 
    	title= "Isosurface ϕ" * string(component) * " = " * string(iso_value),
    	aspect = (skyrmion.lp[1],skyrmion.lp[2],skyrmion.lp[3])
    	)

	volume!(ax, skyrmion.phi[:,:,:,component], algorithm = :iso, isorange = 0.2, isovalue = iso_value )

	return fig

end


function getmesh(a_density, iso_value,x)

	points,faces = isosurface(a_density, MarchingCubes(iso=iso_value), origin= [ x[1][1], x[2][1], x[3][1] ],  widths = [ x[1][end] - x[1][1], x[2][end] - x[2][1], x[3][end] - x[3][1] ]  )#, samples=(24,24,24) )

	Npts = size(points)[1]
	Nfaces = size(faces)[1]

    pointsaspoints = fill(Point3f(1.0,1.0,1.0),Npts);
	facesagain = TriangleFace{Int}[faces[a] for a in 1:Nfaces]

	for i in 1:Npts
	    pointsaspoints[i] = Point3f(points[i][1], points[i][2], points[i][3])
	end

	return GeometryBasics.Mesh(meta(pointsaspoints), facesagain)

end

function make_color_map_juggle(skyrmion, BDmesh)

	x = skyrmion.x

    ϕinterp = [ linear_interpolation((x[1],x[2],x[3]), skyrmion.phi[:,:,:,a] )  for a in 1:3 ]
    Npts = size(coordinates(BDmesh))[1]
    phi_on_mesh = zeros(3, Npts)

	the_color_function = zeros(RGB{Float64}, Npts)

	cube_faces = [ [1.0,0.0,0.0], [-1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,-1.0,0.0], [0.0,0.0,1.0], [0.0,0.0,-1.0] ]
	distances = zeros(6)

	which_is_min = 1

    a=1
	for vertex in coordinates(BDmesh)

	    for b in 1:3
	        phi_on_mesh[b,a] = ϕinterp[b](vertex[1],vertex[2],vertex[3])
	    end

		for c in 1:6
			distances[c] = 0.0
			for b in 1:3
				distances[c] += (phi_on_mesh[b,a] - cube_faces[c][b])^2
			end
		end

		which_is_min = findmin( distances )[2]

		if which_is_min == 1
			the_color_function[a] = Colors.RGB(1.0,0.0,0.0)
		elseif which_is_min == 2
			the_color_function[a] = Colors.RGB(0.0,1.0,0.0)
		elseif which_is_min == 3
			the_color_function[a] = Colors.RGB(0.0,0.0,1.0)
		elseif which_is_min == 4
			the_color_function[a] = Colors.RGB(1.0,1.0,0.0)
		elseif which_is_min == 5
			the_color_function[a] = Colors.RGB(0.0,1.0,1.0)
		else
			the_color_function[a] = Colors.RGB(1.0,0.0,1.0)
		end

	    a += 1
	end

	return the_color_function 
end



function make_color_map(skyrmion, BDmesh)

	x = skyrmion.x

    ϕinterp = [ linear_interpolation((x[1],x[2],x[3]), skyrmion.phi[:,:,:,a] )  for a in 1:3 ]
    Npts = size(coordinates(BDmesh))[1]
    phi_on_mesh = zeros(3, Npts)

    a=1
	for vertex in coordinates(BDmesh)
	    for b in 1:3
	        phi_on_mesh[b,a] = ϕinterp[b](vertex[1],vertex[2],vertex[3])
	    end
	    a += 1
	end

	maxp = maximum(phi_on_mesh[3,:])
	minp = minimum(phi_on_mesh[3,:])

   p3color = (phi_on_mesh[3,:] .- minp)./(maxp - minp)

	return [ Colors.HSL( 360*(atan.(phi_on_mesh[1,a], -phi_on_mesh[2,a]) .+ pi)./(2pi) , 1, p3color[a] ) for a in 1:Npts ];

end


function plot_baryon_density(skyrmion; juggling = false, iso_value = 0.5, kwargs...)
    
	x = skyrmion.x
	lp = skyrmion.lp

	BD = Baryon(skyrmion,density=true)

	BDmesh = getmesh(BD, iso_value, x)
	if juggling == false
		skcolormap = make_color_map(skyrmion, BDmesh)
	else
		skcolormap = make_color_map_juggle(skyrmion, BDmesh)
	end

    fig = Figure()

	the_extrema = [ extrema( [ BDmesh[a][1][b] for a in 1:7840 ] ) for b in 1:3 ]
	the_aspect = ( the_extrema[1][2] - the_extrema[1][1], the_extrema[2][2] - the_extrema[2][1], the_extrema[3][2] - the_extrema[3][1] )

	ax = Axis3(fig[1,1], aspect = the_aspect ;kwargs...)

        Makie.mesh!(ax,BDmesh,
        	color = skcolormap,
        	shading=false,

        	)


	return fig

end


