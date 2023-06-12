function imusingnotebook()
	CairoMakie.activate!()
end

function imusingterminal()
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


function plot_baryon_density(skyrmion; iso_value = 0.5, kwargs...)
    
	x = skyrmion.x
	lp = skyrmion.lp

	BD = BaryonD(skyrmion)

	BDmesh = getmesh(BD, iso_value, x)
	skcolormap = make_color_map(skyrmion, BDmesh)

    fig = Figure()

	ax = Axis3(fig[1,1],aspect=(lp[1],lp[2],lp[3]), ;kwargs...)

        Makie.mesh!(ax,BDmesh,
        	color = skcolormap,
        	shading=false,

        	)


	return fig

end


