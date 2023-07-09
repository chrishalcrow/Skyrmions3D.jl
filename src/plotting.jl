function imusingCairo()
	CairoMakie.activate!()
end

function imusingGL()
	GLMakie.activate!()
end

"""
    plot_field!(skyrmion; component=3, iso_value=0.5 )
    
Plots an isosurface of constant value, `skyrme_field[component] = iso_value`

"""
function plot_field(skyrmion; component=3, iso_value = 0.5)
    
    fig = Figure()

    ax = Axis3(fig[1,1], 
    	title= "Isosurface ϕ" * string(component) * " = " * string(iso_value),
    	aspect = (skyrmion.lp[1],skyrmion.lp[2],skyrmion.lp[3])
    	)

	volume!(ax, skyrmion.phi[:,:,:,component], algorithm = :iso, isorange = 0.2, isovalue = iso_value )

	return fig

end


"""
    plot_baryon_density(skyrmion; iso_value = 0.5, juggling = false, kwargs...)
    
Plots an isosurface of constant baryon density, with value `iso_value`, coloured to reveal the pion field structure, originally described in [].

Can use a _juggling ball_ colouring scheme by setting `juggling = true`.

# Optional argument

Can accept any arguments used in `Axis3` from the `Makie` package. See more: [].

"""
function plot_baryon_density(skyrmion; juggling = false, iso_value = 0.5, kwargs...)
    
	x = skyrmion.x
	lp = skyrmion.lp

	BD = Baryon(skyrmion,density=true)
	bdmax = maximum(BD)
    bdmin = minimum(BD)

    if iso_value > bdmax || iso_value < bdmin
        println("ERROR: Your iso_value is out of range. The baryon density of your skymion has a minimum ", bdmin, " and maximum ", bdmax)
        return
    end

	BDmesh = getmesh(BD, iso_value, x)
	obBD = Observable(BDmesh)

	if juggling == false
		skcolormap = make_color_map(skyrmion, BDmesh)
	else
		skcolormap = make_color_map_juggle(skyrmion, BDmesh)
	end

    fig = Figure()

	the_extrema = [ extrema( [ BDmesh[a][1][b] for a in 1:size(BDmesh)[1] ] ) for b in 1:3 ]
	the_aspect = ( the_extrema[1][2] - the_extrema[1][1], the_extrema[2][2] - the_extrema[2][1], the_extrema[3][2] - the_extrema[3][1] )


	ax = Axis3(fig[1,1], aspect = the_aspect ; kwargs...)

        Makie.mesh!(ax,BDmesh,
        	color = skcolormap,
        	shading=false,

        	)


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


"""
	interactive_flow(my_skyrmion; iso_value=2.0, kwargs...)
    
Initialises a interface, used for dynamics. Requires GLMakie to be active. Once activated, can interactively apply a gradient flow or arrested newton flow to the initial skyrme field `my_skyrmion'.

"""
function interactive_flow(my_skyrmion; iso_value=2.0, kwargs... )

	skd = 0.0 .* similar(my_skyrmion.phi)

	juggling=false
	x = my_skyrmion.x
	lp = my_skyrmion.lp

	BD = Baryon(my_skyrmion, density=true)

	BDmesh = getmesh(BD, iso_value, x)
	skcolormap = make_color_map(my_skyrmion, BDmesh)
	#the_energy = Energy(my_skyrmion)
	#the_baryon = Baryon(my_skyrmion)

	obBD = Observable(BDmesh)
	obCM = Observable(skcolormap)
	#displayed_energy = Observable(the_energy)
	#displayed_baryon = Observable(the_baryon)

	obBD[] = BDmesh
	obCM[] = skcolormap
	#displayed_energy[] = the_energy
	#displayed_baryon[] = the_baryon


	fig = Figure(clear=true)


	g_skyrmion = fig[1,1] =GridLayout()
	g_info = fig[1,2] = GridLayout(tellwidth=false, tellheight=false, valign=:top)

	#g_attributes = g_info[1,1] = GridLayout(tellwidth=false, tellheight=true, halign=:left)
	g_plotting_options = g_info[2,1] = GridLayout(tellwidth=false, tellheight=true, halign=:left)
	g_dynamics = g_info[3,1] = GridLayout(tellwidth=false, tellheight=true, halign=:left)

	#=Label(g_attributes[1,1:2], text="Skyrmion attributes", font = :bold, halign=:left)

	Label(g_attributes[2,1], text="Baryon: ",  halign=:left)
	Label(g_attributes[2,2], text=string(round(to_value(displayed_baryon),digits=4)), halign=:left)

	Label(g_attributes[3,1], text="Energy: ",  halign=:left)
	Label(g_attributes[3,2], text=string(round(to_value(displayed_energy),digits=4)), halign=:left)

	Label(g_attributes[3,2], text=string(round(to_value(displayed_energy),digits=4)), halign=:left, visible=false)
	=#

	Label(g_plotting_options[1,1:2], text="Plotting options", font = :bold, halign=:left)
	Label(g_plotting_options[2,1], text="Isosurface:  ", halign=:left)
	bd_tb = Textbox(g_plotting_options[2,2]; stored_string=string(iso_value),validator = Float64, tellwidth=false, halign=:left)

	Label(g_dynamics[1,1:2], text="Dynamics", font = :bold, halign=:left)
	menu1 = Menu(g_dynamics[2, 1:2], options = ["Gradient flow", "Arrested NF"], prompt="Flow algorithm...")

	Label(g_dynamics[3,1], text="Steps: ",halign=:right)
	flow_tb = Textbox(g_dynamics[3,2]; placeholder="100",validator = Int64, tellwidth=false, halign=:left)

	Label(g_dynamics[4,1], text="dt: ",halign=:right)
	dt_tb = Textbox(g_dynamics[4,2]; placeholder=string(my_skyrmion.ls[1]*my_skyrmion.ls[2]/80.0),validator = Float64, tellwidth=false, halign=:left)

	flow_button = Button(g_dynamics[5,1:2]; label="Flow!", tellwidth=false, halign=:center, font=:bold, width=180)


	the_extrema = [ extrema( [ BDmesh[a][1][b] for a in 1:size(BDmesh)[1] ] ) for b in 1:3 ]
	the_aspect = ( the_extrema[1][2] - the_extrema[1][1], the_extrema[2][2] - the_extrema[2][1], the_extrema[3][2] - the_extrema[3][1] )

	ax = Axis3(g_skyrmion[1,1], aspect = the_aspect, tellwidth=false; kwargs...)#, height=950, width=950 )


	Makie.mesh!(ax,obBD,
				color = obCM,
				shading=false
	)	#colgap!(fig.layout,1,0.0)
	#rowgap!(g_info,0.0)
	#colsize!(fig.layout, 1, Auto(1.5))
	colsize!(fig.layout, 2, Fixed(180))

	rowgap!(g_info,20.0)
	rowgap!(g_dynamics,10.0)
	rowgap!(g_plotting_options,10.0)
	#rowgap!(g_attributes,10.0)

	colgap!(g_dynamics,5.0)

	

	on(bd_tb.stored_string) do s

		iso_val = parse(Float64,to_value(bd_tb.displayed_string))

		BD = Baryon(my_skyrmion,density=true)
		BDmesh = getmesh(BD, parse(Float64,s), x)
		skcolormap = make_color_map(my_skyrmion, BDmesh)

		obBD[] = BDmesh
		obCM[] = skcolormap

		sleep(0.001)
	end

	#on(menu1.selection) do e
		#if menu1.selection == "Arrested NF"
			#println("oh!")
			#dt_tb.displayed_string = 0.5
		#end	
	#end

	on(flow_button.clicks) do clicks

		#flow_tb.stored_string = flow_tb.displayed_string
		#bd_tb.stored_string = bd_tb.displayed_string

		dt = parse(Float64,to_value(dt_tb.displayed_string))
		println("dt = ", dt)

		total_runs = deepcopy(parse(Int64, to_value(flow_tb.displayed_string)))
		iso_val = deepcopy(parse(Float64,to_value(bd_tb.displayed_string)))


		if to_value(menu1.selection) == "Gradient flow"
			flow!(my_skyrmion; steps=total_runs, dt=dt )
		elseif to_value(menu1.selection) == "Arrested NF"
			ANFflow!(my_skyrmion,skd, dt, total_runs )
		end

		BD = Baryon(my_skyrmion,density=true)
		BDmesh = getmesh(BD, iso_val, x)
		skcolormap = make_color_map(my_skyrmion, BDmesh)

		#displayed_energy[] = Energy(my_skyrmion)
		#displayed_baryon[] = sum(BD)
		obBD[] = BDmesh
		obCM[] = skcolormap

		sleep(0.001)

	end

	display(GLMakie.Screen(),fig)

	return

end




