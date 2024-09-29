function activate_GLMakie()
	GLMakie.activate!()
end


"""
	interactive_flow(my_skyrmion; iso_value=2.0, kwargs...)
    
Initialises a interface, used for dynamics. Requires GLMakie to be active. Once activated, can interactively apply a gradient flow or arrested newton flow to the initial skyrme field `my_skyrmion'.

"""
function interactive_flow(my_skyrmion; iso_value=2.0, kwargs... )

	which_flow = 1

	skd = 0.0 .* similar(my_skyrmion.pion_field)

	x = my_skyrmion.x

	active_button_color = :white
	active_button_strokecolor = :red

	not_active_button_color = :lightgray
	not_active_button_strokecolor = :gray

	BD = Baryon(my_skyrmion, density=true)

	BDmesh = getmesh(BD, iso_value, x)
	skcolormap = make_color_map(my_skyrmion, BDmesh)

	obBD = Observable(BDmesh)
	obCM = Observable(skcolormap)

	obBD[] = BDmesh
	obCM[] = skcolormap


	fig = Figure(clear=true)


	g_skyrmion = fig[1,1] =GridLayout()
	g_info = fig[1,2] = GridLayout(tellwidth=true, tellheight=true, valign=:top)

	g_plotting_options = g_info[1,1] = GridLayout(tellwidth=false, tellheight=false, halign=:left, valign=:top)
	g_dynamics = g_info[2,1] = GridLayout(tellwidth=false, tellheight=false, halign=:left, valign=:center)
	g_export = g_info[3,1] = GridLayout(tellwidth=false, tellheight=false, halign=:left, valign=:bottom)


	Label(g_plotting_options[1,1:2], text="Plotting options", font = :bold, halign=:left)
	Label(g_plotting_options[2,1], text="Isosurface:  ", halign=:left)
	bd_tb = Textbox(g_plotting_options[2,2]; stored_string=string(iso_value),validator = Float64, tellwidth=false, halign=:left)

	Label(g_dynamics[1,1:2], text="Dynamics", font = :bold, halign=:left)
	
	gf_button = Button(g_dynamics[2,1]; label="Grad flow", tellwidth=true, halign=:left, font="Courier")
	an_button = Button(g_dynamics[2,2]; label="Ar N flow", tellwidth=true, halign=:left, font="Courier")
	newt_button = Button(g_dynamics[3,1]; label="Newt flow", tellwidth=true, halign=:left, font="Courier")
	#full_button = Button(g_dynamics[3,2]; label="no  no no", tellwidth=true, halign=:left, font="Courier")
	
	gf_button.strokecolor = active_button_strokecolor
	gf_button.buttoncolor = active_button_color

	an_button.strokecolor = not_active_button_strokecolor
	an_button.buttoncolor = not_active_button_color

	newt_button.strokecolor = not_active_button_strokecolor
	newt_button.buttoncolor = not_active_button_color

	#full_button.strokecolor = not_active_button_strokecolor
	#full_button.buttoncolor = not_active_button_color

	Label(g_dynamics[4,1], text="Steps: ",halign=:right)
	flow_tb = Textbox(g_dynamics[4,2]; placeholder="100",validator = Int64, tellwidth=false, halign=:left)

	Label(g_dynamics[5,1], text="dt: ",halign=:right)
	dt_tb = Textbox(g_dynamics[5,2]; placeholder=string(0.0),validator = Float64, tellwidth=false, halign=:left)

	flow_button = Button(g_dynamics[6,1:2]; label="Flow!", tellwidth=false, halign=:center, font=:bold, width=220)


	export_tb = Textbox(g_export[1,1]; stored_string="image_1.png", tellwidth=false, halign=:left)
	export_button = Button(g_export[1,2]; label="Picture", tellwidth=false, halign=:right, font=:bold)



	the_extrema = [ extrema( [ BDmesh[a][1][b] for a in 1:size(BDmesh)[1] ] ) for b in 1:3 ]
	the_aspect = ( the_extrema[1][2] - the_extrema[1][1], the_extrema[2][2] - the_extrema[2][1], the_extrema[3][2] - the_extrema[3][1] )

	ax = Axis3(g_skyrmion[1,1], aspect = the_aspect, tellwidth=true; kwargs...)#, height=950, width=950 )


	Makie.mesh!(ax,obBD,
				color = obCM,
				shading=NoShading
	)

	colsize!(fig.layout, 2, Fixed(220))


	rowgap!(g_dynamics,10.0)
	rowgap!(g_plotting_options,10.0)

	colgap!(g_dynamics,5.0)


	on(bd_tb.stored_string) do s

		BD = Baryon(my_skyrmion,density=true)
		BDmesh = getmesh(BD, parse(Float64,s), x)
		skcolormap = make_color_map(my_skyrmion, BDmesh)

		obBD[] = BDmesh
		obCM[] = skcolormap

		sleep(0.001)
	end

	on(gf_button.clicks) do clicks
		
		which_flow = 1

		gf_button.strokecolor = active_button_strokecolor
		gf_button.buttoncolor = active_button_color

		an_button.strokecolor = not_active_button_strokecolor
		an_button.buttoncolor = not_active_button_color

		newt_button.strokecolor = not_active_button_strokecolor
		newt_button.buttoncolor = not_active_button_color

		#full_button.strokecolor = not_active_button_strokecolor
		#full_button.buttoncolor = not_active_button_color

	end
	on(an_button.clicks) do clicks
		
		which_flow = 2

		an_button.strokecolor = active_button_strokecolor
		an_button.buttoncolor = active_button_color

		gf_button.strokecolor = not_active_button_strokecolor
		gf_button.buttoncolor = not_active_button_color

		newt_button.strokecolor = not_active_button_strokecolor
		newt_button.buttoncolor = not_active_button_color

		#full_button.strokecolor = not_active_button_strokecolor
		#full_button.buttoncolor = not_active_button_color

	end
	on(newt_button.clicks) do clicks
		
		which_flow = 3

		newt_button.strokecolor = active_button_strokecolor
		newt_button.buttoncolor = active_button_color

		gf_button.strokecolor = not_active_button_strokecolor
		gf_button.buttoncolor = not_active_button_color

		an_button.strokecolor = not_active_button_strokecolor
		an_button.buttoncolor = not_active_button_color

		#full_button.strokecolor = not_active_button_strokecolor
		#full_button.buttoncolor = not_active_button_color

	end
	#=on(full_button.clicks) do clicks
		
		which_flow = 4

		full_button.strokecolor = active_button_strokecolor
		full_button.buttoncolor = active_button_color

		gf_button.strokecolor = not_active_button_strokecolor
		gf_button.buttoncolor = not_active_button_color

		an_button.strokecolor = not_active_button_strokecolor
		an_button.buttoncolor = not_active_button_color

		newt_button.strokecolor = not_active_button_strokecolor
		newt_button.buttoncolor = not_active_button_color

	end=#

	on(export_button.clicks) do clicks

		save( to_value(export_tb.displayed_string), fig)

	end


	on(flow_button.clicks) do clicks

		dt = parse(Float64,to_value(dt_tb.displayed_string))

		total_runs = deepcopy(parse(Int64, to_value(flow_tb.displayed_string)))
		iso_val = deepcopy(parse(Float64,to_value(bd_tb.displayed_string)))

		if which_flow == 1
			if dt == 0.0
				gradient_flow!(my_skyrmion; steps=total_runs)
			else
				gradient_flow!(my_skyrmion; steps=total_runs, dt=dt )
			end
		elseif which_flow == 2
			if dt == 0.0
				arrested_newton_flow!(my_skyrmion; steps=total_runs )
			else
				arrested_newton_flow!(my_skyrmion; steps=total_runs, dt = dt )
			end
			
		elseif which_flow == 3
			if dt == 0.0
				newton_flow!(my_skyrmion; ϕd = skd, steps=total_runs )
			else
				newton_flow!(my_skyrmion; ϕd = skd, steps=total_runs, dt = dt  )
			end
		end

		BD = Baryon(my_skyrmion,density=true)
		BDmesh = getmesh(BD, iso_val, x)
		skcolormap = make_color_map(my_skyrmion, BDmesh)

		obBD[] = BDmesh
		obCM[] = skcolormap

		sleep(0.001)

	end

	display(GLMakie.Screen(),fig)

	return

end
