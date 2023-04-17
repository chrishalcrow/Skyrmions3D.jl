function plot_skyrmion(phi)
    
	Makie.volume(phi.ED,algorithm = :iso, isorange = 0.5, isovalue = 7.0)

end