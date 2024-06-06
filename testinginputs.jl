
using Skyrmions3D
using GLMakie
GLMakie.activate!()
Makie.inline!(false)
display(plot( rand(3), rand(3) ) )

nuc = Skyrmion( [40, 40, 40], [0.2, 0.2, 0.2])
set_metric_var!(nuc,metric_var=2.0)
overview(nuc)