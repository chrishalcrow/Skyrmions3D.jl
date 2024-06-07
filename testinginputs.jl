
using Skyrmions3D
using GLMakie
GLMakie.activate!()
Makie.inline!(false)
#display(plot( rand(3), rand(3) ) )

#nuc = Skyrmion( [40, 40, 40], [0.2, 0.2, 0.2])
#set_metric!(nuc,1.0)
#overview(nuc)

nuc = Skyrmion( [40, 40, 40], [0.2, 0.2, 0.2], mpi = 1.0, Fpi=186, ee=4.0, boundary_conditions="periodic")


overview(nuc)

set_mpi!(nuc, 0.5)
set_Fpi!(nuc, 100)
set_ee!(nuc, 6.5)
set_lattice!(nuc, [60,60,60], [0.2,0.2,0.2])
set_periodic!(nuc)

overview(nuc)



p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0;
q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0;
f4(r) = pi*exp( -(r.^3)./12.0 )

make_rational_map!(nuc, p4, q4, f4)

Baryon(nuc)
Energy(nuc)

# Without a profile function, the make_rational_map! function will find an OK approximate

make_rational_map!(nuc, p4, q4; baryon=4)
Energy(nuc)