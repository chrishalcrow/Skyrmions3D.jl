"""
    overview(skyrmion)

Displays an overview of `skyrmion`'s properties.

"""
function overview(sk)

    hbar = 197

    println("This skyrmion is on a ", sk.lp[1],"x",sk.lp[2],"x",sk.lp[3]," grid, with lattice spacing [", sk.ls[1],", ", sk.ls[2],", ", sk.ls[3], "]. ")
    println("The boundary conditions are "*sk.boundary_conditions)
    println()
    println("            m = ", round(sk.mpi, sigdigits=5) ) 
    println("Baryon number = ", round(Baryon(sk), sigdigits=5) ) 
    println("Using metric parameter = ", sk.metric, ", Energy = ", round(12*pi*pi*Energy(sk), sigdigits=5) ) 
    println("   Baryon rms = ", round(rms_baryon(sk), sigdigits=5) ) 
    println()
    println("With physical constants Fpi = ", sk.Fpi, " and e = ", sk.ee,",")
    println("the energy and length units are ", round(sk.Fpi/(4*sk.ee),sigdigits=4), " MeV and ", round(hbar*2/(sk.ee*sk.Fpi),sigdigits=4), " fm." )
    println("So physical E = ", round(12*pi*pi*Energy(sk)*sk.Fpi/(4*sk.ee), sigdigits=5), " MeV." )
    println("   Baryon rms = ", round(rms_baryon(sk)*2*hbar/(sk.ee*sk.Fpi), sigdigits=5), " fm." )


end


"""
    Energy(skyrmion; density=false)

Compute energy of `skyrmion`.

Set 'density = true' to output the energy density and moment to `n` to calculate the nth moment of the energy density.

"""
function Energy(sk; density=false, moment=0)

    ED = zeros(sk.lp[1], sk.lp[2], sk.lp[3])

    get_energy_density!(ED,sk,moment=moment)

    if density == false
        engtot = sum(ED)*sk.ls[1]*sk.ls[2]*sk.ls[3]
        if sk.physical == false
            return engtot/(12.0*pi^2)
        else
            if moment == 0
                return (engtot*sk.Fpi/(4.0*sk.ee), "MeV" )
            else
                return (engtot*sk.Fpi*(197.327*2.0/(sk.ee*sk.Fpi))/(4.0*sk.ee), "MeV fm^"*string(moment))
            end
        end
    else
        return ED
    end

end 

function get_energy_density!(density, sk ;moment=0)

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            dp = getDP(sk ,i, j, k )
            rm = sqrt( sk.x[1][i]^2 + sk.x[2][j]^2 + sk.x[3][k]^2 )^moment

            density[i,j,k] = engpt(dp,sk.pion_field[i,j,k,:],sk.mpi,sk.metric) * rm

        
        end
    end

end

function engpt(dp,p,mpi,alpha)
    L3_1 = p[4] * dp[1,3] - p[3] * dp[1,4] + p[1] * dp[1,2] - p[2] * dp[1,1]
    L3_2 = p[4] * dp[2,3] - p[3] * dp[2,4] + p[1] * dp[2,2] - p[2] * dp[2,1]
    L3_3 = p[4] * dp[3,3] - p[3] * dp[3,4] + p[1] * dp[3,2] - p[2] * dp[3,1]

    LB12 = dp[1,4] * dp[2,3] - dp[1,3] * dp[2,4] + dp[1,1] * dp[2,2] - dp[1,2] * dp[2,1]
    LB13 = dp[1,4] * dp[3,3] - dp[1,3] * dp[3,4] + dp[1,1] * dp[3,2] - dp[1,2] * dp[3,1]
    LB23 = dp[2,4] * dp[3,3] - dp[2,3] * dp[3,4] + dp[2,1] * dp[3,2] - dp[2,2] * dp[3,1]

    e_0 = 2*mpi^2*(1 - p[4])
    e_2 = (dp[1,1]^2 + dp[1,2]^2 + dp[1,3]^2 + dp[1,4]^2 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2 + dp[3,4]^2)
    e_4 = (dp[1,4]^2*dp[2,1]^2 + dp[1,4]^2*dp[2,2]^2 + dp[1,4]^2*dp[2,3]^2 + dp[1,1]^2*(dp[2,2]^2 + dp[2,3]^2) - 2*dp[1,1]*dp[1,4]*dp[2,1]*dp[2,4] + dp[1,1]^2*dp[2,4]^2 + dp[1,4]^2*dp[3,1]^2 + dp[2,2]^2*dp[3,1]^2 + dp[2,3]^2*dp[3,1]^2 + dp[2,4]^2*dp[3,1]^2 - 2*dp[2,1]*dp[2,2]*dp[3,1]*dp[3,2] + dp[1,1]^2*dp[3,2]^2 + dp[1,4]^2*dp[3,2]^2 + dp[2,1]^2*dp[3,2]^2 + dp[2,3]^2*dp[3,2]^2 + dp[2,4]^2*dp[3,2]^2 - 2*dp[2,1]*dp[2,3]*dp[3,1]*dp[3,3] - 2*dp[2,2]*dp[2,3]*dp[3,2]*dp[3,3] + dp[1,1]^2*dp[3,3]^2 + dp[1,4]^2*dp[3,3]^2 + dp[2,1]^2*dp[3,3]^2 + dp[2,2]^2*dp[3,3]^2 + dp[2,4]^2*dp[3,3]^2 - 2*(dp[1,1]*dp[1,4]*dp[3,1] + dp[2,4]*(dp[2,1]*dp[3,1] + dp[2,2]*dp[3,2] + dp[2,3]*dp[3,3]))*dp[3,4] + (dp[1,1]^2 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2)*dp[3,4]^2 + dp[1,3]^2*(dp[2,1]^2 + dp[2,2]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,4]^2) + dp[1,2]^2*(dp[2,1]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,3]^2 + dp[3,4]^2) - 2*dp[1,2]*(dp[1,1]*(dp[2,1]*dp[2,2] + dp[3,1]*dp[3,2]) + dp[1,3]*(dp[2,2]*dp[2,3] + dp[3,2]*dp[3,3]) + dp[1,4]*(dp[2,2]*dp[2,4] + dp[3,2]*dp[3,4])) - 2*dp[1,3]*(dp[1,1]*(dp[2,1]*dp[2,3] + dp[3,1]*dp[3,3]) + dp[1,4]*(dp[2,3]*dp[2,4] + dp[3,3]*dp[3,4])))
    e_0_star = mpi^2*((p[3])^2)
    e_2_star = (L3_1)^2 + (L3_2)^2 + (L3_3)^2
    e_4_star = (LB12)^2 + (LB13)^2 + (LB23)^2


    return e_0 + e_2 + e_4 + (alpha^2 - 1)*e_0_star + (alpha^2 - 1)*e_2_star + (alpha^2 - 1)*e_4_star

end

"""
    Baryon(skyrmion; density=false)

Compute baryon number of `skyrmion`.

Set 'density = true' to output the baryon density and moment to `n` to calculate the nth moment of the baryon density.

"""
function Baryon(sk; density=false, moment = 0, component = 0)

    BD = zeros(sk.lp[1],sk.lp[2],sk.lp[3])
    
    get_baryon_density!(BD,sk,moment=moment, component = component)

    if density == false
        bartot = sum(BD)*sk.ls[1]*sk.ls[2]*sk.ls[3]/(2.0*pi^2)
        return bartot
    else
        return BD
    end
    
end 

function get_baryon_density!(baryon_density, sk ;moment=0, component = 0)

    vc = zeros(3)

    if component == 0
        vc = [1,1,1]
    elseif component == 1
        vc = [1,0,0]
    elseif  component == 2
        vc = [0,1,0]
    elseif  component == 3
        vc = [0,0,1]
    end    

    Threads.@threads for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            pp = getX(sk,i,j,k)
            dp = getDP(sk ,i, j, k )
            rm = sqrt( vc[1]*sk.x[1][i]^2 + vc[2]*sk.x[2][j]^2 + vc[3]*sk.x[3][k]^2 )^moment
            
            baryon_density[i,j,k] = barypt(dp,pp) * rm

        end
    end

end

function barypt(dp,pp)
    return pp[4]*dp[1,3]*dp[2,2]*dp[3,1] - pp[3]*dp[1,4]*dp[2,2]*dp[3,1] - pp[4]*dp[1,2]*dp[2,3]*dp[3,1] + pp[2]*dp[1,4]*dp[2,3]*dp[3,1] + pp[3]*dp[1,2]*dp[2,4]*dp[3,1] - pp[2]*dp[1,3]*dp[2,4]*dp[3,1] - pp[4]*dp[1,3]*dp[2,1]*dp[3,2] + pp[3]*dp[1,4]*dp[2,1]*dp[3,2] + pp[4]*dp[1,1]*dp[2,3]*dp[3,2] - pp[1]*dp[1,4]*dp[2,3]*dp[3,2] - pp[3]*dp[1,1]*dp[2,4]*dp[3,2] + pp[1]*dp[1,3]*dp[2,4]*dp[3,2] + pp[4]*dp[1,2]*dp[2,1]*dp[3,3] - pp[2]*dp[1,4]*dp[2,1]*dp[3,3] - pp[4]*dp[1,1]*dp[2,2]*dp[3,3] + pp[1]*dp[1,4]*dp[2,2]*dp[3,3] + pp[2]*dp[1,1]*dp[2,4]*dp[3,3] - pp[1]*dp[1,2]*dp[2,4]*dp[3,3] - pp[3]*dp[1,2]*dp[2,1]*dp[3,4] + pp[2]*dp[1,3]*dp[2,1]*dp[3,4] + pp[3]*dp[1,1]*dp[2,2]*dp[3,4] - pp[1]*dp[1,3]*dp[2,2]*dp[3,4] - pp[2]*dp[1,1]*dp[2,3]*dp[3,4] + pp[1]*dp[1,2]*dp[2,3]*dp[3,4]
end


"""
    center_of_mass(skyrmion; density=false, moment=0)

Compute the center of mass of `skyrmion`, based on the energy density.

"""
function center_of_mass(sk)

    com = zeros(3)

    the_engpt = 0.0
    toteng = 0.0

    @inbounds for i in 3:sk.lp[1]-2, j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
    
        dp = getDP(sk ,i, j, k )

        the_engpt = engpt(dp,sk.pion_field[i,j,k,:],sk.mpi,sk.metric)

        com[1] += the_engpt*sk.x[1][i]
        com[2] += the_engpt*sk.x[2][j]
        com[3] += the_engpt*sk.x[3][k]

        toteng += the_engpt
   
    end

    if toteng ≈ 0.0
        return [0.0,0.0,0.0]
    else
        return com/toteng
    end

end 


"""
    rms_baryon(skyrmion)

Compute root mean square charge radius of a skyrmion, using the baryon density.
"""
function rms_baryon(sk; component = 0)

    B0 = Baryon(sk; moment = 0)
    if B0 == 0.0
        return 0.0
    end

    rms = sqrt(Baryon(sk; moment = 2, component=component)/B0)

    if sk.physical == false
        return rms
    else
        return ( rms*197.327*(2.0/(sk.ee*sk.Fpi)), "fm" )
    end

end

function sphericity(my_skyrmion)

    (λ1,λ2,λ3) = sort(eigvals(compute_current(my_skyrmion, label="vMOI")))

    Δ1 = (λ2 - λ1)/(λ2 + λ1)
    Δ2 = (λ3 - λ2)/(λ3 + λ2)

    return sqrt( Δ2^2 + Δ1^2 )

end

"""
    compute_current(skyrmion; label="uMOI", indices=[0,0], density = false, moment=0)

Compute a variety of currents in the Skyrme model, based on a `skyrmion`. 

You can calculate specific indices using e.g. `indices = [1,2]`. If `indices = [0,0]`, the function will calculate all indices. If `density = false`, the function will return the integrated densities, while `density = true` it will return the densities. 

The possible currents are (currently):
- `uMOI`: the isospin moment of inertia.
- `wMOI`: the mixed moment of inertia.
- `vMOI`: the spin moment of inertia.
- `uAxial`: the u-axial moment of inertia.
- `wAxial`: the w-axial moment of inertia.
- `NoetherIso`: the Noether vector current.
- `NoetherAxial`: the Noether axial current.
- `stress`: the stress tensor.

"""
function compute_current(sk; label="NULL", indices=[0,0], density=false, moment=0  )

    if label != "energy" && label != "uMOI" && label != "wMOI" && label != "vMOI" && label != "uAxial" && label != "wAxial" && label != "NoetherIso" && label != "NoetherAxial" && label != "stress"
        println("input label '", label, "' unknown. Type '?compute_current' for help.")
        return
    end

    current_density = zeros(3,3,sk.lp[1],sk.lp[2],sk.lp[3])

    aindices = 1:3
    bindices = 1:3


    if indices != [0,0]
        aindices = indices[1]:indices[1]
        bindices = indices[2]:indices[2]
    end
    

   for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]

            KD = diagm([1.0,1.0,1.0])

            xxx = SVector{3,Float64}(sk.x[1][i], sk.x[2][j], sk.x[3][k] )
            rm = sqrt( xxx[1]^2 + xxx[2]^2 + xxx[3]^2 )^moment

            p = getX(sk,i,j,k)
            dp = getDP(sk,i,j,k)

            if label == "uMOI"

                Tiam, Lia = getTiam(p), getLka(p,dp)

                for a in aindices, b in bindices
                    current_density[a,b,i,j,k] = -trace_su2_ij(Tiam,Tiam,a,b)*rm
                    for c in 1:3
                        current_density[a,b,i,j,k] -= 0.25*trace_su2_ijkl(Tiam,Lia,Tiam,Lia,a,c,b,c)*rm
                    end
                end

            elseif label == "wMOI"

                Tiam, Lia = getTiam(p), getLka(p,dp)
                Gia = -1.0.*getGia(Lia,xxx)

                for a in aindices, b in bindices
                    current_density[a,b,i,j,k] = -trace_su2_ij(Gia,Tiam,a,b)*rm
                    for c in 1:3
                        current_density[a,b,i,j,k] -= 0.25*trace_su2_ijkl(Gia,Lia,Tiam,Lia,a,c,b,c)*rm
                    end
                end

            elseif label == "vMOI"

                Lia = getLka(p,dp)
                Gia = -1.0.*getGia(Lia,xxx)

                for a in aindices, b in bindices
                    current_density[a,b,i,j,k] = -trace_su2_ij(Gia,Gia,a,b)*rm
                    for c in 1:3
                        current_density[a,b,i,j,k] -= 0.25*trace_su2_ijkl(Gia,Lia,Gia,Lia,a,c,b,c)*rm
                    end
                end

            elseif label == "uAxial"

                Tiam, Tiap, Lia = getTiam(p), getTiap(p), getRka(p,dp)

                for a in aindices, b in bindices
                    current_density[a,b,i,j,k] = -trace_su2_ij(Tiam,Tiap,a,b)*rm
                    for c in 1:3
                        current_density[a,b,i,j,k] -= 0.25*trace_su2_ijkl(Tiam,Lia,Tiap,Lia,a,c,b,c)*rm
                    end
                end

            elseif label == "wAxial"

                Tiap, Lia = getTiap(p), getLka(p,dp)
                Gia = -1.0.*getGia(Lia,xxx)

                for a in aindices, b in bindices
                    current_density[a,b,i,j,k] = trace_su2_ij(Tiap,Gia,a,b) *rm
                    for c in 1:3
                        current_density[a,b,i,j,k] += 0.25*trace_su2_ijkl(Tiap,Lia,Gia,Lia,a,c,b,c)*rm
                    end
                end
                

            elseif label == "stress"

                Lia = getLka(p,dp)

                for a in aindices, b in bindices

                    current_density[a,b,i,j,k] -= 2.0*(sk.mpi^2)*(1.0 - p[4])*KD[a,b]*rm
                    
                    current_density[a,b,i,j,k] -= trace_su2_ij(Lia,Lia,a,b)*rm

                    for c in 1:3

                        current_density[a,b,i,j,k] += 0.5*trace_su2_ij(Lia,Lia,c,c)*KD[a,b]*rm
                        
                        current_density[a,b,i,j,k] -= 1/4.0*trace_su2_ijkl(Lia,Lia,Lia,Lia,a,c,b,c)*rm

                        for d in 1:3
                            current_density[a,b,i,j,k] += 1/16.0*trace_su2_ijkl(Lia,Lia,Lia,Lia,c,d,c,d)*KD[a,b]*rm
                        end
                    end
                    
                end

            elseif label == "NoetherIso"

                Lia, Tiam = getLka(p,dp), getTiam(p)

                for a in aindices, b in bindices

                    current_density[a,b,i,j,k] -= trace_su2_ij(Lia,Tiam,b,a)*rm

                    for c in 1:3
                        current_density[a,b,i,j,k] -= 0.25*trace_su2_ijkl(Lia,Lia,Lia,Tiam,c,b,a,c)*rm
                    end

                end
            elseif label == "NoetherAxial"

                Lia, Tiap = getLka(p,dp), getTiap(p)

                for a in aindices, b in bindices

                    current_density[a,b,i,j,k] += trace_su2_ij(Lia,Tiap,b,a)*rm

                    for c in 1:3
                        current_density[a,b,i,j,k] += 0.25*trace_su2_ijkl(Lia,Lia,Lia,Tiap,c,b,a,c)*rm
                    end

                end
            
            end

        
        end

    end

    if ( label == "uMOI" || label == "wMOI" || label == "vMOI" ) && sk.physical == true
        println("Using physical units, result is in MeV fm^2.")
        current_density .*= 197^2/(sk.Fpi*sk.ee^3)
    end

    if density == false
        tot_den = zeros(3,3)
        for a in aindices, b in bindices
            tot_den[a,b] = sum(current_density[a,b,:,:,:])*sk.ls[1]*sk.ls[2]*sk.ls[3]
        end
        if indices == [0,0]
            return tot_den
        else
            return tot_den[aindices,bindices]
        end
    else
        if indices == [0,0]
            return current_density
        else
            return current_density[aindices,bindices,:,:,:]
        end
    end
    

end


function getGia(Lia,x)

    return SMatrix{3,3,Float64, 9}(
        Lia[3,1]*x[2] - Lia[2,1]*x[3],
        -Lia[3,1]*x[1] + Lia[1,1]*x[3],
        Lia[2,1]*x[1] - Lia[1,1]*x[2],
        Lia[3,2]*x[2] - Lia[2,2]*x[3],
        -Lia[3,2]*x[1] + Lia[1,2]*x[3],
        Lia[2,2]*x[1] - Lia[1,2]*x[2],
        Lia[3,3]*x[2] - Lia[2,3]*x[3],
        -Lia[3,3]*x[1] + Lia[1,3]*x[3],
        Lia[2,3]*x[1] - Lia[1,3]*x[2]
    )

end

function getTiam(p)

    return SMatrix{3,3,Float64, 9}(
        -p[2]^2 - p[3]^2,
        p[1]*p[2] - p[3]*p[4],
        p[1]*p[3] + p[2]*p[4],
        p[1]*p[2] + p[3]*p[4],
        -p[1]^2 - p[3]^2,
        p[2]*p[3] - p[1]*p[4],
        p[1]*p[3] - p[2]*p[4],
        p[2]*p[3] + p[1]*p[4],
        -p[1]^2 - p[2]^2
    ) 
    
end

function getTiap(p)
    
    return SMatrix{3,3,Float64, 9}(
        p[1]^2 + p[4]^2,
        p[1]*p[2] - p[3]*p[4],
        p[1]*p[3] + p[2]*p[4],
        p[1]*p[2] + p[3]*p[4],
        p[2]^2 + p[4]^2,
        p[2]*p[3] - p[1]*p[4],
        p[1]*p[3] - p[2]*p[4],
        p[2]*p[3] + p[1]*p[4],
        p[3]^2 + p[4]^2
    ) 
    
end

function getRka(p,dp)
    
    return SMatrix{3,3,Float64, 9}(

        -(dp[1,4]*p[1]) - dp[1,3]*p[2] + dp[1,2]*p[3] + dp[1,1]*p[4],
        -(dp[2,4]*p[1]) - dp[2,3]*p[2] + dp[2,2]*p[3] + dp[2,1]*p[4],
        -(dp[3,4]*p[1]) - dp[3,3]*p[2] + dp[3,2]*p[3] + dp[3,1]*p[4],
        dp[1,3]*p[1] - dp[1,4]*p[2] - dp[1,1]*p[3] + dp[1,2]*p[4],
        dp[2,3]*p[1] - dp[2,4]*p[2] - dp[2,1]*p[3] + dp[2,2]*p[4],
        dp[3,3]*p[1] - dp[3,4]*p[2] - dp[3,1]*p[3] + dp[3,2]*p[4],
        -(dp[1,2]*p[1]) + dp[1,1]*p[2] - dp[1,4]*p[3] + dp[1,3]*p[4],
        -(dp[2,2]*p[1]) + dp[2,1]*p[2] - dp[2,4]*p[3] + dp[2,3]*p[4],
        -(dp[3,2]*p[1]) + dp[3,1]*p[2] - dp[3,4]*p[3] + dp[3,3]*p[4],
    ) 
    
end

function getLka(p,dp)

    return SMatrix{3,3,Float64, 9}(

    -(dp[1,4]*p[1]) + dp[1,3]*p[2] - dp[1,2]*p[3] + dp[1,1]*p[4],
    -(dp[2,4]*p[1]) + dp[2,3]*p[2] - dp[2,2]*p[3] + dp[2,1]*p[4],
    -(dp[3,4]*p[1]) + dp[3,3]*p[2] - dp[3,2]*p[3] + dp[3,1]*p[4],
    -(dp[1,3]*p[1]) - dp[1,4]*p[2] + dp[1,1]*p[3] + dp[1,2]*p[4],
    -(dp[2,3]*p[1]) - dp[2,4]*p[2] + dp[2,1]*p[3] + dp[2,2]*p[4],
    -(dp[3,3]*p[1]) - dp[3,4]*p[2] + dp[3,1]*p[3] + dp[3,2]*p[4],
    dp[1,2]*p[1] - dp[1,1]*p[2] - dp[1,4]*p[3] + dp[1,3]*p[4],
    dp[2,2]*p[1] - dp[2,1]*p[2] - dp[2,4]*p[3] + dp[2,3]*p[4],
    dp[3,2]*p[1] - dp[3,1]*p[2] - dp[3,4]*p[3] + dp[3,3]*p[4]

    )

end

function trace_su2_ij(Lia,Lib,i,j)
    return -2.0*(Lia[1,i]*Lib[1,j] + Lia[2,i]*Lib[2,j] + Lia[3,i]*Lib[3,j])
end

function trace_su2_ijkl(L1,L2,L3,L4,i,j,k,l)
    return -8.0*(L1[i,1]*(L2[j,2]*(-(L3[k,2]*L4[l,1]) + L3[k,1]*L4[l,2]) + L2[j,3]*(-(L3[k,3]*L4[l,1]) + L3[k,1]*L4[l,3])) + L1[i,3]*(L2[j,1]*(L3[k,3]*L4[l,1] - L3[k,1]*L4[l,3]) + L2[j,2]*(L3[k,3]*L4[l,2] - L3[k,2]*L4[l,3])) + L1[i,2]*(L2[j,1]*(L3[k,2]*L4[l,1] - L3[k,1]*L4[l,2]) + L2[j,3]*(-(L3[k,3]*L4[l,2]) + L3[k,2]*L4[l,3])))
end

function Berger_Isospin(sk)

    ID = zeros(sk.lp[1], sk.lp[2], sk.lp[3])

    get_BT_density!(ID,sk)

    total = sum(ID)*sk.ls[1]*sk.ls[2]*sk.ls[3]

    return total

end    

function get_BT_density!(density, sk)

    for k in sk.sum_grid[3]
        @inbounds for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
            dp = getDP(sk ,i, j, k )

            density[i,j,k] = BT_engpt(sk,dp,i,j,k,sk.metric) 

        end
    end

end

function left_t(p,dp)
    
    p0 = p[4]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]
    phi = [p1, p2, p3]

    c = zeros(3,3)
    dp_s = dp[:, 1:3]  
    dp_t = dp[:, 4]    

    for a in 1:3
        v = dp_s[a,:]
        v0 = dp_t[a]
        c[a,:] = (p0 * v) - (v0 * phi) + [(phi[2] * v[3] - phi[3] * v[2]),
        (phi[3] * v[1] - phi[1] * v[3]),
        (phi[1] * v[2] - phi[2] * v[1])]
    end

    return c
end

function get_left_currents(p,dp)

    lc = left_t(p,dp)

    L_1 = lc[1,:]
    L_2 = lc[2,:]
    L_3 = lc[3,:]

    return (L_1,L_2,L_3)
end

function b_metric_su2(alpha,u,v)
   
    met_su2 = u[1]*v[1]+u[2]*v[2]+(alpha^2)*(u[3]*v[3])

    return met_su2
end

function lie_bracket(x, y)
    return [-2 * (x[2] * y[3] - x[3] * y[2]),
            -2 * (x[3] * y[1] - x[1] * y[3]),
            -2 * (x[1] * y[2] - x[2] * y[1])]
end

function T2_energy(alpha,L_0)

    return  b_metric_su2(alpha,L_0,L_0)

end

function T4_energy(alpha,L0,L1,L2,L3)
    t01 = lie_bracket(L0,L1)
    t02 = lie_bracket(L0,L2)
    t03 = lie_bracket(L0,L3)

    v = t01 + t02 + t03
    
    return (1/4) * b_metric_su2(alpha,v,v)
end    

function BT_engpt(sk,dp,i,j,k,alpha)

    p = sk.pion_field[i,j,k,:]

    LC = get_left_currents(p,dp)

    L_0 = (2 * (p[4]*p[2] + p[3]*p[1]), 2 * (-p[4]*p[1] + p[3]*p[2]), 2 * (-(p[1])^2 - (p[2])^2))
    L_1 = LC[1]
    L_2 = LC[2]
    L_3 = LC[3]

    T2 = T2_energy(alpha,L_0)
    T4 = T4_energy(alpha,L_0,L_1,L_2,L_3)
    
    return 2*(T2+T4)

end

