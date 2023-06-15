"""
    Energy(skyrmion; density=false)

Compute energy of `skyrmion`.

Set 'density = true' to output the energy density and moment to `n` to calculate the nth moment of the energy density.

"""
function Energy(sk,pion_mass; density=false, moment=0)

    sk.mpi = pion_mass

    ED = zeros(sk.lp[1], sk.lp[2], sk.lp[3])

    engtot = 0.0
    dp = zeros(3,4)

    @simd for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
    
            getDX!(dp, sk ,i, j, k )

            if moment == 0
                @inbounds ED[i,j,k] = engpt(dp,sk.phi[i,j,k,4],sk.mpi)
            else
                r = sqrt( sk.x[1][i]^2 + sk.x[2][j]^2 + sk.x[3][k]^2 )
                @inbounds ED[i,j,k] = engpt(dp,sk.phi[i,j,k,4],sk.mpi) * r^moment
            end

        end
        
    end

    if density == false
        engtot = sum(ED)*sk.ls[1]*sk.ls[2]*sk.ls[3]
        if sk.physical == false
            return engtot/(12.0*pi^2)
        else
            return (engtot*sk.Fpi/(4.0*sk.ee), "MeV" )
        end
    else
        return ED
    end

end 


function engpt(dp,p4,mpi)
    return 2.0*mpi*(1.0 - p4) + (dp[1,1]^2 + dp[1,2]^2 + dp[1,3]^2 + dp[1,4]^2 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2 + dp[3,4]^2) + (dp[1,4]^2*dp[2,1]^2 + dp[1,4]^2*dp[2,2]^2 + dp[1,4]^2*dp[2,3]^2 + dp[1,1]^2*(dp[2,2]^2 + dp[2,3]^2) - 2*dp[1,1]*dp[1,4]*dp[2,1]*dp[2,4] + dp[1,1]^2*dp[2,4]^2 + dp[1,4]^2*dp[3,1]^2 + dp[2,2]^2*dp[3,1]^2 + dp[2,3]^2*dp[3,1]^2 + dp[2,4]^2*dp[3,1]^2 - 2*dp[2,1]*dp[2,2]*dp[3,1]*dp[3,2] + dp[1,1]^2*dp[3,2]^2 + dp[1,4]^2*dp[3,2]^2 + dp[2,1]^2*dp[3,2]^2 + dp[2,3]^2*dp[3,2]^2 + dp[2,4]^2*dp[3,2]^2 - 2*dp[2,1]*dp[2,3]*dp[3,1]*dp[3,3] - 2*dp[2,2]*dp[2,3]*dp[3,2]*dp[3,3] + dp[1,1]^2*dp[3,3]^2 + dp[1,4]^2*dp[3,3]^2 + dp[2,1]^2*dp[3,3]^2 + dp[2,2]^2*dp[3,3]^2 + dp[2,4]^2*dp[3,3]^2 - 2*(dp[1,1]*dp[1,4]*dp[3,1] + dp[2,4]*(dp[2,1]*dp[3,1] + dp[2,2]*dp[3,2] + dp[2,3]*dp[3,3]))*dp[3,4] + (dp[1,1]^2 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2)*dp[3,4]^2 + dp[1,3]^2*(dp[2,1]^2 + dp[2,2]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,4]^2) + dp[1,2]^2*(dp[2,1]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,3]^2 + dp[3,4]^2) - 2*dp[1,2]*(dp[1,1]*(dp[2,1]*dp[2,2] + dp[3,1]*dp[3,2]) + dp[1,3]*(dp[2,2]*dp[2,3] + dp[3,2]*dp[3,3]) + dp[1,4]*(dp[2,2]*dp[2,4] + dp[3,2]*dp[3,4])) - 2*dp[1,3]*(dp[1,1]*(dp[2,1]*dp[2,3] + dp[3,1]*dp[3,3]) + dp[1,4]*(dp[2,3]*dp[2,4] + dp[3,3]*dp[3,4])))
end

"""
    Baryon(skyrmion; density=false)

Compute baryon number of `skyrmion`.

Set 'density = true' to output the baryon density and moment to `n` to calculate the nth moment of the baryon density.

"""
function Baryon(sk; density=false, moment = 0)

    BD = zeros(sk.lp[1],sk.lp[2],sk.lp[3])

    bartot = 0.0
    dp = zeros(3,4)
    pp = zeros(4)
    
    for i in 3:sk.lp[1]-2, j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
        getDX!(dp, sk ,i, j, k )
        getX!(pp,sk,i,j,k)

        if moment == 0
            @inbounds BD[i,j,k] = barypt(dp,pp)
        else
            r = sqrt( sk.x[1][i]^2 + sk.x[2][j]^2 + sk.x[3][k]^2 )
            @inbounds BD[i,j,k] = barypt(dp,pp) * r^moment
        end
            
    end
    
    if density == false
        bartot = sum(BD)*sk.ls[1]*sk.ls[2]*sk.ls[3]/(2.0*pi^2)
        return bartot
    else
        return BD
    end
    
end 

function barypt(dp,pp)
    return pp[4]*dp[1,3]*dp[2,2]*dp[3,1] - pp[3]*dp[1,4]*dp[2,2]*dp[3,1] - pp[4]*dp[1,2]*dp[2,3]*dp[3,1] + pp[2]*dp[1,4]*dp[2,3]*dp[3,1] + pp[3]*dp[1,2]*dp[2,4]*dp[3,1] - pp[2]*dp[1,3]*dp[2,4]*dp[3,1] - pp[4]*dp[1,3]*dp[2,1]*dp[3,2] + pp[3]*dp[1,4]*dp[2,1]*dp[3,2] + pp[4]*dp[1,1]*dp[2,3]*dp[3,2] - pp[1]*dp[1,4]*dp[2,3]*dp[3,2] - pp[3]*dp[1,1]*dp[2,4]*dp[3,2] + pp[1]*dp[1,3]*dp[2,4]*dp[3,2] + pp[4]*dp[1,2]*dp[2,1]*dp[3,3] - pp[2]*dp[1,4]*dp[2,1]*dp[3,3] - pp[4]*dp[1,1]*dp[2,2]*dp[3,3] + pp[1]*dp[1,4]*dp[2,2]*dp[3,3] + pp[2]*dp[1,1]*dp[2,4]*dp[3,3] - pp[1]*dp[1,2]*dp[2,4]*dp[3,3] - pp[3]*dp[1,2]*dp[2,1]*dp[3,4] + pp[2]*dp[1,3]*dp[2,1]*dp[3,4] + pp[3]*dp[1,1]*dp[2,2]*dp[3,4] - pp[1]*dp[1,3]*dp[2,2]*dp[3,4] - pp[2]*dp[1,1]*dp[2,3]*dp[3,4] + pp[1]*dp[1,2]*dp[2,3]*dp[3,4]
end


"""
    getMOI(skyrmion; density=false, moment=0)

Compute the moments of inertia of `skyrmion`.

Set 'density = true' to output the MOI densities and moment to `n` to calculate the nth moment of the MOI.

"""
function getMOI(sk; density = false, moment=0)
	
    x = sk.x

    MOI_D = zeros(6,6,sk.lp[1],sk.lp[2],sk.lp[3])

    epsilon = make_levi_civita()

	VV = zeros(6,6)
	
    mm = zeros(3,4)
    po = zeros(4)

    GG = zeros(6,3)
    RR = zeros(3,3)
	
	xxx = zeros(3)

	for i = 3:sk.lp[1]-2,  j = 3:sk.lp[2]-2, k = 3:sk.lp[3]-2

	    xxx = [x[1][i], x[2][j], x[3][k]]
		
        getDX!(mm, sk ,i, j, k)
		getX!(po, sk, i, j, k)

	    for a in 1:3, b in 1:3

	        RR[a,b] = po[4]*mm[a,b] - mm[a,4]*po[b]
	        GG[a+3,b] = -po[a]*po[b]
	        GG[a+3,a] += po[b]*po[b]

	        for c in 1:3
	            GG[a+3,b] -= po[4]*po[c]*epsilon[a,c,b]
	            for d in 1:3
	                RR[a,b] += epsilon[b,c,d]*mm[a,c]*po[d]
	            end
	        end
	    end

	    for a in 1:3, d in 1:3
            GG[a,d] = 0.0
            for b in 1:3, c in 1:3
	    	    GG[a,d] += epsilon[a,b,c]*xxx[b]*RR[c,d]
            end
	    end

        for a in 1:6, b in 1:6, c in 1:3

            if moment == 0
                MOI_D[a,b,i,j,k] += GG[a,c]*GG[b,c]
            else
                r = sqrt( sk.x[1][i]^2 + sk.x[2][j]^2 + sk.x[3][k]^2 )
                MOI_D[a,b,i,j,k] += GG[a,c]*GG[b,c] * r^moment
            end

            for d = 1:3, e = 1:3

                if moment == 0
                    MOI_D[a,b,i,j,k] -= 2.0*(GG[a,c]*GG[b,d]*RR[e,c]*RR[e,d] - GG[a,c]*GG[b,c]*RR[d,e]*RR[d,e])
                else
                    r = sqrt( sk.x[1][i]^2 + sk.x[2][j]^2 + sk.x[3][k]^2 )
                    MOI_D[a,b,i,j,k] -= 2.0*(GG[a,c]*GG[b,d]*RR[e,c]*RR[e,d] - GG[a,c]*GG[b,c]*RR[d,e]*RR[d,e])* r^moment
                end

            end
		end
		
	end

    if density == false
        for a in 1:6, b in 1:6
            for i = 3:sk.lp[1]-2,  j = 3:sk.lp[2]-2, k = 3:sk.lp[3]-2
                VV[a,b] += MOI_D[a,b,i,j,k]*sk.ls[1]*sk.ls[2]*sk.ls[3]
            end
        end
        if sk.physical == false
            return VV
        else
            return (VV*sk.Fpi/(4.0*sk.ee)*(2.0/(sk.ee*sk.Fpi))^2, "MeV fm^2")
        end
    else
        return MOI_D
    end
	
	return VV
	
end

function make_levi_civita()

    epsilon = zeros(3,3,3)
	
	epsilon[1,2,3] = 1.0
	epsilon[1,3,2] = -1.0
	epsilon[2,3,1] = 1.0
	epsilon[2,1,3] = -1.0
	epsilon[3,1,2] = 1.0
	epsilon[3,2,1] = -1.0

    return epsilon

end

"""
    center_of_mass(skyrmion; density=false, moment=0)

Compute the center of mass of `skyrmion`, based on the energy density.


"""
function center_of_mass(sk)

    com = zeros(3)

    the_engpt = 0.0
    toteng = 0.0

    dp = zeros(3,4)

    @simd for i in 3:sk.lp[1]-2
        for j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
    
            getDX!(dp, sk ,i, j, k )

            the_engpt = engpt(dp,sk.phi[i,j,k,4],sk.mpi)

            com[1] += the_engpt*sk.x[1][i]
            com[2] += the_engpt*sk.x[2][j]
            com[3] += the_engpt*sk.x[3][k]

            toteng += the_engpt
   
        end
        
    end

    return com/toteng

end 


"""
    rms_baryon(skyrmion)

Compute root mean square charge radius of a skyrmion, using the baryon density.
"""
function rms_baryon(sk)

    if sk.physical == false
        return sqrt(Baryon(sk; moment = 2))
    else
        return ( sqrt(Baryon(sk; moment = 2))*(2.0/(sk.ee*sk.Fpi)), "fm" )
    end

end


