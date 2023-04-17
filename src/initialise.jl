function makeRM(prof, pfn, qfn,x,lp)
    
    phi = zeros(4,lp[1],lp[2],lp[3])

    r = zeros(lp[1],lp[2],lp[3])

    zRM = complex(zeros(lp[1],lp[2],lp[3]));
    pRM = complex(zeros(lp[1],lp[2],lp[3]));
    qRM = complex(zeros(lp[1],lp[2],lp[3]));
    den = zeros(lp[1],lp[2],lp[3]);

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
        r[i,j,k] = sqrt( x[1][i]^2 + x[2][j]^2 + x[3][k]^2 )
        zRM[i,j,k] = (x[1][i]/r[i,j,k])/(1.0 + x[3][k]/r[i,j,k]) + 1.0im*(x[2][j]/r[i,j,k])/(1.0 + x[3][k]/r[i,j,k])
    end

    pRM = pfn.(zRM)
    qRM = qfn.(zRM)
    den = real.( qRM.*conj.(qRM) .+ pRM.*conj.(pRM) )

    phi[1,:,:,:] = ((sin.(prof.(r)))./den).*real( pRM.*conj.(qRM) .+ qRM.*conj.(pRM) )
    phi[2,:,:,:] = ((sin.(prof.(r)))./den).*imag( pRM.*conj.(qRM) .- qRM.*conj.(pRM) )
    phi[3,:,:,:] = ((sin.(prof.(r)))./den).*real( qRM.*conj.(qRM) .- pRM.*conj.(pRM) )
    phi[4,:,:,:] = cos.(prof.(r))
    
    return phi

end

function multicubes!(phi3,phi1,phi2,phi,Xs)
	
	
	for a in 1:8

	    SkyrShift!(phi1,phi,Xs[a][1],Xs[a][2],Xs[a][3])
	    #SkyrIso!(ϕ1,αs[a])
	    MakeProduct!(phi3,phi1,phi2)

	    phi2.phi .= phi3.phi
    
	end

end
	



function makeVac(lp)
    
    phi = zeros(4,lp[1],lp[2],lp[3])
    phi[4,:,:,:] = ones(lp[1],lp[2],lp[3])
    
    return phi

end

function SkyrIso!(ϕ,α)

    newphi = makeVac()

    for i in 1:ϕ.lp[1], j in 1:ϕ.lp[2], k in 1:ϕ.lp[3]
        newphi[1,i,j,k] = cos(α)*ϕ.phi[1,i,j,k] + sin(α)*ϕ.phi[2,i,j,k]
        newphi[2,i,j,k] = -sin(α)*ϕ.phi[1,i,j,k] + cos(α)*ϕ.phi[2,i,j,k]
        newphi[3,i,j,k] = ϕ.phi[3,i,j,k]
        newphi[4,i,j,k] = ϕ.phi[4,i,j,k]
    end
    
    ϕ.phi .= newphi
    
end



function MakeProduct!(ϕ3,ϕ1,ϕ2)
    
    for i in 1:ϕ3.lp[1], j in 1:ϕ3.lp[2], k in 1:ϕ3.lp[3]
        
        ϕ3.phi[4,i,j,k] = ϕ1.phi[4,i,j,k]*ϕ2.phi[4,i,j,k]
        
        for a in 1:3
            ϕ3.phi[a,i,j,k] = ϕ1.phi[4,i,j,k]*ϕ2.phi[a,i,j,k] + ϕ2.phi[4,i,j,k]*ϕ1.phi[a,i,j,k]
            ϕ3.phi[4,i,j,k] -= ϕ1.phi[a,i,j,k]*ϕ2.phi[a,i,j,k]
        end
    end
    
    normer(ϕ3)
    
end
        




function SkyrShift!(ϕnew,ϕ, xs, ys, zs)

    ϕnew.phi = zeros(4,ϕ.lp[1],ϕ.lp[2],ϕ.lp[3])
    ϕnew.phi[4,:,:,:] = ones(ϕ.lp[1],ϕ.lp[2],ϕ.lp[3])

    xdn = 1
    xup = ϕnew.lp[1]

    ydn = 1
    yup = ϕnew.lp[2]

    zdn = 1
    zup = ϕnew.lp[3]

    if xs > 0
        xdn = 1 + xs
    elseif xs < 0
        xup = ϕnew.lp[1] + xs
    end

    if ys > 0
        ydn = 1 + ys
    elseif ys < 0
        yup = ϕnew.lp[2] + ys
    end

    if zs > 0
        zdn = 1 + zs
    elseif zs < 0
        zup = ϕnew.lp[3] + zs
    end

    for i in xdn:xup, j in ydn:yup, k in zdn:zup, a in 1:4
        ϕnew.phi[a,i,j,k] = ϕ.phi[a,i-xs,j-ys,k-zs]
    end
    
    

end