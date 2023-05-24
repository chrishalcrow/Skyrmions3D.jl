
function makeRM!(skyrmion, prof, pfn, qfn)
    
    lp, x = skyrmion.lp, skyrmion.x

    r = zeros(lp[1],lp[2],lp[3])

    zRM = complex(0.0,0.0)  
    pRM = complex(0.0,0.0)
    qRM = complex(0.0,0.0) 
    den = complex(0.0,0.0)

    sine_of_prof_r = 0.0

    @inbounds for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        r = sqrt( x[1][i]^2 + x[2][j]^2 + x[3][k]^2 )
        sine_of_prof_r = sin(prof(r))

        zRM = (x[1][i] + 1.0im*x[2][j])/(r + x[3][k])

        pRM = pfn(zRM)
        qRM = qfn(zRM)
        den = real( qRM*conj(qRM) + pRM*conj(pRM) )

        skyrmion.phi[i,j,k,1] = (sine_of_prof_r/den)*real( pRM*conj(qRM) + qRM*conj(pRM) )
        skyrmion.phi[i,j,k,2] = (sine_of_prof_r/den)*imag( pRM*conj(qRM) - qRM*conj(pRM) )
        skyrmion.phi[i,j,k,3] = (sine_of_prof_r/den)*real( qRM*conj(qRM) - pRM*conj(pRM) )
        skyrmion.phi[i,j,k,4] = cos(prof(r))

    end
    
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






# ADHM STUFF

function q1(l,p,a)
    @inbounds l[a,1]*p[a,1] - l[a,2]*p[a,2] - l[a,3]*p[a,3] - l[a,4]*p[a,4]
end
function q2(l,p,a)
    @inbounds l[a,2]*p[a,1] + l[a,1]*p[a,2] - l[a,4]*p[a,3] + l[a,3]*p[a,4]
end
function q3(l,p,a)
    @inbounds l[a,3]*p[a,1] + l[a,4]*p[a,2] + l[a,1]*p[a,3] - l[a,2]*p[a,4]
end
function q4(l,p,a)
    @inbounds l[a,4]*p[a,1] - l[a,3]*p[a,2] + l[a,2]*p[a,3] + l[a,1]*p[a,4]
end


function q1(l,p,a,b)
    @inbounds l[a,b,1]*p[b,1] - l[a,b,2]*p[b,2] - l[a,b,3]*p[b,3] - l[a,b,4]*p[b,4]
end
function q2(l,p,a,b)
    @inbounds l[a,b,2]*p[b,1] + l[a,b,1]*p[b,2] - l[a,b,4]*p[b,3] + l[a,b,3]*p[b,4]
end
function q3(l,p,a,b)
    @inbounds l[a,b,3]*p[b,1] + l[a,b,4]*p[b,2] + l[a,b,1]*p[b,3] - l[a,b,2]*p[b,4]
end
function q4(l,p,a,b)
    @inbounds l[a,b,4]*p[b,1] - l[a,b,3]*p[b,2] + l[a,b,2]*p[b,3] + l[a,b,1]*p[b,4]
end


 function conj!(Np,B)
    for a in 1:B, c in 2:4
        Np[a,c] *= -1.0
    end
end

function conj!(Np,tint,B)
    for a in 1:B, c in 2:4
        Np[tint,a,c] *= -1.0
    end
end

function makeQmultmatrix(l)

    return [ l[1] -l[2] -l[3] -l[4]
             l[2]  l[1]  l[4] -l[3]  
             l[3]  l[4]  l[1]  l[2] 
             l[4]  l[3] -l[2]  l[1] ]

end

function makeQmultmatrix2(q)


    return [ q[1] -q[2] -q[3] -q[4] 
             q[2]  q[1] -q[4] q[3]
             q[3]  q[4]  q[1]  -q[2]
             q[4] -q[3]  q[2] q[1] ]

end



function ADHMpt2(L,M,y,B)

    U = [1.0 ; 0.0; 0; 0]




    Mmdysp = zeros(B,B,4)
    p = zeros(B,4)


    Ω1M = zeros(4,4)
    Ω2M = zeros(4,4)
    Ω3M = zeros(4,4)

    Ω1 = zeros(4)
    Ω1i = zeros(4)
    Ω2 = zeros(4)
    Ω3 = zeros(4)


    tsteps = 42;
    lstime = pi/tsteps;

    allN = zeros(tsteps+1,B+1,4)

   

    

    

    for tint in 1:tsteps+1

        #Nfy!(allN,tint,L,M,[t,y[1],y[2],y[3]], Mmdysp ,p, B,1.0)
        #t += δt

        ct = cos(lstime*(tint-1));
        x0 = sin(lstime*(tint-1));

        x1t = y[1]*ct; x2t = y[2]*ct; x3t = y[3]*ct;
        Nfy!(allN,tint,L.*ct,M.*ct,[x0,x1t,x2t,x3t], Mmdysp ,p, B,ct)

        
        
        #if tint == 2; println(t, allN[tint,:,:]); end;
        
        #println(t)

        #Nfy!(allN,tint,L.*ct,M.*ct,[x0,x1t,x2t,x3t], Mmdysp ,p, B,ct)

        
    end

    for tint in 1:2:tsteps-1


        
        #= FIRST ORDER 
        getΩ!(Ω1,allN[tint+1,:,:],allN[tint,:,:],B)
        Ω1M = makeQmultmatrix2(Ω1)
        U = Ω1M*U =#

        #= SECOND ORDER =#
        getΩ!(Ω1,allN[tint+2,:,:],allN[tint+1,:,:],B)
        getΩ!(Ω2,allN[tint+1,:,:],allN[tint,:,:],B)
        getΩ!(Ω3,allN[tint+2,:,:],allN[tint,:,:],B)
        
        Ω1M = makeQmultmatrix2(Ω1)
        Ω2M = makeQmultmatrix2(Ω2)
        Ω3M = makeQmultmatrix2(Ω3)

        U = (4.0.*Ω1M*Ω2M - Ω3M)*U./3.0
        

    end

    normer = sqrt( U[1]^2 + U[2]^2 + U[3]^2 + U[4]^2 )
    
    return [U[2],U[3],U[4],U[1]]./normer

end

function getΩ!(Ω1,vp,v,B)

    Ω1[1] = 0.0; 
    Ω1[2] = 0.0; 
    Ω1[3] = 0.0; 
    Ω1[4] = 0.0;

    for a in 1:B+1
        Ω1[1] += v[a,1]*vp[a,1] + v[a,2]*vp[a,2] + v[a,3]*vp[a,3] + v[a,4]*vp[a,4]
        Ω1[2] += v[a,2]*vp[a,1] - v[a,1]*vp[a,2] - v[a,4]*vp[a,3] + v[a,3]*vp[a,4]
        Ω1[3] += v[a,3]*vp[a,1] + v[a,4]*vp[a,2] - v[a,1]*vp[a,3] - v[a,2]*vp[a,4]
        Ω1[4] += v[a,4]*vp[a,1] - v[a,3]*vp[a,2] + v[a,2]*vp[a,3] - v[a,1]*vp[a,4]
    end

end




function Nfy!(Nα,L,M,y,Mmdysp,p,B)

    Rnm = zeros(B,B)#zeros( MMatrix{2,2,Float64} )
    
    for a in 1:B, b in 1:B, c in 1:4
        Mmdysp[a,b,c] = 0.0
    end

    for a in 1:B, c in 1:4
        p[a,c] = 0.0
    end
    
    @inbounds for b in 1:B, a in 1:4
        Mmdysp[b,b,a] -= y[a]
        for c in 1:B
            Mmdysp[b,c,a] += M[b,c,a]
        end
    end

    makeRnm!(Rnm, L, Mmdysp, B)
    
    
    iRnm = inv(Rnm)


    @inbounds for a in 1:B

        for b in 1:B
            
            p[a,1] += iRnm[a,b]*L[b,1]

            for c in 2:4
                p[a,c] -= iRnm[a,b]*L[b,c]
            end

        end
        
        for c in 2:4
            Nα[a,c] = 0.0
        end

    end

    Nα[1,1] = 1.0



    for c in 1:4
        Nα[B+1,c] = 0.0
    end

    normer = 0.0

    @inbounds for a in 1:B

        Nα[1,1] -= q1(L,p,a)
        Nα[1,2] -= q2(L,p,a)
        Nα[1,3] -= q3(L,p,a)
        Nα[1,4] -= q4(L,p,a)

        for b in 1:B

            Nα[b+1,1] -= q1(Mmdysp,p,b,a)
            Nα[b+1,2] -= q2(Mmdysp,p,b,a)
            Nα[b+1,3] -= q3(Mmdysp,p,b,a)
            Nα[b+1,4] -= q4(Mmdysp,p,b,a)

        end

    end

    @inbounds for a in 1:B+1, b in 1:4 
        normer += Nα[a,b]^2
    end
    normer = sqrt(normer)

    @inbounds for a in 1:B+1, b in 1:4 
        Nα[a,b] /= normer
    end
    

end






function Nfy!(Nα,tint,L,M,y,Mmdysp,p,B,ct)

    Rnm = zeros(B,B)#zeros( MMatrix{2,2,Float64} )
    
    for a in 1:B, b in 1:B, c in 1:4
        Mmdysp[a,b,c] = 0.0
    end

    for a in 1:B, c in 1:4
        p[a,c] = 0.0
    end
    
    @inbounds for b in 1:B, a in 1:4
        Mmdysp[b,b,a] -= y[a]
        for c in 1:B
            Mmdysp[b,c,a] += M[b,c,a]
        end
    end

    makeRnm!(Rnm, L, Mmdysp, B)
    
    
    iRnm = inv(Rnm)


    @inbounds for a in 1:B

        for b in 1:B
            
            p[a,1] += iRnm[a,b]*L[b,1]

            for c in 2:4
                p[a,c] -= iRnm[a,b]*L[b,c]
            end

        end
        
        for c in 2:4
            Nα[tint,a,c] = 0.0
        end

    end

    Nα[tint,1,1] = 1.0



    for c in 1:4
        Nα[tint,B+1,c] = 0.0
    end

    

    @inbounds for a in 1:B

        Nα[tint,1,1] -= q1(L,p,a)
        Nα[tint,1,2] -= q2(L,p,a)
        Nα[tint,1,3] -= q3(L,p,a)
        Nα[tint,1,4] -= q4(L,p,a)

        for b in 1:B

            Nα[tint,b+1,1] -= q1(Mmdysp,p,b,a)
            Nα[tint,b+1,2] -= q2(Mmdysp,p,b,a)
            Nα[tint,b+1,3] -= q3(Mmdysp,p,b,a)
            Nα[tint,b+1,4] -= q4(Mmdysp,p,b,a)

        end

    end

    normer = 0.0

    @inbounds for a in 1:B+1, b in 1:4 
        normer += Nα[tint,a,b]^2
    end
    normer = sqrt(normer)

    @inbounds for a in 1:B+1, b in 1:4 
        Nα[tint,a,b] /= normer
    end
    

end




function makeRnm!(Rnm,L,M,B)
    @inbounds for a in 1:B, b in 1:B
        @fastmath Rnm[a,b] = L[a,1]*L[b,1] + L[a,2]*L[b,2] + L[a,3]*L[b,3] + L[a,4]*L[b,4]
        for c in 1:B, d in 1:4
            @fastmath Rnm[a,b] += (M[a,c,d]*M[c,b,d])
        end
    end
    
    #return Rnm
        
end
























