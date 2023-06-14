function R_from_axis_angle(th, n)

    if th == 0.0
        return [ 1.0 0 0 ; 0 1.0 0 ; 0 0 1.0 ]
    end

    if n == [0.0,0.0,0.0]
        println("ERROR: your vector is zero.")
    end

    n1 = n[1]
    n2 = n[2]
    n3 = n[3]

    normer = sqrt( n1^2 + n2^2 + n3^2 )

    n1 /= normer
    n2 /= normer
    n3 /= normer

     return [ n1^2 + (n2^2 + n3^2)*cos(th) 2*sin(th/2.0)*(-(n3*cos(th/2.0 + n1*n3*sin(th/2.0)))) 2*sin(th/2.0)*(n2*cos(th/2.0) + n1*n3*sin(th/2.0)) ;  2*sin(th/2.0)*(n3*cos(th/2.0) + n1*n2*sin(th/2.0)) n2^2 + (n1^2 + n3^2)*cos(th) n2*n3 - n2*n3*cos(th) - n1*sin(th) ; 2*n1*n3*sin(th/2.0)^2 - n2*sin(th) 2*sin(th/2.0)*(n1*cos(th/2.0) + n2*n3*sin(th/2.0)) n3^2 + (n1^2 + n2^2)*cos(th) ]
        

end


"""
    makeRM!(skyrmion, prof, pfn, qfn; kwargs... )
    
Writes a rational map skyrmion in to `skyrmion`. The rational map is given by the polynomials R(z) = p(z)/q(z) and the profile f(r).

# Optional arguments
-  `X=[0.0,0.0,0.0]`: translate the initial skyrmion by `X`
-  `iTH = 0.0`: isorotate by initial skyrmion by `iTH`
-  `i_n = 0.0`: isorotate initial skyrmion around `i_n`
-  `jTH = 0.0`: isorotate by initial skyrmion by `jTH`
-  `j_n = 0.0`: isorotate initial skyrmion around `j_n`

"""
function makeRM!(skyrmion, prof, pfn, qfn; X=[0.0,0.0,0.0], iTH=0.0, i_n = [0.0,0.0,1.0], jTH = 0.0, j_n = [0.0,0.0,0.0] )
    
    lp, x = skyrmion.lp, skyrmion.x

    r = zeros(lp[1],lp[2],lp[3])

    zRM = complex(0.0,0.0)  
    pRM = complex(0.0,0.0)
    qRM = complex(0.0,0.0) 
    den = complex(0.0,0.0)

    RI = R_from_axis_angle(iTH, i_n)
    RJ = R_from_axis_angle(jTH, j_n)

    sine_of_prof_r = 0.0

    Xt = zeros(3)

    @inbounds for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        Xt[1] = x[1][i]-X[1];
        Xt[2] = x[2][j]-X[2];
        Xt[3] = x[3][k]-X[3];

        Xt = RJ*Xt;

        r = sqrt( Xt[1]^2 + Xt[2]^2 + Xt[3]^2 )

        sine_of_prof_r = sin(prof(r))

        zRM = ( Xt[1] + 1.0im*Xt[2] )/(r + Xt[3])

        pRM = pfn(zRM)
        qRM = qfn(zRM)

        den = real( qRM*conj(qRM) + pRM*conj(pRM) )

        skyrmion.phi[i,j,k,1] = (sine_of_prof_r/den)*real( pRM*conj(qRM) + qRM*conj(pRM) )
        skyrmion.phi[i,j,k,2] = (sine_of_prof_r/den)*imag( pRM*conj(qRM) - qRM*conj(pRM) )
        skyrmion.phi[i,j,k,3] = (sine_of_prof_r/den)*real( qRM*conj(qRM) - pRM*conj(pRM) )
        skyrmion.phi[i,j,k,4] = cos(prof(r))

        if iTH != 0.0
            skyrmion.phi[i,j,k,1:3] = RI*skyrmion.phi[i,j,k,1:3]
        end

    end
    
end







# ADHM stuff



"""
    makeADHM!(skyrmion, L, M )
    
Writes an ADHM skyrmion in to `skyrmion`. The ADHM data is given by L and M. L and M can be given by `Bx4` and `BxBx4` arrays or as `B` and `BxB` arrays of Quaternions, from the `Quaternionic` package.

# Example

"""
function makeADHM!(an_ADHM_skyrmion, L, M)

    B = size(L)[1]

    L_final = zeros(B,4)
    M_final = zeros(B,B,4)


   if typeof(L[end]) == QuaternionF64

        for a in 1:B
            L_final[a,1] = L[a].w
            L_final[a,2] = L[a].x
            L_final[a,3] = L[a].y
            L_final[a,4] = L[a].z
        end
    else
        for a in 1:B, c in 1:4
            L_final[a,c] = L[a,c]
        end
    end

    if typeof(M[end]) == QuaternionF64

        for a in 1:B, b in 1:B
            M_final[a,b,1] = M[a,b].w
            M_final[a,b,2] = M[a,b].x
            M_final[a,b,3] = M[a,b].y
            M_final[a,b,4] = M[a,b].z
        end
    else

        for a in 1:B, b in 1:B, c in 1:4
            M_final[a,b,c] = M[a,b,c]
        end
    end



    tsteps = 42
    lstime = pi/tsteps

    ctL = [ cos(lstime*(tint-1)) for tint in 1:tsteps+1 ]
    stL = [ sin(lstime*(tint-1)) for tint in 1:tsteps+1 ]




    x = an_ADHM_skyrmion.x
    lp = an_ADHM_skyrmion.lp

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
        an_ADHM_skyrmion.phi[i,j,k,:] = ADHMpt2(L_final,M_final,[x[1][i],x[2][j],x[3][k]], B, tsteps,ctL,stL)
    end

end

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

function makeQmultmatrix2!(qM,q)

    qM[1,1] = q[1];  qM[1,2] = -q[2]; qM[1,3] = -q[3]; qM[1,4] = -q[4]; 
    qM[2,1] = q[2];  qM[2,2] = q[1]; qM[2,3] = -q[4]; qM[2,4] = q[3]; 
    qM[3,1] = q[3];  qM[3,2] = q[4]; qM[3,3] = q[1]; qM[3,4] = -q[2]; 
    qM[4,1] = q[4];  qM[4,2] = -q[3]; qM[4,3] = q[2]; qM[4,4] = q[1]; 

end

function third_order_update!(Q,q1,q2,q3)

    Q[1,1] = (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0
    Q[1,2] = (-4.0*q1[2]*q2[1] - 4.0*q1[1]*q2[2] + 4.0*q1[4]*q2[3] - 4.0*q1[3]*q2[4] + q3[2])/3.0
    Q[1,3] = (-4.0*q1[3]*q2[1] - 4.0*q1[4]*q2[2] - 4.0*q1[1]*q2[3] + 4.0*q1[2]*q2[4] + q3[3])/3.0
    Q[1,4] = (-4.0*q1[4]*q2[1] + 4.0*q1[3]*q2[2] - 4.0*q1[2]*q2[3] - 4.0*q1[1]*q2[4] + q3[4])/3.0
    Q[2,1] = (4.0*q1[2]*q2[1] + 4.0*q1[1]*q2[2] - 4.0*q1[4]*q2[3] + 4.0*q1[3]*q2[4] - q3[2])/3.0
    Q[2,2] = (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0
    Q[2,3] = (-4.0*q1[4]*q2[1] + 4.0*q1[3]*q2[2] - 4.0*q1[2]*q2[3] - 4.0*q1[1]*q2[4] + q3[4])/3.0
    Q[2,4] = (4.0*q1[3]*q2[1] + 4.0*q1[4]*q2[2] + 4.0*q1[1]*q2[3] - 4.0*q1[2]*q2[4] - q3[3])/3.0
    Q[3,1] = (4.0*q1[3]*q2[1] + 4.0*q1[4]*q2[2] + 4.0*q1[1]*q2[3] - 4.0*q1[2]*q2[4] - q3[3])/3.0
    Q[3,2] = (4.0*q1[4]*q2[1] - 4.0*q1[3]*q2[2] + 4.0*q1[2]*q2[3] + 4.0*q1[1]*q2[4] - q3[4])/3.0
    Q[3,3] = (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0
    Q[3,4] = (-4.0*q1[2]*q2[1] - 4.0*q1[1]*q2[2] + 4.0*q1[4]*q2[3] - 4.0*q1[3]*q2[4] + q3[2])/3.0
    Q[4,1] = (4.0*q1[4]*q2[1] - 4.0*q1[3]*q2[2] + 4.0*q1[2]*q2[3] + 4.0*q1[1]*q2[4] - q3[4])/3.0
    Q[4,2] = (-4.0*q1[3]*q2[1] - 4.0*q1[4]*q2[2] - 4.0*q1[1]*q2[3] + 4.0*q1[2]*q2[4] + q3[3])/3.0
    Q[4,3] = (4.0*q1[2]*q2[1] + 4.0*q1[1]*q2[2] - 4.0*q1[4]*q2[3] + 4.0*q1[3]*q2[4] - q3[2])/3.0
    Q[4,4] = (4.0*q1[1]*q2[1] - 4.0*q1[2]*q2[2] - 4.0*q1[3]*q2[3] - 4.0*q1[4]*q2[4] - q3[1])/3.0

end

function B3_tet_data(lam)

    L = zeros(3,4)
    M = zeros(3,3,4)

    L[1,2] = lam
    L[2,3] = lam
    L[3,4] = lam

    M[1,2,4] = lam;    M[1,3,3] = lam;
    M[2,1,4] = lam;    M[2,3,2] = lam;
    M[3,1,3] = lam;    M[3,2,2] = lam;

    return L, M

end

function B4_cube_data(lam)

    L = zeros(4,4)
    M = zeros(4,4,4)

    L[1,1] = lam
    L[2,2] = lam
    L[3,3] = lam
    L[4,4] = lam

    lam /= sqrt(2.0)

    M[1,2,3] = -lam; M[1,2,4] = -lam;
    M[1,3,2] = -lam; M[1,3,4] = -lam;
    M[1,4,2] = -lam; M[1,4,3] = -lam;

    M[2,3,2] = -lam; M[2,3,3] =  lam;
    M[2,4,2] =  lam; M[2,4,4] = -lam;

    M[3,4,3] = -lam; M[3,4,4] =  lam;

    for a in 1:4, b in 1:a, c in 1:4
        M[a,b,c] = M[b,a,c]
    end

    return L, M

end


function ADHMpt2(L,M,y,B,tsteps,ctL,stL)

    if B == 2
        Rnm = zeros( MMatrix{2,2,Float64} )
        iRnm = zeros( MMatrix{2,2,Float64} )
    elseif B == 3
        Rnm = zeros( MMatrix{3,3,Float64} )
        iRnm = zeros( MMatrix{3,3,Float64} )
    elseif B == 4
        Rnm = zeros( MMatrix{4,4,Float64} )
        iRnm = zeros( MMatrix{4,4,Float64} )
    elseif B == 5
        Rnm = zeros( MMatrix{5,5,Float64} )
        iRnm = zeros( MMatrix{5,5,Float64} )
    elseif B == 6
        Rnm = zeros( MMatrix{6,6,Float64} )
        iRnm = zeros( MMatrix{6,6,Float64} )
    elseif B == 7
        Rnm = zeros( MMatrix{7,7,Float64} )
        iRnm = zeros( MMatrix{7,7,Float64} )
    elseif B == 8
        Rnm = zeros( MMatrix{8,8,Float64} )
        iRnm = zeros( MMatrix{8,8,Float64} )
    else
        Rnm = zeros(B, B)
        iRnm = zeros(B,B)
    end

    Ln = zeros(B,4)
    Mn = zeros(B,B,4)

    U = zeros( MVector{4,Float64} )
    U[1] = 1.0


    Mmdysp = zeros(B,B,4)
    p = zeros(B,4)


    Ω1M = zeros( MMatrix{4,4,Float64} )

    Ω1 = zeros( MVector{4,Float64} )
    Ω2 = zeros( MVector{4,Float64} )
    Ω3 = zeros( MVector{4,Float64} )

    lstime = pi/tsteps;

    allN = zeros(tsteps+1,B+1,4)

    for tint in 1:tsteps+1

        ct = ctL[tint]

        @inbounds x1t = y[1]*ct; x2t = y[2]*ct; x3t = y[3]*ct;

        @inbounds for c in 1:4
            for a in 1:B
                Ln[a,c] = ct*L[a,c]
                for b in 1:B
                    Mn[a,b,c] = ct*M[a,b,c]
                end
            end
        end

        
        Nfy!(allN,tint,Ln,Mn,[stL[tint],x1t,x2t,x3t], Mmdysp ,p, B,ct,Rnm,iRnm)
        
    end

    for tint in 1:2:tsteps-1

        #= FIRST ORDER 
        getΩ!(Ω1,allN[tint+1,:,:],allN[tint,:,:],B)
        Ω1M = makeQmultmatrix2(Ω1)
        U = Ω1M*U =#

        getΩf!(Ω1,allN,tint+2,tint+1,B)
        getΩf!(Ω2,allN,tint+1,tint,B)
        getΩf!(Ω3,allN,tint+2,tint,B)

        third_order_update!(Ω1M,Ω1,Ω2,Ω3)

        U = Ω1M*U
        
    end

    normer = sqrt( U[1]^2 + U[2]^2 + U[3]^2 + U[4]^2 )
    @simd for a in 1:4
        @inbounds U[a] /= normer
    end
    
    return [U[2],U[3],U[4],U[1]]

end

function getΩf!(Ω1,allv,t1,t2,B)

    Ω1[1] = 0.0; 
    Ω1[2] = 0.0; 
    Ω1[3] = 0.0; 
    Ω1[4] = 0.0;

    @simd for a in 1:B+1
        @inbounds Ω1[1] += allv[t2,a,1]*allv[t1,a,1] + allv[t2,a,2]*allv[t1,a,2] + allv[t2,a,3]*allv[t1,a,3] + allv[t2,a,4]*allv[t1,a,4]
        @inbounds Ω1[2] += allv[t2,a,2]*allv[t1,a,1] - allv[t2,a,1]*allv[t1,a,2] - allv[t2,a,4]*allv[t1,a,3] + allv[t2,a,3]*allv[t1,a,4]
        @inbounds Ω1[3] += allv[t2,a,3]*allv[t1,a,1] + allv[t2,a,4]*allv[t1,a,2] - allv[t2,a,1]*allv[t1,a,3] - allv[t2,a,2]*allv[t1,a,4]
        @inbounds Ω1[4] += allv[t2,a,4]*allv[t1,a,1] - allv[t2,a,3]*allv[t1,a,2] + allv[t2,a,2]*allv[t1,a,3] - allv[t2,a,1]*allv[t1,a,4]
    end

end

function getΩ!(Ω1,vp,v,B)

    Ω1[1] = 0.0; 
    Ω1[2] = 0.0; 
    Ω1[3] = 0.0; 
    Ω1[4] = 0.0;

    @inbounds for a in 1:B+1
        Ω1[1] += v[a,1]*vp[a,1] + v[a,2]*vp[a,2] + v[a,3]*vp[a,3] + v[a,4]*vp[a,4]
        Ω1[2] += v[a,2]*vp[a,1] - v[a,1]*vp[a,2] - v[a,4]*vp[a,3] + v[a,3]*vp[a,4]
        Ω1[3] += v[a,3]*vp[a,1] + v[a,4]*vp[a,2] - v[a,1]*vp[a,3] - v[a,2]*vp[a,4]
        Ω1[4] += v[a,4]*vp[a,1] - v[a,3]*vp[a,2] + v[a,2]*vp[a,3] - v[a,1]*vp[a,4]
    end

end





function Nfy!(Nα,tint,L,M,y,Mmdysp,p,B,ct,Rnm,iRnm)


    for a in 1:B, c in 1:4
        p[a,c] = 0.0
    end
    
        @inbounds for b in 1:B, a in 1:4
        for c in 1:B
            Mmdysp[b,c,a] = M[b,c,a]
        end
        Mmdysp[b,b,a] -= y[a]
        
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
    @simd for a in 1:B
        for  b in 1:B
            @inbounds @fastmath Rnm[a,b] = L[a,1]*L[b,1] + L[a,2]*L[b,2] + L[a,3]*L[b,3] + L[a,4]*L[b,4]
            for c in 1:B, d in 1:4
                @inbounds @fastmath Rnm[a,b] += (M[a,c,d]*M[c,b,d])
            end
        end
    end    
end
























