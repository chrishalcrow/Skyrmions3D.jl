
"""
    make_RM_product!(skyrmion, X_list) 

Makes a product approximation of many rational map skyrmions, determined through the  list `X_list`. The final field is written into `skyrmion`.

The formatting of the list is as follow:
`X_list = [ data_1, data_2, data_3, ... ]`
where
`data_1 = [ p(z), q(z), f(r), X, θiso, n_iso, θrot, n_rot ]`

See also [`product`]

# Example of list
```
p1(z) = z; q1(z) = 1; f1(r) = 4*atan(exp(-r));
p2(z) = z^2; q2(z) = 1; f2(r) = 4*atan(exp(-0.7*r));
X_list = [ [ p1, q1, f1, [0.0,0.0,1.5], 0.0, [0.0,0.0,1.0], 0.0, [0.0,0.0,1.0] ], [ p2, q2, f2, [0.0,0.0,-1.5], pi, [1.0,0.0,0.0], 0.0, [0.0,0.0,1.0] ] ]
```

# Technical details

The product is taken pairwise in order. E.g. for a list of 3 skyrmions, we first calculate the symmetrised product of the first and second skyrmions then calculate the symmtrised product with the third skyrmion. Hence the final solution is not symmetric under permutations.

"""
function make_RM_product!(sk, Xs)

    x = sk.x
    lp = sk.lp
    ls = sk.ls

    temp_sk = Skyrmion(lp,ls)

    a=1
    make_rational_map!(sk, Xs[a][1],Xs[a][2],Xs[a][3], X = Xs[a][4], iTH = Xs[a][5], i_n = Xs[a][6], jTH = Xs[a][7], j_n = Xs[a][8]  )
    
    for a in 2:size(Xs)[1]

        make_rational_map!(temp_sk, Xs[a][1],Xs[a][2],Xs[a][3], X = Xs[a][4], iTH = Xs[a][5], i_n = Xs[a][6], jTH = Xs[a][7], j_n = Xs[a][8]  )
        product_approx!(sk, temp_sk)
    end

end

"""
    make_rational_map!(skyrmion, prof, pfn, qfn; kwargs... )
    
Writes a rational map skyrmion in to `skyrmion`. The rational map is given by the polynomials R(z) = p(z)/q(z) and the profile f(r).

If no `f` is given, the function will find an OK approximation for the profile.

# Optional arguments
-  `X=[0.0,0.0,0.0]`: translate the initial skyrmion by `X`
-  `iTH = 0.0`: isorotate by initial skyrmion by `iTH`
-  `i_n = 0.0`: isorotate initial skyrmion around `i_n`
-  `jTH = 0.0`: isorotate by initial skyrmion by `jTH`
-  `j_n = 0.0`: isorotate initial skyrmion around `j_n`

"""
function make_rational_map!(skyrmion, pfn, qfn, prof; X=[0.0,0.0,0.0], iTH=0.0, i_n = [0.0,0.0,1.0], jTH = 0.0, j_n = [0.0,0.0,0.0] )
    
    lp, x = skyrmion.lp, skyrmion.x

    RI = R_from_axis_angle(iTH, i_n)
    RJ = R_from_axis_angle(jTH, j_n)

    Threads.@threads for k in 1:lp[3]
        @inbounds for j in 1:lp[2], i in 1:lp[1]

            Xto = SVector{3,Float64}( x[1][i]-X[1], x[2][j]-X[2], x[3][k]-X[3] )
            Xt = RJ*Xto;

            r = sqrt( Xt[1]^2 + Xt[2]^2 + Xt[3]^2 )

            sine_of_prof_r = sin(prof(r))

            zRM = ( Xt[1] + 1.0im*Xt[2] )/(r + Xt[3])

            pRM = pfn(zRM)
            qRM = qfn(zRM)

            den = real( qRM*conj(qRM) + pRM*conj(pRM) )

            skyrmion.pion_field[i,j,k,1] = (sine_of_prof_r/den)*real( pRM*conj(qRM) + qRM*conj(pRM) )
            skyrmion.pion_field[i,j,k,2] = (sine_of_prof_r/den)*imag( pRM*conj(qRM) - qRM*conj(pRM) )
            skyrmion.pion_field[i,j,k,3] = (sine_of_prof_r/den)*real( qRM*conj(qRM) - pRM*conj(pRM) )
            skyrmion.pion_field[i,j,k,4] = cos(prof(r))

            if iTH != 0.0
                skyrmion.pion_field[i,j,k,1:3] = RI*skyrmion.pion_field[i,j,k,1:3]
            end

        end
    end

    if skyrmion.periodic == false
        set_dirichlet!(skyrmion)
    end

    #println("hello.")
    
end

function make_rational_map!(skyrmion, pfn, qfn; baryon=0.0, X=[0.0,0.0,0.0], iTH=0.0, i_n = [0.0,0.0,1.0], jTH = 0.0, j_n = [0.0,0.0,0.0] )

    if baryon == 0.0
        baryon1 = abs( (log(pfn(10000)) - log(pfn(1)))/log(10000) )
        baryon2 = abs( (log(qfn(10000)) - log(qfn(1)))/log(10000) )
        baryon = max( round(baryon1), round(baryon2) )
        println("I think your baryon number is ", baryon, ". If it is not, include '; baryon=B' in your argument.")
    end
    
    R(z) = pfn(z)/qfn(z)
    k1,k2=getOKprofile(1.0,1.0,baryon,getI(R),skyrmion.mpi)
    #println(k1)
    #println(k2)
    prof(r) = pi/(1 - tanh(-k2*k1))*( -tanh(k2*(r - k1)) + 1.0  );
    make_rational_map!(skyrmion, pfn, qfn, prof; X=[0.0,0.0,0.0], iTH=0.0, i_n = [0.0,0.0,1.0], jTH = 0.0, j_n = [0.0,0.0,0.0] )
    

end


function getI(R)

    I_tot = 0.0
    dz_r = 0.05
    dz_i = 0.05
    dx=0.0001;
    for z_real in -10-dz_r/2:dz_r:10+dz_r/2, z_imag in -10-dz_i/2:dz_i:10+dz_i/2
        z = z_real + 1.0im*z_imag
        Rp = ((R(z+dx)-R(z-dx))/(2dx) + (R(z+dx*1.0im)-R(z-dx*1.0im))/(2im*dx))/2.0 
        I_tot += real( ( Rp*conj(Rp) )^2*( 1 + z*conj(z) )^2/(1.0 + R(z)*conj(R(z)))^4 )
    end
    
    return I_tot*dz_r*dz_i/pi

end


function getOKprofile(k1,k2,B,I,m)

    #=dk1=0.001;
    dk2=0.001;

    for _ in 1:10

        dE = [(energy_test(k1+dk1,k2,(B,I,m)) - energy_test(k1-dk1,k2,(B,I,m)))/(2*dk1)   (energy_test(k1,k2+dk2,(B,I,m)) - energy_test(k1,k2-dk2,(B,I,m)))/(2*dk1) ]

        ddE = [ (energy_test(k1+dk1,k2,(B,I,m)) - 2.0*energy_test(k1,k2,(B,I,m)) + energy_test(k1-dk1,k2,(B,I,m)))/(dk1^2) ( energy_test(k1+dk1,k2+dk2,(B,I,m)) + energy_test(k1-dk1,k2-dk2,(B,I,m)) - energy_test(k1+dk1,k2-dk2,(B,I,m)) - energy_test(k1-dk1,k2+dk2,(B,I,m)) )/(4*dk1*dk2)  ;
                ( energy_test(k1+dk1,k2+dk2,(B,I,m)) + energy_test(k1-dk1,k2-dk2,(B,I,m)) - energy_test(k1+dk1,k2-dk2,(B,I,m)) - energy_test(k1-dk1,k2+dk2,(B,I,m)) )/(4*dk1*dk2)  (energy_test(k1,k2+dk2,(B,I,m)) - 2.0*energy_test(k1,k2,(B,I,m)) + energy_test(k1,k2-dk2,(B,I,m)))/(dk2^2) ]

        change = inv(ddE)*dE'

        k1 -= change[1]
        k2 -= change[2]
    
    end

    return k1, k2=#

    dk1=0.001;

    for _ in 1:10

        dE = (energy_test(k1+dk1,k2,(B,I,m)) - energy_test(k1-dk1,k2,(B,I,m)))/(2*dk1)   

        ddE =  (energy_test(k1+dk1,k2,(B,I,m)) - 2.0*energy_test(k1,k2,(B,I,m)) + energy_test(k1-dk1,k2,(B,I,m)))/(dk1^2)  

        change = dE/ddE

        k1 -= change[1]
        #k2 -= change[2]
    
    end

    return k1, 1.0



end


function energy_test(k1,k2,(B,I,m);lp=500,ls=0.05,test_prof = profile(lp, ls))
    
    test_prof.field .= pi/(1.0 - tanh(-k2*k1)).*( -tanh.(k2*(test_prof.r_grid .- k1))  .+ 1.0  );
    return energy(test_prof,(B,I,m))

end

function energy(p,(B,I,m))
    
    ED = zeros(p.lp)

    dp = getdpB(p,2)
    ED[2] = Ept([p.field[2], dp, p.r_grid[2], B,I,m])
    for i in 3:p.lp-2
        
        dp = getdpD(p,i)
        ED[i] = Ept([p.field[i], dp, p.r_grid[i], B,I,m])

    end

    return sum(ED)*p.ls/(3pi)

end

# From symbolics code.
function Ept(ˍ₋arg1,)
    #= /Users/chris/.julia/packages/SymbolicUtils/H684H/src/code.jl:350 =#
    #= /Users/chris/.julia/packages/SymbolicUtils/H684H/src/code.jl:351 =#
    #= /Users/chris/.julia/packages/SymbolicUtils/H684H/src/code.jl:352 =#
    begin
        (/)((+)((+)((+)((+)((+)((*)(ˍ₋arg1[5], (^)((sin)(ˍ₋arg1[1]), 4)), (*)((^)(ˍ₋arg1[3], 4), (^)(ˍ₋arg1[2], 2))), (*)((*)(2, (^)(ˍ₋arg1[6], 2)), (^)(ˍ₋arg1[3], 4))), (*)((*)((*)(2, ˍ₋arg1[4]), (^)(ˍ₋arg1[3], 2)), (^)((sin)(ˍ₋arg1[1]), 2))), (*)((*)((*)(-2, (^)(ˍ₋arg1[6], 2)), (^)(ˍ₋arg1[3], 4)), (cos)(ˍ₋arg1[1]))), (*)((*)((*)((*)(2, ˍ₋arg1[4]), (^)(ˍ₋arg1[3], 2)), (^)(ˍ₋arg1[2], 2)), (^)((sin)(ˍ₋arg1[1]), 2))), (^)(ˍ₋arg1[3], 2))
    end
end


function getdpD(p,i)
    return (-p.field[i+2] + 8.0*p.field[i+1] - 8.0*p.field[i-1] + p.field[i-2])/(12.0*p.ls)
end

function getdpB(p,i)
    return (p.field[i+1] - p.field[i-1])/(2.0*p.ls)
end

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

function make_ADHM!(an_ADHM_skyrmion, LM)
    B = size(LM)[2]
    make_ADHM!(an_ADHM_skyrmion, LM[1,1:B], LM[2:B+1,1:B])
end

"""
    make_ADHM!(skyrmion, L, M )
    
Writes an ADHM skyrmion in to `skyrmion`. The ADHM data is given by L and M. L and M can be given by `Bx4` and `BxBx4` arrays or as `B` and `BxB` arrays of Quaternions, from the `GLMakie` package.

# Example of data
```
B=2

L = [ Quaternion(0.0,0.0,0.0,0.0) for a in 1:B ]
M = [ Quaternion(0.0,0.0,0.0,0.0) for a in 1:B, b in 1:B ]

L[1] = Quaternion(0.0, 0.0, 0.0, sqrt(2.0))
L[2] = Quaternion(0.0, 0.0, sqrt(2.0), 0.0)

M[1,1] = Quaternion(1.0, 0.0, 0.0, 0.0)
M[1,2] = Quaternion(0.0, 1.0, 0.0, 0.0)
M[2,1] = Quaternion(0.0, 1.0, 0.0, 0.0)
M[2,2] = Quaternion(-1.0, 0.0, 0.0, 0.0)
```

"""
function make_ADHM!(an_ADHM_skyrmion, L, M)

    B = size(L)[1]

    L_final = zeros(B,4)
    M_final = zeros(B,B,4)

    if typeof(L[end]) == Quaternion{Float64}

        for a in 1:B
            L_final[a,1] = L[a][4]
            L_final[a,2] = L[a][1]
            L_final[a,3] = L[a][2]
            L_final[a,4] = L[a][3]
        end
    else
        for a in 1:B, c in 1:4
            L_final[a,c] = L[a,c]
        end
    end

    if typeof(M[end]) == Quaternion{Float64}

        for a in 1:B, b in 1:B
            M_final[a,b,1] = M[a,b][4]
            M_final[a,b,2] = M[a,b][1]
            M_final[a,b,3] = M[a,b][2]
            M_final[a,b,4] = M[a,b][3]
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

    Threads.@threads for k in 1:lp[3]
        for j in 1:lp[2], i in 1:lp[1]
            @inbounds an_ADHM_skyrmion.pion_field[i,j,k,:] = ADHMpt2(L_final,M_final,[x[1][i],x[2][j],x[3][k]], B, tsteps,ctL,stL)
        end
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





function get_close_ADHM_data(uADHM, iADHM; included_indices = [0,0])


    B = size(uADHM)[2]
    println("Baryon number is ", B)

    if included_indices == [0,0]
        included_indices = vcat( [ [1,a] for a in 1:B ], [ [(b+1),b] for b in 1:B ] )
    end

    println(included_indices)

    ux = zeros(4*B*(B+1))
    u0 = zeros(4*B*(B+1))

    coeffs = zeros(4*B*(B+1))


    count=1
    for i in 1:B+1, j in 1:B, k in 1:4
        ux[count] = uADHM[i,j][k]
        u0[count] = iADHM[i,j][k]
        coeffs[count] = 0.0;
        count += 1
    end
    for i in 1:size(included_indices)[1]
        a = included_indices[i][1];
        b = included_indices[i][2];
        for c in 1:4
            coeffs[Int(B*4*(a-1) + 4*(b-1) + c)] = 1.0
        end
    end

    println(coeffs)

    newp = vcat(ux, coeffs)

    optproblem = OptimizationFunction(to_minimise, Optimization.AutoForwardDiff(), cons=reality)
    prob = OptimizationProblem(optproblem, u0 , newp , lcons = zeros(Float64,Int(7*B*(B-1)/2)), ucons = zeros(Float64,Int(7*(B-1)*(B)/2)) )
    sol = solve(prob,IPNewton())

    newLM = [ Quaternion(0.0,0.0,0.0,0.0) for a in 1:B+1, b in 1:B ]
    for i in 1:B+1, j in 1:B
        newLM[i,j] = Quaternion( sol[ B*4*(i-1) + 4*(j-1) + 1],
        sol[ B*4*(i-1) + 4*(j-1) + 2],
        sol[ B*4*(i-1) + 4*(j-1) + 3],
        sol[ B*4*(i-1) + 4*(j-1) + 4] )
    end

    return newLM

end




function to_minimise(x,p)

    B = Int(1/2*(-1+sqrt(1.0+size(x)[1])))

    ux = p[1:4*B*(B+1)];
    coeffs = p[ 4*B*(B+1)+1: end]

    tot = 0.0

    # keep all terms
    #for i in 1:Int(round(size(x)[1]/4))
    #        tot += coeffs[Int(i)]*(x[Int(i)] - ux[Int(i)])^2
    #end
    for i in 1:Int(round(size(x)[1]))
            tot += coeffs[i]*(x[Int(i)] - ux[Int(i)])^2
    end

    #if included_indices == [[0 0]]
    #included_indices = vcat( [ [1,a] for a in 1:B ], [ [(b+1),b] for b in 1:B ] )
    #end

    # only L and diag...

    #(ux, included_indices) = p
    #included_indices = p[2]

    #println(size(p[2]))

    #=for i in 1:size(included_indices)[1]
        a = included_indices[i][1]
        b = included_indices[i][2]
        for c in 1:4
            tot += coeffs[B*4*(a-1) + 4*(b-1) + c]*(x[B*4*(a-1) + 4*(b-1) + c] - ux[B*4*(a-1) + 4*(b-1) + c])^2
        end
    end=#

    #=for i in 1:4*(B+1)
        tot += (x[i] - ux[i])^2
    end
    for a in 2:B
        for i in a*4*B + 4*(a-1) + 1 : a*4*B + 4*(a-1) + 4
            tot += (x[i] - ux[i])^2
        end
    end
=#

    

    return tot

end






function reality(res,x,p)

    B = Int(1/2*(-1+sqrt(1.0+size(x)[1])))

    count = 0
    for i in 1:B-1, k in i+1:B

        # imposes reality on the upper triangular part of the ADHM data 

        for a in 1:3
            res[7*count + a] = 0.0
        end

        for j in 1:B+1

            

                res[7*count+1] += x[ B*4*(j-1) + 4*(i-1) + 1]*x[ B*4*(j-1) + 4*(k-1) + 4] - x[ B*4*(j-1) + 4*(i-1) + 4]*x[ B*4*(j-1) + 4*(k-1) + 1] + x[ B*4*(j-1) + 4*(i-1) + 2]*x[ B*4*(j-1) + 4*(k-1) + 3] - x[ B*4*(j-1) + 4*(i-1) + 3]*x[ B*4*(j-1) + 4*(k-1) + 2]

                res[7*count+2] += x[ B*4*(j-1) + 4*(i-1) + 2]*x[ B*4*(j-1) + 4*(k-1) + 4] - x[ B*4*(j-1) + 4*(i-1) + 4]*x[ B*4*(j-1) + 4*(k-1) + 2] + x[ B*4*(j-1) + 4*(i-1) + 3]*x[ B*4*(j-1) + 4*(k-1) + 1] - x[ B*4*(j-1) + 4*(i-1) + 1]*x[ B*4*(j-1) + 4*(k-1) + 3]

                res[7*count+3] += x[ B*4*(j-1) + 4*(i-1) + 3]*x[ B*4*(j-1) + 4*(k-1) + 4] - x[ B*4*(j-1) + 4*(i-1) + 4]*x[ B*4*(j-1) + 4*(k-1) + 3] + x[ B*4*(j-1) + 4*(i-1) + 1]*x[ B*4*(j-1) + 4*(k-1) + 2] - x[ B*4*(j-1) + 4*(i-1) + 2]*x[ B*4*(j-1) + 4*(k-1) + 1]

        end

        # makes the ADHM data symmetric
        for a in 1:4
            res[7*count + 3 + a] = x[ B*4*(i+1 - 1) + 4*(k-1) + a] - x[ B*4*(k+1-1) + 4*(i-1) + a]
        end


        count += 1

    end

end
