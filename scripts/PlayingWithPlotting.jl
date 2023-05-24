
1+1

using Revise

using Skyrmions




an_ADHM_skyrmion = Skyrmion([30,30,20],[0.2,0.2,0.2]);


λ=1.0;
B=2

L = zeros(B,4)
M = zeros(B,B,4)

L[1,:] = sqrt(2.0).*[1.0, 0.0, 0.0, 0.0]
L[2,:] = sqrt(2.0).*[0.0, 0.0, 0.0, 1.0]

M[1,1,:] = [0.0, 1.0, 0.0, 0.0]
M[1,2,:] = [0.0, 0.0, 1.0, 0.0]
M[2,1,:] = [0.0, 0.0, 1.0, 0.0]
M[2,2,:] = [0.0, -1.0, 0.0, 0.0]

makeADHM!(an_ADHM_skyrmion, L, M)

an_ADHM_skyrmion.lp[1]

x = an_ADHM_skyrmion.x
lp = an_ADHM_skyrmion.lp

ADHMpt2(L,M,[0.5,-0.3,0.2],B)

for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
    an_ADHM_skyrmion.phi[i,j,k,:] = ADHMpt2(L,M,[x[1][i],x[2][j],x[3][k]],B)
end

display(plot_baryon_density(an_ADHM_skyrmion,iso_value=2.0))


display(plot_field(an_ADHM_skyrmion, component=4, iso_value = 0.5))


ADHMpt2(L,M,[1.0,2.0,3.0],B)





λ=1.0;
B=1

L = zeros(B,4)
M = zeros(B,B,4)

L[1,:] = [1.0, 0.0, 0.0, 0.0]
M[1,1,:] = [0.0, 0.0, 0.0, 0.0]

ADHMpt2(L,M,[0.0,0.0,1.0],B)












an_ADHM_skyrmion.phi[2,3,4,:] = ADHMpt2(L,M,[1.0,2.0,3.0],B)

an_ADHM_skyrmion.x[1][5]


using Quaternions
using LinearAlgebra


B = 2

ΔQ = Matrix{QuaternionF64}(undef,B+1,B)
LQ = Matrix{QuaternionF64}(undef,1,B)
MQ = Matrix{QuaternionF64}(undef,B,B)
MmdyspQ = Matrix{QuaternionF64}(undef,B,B)
NαQ = Vector{QuaternionF64}(undef,B+1)
ysQ = Vector{QuaternionF64}(undef,B)


λ=1.0;

LQ[1,1] = sqrt(2.0).*quat(λ, 0.0, 0.0, 0.0);    LQ[1,2] = sqrt(2.0).*quat(0.0, 0.0, 0.0, λ)

MQ[1,1] = λ*quat(0.0, 1.0, 0.0, 0.0);     MQ[1,2] = λ*quat(0.0, 0.0, 1.0, 0.0)
MQ[2,1] = λ*quat(0.0, 0.0, 1.0, 0.0);     MQ[2,2] = λ*quat(0.0, -1.0, 0.0, 0.0)



an_ADHM_skyrmion = Skyrmion(30,0.2);

lp, x = an_ADHM_skyrmion.lp, an_ADHM_skyrmion.x

for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

    Skpt = ADHMpt(LQ,MQ,[x[1][i],x[2][j],x[3][k]], ysQ, MmdyspQ,B)
    an_ADHM_skyrmion.phi[i,j,k,:] = [imag_part(Skpt)[1], imag_part(Skpt)[2], imag_part(Skpt)[3], real(Skpt) ]
    
end



display(plot_baryon_density(an_ADHM_skyrmion,iso_value=0.8))



function ADHMpt(L,M,y,ysQ,MmdyspQ,B)

    U = quat(1.0,0.0,0,0)
    δt = 0.2

    t = -10.0

    for tint in 1:100

        v0 = NfyQ(L,M,B,quat(t,y[1],y[2],y[3]), ysQ, MmdyspQ )
        vp = NfyQ(L,M,B,quat(t + δt,y[1],y[2],y[3]), ysQ, MmdyspQ       )

        

        U *= (vp')*v0

        t += δt
    end
    return U
end


@btime ADHMpt(LQ,MQ,[0.5,-0.3,0.2],ysQ,MmdyspQ,2)

using BenchmarkTools

ADHMpt2(L,M,[0.5,-0.3,0.2],B)


NfyQ(LQ,MQ,B,quat(0.0,0.5,-0.3,0.2),ysQ,MmdyspQ)




function NfyQ(L,M,B,y,ys,Mmdysp)
    
    Nα = Vector{QuaternionF64}(undef,B+1)
    
    for b in 1:B; ys[b] = y; end
    
    Mmdysp = M - diagm(ys)
    Rnm = makeRnmQ(L,Mmdysp)
    
    p = conj( inv(Rnm)*transpose(L) )
    

    Nα[1] = quat(1.0,0,0,0) - (L*p)[1]
    Nα[2:B+1] = -Mmdysp*p

    norma = conj(Nα[1])*(Nα[1]) + conj(Nα[2])*(Nα[2]) + conj(Nα[3])*(Nα[3])
    
    return Nα ./ sqrt(norma)
    
end




function makeRnmQ(L,Mmdysp)
    return real( (L')*L +  Mmdysp'*Mmdysp ) 
end



















Tempsk = Skyrmion(40,0.2);
B4sk = Skyrmion(40,0.2);

p4(z) = z^4 + 2.0*sqrt(3.0)*im*z^2 + 1.0
q4(z) = z^4 - 2.0*sqrt(3.0)*im*z^2 + 1.0
f4(r) = pi*exp( -(r.^3)./12.0 )


makeRM!(B4sk,f4,p4,q4)

makeRM!(B4sk,f4,p4,q4,[0.0,0.0,0.0])

display(plot_baryon_density(B4sk, iso_value=0.5))


shift!(Tempsk,B4sk,[0.0,0.0,3.0])

isorotate!(Tempsk,B4sk,6pi/3,[0,0,1])

display(plot_baryon_density(Tempsk, iso_value=0.5))

display(plot_field(Tempsk, iso_value=0.01))

display(plot_field(Tempsk, component=4, iso_value=0.3))


display(plot_field(B4sk, component=4, iso_value=0.3))