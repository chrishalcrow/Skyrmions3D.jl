
function engpt(dp,p4,mpi)
    return 2.0*mpi*(1.0 - p4) + (dp[1,1]^2 + dp[1,2]^2 + dp[1,3]^2 + dp[1,4]^2 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,3]^2 + dp[3,4]^2) + (dp[1,4]^2*dp[2,1]^2 + dp[1,4]^2*dp[2,2]^2 + dp[1,4]^2*dp[2,3]^2 + dp[1,1]^2*(dp[2,2]^2 + dp[2,3]^2) - 2*dp[1,1]*dp[1,4]*dp[2,1]*dp[2,4] + dp[1,1]^2*dp[2,4]^2 + dp[1,4]^2*dp[3,1]^2 + dp[2,2]^2*dp[3,1]^2 + dp[2,3]^2*dp[3,1]^2 + dp[2,4]^2*dp[3,1]^2 - 2*dp[2,1]*dp[2,2]*dp[3,1]*dp[3,2] + dp[1,1]^2*dp[3,2]^2 + dp[1,4]^2*dp[3,2]^2 + dp[2,1]^2*dp[3,2]^2 + dp[2,3]^2*dp[3,2]^2 + dp[2,4]^2*dp[3,2]^2 - 2*dp[2,1]*dp[2,3]*dp[3,1]*dp[3,3] - 2*dp[2,2]*dp[2,3]*dp[3,2]*dp[3,3] + dp[1,1]^2*dp[3,3]^2 + dp[1,4]^2*dp[3,3]^2 + dp[2,1]^2*dp[3,3]^2 + dp[2,2]^2*dp[3,3]^2 + dp[2,4]^2*dp[3,3]^2 - 2*(dp[1,1]*dp[1,4]*dp[3,1] + dp[2,4]*(dp[2,1]*dp[3,1] + dp[2,2]*dp[3,2] + dp[2,3]*dp[3,3]))*dp[3,4] + (dp[1,1]^2 + dp[2,1]^2 + dp[2,2]^2 + dp[2,3]^2)*dp[3,4]^2 + dp[1,3]^2*(dp[2,1]^2 + dp[2,2]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,2]^2 + dp[3,4]^2) + dp[1,2]^2*(dp[2,1]^2 + dp[2,3]^2 + dp[2,4]^2 + dp[3,1]^2 + dp[3,3]^2 + dp[3,4]^2) - 2*dp[1,2]*(dp[1,1]*(dp[2,1]*dp[2,2] + dp[3,1]*dp[3,2]) + dp[1,3]*(dp[2,2]*dp[2,3] + dp[3,2]*dp[3,3]) + dp[1,4]*(dp[2,2]*dp[2,4] + dp[3,2]*dp[3,4])) - 2*dp[1,3]*(dp[1,1]*(dp[2,1]*dp[2,3] + dp[3,1]*dp[3,3]) + dp[1,4]*(dp[2,3]*dp[2,4] + dp[3,3]*dp[3,4])))
end

function Energy(sk, mpi)

    ED = zeros(sk.lp[1], sk.lp[2], sk.lp[3])

    engtot = 0.0
    dp = zeros(3,4)

    for i in 3:sk.lp[1]-2, j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
    
        getDX!(dp, sk ,i, j, k )
        @inbounds ED[i,j,k] = engpt(dp,sk.phi[i,j,k,4],mpi)
        
    end

    engtot = sum(ED)*sk.ls[1]*sk.ls[2]*sk.ls[3]/(12.0*pi^2)

    return engtot

end 


function EnergyD(sk, mpi)

	engD = zeros(sk.lp[1], sk.lp[2], sk.lp[3] )
    dp = zeros(3,4)

    @inbounds for i in 3:sk.lp[1]-2, j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
    
        getDX!(dp, sk ,i, j, k )
        engD[i,j,k] = engpt(dp,sk.phi[i,j,k,4],mpi)
        
    end

    return engD

end 

function BaryonD(sk)

	barD = zeros(sk.lp[1], sk.lp[2], sk.lp[3] )

    dp = zeros(3,4)
    pp = zeros(4)

    @inbounds for i in 3:sk.lp[1]-2, j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
    
        getDX!(dp, sk ,i, j, k )
        getX!(pp,sk,i,j,k)

        barD[i,j,k] = barypt(dp,pp)
        
    end
    
    return barD

end 


function barypt(dp,pp)
    return pp[4]*dp[1,3]*dp[2,2]*dp[3,1] - pp[3]*dp[1,4]*dp[2,2]*dp[3,1] - pp[4]*dp[1,2]*dp[2,3]*dp[3,1] + pp[2]*dp[1,4]*dp[2,3]*dp[3,1] + pp[3]*dp[1,2]*dp[2,4]*dp[3,1] - pp[2]*dp[1,3]*dp[2,4]*dp[3,1] - pp[4]*dp[1,3]*dp[2,1]*dp[3,2] + pp[3]*dp[1,4]*dp[2,1]*dp[3,2] + pp[4]*dp[1,1]*dp[2,3]*dp[3,2] - pp[1]*dp[1,4]*dp[2,3]*dp[3,2] - pp[3]*dp[1,1]*dp[2,4]*dp[3,2] + pp[1]*dp[1,3]*dp[2,4]*dp[3,2] + pp[4]*dp[1,2]*dp[2,1]*dp[3,3] - pp[2]*dp[1,4]*dp[2,1]*dp[3,3] - pp[4]*dp[1,1]*dp[2,2]*dp[3,3] + pp[1]*dp[1,4]*dp[2,2]*dp[3,3] + pp[2]*dp[1,1]*dp[2,4]*dp[3,3] - pp[1]*dp[1,2]*dp[2,4]*dp[3,3] - pp[3]*dp[1,2]*dp[2,1]*dp[3,4] + pp[2]*dp[1,3]*dp[2,1]*dp[3,4] + pp[3]*dp[1,1]*dp[2,2]*dp[3,4] - pp[1]*dp[1,3]*dp[2,2]*dp[3,4] - pp[2]*dp[1,1]*dp[2,3]*dp[3,4] + pp[1]*dp[1,2]*dp[2,3]*dp[3,4]
end

function Baryon(sk)

    BD = zeros(sk.lp[1],sk.lp[2],sk.lp[3])

    bartot = 0.0
    dp = zeros(3,4)
    pp = zeros(4)
    
    for i in 3:sk.lp[1]-2, j in 3:sk.lp[2]-2, k in 3:sk.lp[3]-2
        
        getDX!(dp, sk ,i, j, k )
        getX!(pp,sk,i,j,k)
        
        BD[i,j,k] = barypt(dp,pp)
            
    end
    
    bartot = sum(BD)*sk.ls[1]*sk.ls[2]*sk.ls[3]/(2.0*pi^2)
    
    return bartot
    
end 



function getMOI(sk,x)
	
	epsilon = zeros(3,3,3)
	
	epsilon[1,2,3] = 1.0
	epsilon[1,3,2] = -1.0
	epsilon[2,3,1] = 1.0
	epsilon[2,1,3] = -1.0
	epsilon[3,1,2] = 1.0
	epsilon[3,2,1] = -1.0
	
	VV = zeros(6,6)
	
    mm = zeros(3,4)
    po = zeros(4)

    TT = zeros(3,3)
    GG = zeros(6,3)
    RR = zeros(3,3)
	
	xxx = zeros(3)

	for i = 3:sk.lp[1]-2,  j = 3:sk.lp[2]-2, k = 3:sk.lp[3]-2
	
	    mm = zeros(3,4)
	    po = zeros(4)

	    TT = zeros(3,3)
	    GG = zeros(6,3)
	    RR = zeros(3,3)

	    xxx = [x[1][i], x[2][j], x[3][k]]
		
        getDX!(mm, sk ,i, j, k)
		getX!(po, sk, i, j, k)

	    for a = 1:3, b = 1:3

	        RR[a,b] += po[4]*mm[a,b] - mm[a,4]*po[b]
	        GG[a+3,b] -= po[a]*po[b]
	        GG[a+3,a] += po[b]*po[b]

	        for c = 1:3
	            GG[a+3,b] -= po[4]*po[c]*epsilon[a,c,b]
	            for d = 1:3
	                RR[a,b] += epsilon[b,c,d]*mm[a,c]*po[d]
	            end
	        end
	    end

	    for a = 1:3, b = 1:3, c = 1:3, d = 1:3
	    	GG[a,d] += epsilon[a,b,c]*xxx[b]*RR[c,d]
	    end

        for a = 1:6, b = 1:6, c = 1:3

        	VV[a,b] += GG[a,c]*GG[b,c]

            for d = 1:3,  e = 1:3
            	VV[a,b] -= 2.0*(GG[a,c]*GG[b,d]*RR[e,c]*RR[e,d] - GG[a,c]*GG[b,c]*RR[d,e]*RR[d,e])
            end
		end
		
	end

	for a = 1:6, b = 1:6
        VV[a,b] *= sk.ls[1]*sk.ls[2]*sk.ls[3]
    end
	
	print(VV)
	
end







