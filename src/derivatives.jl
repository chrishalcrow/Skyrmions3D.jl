
function dxD(phi, a, i, j, k, lsx)
    @inbounds (-phi[i+2,j,k,a] + 8.0*phi[i+1,j,k,a] - 8.0*phi[i-1,j,k,a] + phi[i-2,j,k,a])/(12.0*lsx)
end
function dyD(phi, a, i, j, k, lsy)
    @inbounds (-phi[i,j+2,k,a] + 8.0*phi[i,j+1,k,a] - 8.0*phi[i,j-1,k,a] + phi[i,j-2,k,a])/(12.0*lsy)
end
function dzD(phi, a, i, j, k, lsz)
    @inbounds (-phi[i,j,k+2,a] + 8.0*phi[i,j,k+1,a] - 8.0*phi[i,j,k-1,a] + phi[i,j,k-2,a])/(12.0*lsz)
end

function d2xD(phi, a, i, j, k, lsx)
    @inbounds (-phi[i+2,j,k,a] + 16.0*phi[i+1,j,k,a] - 30.0*phi[i,j,k,a] + 16.0*phi[i-1,j,k,a] - phi[i-2,j,k,a])/(12.0*lsx^2)
end
function d2yD(phi, a, i, j, k, lsx)
    @inbounds (-phi[i,j+2,k,a] + 16.0*phi[i,j+1,k,a] - 30.0*phi[i,j,k,a] + 16.0*phi[i,j-1,k,a] - phi[i,j-2,k,a])/(12.0*lsx^2)
end
function d2zD(phi, a, i, j, k, lsx)
    @inbounds (-phi[i,j,k+2,a] + 16.0*phi[i,j,k+1,a] - 30.0*phi[i,j,k,a] + 16.0*phi[i,j,k-1,a] - phi[i,j,k-2,a])/(12.0*lsx^2)
end

function dxdydiffD(phi, a, i, j, k, lsx, lsy)
    @inbounds (-phi[i+2,j+2,k,a] + 16.0*phi[i+1,j+1,k,a] - 30.0*phi[i,j,k,a] + 16.0*phi[i-1,j-1,k,a] - phi[i-2,j-2,k,a])/(12.0*lsx*lsy)
end
function dxdzdiffD(phi, a, i, j, k, lsx, lsy)
    @inbounds (-phi[i+2,j,k+2,a] + 16.0*phi[i+1,j,k+1,a] - 30.0*phi[i,j,k,a] + 16.0*phi[i-1,j,k-1,a] - phi[i-2,j,k-2,a])/(12.0*lsx*lsy)
end
function dydzdiffD(phi, a, i, j, k, lsx, lsy)
    @inbounds (-phi[i,j+2,k+2,a] + 16.0*phi[i,j+1,k+1,a] - 30.0*phi[i,j,k,a] + 16.0*phi[i,j-1,k-1,a] - phi[i,j-2,k-2,a])/(12.0*lsx*lsy)
end




#=
function d2xD(phi, a, i, j, k, lsx)
    @inbounds (-phi[a,i+2,j,k] + 16.0*phi[a,i+1,j,k] - 30.0*phi[a,i,j,k] + 16.0*phi[a,i-1,j,k] - phi[a,i-2,j,k])/(12.0*lsx^2)
end
function d2yD(phi, a, i, j, k, lsx)
    @inbounds (-phi[a,i,j+2,k] + 16.0*phi[a,i,j+1,k] - 30.0*phi[a,i,j,k] + 16.0*phi[a,i,j-1,k] - phi[a,i,j-2,k])/(12.0*lsx^2)
end
function d2zD(phi, a, i, j, k, lsx)
    @inbounds (-phi[a,i,j,k+2] + 16.0*phi[a,i,j,k+1] - 30.0*phi[a,i,j,k] + 16.0*phi[a,i,j,k-1] - phi[a,i,j,k-2])/(12.0*lsx^2)
end

function dxdydiffD(phi, a, i, j, k, lsx, lsy)
    @inbounds (-phi[a,i+2,j+2,k] + 16.0*phi[a,i+1,j+1,k] - 30.0*phi[a,i,j,k] + 16.0*phi[a,i-1,j-1,k] - phi[a,i-2,j-2,k])/(12.0*lsx*lsy)
end
function dxdzdiffD(phi, a, i, j, k, lsx, lsy)
    @inbounds (-phi[a,i+2,j,k+2] + 16.0*phi[a,i+1,j,k+1] - 30.0*phi[a,i,j,k] + 16.0*phi[a,i-1,j,k-1] - phi[a,i-2,j,k-2])/(12.0*lsx*lsy)
end
function dydzdiffD(phi, a, i, j, k, lsx, lsy)
    @inbounds (-phi[a,i,j+2,k+2] + 16.0*phi[a,i,j+1,k+1] - 30.0*phi[a,i,j,k] + 16.0*phi[a,i,j-1,k-1] - phi[a,i,j-2,k-2])/(12.0*lsx*lsy)
end
=#


function getDXf!(dp, phi,i,j,k,ls)
    @simd for a in 1:4
        @inbounds dp[1,a] = dxD(phi,a,i,j,k,ls[1])
        @inbounds dp[2,a] = dyD(phi,a,i,j,k,ls[2])
        @inbounds dp[3,a] = dzD(phi,a,i,j,k,ls[3]) 
    end
end

function getDDXf!(ddp, phi ,i,j,k,ls)
    
    for a in 1:4
            
        @inbounds ddp[1,a] = d2xD(phi,a,i,j,k,ls[1])
        @inbounds ddp[2,a] = d2yD(phi,a,i,j,k,ls[2])
        @inbounds ddp[3,a] = d2zD(phi,a,i,j,k,ls[3]) 

        @inbounds ddp[6,a] = (dxdydiffD(phi, a, i, j, k, ls[1], ls[2]) - ddp[1,a] - ddp[2,a])/2.0
        @inbounds ddp[5,a] = (dxdzdiffD(phi, a, i, j, k, ls[1], ls[3]) - ddp[1,a] - ddp[3,a])/2.0
        @inbounds ddp[4,a] = (dydzdiffD(phi, a, i, j, k, ls[2], ls[3]) - ddp[2,a] - ddp[3,a])/2.0

    end
end





function getDX_SA!(dp, ϕ,i,j,k)
    @simd for a in 1:4
        @inbounds dp[1,a] = dxD(ϕ.phi,a,i,j,k,ϕ.ls[1])
        @inbounds dp[2,a] = dyD(ϕ.phi,a,i,j,k,ϕ.ls[2])
        @inbounds dp[3,a] = dzD(ϕ.phi,a,i,j,k,ϕ.ls[3]) 
    end
end

function getDDX_SA!(ddp, sk ,i,j,k)
    
    @simd for a in 1:4
            
        @inbounds ddp[1,1,a] = d2xD(sk.phi,a,i,j,k,sk.ls[1])
        @inbounds ddp[2,2,a] = d2yD(sk.phi,a,i,j,k,sk.ls[2])
        @inbounds ddp[3,3,a] = d2zD(sk.phi,a,i,j,k,sk.ls[3]) 

        @inbounds ddp[1,2,a] = (dxdydiffD(sk.phi, a, i, j, k, sk.ls[1], sk.ls[2]) - ddp[1,a] - ddp[2,a])/2.0
        @inbounds ddp[1,3,a] = (dxdzdiffD(sk.phi, a, i, j, k, sk.ls[1], sk.ls[3]) - ddp[1,a] - ddp[3,a])/2.0
        @inbounds ddp[2,3,a] = (dydzdiffD(sk.phi, a, i, j, k, sk.ls[2], sk.ls[3]) - ddp[2,a] - ddp[3,a])/2.0

        @inbounds ddp[2,1,a] = ddp[1,2,a]
        @inbounds ddp[3,2,a] = ddp[2,3,a]
        @inbounds ddp[3,1,a] = ddp[1,3,a]

    end
end

function getDX(ϕ,i,j,k)

    return SMatrix{3,4,Float64, 12}(
        dxD(ϕ.phi,1,i,j,k,ϕ.ls[1]),
        dyD(ϕ.phi,1,i,j,k,ϕ.ls[2]),
        dzD(ϕ.phi,1,i,j,k,ϕ.ls[3]),

        dxD(ϕ.phi,2,i,j,k,ϕ.ls[1]),
        dyD(ϕ.phi,2,i,j,k,ϕ.ls[2]),
        dzD(ϕ.phi,2,i,j,k,ϕ.ls[3]),

        dxD(ϕ.phi,3,i,j,k,ϕ.ls[1]),
        dyD(ϕ.phi,3,i,j,k,ϕ.ls[2]),
        dzD(ϕ.phi,3,i,j,k,ϕ.ls[3]),

        dxD(ϕ.phi,4,i,j,k,ϕ.ls[1]),
        dyD(ϕ.phi,4,i,j,k,ϕ.ls[2]),
        dzD(ϕ.phi,4,i,j,k,ϕ.ls[3])
    )
    

end

function getDX!(dp, ϕ,i,j,k)
    @simd for a in 1:4
        @inbounds dp[1,a] = dxD(ϕ.phi,a,i,j,k,ϕ.ls[1])
        @inbounds dp[2,a] = dyD(ϕ.phi,a,i,j,k,ϕ.ls[2])
        @inbounds dp[3,a] = dzD(ϕ.phi,a,i,j,k,ϕ.ls[3]) 
    end
end

function getDDX!(ddp, sk ,i,j,k)
    
    @simd for a in 1:4
            
        @inbounds ddp[1,a] = d2xD(sk.phi,a,i,j,k,sk.ls[1])
        @inbounds ddp[2,a] = d2yD(sk.phi,a,i,j,k,sk.ls[2])
        @inbounds ddp[3,a] = d2zD(sk.phi,a,i,j,k,sk.ls[3]) 

        @inbounds ddp[6,a] = (dxdydiffD(sk.phi, a, i, j, k, sk.ls[1], sk.ls[2]) - ddp[1,a] - ddp[2,a])/2.0
        @inbounds ddp[5,a] = (dxdzdiffD(sk.phi, a, i, j, k, sk.ls[1], sk.ls[3]) - ddp[1,a] - ddp[3,a])/2.0
        @inbounds ddp[4,a] = (dydzdiffD(sk.phi, a, i, j, k, sk.ls[2], sk.ls[3]) - ddp[2,a] - ddp[3,a])/2.0

    end
end

function getX!(pp,sk,i,j,k)
    @inbounds for a in 1:4
        pp[a] = sk.phi[i,j,k,a]
    end
end
