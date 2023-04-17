
function dxD(phi, a, i, j, k, lsx)
    @inbounds (-phi[a,i+2,j,k] + 8.0*phi[a,i+1,j,k] - 8.0*phi[a,i-1,j,k] + phi[a,i-2,j,k])/(12.0*lsx)
end
function dyD(phi, a, i, j, k, lsy)
    @inbounds (-phi[a,i,j+2,k] + 8.0*phi[a,i,j+1,k] - 8.0*phi[a,i,j-1,k] + phi[a,i,j-2,k])/(12.0*lsy)
end
function dzD(phi, a, i, j, k, lsz)
    @inbounds (-phi[a,i,j,k+2] + 8.0*phi[a,i,j,k+1] - 8.0*phi[a,i,j,k-1] + phi[a,i,j,k-2])/(12.0*lsz)
end

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


function getDXf!(dp, phi,i,j,k,ls)
    @simd for a in 1:4
        @inbounds dp[1,a] = dxD(phi,a,i,j,k,ls[1])
        @inbounds dp[2,a] = dyD(phi,a,i,j,k,ls[2])
        @inbounds dp[3,a] = dzD(phi,a,i,j,k,ls[3]) 
    end
end

function getDDXf!(ddp, phi ,i,j,k,ls)
    
    @simd for a in 1:4
            
        @inbounds ddp[1,a] = d2xD(phi,a,i,j,k,ls[1])
        @inbounds ddp[2,a] = d2yD(phi,a,i,j,k,ls[2])
        @inbounds ddp[3,a] = d2zD(phi,a,i,j,k,ls[3]) 

        @inbounds ddp[6,a] = (dxdydiffD(phi, a, i, j, k, ls[1], ls[2]) - ddp[1,a] - ddp[2,a])/2.0
        @inbounds ddp[5,a] = (dxdzdiffD(phi, a, i, j, k, ls[1], ls[3]) - ddp[1,a] - ddp[3,a])/2.0
        @inbounds ddp[4,a] = (dydzdiffD(phi, a, i, j, k, ls[2], ls[3]) - ddp[2,a] - ddp[3,a])/2.0

    end
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
    for a in 1:4
        @inbounds pp[a] = sk.phi[a,i,j,k]
    end
end
