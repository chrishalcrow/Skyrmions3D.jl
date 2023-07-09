# functions to get local field and their derivatives. Initialised as immutable static arrays, so that they are stored in the stack. Means we can localise, for multithreading, easily.
function getX(sk,i,j,k)

    return SVector{4,Float64}(
        sk.phi[i,j,k,1],
        sk.phi[i,j,k,2],
        sk.phi[i,j,k,3],
        sk.phi[i,j,k,4]
    )
    
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

function getDDX(ϕ,i,j,k)

    return SMatrix{6,4,Float64, 24}(

        d2xD(ϕ.phi,1,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.phi,1,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.phi,1,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.phi,1,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.phi,1,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.phi,1,i,j,k,ϕ.ls[1],ϕ.ls[2]),

        d2xD(ϕ.phi,2,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.phi,2,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.phi,2,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.phi,2,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.phi,2,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.phi,2,i,j,k,ϕ.ls[1],ϕ.ls[2]),

        d2xD(ϕ.phi,3,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.phi,3,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.phi,3,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.phi,3,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.phi,3,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.phi,3,i,j,k,ϕ.ls[1],ϕ.ls[2]),

        d2xD(ϕ.phi,4,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.phi,4,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.phi,4,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.phi,4,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.phi,4,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.phi,4,i,j,k,ϕ.ls[1],ϕ.ls[2])
    )
end


function getDXp(ϕ,i,j,k)

    return SMatrix{3,4,Float64, 12}(

        dxDp(ϕ.phi,1,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.phi,1,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.phi,1,i,j,k,ϕ.ls[3], ϕ.index_grid_z),

        dxDp(ϕ.phi,2,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.phi,2,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.phi,2,i,j,k,ϕ.ls[3], ϕ.index_grid_z),

        dxDp(ϕ.phi,3,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.phi,3,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.phi,3,i,j,k,ϕ.ls[3], ϕ.index_grid_z),

        dxDp(ϕ.phi,4,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.phi,4,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.phi,4,i,j,k,ϕ.ls[3], ϕ.index_grid_z)
    )
    
end

function getDDXp(ϕ,i,j,k)

    return SMatrix{6,4,Float64, 24}(

        d2xDp(ϕ.phi,1,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.phi,1,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.phi,1,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.phi,1,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.phi,1,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.phi,1,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y),

        d2xDp(ϕ.phi,2,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.phi,2,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.phi,2,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.phi,2,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.phi,2,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.phi,2,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y),

        d2xDp(ϕ.phi,3,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.phi,3,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.phi,3,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.phi,3,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.phi,3,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.phi,3,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y),

        d2xDp(ϕ.phi,4,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.phi,4,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.phi,4,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.phi,4,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.phi,4,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.phi,4,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y)
    )
end



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

function dxdyD(phi, a, i, j, k, lsx, lsy)
    @inbounds 0.5*( dxdydiffD(phi, a, i, j, k, lsx, lsy) - d2xD(phi, a, i, j, k, lsx) - d2yD(phi, a, i, j, k, lsx) )
end
function dxdzD(phi, a, i, j, k, lsx, lsz)
    @inbounds 0.5*( dxdzdiffD(phi, a, i, j, k, lsx, lsz) - d2xD(phi, a, i, j, k, lsx) - d2zD(phi, a, i, j, k, lsz) )
end
function dydzD(phi, a, i, j, k, lsy, lsz)
    @inbounds 0.5*( dydzdiffD(phi, a, i, j, k, lsy, lsz) - d2yD(phi, a, i, j, k, lsy) - d2zD(phi, a, i, j, k, lsz) )
end


# periodic stuff

function dxDp(phi, a, i, j, k, lsx, ig)
    @inbounds (-phi[ig[i+4],j,k,a] + 8.0*phi[ig[i+3],j,k,a] - 8.0*phi[ig[i+1],j,k,a] + phi[ig[i],j,k,a])/(12.0*lsx)
end
function dyDp(phi, a, i, j, k, lsy, ig)
    @inbounds (-phi[i,ig[j+4],k,a] + 8.0*phi[i,ig[j+3],k,a] - 8.0*phi[i,ig[j+1],k,a] + phi[i,ig[j],k,a])/(12.0*lsy)
end
function dzDp(phi, a, i, j, k, lsz, ig)
    @inbounds (-phi[i,j,ig[k+4],a] + 8.0*phi[i,j,ig[k+3],a] - 8.0*phi[i,j,ig[k+1],a] + phi[i,j,ig[k],a])/(12.0*lsz)
end


function d2xDp(phi, a, i, j, k, lsx, ig)
    @inbounds (-phi[ig[i+4],j,k,a] + 16.0*phi[ig[i+3],j,k,a] - 30.0*phi[ig[i+2],j,k,a] + 16.0*phi[ig[i+1],j,k,a] - phi[ig[i],j,k,a])/(12.0*lsx^2)
end
function d2yDp(phi, a, i, j, k, lsx, ig)
    @inbounds (-phi[i,ig[j+4],k,a] + 16.0*phi[i,ig[j+3],k,a] - 30.0*phi[i,ig[j+2],k,a] + 16.0*phi[i,ig[j+1],k,a] - phi[i,ig[j],k,a])/(12.0*lsx^2)
end
function d2zDp(phi, a, i, j, k, lsx, ig)
    @inbounds (-phi[i,j,ig[k+4],a] + 16.0*phi[i,j,ig[k+3],a] - 30.0*phi[i,j,ig[k+2],a] + 16.0*phi[i,j,ig[k+1],a] - phi[i,j,ig[k],a])/(12.0*lsx^2)
end

function dxdydiffDp(phi, a, i, j, k, lsx, lsy, igx, igy)
    @inbounds (-phi[igx[i+4],igy[j+4],k,a] + 16.0*phi[igx[i+3],igy[j+3],k,a] - 30.0*phi[igx[i+2],igy[j+2],k,a] + 16.0*phi[igx[i+1],igy[j+1],k,a] - phi[igx[i],igy[j],k,a])/(12.0*lsx*lsy)
end
function dxdzdiffDp(phi, a, i, j, k, lsx, lsy, igx, igz)
    @inbounds (-phi[igx[i+4],j,igz[k+4],a] + 16.0*phi[igx[i+3],j,igz[k+3],a] - 30.0*phi[igx[i+2],j,igz[k+2],a] + 16.0*phi[igx[i+1],j,igz[k+1],a] - phi[igx[i],j,igz[k],a])/(12.0*lsx*lsy)
end
function dydzdiffDp(phi, a, i, j, k, lsx, lsy, igy, igz)
    @inbounds (-phi[i,igy[j+4],igz[k+4],a] + 16.0*phi[i,igy[j+3],igz[k+3],a] - 30.0*phi[i,igy[j+2],igz[k+2],a] + 16.0*phi[i,igy[j+1],igz[k+1],a] - phi[i,igy[j],igz[k],a])/(12.0*lsx*lsy)
end

function dxdyDp(phi, a, i, j, k, lsx, lsy, igx, igy)
    @inbounds 0.5*( dxdydiffDp(phi, a, i, j, k, lsx, lsy, igx, igy) - d2xDp(phi, a, i, j, k, lsx, igx) - d2yDp(phi, a, i, j, k, lsx, igy) )
end
function dxdzDp(phi, a, i, j, k, lsx, lsz, igx, igz)
    @inbounds 0.5*( dxdzdiffDp(phi, a, i, j, k, lsx, lsz, igx, igz) - d2xDp(phi, a, i, j, k, lsx, igx) - d2zDp(phi, a, i, j, k, lsz, igz) )
end
function dydzDp(phi, a, i, j, k, lsy, lsz, igy, igz)
    @inbounds 0.5*( dydzdiffDp(phi, a, i, j, k, lsy, lsz, igy, igz) - d2yDp(phi, a, i, j, k, lsy, igy) - d2zDp(phi, a, i, j, k, lsz, igz) )
end



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



