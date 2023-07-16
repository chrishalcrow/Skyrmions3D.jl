# functions to get local field and their derivatives. Initialised as immutable static arrays, so that they are stored in the stack. Means we can localise, for multithreading, easily.
function getX(sk,i,j,k)

    return SVector{4,Float64}(
        sk.pion_field[i,j,k,1],
        sk.pion_field[i,j,k,2],
        sk.pion_field[i,j,k,3],
        sk.pion_field[i,j,k,4]
    )
    
end

function getDX(ϕ,i,j,k)

    return SMatrix{3,4,Float64, 12}(
        dxD(ϕ.pion_field,1,i,j,k,ϕ.ls[1]),
        dyD(ϕ.pion_field,1,i,j,k,ϕ.ls[2]),
        dzD(ϕ.pion_field,1,i,j,k,ϕ.ls[3]),

        dxD(ϕ.pion_field,2,i,j,k,ϕ.ls[1]),
        dyD(ϕ.pion_field,2,i,j,k,ϕ.ls[2]),
        dzD(ϕ.pion_field,2,i,j,k,ϕ.ls[3]),

        dxD(ϕ.pion_field,3,i,j,k,ϕ.ls[1]),
        dyD(ϕ.pion_field,3,i,j,k,ϕ.ls[2]),
        dzD(ϕ.pion_field,3,i,j,k,ϕ.ls[3]),

        dxD(ϕ.pion_field,4,i,j,k,ϕ.ls[1]),
        dyD(ϕ.pion_field,4,i,j,k,ϕ.ls[2]),
        dzD(ϕ.pion_field,4,i,j,k,ϕ.ls[3])
    )
    

end

function getDDX(ϕ,i,j,k)

    return SMatrix{6,4,Float64, 24}(

        d2xD(ϕ.pion_field,1,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.pion_field,1,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.pion_field,1,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.pion_field,1,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.pion_field,1,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.pion_field,1,i,j,k,ϕ.ls[1],ϕ.ls[2]),

        d2xD(ϕ.pion_field,2,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.pion_field,2,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.pion_field,2,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.pion_field,2,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.pion_field,2,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.pion_field,2,i,j,k,ϕ.ls[1],ϕ.ls[2]),

        d2xD(ϕ.pion_field,3,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.pion_field,3,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.pion_field,3,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.pion_field,3,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.pion_field,3,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.pion_field,3,i,j,k,ϕ.ls[1],ϕ.ls[2]),

        d2xD(ϕ.pion_field,4,i,j,k,ϕ.ls[1]),
        d2yD(ϕ.pion_field,4,i,j,k,ϕ.ls[2]),
        d2zD(ϕ.pion_field,4,i,j,k,ϕ.ls[3]),
        dydzD(ϕ.pion_field,4,i,j,k,ϕ.ls[2],ϕ.ls[3]),
        dxdzD(ϕ.pion_field,4,i,j,k,ϕ.ls[1],ϕ.ls[3]),
        dxdyD(ϕ.pion_field,4,i,j,k,ϕ.ls[1],ϕ.ls[2])
    )
end


function getDXp(ϕ,i,j,k)

    return SMatrix{3,4,Float64, 12}(

        dxDp(ϕ.pion_field,1,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field,1,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field,1,i,j,k,ϕ.ls[3], ϕ.index_grid_z),

        dxDp(ϕ.pion_field,2,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field,2,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field,2,i,j,k,ϕ.ls[3], ϕ.index_grid_z),

        dxDp(ϕ.pion_field,3,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field,3,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field,3,i,j,k,ϕ.ls[3], ϕ.index_grid_z),

        dxDp(ϕ.pion_field,4,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field,4,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field,4,i,j,k,ϕ.ls[3], ϕ.index_grid_z)
    )
    
end

function getDDXp(ϕ,i,j,k)

    return SMatrix{6,4,Float64, 24}(

        d2xDp(ϕ.pion_field,1,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.pion_field,1,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.pion_field,1,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.pion_field,1,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.pion_field,1,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.pion_field,1,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y),

        d2xDp(ϕ.pion_field,2,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.pion_field,2,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.pion_field,2,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.pion_field,2,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.pion_field,2,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.pion_field,2,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y),

        d2xDp(ϕ.pion_field,3,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.pion_field,3,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.pion_field,3,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.pion_field,3,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.pion_field,3,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.pion_field,3,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y),

        d2xDp(ϕ.pion_field,4,i,j,k,ϕ.ls[1], ϕ.index_grid_x),
        d2yDp(ϕ.pion_field,4,i,j,k,ϕ.ls[2], ϕ.index_grid_y),
        d2zDp(ϕ.pion_field,4,i,j,k,ϕ.ls[3], ϕ.index_grid_z),
        dydzDp(ϕ.pion_field,4,i,j,k,ϕ.ls[2],ϕ.ls[3], ϕ.index_grid_y, ϕ.index_grid_z),
        dxdzDp(ϕ.pion_field,4,i,j,k,ϕ.ls[1],ϕ.ls[3], ϕ.index_grid_x, ϕ.index_grid_z),
        dxdyDp(ϕ.pion_field,4,i,j,k,ϕ.ls[1],ϕ.ls[2], ϕ.index_grid_x, ϕ.index_grid_y)
    )
end



# The same functions, but going directly to the pion field ()
function getX(pion_field::Array{Float64, 4},i,j,k)

    return SVector{4,Float64}(
        pion_field[i,j,k,1],
        pion_field[i,j,k,2],
        pion_field[i,j,k,3],
        pion_field[i,j,k,4]
    )
    
end

function getDX(pion_field::Array{Float64, 4},i,j,k,ls)

    return SMatrix{3,4,Float64, 12}(
        dxD(pion_field,1,i,j,k,ls[1]),
        dyD(pion_field,1,i,j,k,ls[2]),
        dzD(pion_field,1,i,j,k,ls[3]),

        dxD(pion_field,2,i,j,k,ls[1]),
        dyD(pion_field,2,i,j,k,ls[2]),
        dzD(pion_field,2,i,j,k,ls[3]),

        dxD(pion_field,3,i,j,k,ls[1]),
        dyD(pion_field,3,i,j,k,ls[2]),
        dzD(pion_field,3,i,j,k,ls[3]),

        dxD(pion_field,4,i,j,k,ls[1]),
        dyD(pion_field,4,i,j,k,ls[2]),
        dzD(pion_field,4,i,j,k,ls[3])
    )
    

end

function getDDX(pion_field::Array{Float64, 4},i,j,k,ls)

    return SMatrix{6,4,Float64, 24}(

        d2xD(pion_field,1,i,j,k,ls[1]),
        d2yD(pion_field,1,i,j,k,ls[2]),
        d2zD(pion_field,1,i,j,k,ls[3]),
        dydzD(pion_field,1,i,j,k,ls[2],ls[3]),
        dxdzD(pion_field,1,i,j,k,ls[1],ls[3]),
        dxdyD(pion_field,1,i,j,k,ls[1],ls[2]),

        d2xD(pion_field,2,i,j,k,ls[1]),
        d2yD(pion_field,2,i,j,k,ls[2]),
        d2zD(pion_field,2,i,j,k,ls[3]),
        dydzD(pion_field,2,i,j,k,ls[2],ls[3]),
        dxdzD(pion_field,2,i,j,k,ls[1],ls[3]),
        dxdyD(pion_field,2,i,j,k,ls[1],ls[2]),

        d2xD(pion_field,3,i,j,k,ls[1]),
        d2yD(pion_field,3,i,j,k,ls[2]),
        d2zD(pion_field,3,i,j,k,ls[3]),
        dydzD(pion_field,3,i,j,k,ls[2],ls[3]),
        dxdzD(pion_field,3,i,j,k,ls[1],ls[3]),
        dxdyD(pion_field,3,i,j,k,ls[1],ls[2]),

        d2xD(pion_field,4,i,j,k,ls[1]),
        d2yD(pion_field,4,i,j,k,ls[2]),
        d2zD(pion_field,4,i,j,k,ls[3]),
        dydzD(pion_field,4,i,j,k,ls[2],ls[3]),
        dxdzD(pion_field,4,i,j,k,ls[1],ls[3]),
        dxdyD(pion_field,4,i,j,k,ls[1],ls[2])
    )
end

# The actual derivatives

function dxD(pion_field, a, i, j, k, lsx)
    @inbounds (-pion_field[i+2,j,k,a] + 8.0*pion_field[i+1,j,k,a] - 8.0*pion_field[i-1,j,k,a] + pion_field[i-2,j,k,a])/(12.0*lsx)
end
function dyD(pion_field, a, i, j, k, lsy)
    @inbounds (-pion_field[i,j+2,k,a] + 8.0*pion_field[i,j+1,k,a] - 8.0*pion_field[i,j-1,k,a] + pion_field[i,j-2,k,a])/(12.0*lsy)
end
function dzD(pion_field, a, i, j, k, lsz)
    @inbounds (-pion_field[i,j,k+2,a] + 8.0*pion_field[i,j,k+1,a] - 8.0*pion_field[i,j,k-1,a] + pion_field[i,j,k-2,a])/(12.0*lsz)
end

function d2xD(pion_field, a, i, j, k, lsx)
    @inbounds (-pion_field[i+2,j,k,a] + 16.0*pion_field[i+1,j,k,a] - 30.0*pion_field[i,j,k,a] + 16.0*pion_field[i-1,j,k,a] - pion_field[i-2,j,k,a])/(12.0*lsx^2)
end
function d2yD(pion_field, a, i, j, k, lsx)
    @inbounds (-pion_field[i,j+2,k,a] + 16.0*pion_field[i,j+1,k,a] - 30.0*pion_field[i,j,k,a] + 16.0*pion_field[i,j-1,k,a] - pion_field[i,j-2,k,a])/(12.0*lsx^2)
end
function d2zD(pion_field, a, i, j, k, lsx)
    @inbounds (-pion_field[i,j,k+2,a] + 16.0*pion_field[i,j,k+1,a] - 30.0*pion_field[i,j,k,a] + 16.0*pion_field[i,j,k-1,a] - pion_field[i,j,k-2,a])/(12.0*lsx^2)
end

function dxdydiffD(pion_field, a, i, j, k, lsx, lsy)
    @inbounds (-pion_field[i+2,j+2,k,a] + 16.0*pion_field[i+1,j+1,k,a] - 30.0*pion_field[i,j,k,a] + 16.0*pion_field[i-1,j-1,k,a] - pion_field[i-2,j-2,k,a])/(12.0*lsx*lsy)
end
function dxdzdiffD(pion_field, a, i, j, k, lsx, lsy)
    @inbounds (-pion_field[i+2,j,k+2,a] + 16.0*pion_field[i+1,j,k+1,a] - 30.0*pion_field[i,j,k,a] + 16.0*pion_field[i-1,j,k-1,a] - pion_field[i-2,j,k-2,a])/(12.0*lsx*lsy)
end
function dydzdiffD(pion_field, a, i, j, k, lsx, lsy)
    @inbounds (-pion_field[i,j+2,k+2,a] + 16.0*pion_field[i,j+1,k+1,a] - 30.0*pion_field[i,j,k,a] + 16.0*pion_field[i,j-1,k-1,a] - pion_field[i,j-2,k-2,a])/(12.0*lsx*lsy)
end

function dxdyD(pion_field, a, i, j, k, lsx, lsy)
    @inbounds 0.5*( dxdydiffD(pion_field, a, i, j, k, lsx, lsy) - d2xD(pion_field, a, i, j, k, lsx) - d2yD(pion_field, a, i, j, k, lsx) )
end
function dxdzD(pion_field, a, i, j, k, lsx, lsz)
    @inbounds 0.5*( dxdzdiffD(pion_field, a, i, j, k, lsx, lsz) - d2xD(pion_field, a, i, j, k, lsx) - d2zD(pion_field, a, i, j, k, lsz) )
end
function dydzD(pion_field, a, i, j, k, lsy, lsz)
    @inbounds 0.5*( dydzdiffD(pion_field, a, i, j, k, lsy, lsz) - d2yD(pion_field, a, i, j, k, lsy) - d2zD(pion_field, a, i, j, k, lsz) )
end


# periodic stuff

function dxDp(pion_field, a, i, j, k, lsx, ig)
    @inbounds (-pion_field[ig[i+4],j,k,a] + 8.0*pion_field[ig[i+3],j,k,a] - 8.0*pion_field[ig[i+1],j,k,a] + pion_field[ig[i],j,k,a])/(12.0*lsx)
end
function dyDp(pion_field, a, i, j, k, lsy, ig)
    @inbounds (-pion_field[i,ig[j+4],k,a] + 8.0*pion_field[i,ig[j+3],k,a] - 8.0*pion_field[i,ig[j+1],k,a] + pion_field[i,ig[j],k,a])/(12.0*lsy)
end
function dzDp(pion_field, a, i, j, k, lsz, ig)
    @inbounds (-pion_field[i,j,ig[k+4],a] + 8.0*pion_field[i,j,ig[k+3],a] - 8.0*pion_field[i,j,ig[k+1],a] + pion_field[i,j,ig[k],a])/(12.0*lsz)
end


function d2xDp(pion_field, a, i, j, k, lsx, ig)
    @inbounds (-pion_field[ig[i+4],j,k,a] + 16.0*pion_field[ig[i+3],j,k,a] - 30.0*pion_field[ig[i+2],j,k,a] + 16.0*pion_field[ig[i+1],j,k,a] - pion_field[ig[i],j,k,a])/(12.0*lsx^2)
end
function d2yDp(pion_field, a, i, j, k, lsx, ig)
    @inbounds (-pion_field[i,ig[j+4],k,a] + 16.0*pion_field[i,ig[j+3],k,a] - 30.0*pion_field[i,ig[j+2],k,a] + 16.0*pion_field[i,ig[j+1],k,a] - pion_field[i,ig[j],k,a])/(12.0*lsx^2)
end
function d2zDp(pion_field, a, i, j, k, lsx, ig)
    @inbounds (-pion_field[i,j,ig[k+4],a] + 16.0*pion_field[i,j,ig[k+3],a] - 30.0*pion_field[i,j,ig[k+2],a] + 16.0*pion_field[i,j,ig[k+1],a] - pion_field[i,j,ig[k],a])/(12.0*lsx^2)
end

function dxdydiffDp(pion_field, a, i, j, k, lsx, lsy, igx, igy)
    @inbounds (-pion_field[igx[i+4],igy[j+4],k,a] + 16.0*pion_field[igx[i+3],igy[j+3],k,a] - 30.0*pion_field[igx[i+2],igy[j+2],k,a] + 16.0*pion_field[igx[i+1],igy[j+1],k,a] - pion_field[igx[i],igy[j],k,a])/(12.0*lsx*lsy)
end
function dxdzdiffDp(pion_field, a, i, j, k, lsx, lsy, igx, igz)
    @inbounds (-pion_field[igx[i+4],j,igz[k+4],a] + 16.0*pion_field[igx[i+3],j,igz[k+3],a] - 30.0*pion_field[igx[i+2],j,igz[k+2],a] + 16.0*pion_field[igx[i+1],j,igz[k+1],a] - pion_field[igx[i],j,igz[k],a])/(12.0*lsx*lsy)
end
function dydzdiffDp(pion_field, a, i, j, k, lsx, lsy, igy, igz)
    @inbounds (-pion_field[i,igy[j+4],igz[k+4],a] + 16.0*pion_field[i,igy[j+3],igz[k+3],a] - 30.0*pion_field[i,igy[j+2],igz[k+2],a] + 16.0*pion_field[i,igy[j+1],igz[k+1],a] - pion_field[i,igy[j],igz[k],a])/(12.0*lsx*lsy)
end

function dxdyDp(pion_field, a, i, j, k, lsx, lsy, igx, igy)
    @inbounds 0.5*( dxdydiffDp(pion_field, a, i, j, k, lsx, lsy, igx, igy) - d2xDp(pion_field, a, i, j, k, lsx, igx) - d2yDp(pion_field, a, i, j, k, lsx, igy) )
end
function dxdzDp(pion_field, a, i, j, k, lsx, lsz, igx, igz)
    @inbounds 0.5*( dxdzdiffDp(pion_field, a, i, j, k, lsx, lsz, igx, igz) - d2xDp(pion_field, a, i, j, k, lsx, igx) - d2zDp(pion_field, a, i, j, k, lsz, igz) )
end
function dydzDp(pion_field, a, i, j, k, lsy, lsz, igy, igz)
    @inbounds 0.5*( dydzdiffDp(pion_field, a, i, j, k, lsy, lsz, igy, igz) - d2yDp(pion_field, a, i, j, k, lsy, igy) - d2zDp(pion_field, a, i, j, k, lsz, igz) )
end