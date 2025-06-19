# functions to get local field and their derivatives. Initialised as immutable static arrays, so that they are stored in the stack. Means we can localise, for multithreading, easily.

function getX(sk, i, j, k)

    return SVector{4,Float64}(
        sk.pion_field[i, j, k, 1],
        sk.pion_field[i, j, k, 2],
        sk.pion_field[i, j, k, 3],
        sk.pion_field[i, j, k, 4],
    )

end

function getDP(sk, i, j, k)
    if sk.dirichlet == true
        return getDX(sk, i, j, k)
    else
        return getDXp(sk, i, j, k)
    end
end

function getDX(ϕ, i, j, k)

    return SMatrix{3,4,Float64,12}(
        dxD(ϕ.pion_field, 1, i, j, k, ϕ.ls[1]),
        dyD(ϕ.pion_field, 1, i, j, k, ϕ.ls[2]),
        dzD(ϕ.pion_field, 1, i, j, k, ϕ.ls[3]),
        dxD(ϕ.pion_field, 2, i, j, k, ϕ.ls[1]),
        dyD(ϕ.pion_field, 2, i, j, k, ϕ.ls[2]),
        dzD(ϕ.pion_field, 2, i, j, k, ϕ.ls[3]),
        dxD(ϕ.pion_field, 3, i, j, k, ϕ.ls[1]),
        dyD(ϕ.pion_field, 3, i, j, k, ϕ.ls[2]),
        dzD(ϕ.pion_field, 3, i, j, k, ϕ.ls[3]),
        dxD(ϕ.pion_field, 4, i, j, k, ϕ.ls[1]),
        dyD(ϕ.pion_field, 4, i, j, k, ϕ.ls[2]),
        dzD(ϕ.pion_field, 4, i, j, k, ϕ.ls[3]),
    )

end

function getDXp(ϕ, i, j, k)

    return SMatrix{3,4,Float64,12}(
        dxDp(ϕ.pion_field, 1, i, j, k, ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field, 1, i, j, k, ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field, 1, i, j, k, ϕ.ls[3], ϕ.index_grid_z),
        dxDp(ϕ.pion_field, 2, i, j, k, ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field, 2, i, j, k, ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field, 2, i, j, k, ϕ.ls[3], ϕ.index_grid_z),
        dxDp(ϕ.pion_field, 3, i, j, k, ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field, 3, i, j, k, ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field, 3, i, j, k, ϕ.ls[3], ϕ.index_grid_z),
        dxDp(ϕ.pion_field, 4, i, j, k, ϕ.ls[1], ϕ.index_grid_x),
        dyDp(ϕ.pion_field, 4, i, j, k, ϕ.ls[2], ϕ.index_grid_y),
        dzDp(ϕ.pion_field, 4, i, j, k, ϕ.ls[3], ϕ.index_grid_z),
    )

end



# The actual derivatives

function dxD(pion_field, a, i, j, k, lsx)
    @fastmath @inbounds (
        -pion_field[i+2, j, k, a] + 8.0*pion_field[i+1, j, k, a] -
        8.0*pion_field[i-1, j, k, a] + pion_field[i-2, j, k, a]
    )/(12.0*lsx)
end
function dyD(pion_field, a, i, j, k, lsy)
    @fastmath @inbounds (
        -pion_field[i, j+2, k, a] + 8.0*pion_field[i, j+1, k, a] -
        8.0*pion_field[i, j-1, k, a] + pion_field[i, j-2, k, a]
    )/(12.0*lsy)
end
function dzD(pion_field, a, i, j, k, lsz)
    @fastmath @inbounds (
        -pion_field[i, j, k+2, a] + 8.0*pion_field[i, j, k+1, a] -
        8.0*pion_field[i, j, k-1, a] + pion_field[i, j, k-2, a]
    )/(12.0*lsz)
end


# periodic stuff

function dxDp(pion_field, a, i, j, k, lsx, ig)
    @fastmath @inbounds (
        -pion_field[ig[i+4], j, k, a] + 8.0*pion_field[ig[i+3], j, k, a] -
        8.0*pion_field[ig[i+1], j, k, a] + pion_field[ig[i], j, k, a]
    )/(12.0*lsx)
end
function dyDp(pion_field, a, i, j, k, lsy, ig)
    @fastmath @inbounds (
        -pion_field[i, ig[j+4], k, a] + 8.0*pion_field[i, ig[j+3], k, a] -
        8.0*pion_field[i, ig[j+1], k, a] + pion_field[i, ig[j], k, a]
    )/(12.0*lsy)
end
function dzDp(pion_field, a, i, j, k, lsz, ig)
    @fastmath @inbounds (
        -pion_field[i, j, ig[k+4], a] + 8.0*pion_field[i, j, ig[k+3], a] -
        8.0*pion_field[i, j, ig[k+1], a] + pion_field[i, j, ig[k], a]
    )/(12.0*lsz)
end

# Functions below used for functional gradients. Some efficient choices here, which make the code slightly awkward:
#  1. You can use information about dx^2 and dy^2 to calculate dxdy faster. So we first calc dx^2 dy^2, then dxdy.
#  2. We store each stencil seperately, to avoid vector operations. This might be unnecessary.


function getders_local_np(sk, i, j, k)

    ls = SVector{3,Float64}(sk.ls[1], sk.ls[2], sk.ls[3])

    px1 = dx_stencil(sk, i, j, k, 1)
    py1 = dy_stencil(sk, i, j, k, 1)
    pz1 = dz_stencil(sk, i, j, k, 1)

    px2 = dx_stencil(sk, i, j, k, 2)
    py2 = dy_stencil(sk, i, j, k, 2)
    pz2 = dz_stencil(sk, i, j, k, 2)

    px3 = dx_stencil(sk, i, j, k, 3)
    py3 = dy_stencil(sk, i, j, k, 3)
    pz3 = dz_stencil(sk, i, j, k, 3)

    px4 = dx_stencil(sk, i, j, k, 4)
    py4 = dy_stencil(sk, i, j, k, 4)
    pz4 = dz_stencil(sk, i, j, k, 4)

    pxy1 = dxy_stencil(sk, i, j, k, 1)
    pxz1 = dxz_stencil(sk, i, j, k, 1)
    pyz1 = dyz_stencil(sk, i, j, k, 1)

    pxy2 = dxy_stencil(sk, i, j, k, 2)
    pxz2 = dxz_stencil(sk, i, j, k, 2)
    pyz2 = dyz_stencil(sk, i, j, k, 2)

    pxy3 = dxy_stencil(sk, i, j, k, 3)
    pxz3 = dxz_stencil(sk, i, j, k, 3)
    pyz3 = dyz_stencil(sk, i, j, k, 3)

    pxy4 = dxy_stencil(sk, i, j, k, 4)
    pxz4 = dxz_stencil(sk, i, j, k, 4)
    pyz4 = dyz_stencil(sk, i, j, k, 4)

    dxV = SMatrix{3,4,Float64,12}(
        dx(px1, ls[1]),
        dx(py1, ls[2]),
        dx(pz1, ls[3]),
        dx(px2, ls[1]),
        dx(py2, ls[2]),
        dx(pz2, ls[3]),
        dx(px3, ls[1]),
        dx(py3, ls[2]),
        dx(pz3, ls[3]),
        dx(px4, ls[1]),
        dx(py4, ls[2]),
        dx(pz4, ls[3]),
    )

    d2xV1 = SMatrix{3,4,Float64,12}(
        ddx(px1, ls[1]),
        ddx(py1, ls[2]),
        ddx(pz1, ls[3]),
        ddx(px2, ls[1]),
        ddx(py2, ls[2]),
        ddx(pz2, ls[3]),
        ddx(px3, ls[1]),
        ddx(py3, ls[2]),
        ddx(pz3, ls[3]),
        ddx(px4, ls[1]),
        ddx(py4, ls[2]),
        ddx(pz4, ls[3]),
    )

    d2xV2 = SMatrix{3,4,Float64,12}(
        @fastmath 0.5*(dm(pyz1, ls[2], ls[3]) - d2xV1[2, 1] - d2xV1[3, 1]),
        0.5*(dm(pxz1, ls[1], ls[3]) - d2xV1[1, 1] - d2xV1[3, 1]),
        0.5*(dm(pxy1, ls[2], ls[1]) - d2xV1[2, 1] - d2xV1[1, 1]),
        0.5 * (dm(pyz2, ls[2], ls[3]) - d2xV1[2, 2] - d2xV1[3, 2]),
        0.5*(dm(pxz2, ls[1], ls[3]) - d2xV1[1, 2] - d2xV1[3, 2]),
        0.5*(dm(pxy2, ls[2], ls[1]) - d2xV1[2, 2] - d2xV1[1, 2]),
        0.5 * (dm(pyz3, ls[2], ls[3]) - d2xV1[2, 3] - d2xV1[3, 3]),
        0.5*(dm(pxz3, ls[1], ls[3]) - d2xV1[1, 3] - d2xV1[3, 3]),
        0.5*(dm(pxy3, ls[2], ls[1]) - d2xV1[2, 3] - d2xV1[1, 3]),
        0.5 * (dm(pyz4, ls[2], ls[3]) - d2xV1[2, 4] - d2xV1[3, 4]),
        0.5*(dm(pxz4, ls[1], ls[3]) - d2xV1[1, 4] - d2xV1[3, 4]),
        0.5*(dm(pxy4, ls[2], ls[1]) - d2xV1[2, 4] - d2xV1[1, 4])
    )

    return SVector{4,Float64}(
        sk.pion_field[i, j, k, 1],
        sk.pion_field[i, j, k, 2],
        sk.pion_field[i, j, k, 3],
        sk.pion_field[i, j, k, 4],
    ),
    dxV,
    d2xV1,
    d2xV2

end

function getders_local_p(sk, i, j, k)

    ls = SVector{3,Float64}(sk.ls[1], sk.ls[2], sk.ls[3])

    px1 = dx_stencilp(sk, i, j, k, 1, sk.index_grid_x)
    py1 = dy_stencilp(sk, i, j, k, 1, sk.index_grid_y)
    pz1 = dz_stencilp(sk, i, j, k, 1, sk.index_grid_z)

    px2 = dx_stencilp(sk, i, j, k, 2, sk.index_grid_x)
    py2 = dy_stencilp(sk, i, j, k, 2, sk.index_grid_y)
    pz2 = dz_stencilp(sk, i, j, k, 2, sk.index_grid_z)

    px3 = dx_stencilp(sk, i, j, k, 3, sk.index_grid_x)
    py3 = dy_stencilp(sk, i, j, k, 3, sk.index_grid_y)
    pz3 = dz_stencilp(sk, i, j, k, 3, sk.index_grid_z)

    px4 = dx_stencilp(sk, i, j, k, 4, sk.index_grid_x)
    py4 = dy_stencilp(sk, i, j, k, 4, sk.index_grid_y)
    pz4 = dz_stencilp(sk, i, j, k, 4, sk.index_grid_z)

    pxy1 = dxy_stencilp(sk, i, j, k, 1, sk.index_grid_x, sk.index_grid_y)
    pxz1 = dxz_stencilp(sk, i, j, k, 1, sk.index_grid_x, sk.index_grid_z)
    pyz1 = dyz_stencilp(sk, i, j, k, 1, sk.index_grid_y, sk.index_grid_z)

    pxy2 = dxy_stencilp(sk, i, j, k, 2, sk.index_grid_x, sk.index_grid_y)
    pxz2 = dxz_stencilp(sk, i, j, k, 2, sk.index_grid_x, sk.index_grid_z)
    pyz2 = dyz_stencilp(sk, i, j, k, 2, sk.index_grid_y, sk.index_grid_z)

    pxy3 = dxy_stencilp(sk, i, j, k, 3, sk.index_grid_x, sk.index_grid_y)
    pxz3 = dxz_stencilp(sk, i, j, k, 3, sk.index_grid_x, sk.index_grid_z)
    pyz3 = dyz_stencilp(sk, i, j, k, 3, sk.index_grid_y, sk.index_grid_z)

    pxy4 = dxy_stencilp(sk, i, j, k, 4, sk.index_grid_x, sk.index_grid_y)
    pxz4 = dxz_stencilp(sk, i, j, k, 4, sk.index_grid_x, sk.index_grid_z)
    pyz4 = dyz_stencilp(sk, i, j, k, 4, sk.index_grid_y, sk.index_grid_z)


    dxV = SMatrix{3,4,Float64,12}(
        dx(px1, ls[1]),
        dx(py1, ls[2]),
        dx(pz1, ls[3]),
        dx(px2, ls[1]),
        dx(py2, ls[2]),
        dx(pz2, ls[3]),
        dx(px3, ls[1]),
        dx(py3, ls[2]),
        dx(pz3, ls[3]),
        dx(px4, ls[1]),
        dx(py4, ls[2]),
        dx(pz4, ls[3]),
    )

    d2xV1 = SMatrix{3,4,Float64,12}(
        ddx(px1, ls[1]),
        ddx(py1, ls[2]),
        ddx(pz1, ls[3]),
        ddx(px2, ls[1]),
        ddx(py2, ls[2]),
        ddx(pz2, ls[3]),
        ddx(px3, ls[1]),
        ddx(py3, ls[2]),
        ddx(pz3, ls[3]),
        ddx(px4, ls[1]),
        ddx(py4, ls[2]),
        ddx(pz4, ls[3]),
    )

    d2xV2 = @fastmath SMatrix{3,4,Float64,12}(
        0.5*(dm(pyz1, ls[2], ls[3]) - d2xV1[2, 1] - d2xV1[3, 1]),
        0.5*(dm(pxz1, ls[1], ls[3]) - d2xV1[1, 1] - d2xV1[3, 1]),
        0.5*(dm(pxy1, ls[2], ls[1]) - d2xV1[2, 1] - d2xV1[1, 1]),
        0.5 * (dm(pyz2, ls[2], ls[3]) - d2xV1[2, 2] - d2xV1[3, 2]),
        0.5*(dm(pxz2, ls[1], ls[3]) - d2xV1[1, 2] - d2xV1[3, 2]),
        0.5*(dm(pxy2, ls[2], ls[1]) - d2xV1[2, 2] - d2xV1[1, 2]),
        0.5 * (dm(pyz3, ls[2], ls[3]) - d2xV1[2, 3] - d2xV1[3, 3]),
        0.5*(dm(pxz3, ls[1], ls[3]) - d2xV1[1, 3] - d2xV1[3, 3]),
        0.5*(dm(pxy3, ls[2], ls[1]) - d2xV1[2, 3] - d2xV1[1, 3]),
        0.5 * (dm(pyz4, ls[2], ls[3]) - d2xV1[2, 4] - d2xV1[3, 4]),
        0.5*(dm(pxz4, ls[1], ls[3]) - d2xV1[1, 4] - d2xV1[3, 4]),
        0.5*(dm(pxy4, ls[2], ls[1]) - d2xV1[2, 4] - d2xV1[1, 4]),
    )

    return SVector{4,Float64}(
        sk.pion_field[i, j, k, 1],
        sk.pion_field[i, j, k, 2],
        sk.pion_field[i, j, k, 3],
        sk.pion_field[i, j, k, 4],
    ),
    dxV,
    d2xV1,
    d2xV2

end




function dx(pion_field, lsx)
    @fastmath @inbounds (
        -pion_field[5] + 8.0*pion_field[4] - 8.0*pion_field[2] + pion_field[1]
    )/(12.0*lsx)
end

function ddx(pion_field, lsx)
    @fastmath @inbounds (
        -pion_field[5] + 16.0*pion_field[4] - 30.0*pion_field[3] + 16.0*pion_field[2] -
        pion_field[1]
    )/(12.0*lsx^2)
end

function dm(pion_field, lsx, lsy)
    @fastmath @inbounds (
        -pion_field[5] + 16.0*pion_field[4] - 30.0*pion_field[3] + 16.0*pion_field[2] -
        pion_field[1]
    )/(12.0*lsx*lsy)
end

function dx_stencil(sk, i, j, k, a)
    return SVector{5,Float64}(
        sk.pion_field[i-2, j, k, a],
        sk.pion_field[i-1, j, k, a],
        sk.pion_field[i, j, k, a],
        sk.pion_field[i+1, j, k, a],
        sk.pion_field[i+2, j, k, a],
    )
end
function dy_stencil(sk, i, j, k, a)
    return SVector{5,Float64}(
        sk.pion_field[i, j-2, k, a],
        sk.pion_field[i, j-1, k, a],
        sk.pion_field[i, j, k, a],
        sk.pion_field[i, j+1, k, a],
        sk.pion_field[i, j+2, k, a],
    )
end
function dz_stencil(sk, i, j, k, a)
    return SVector{5,Float64}(
        sk.pion_field[i, j, k-2, a],
        sk.pion_field[i, j, k-1, a],
        sk.pion_field[i, j, k, a],
        sk.pion_field[i, j, k+1, a],
        sk.pion_field[i, j, k+2, a],
    )
end
function dxy_stencil(sk, i, j, k, a)
    return SVector{5,Float64}(
        sk.pion_field[i-2, j-2, k, a],
        sk.pion_field[i-1, j-1, k, a],
        sk.pion_field[i, j, k, a],
        sk.pion_field[i+1, j+1, k, a],
        sk.pion_field[i+2, j+2, k, a],
    )
end
function dxz_stencil(sk, i, j, k, a)
    return SVector{5,Float64}(
        sk.pion_field[i-2, j, k-2, a],
        sk.pion_field[i-1, j, k-1, a],
        sk.pion_field[i, j, k, a],
        sk.pion_field[i+1, j, k+1, a],
        sk.pion_field[i+2, j, k+2, a],
    )
end
function dyz_stencil(sk, i, j, k, a)
    return SVector{5,Float64}(
        sk.pion_field[i, j-2, k-2, a],
        sk.pion_field[i, j-1, k-1, a],
        sk.pion_field[i, j, k, a],
        sk.pion_field[i, j+1, k+1, a],
        sk.pion_field[i, j+2, k+2, a],
    )
end



function dx_stencilp(sk, i, j, k, a, ig)
    return SVector{5,Float64}(
        sk.pion_field[ig[i], j, k, a],
        sk.pion_field[ig[i+1], j, k, a],
        sk.pion_field[ig[i+2], j, k, a],
        sk.pion_field[ig[i+3], j, k, a],
        sk.pion_field[ig[i+4], j, k, a],
    )
end
function dy_stencilp(sk, i, j, k, a, ig)
    return SVector{5,Float64}(
        sk.pion_field[i, ig[j], k, a],
        sk.pion_field[i, ig[j+1], k, a],
        sk.pion_field[i, ig[j+2], k, a],
        sk.pion_field[i, ig[j+3], k, a],
        sk.pion_field[i, ig[j+4], k, a],
    )
end
function dz_stencilp(sk, i, j, k, a, ig)
    return SVector{5,Float64}(
        sk.pion_field[i, j, ig[k], a],
        sk.pion_field[i, j, ig[k+1], a],
        sk.pion_field[i, j, ig[k+2], a],
        sk.pion_field[i, j, ig[k+3], a],
        sk.pion_field[i, j, ig[k+4], a],
    )
end
function dxy_stencilp(sk, i, j, k, a, ig1, ig2)
    return SVector{5,Float64}(
        sk.pion_field[ig1[i], ig2[j], k, a],
        sk.pion_field[ig1[i+1], ig2[j+1], k, a],
        sk.pion_field[ig1[i+2], ig2[j+2], k, a],
        sk.pion_field[ig1[i+3], ig2[j+3], k, a],
        sk.pion_field[ig1[i+4], ig2[j+4], k, a],
    )
end
function dxz_stencilp(sk, i, j, k, a, ig1, ig2)
    return SVector{5,Float64}(
        sk.pion_field[ig1[i], j, ig2[k], a],
        sk.pion_field[ig1[i+1], j, ig2[k+1], a],
        sk.pion_field[ig1[i+2], j, ig2[k+2], a],
        sk.pion_field[ig1[i+3], j, ig2[k+3], a],
        sk.pion_field[ig1[i+4], j, ig2[k+4], a],
    )
end
function dyz_stencilp(sk, i, j, k, a, ig1, ig2)
    return SVector{5,Float64}(
        sk.pion_field[i, ig1[j], ig2[k], a],
        sk.pion_field[i, ig1[j+1], ig2[k+1], a],
        sk.pion_field[i, ig1[j+2], ig2[k+2], a],
        sk.pion_field[i, ig1[j+3], ig2[k+3], a],
        sk.pion_field[i, ig1[j+4], ig2[k+4], a],
    )
end
