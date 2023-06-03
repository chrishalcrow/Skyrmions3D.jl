function make_RM_product!(sk,f,p,q,Xs)

    x = sk.x
    lp = sk.lp
    ls = sk.ls

    temp_sk = Skyrmion(lp,ls)

    makeRM!(sk, f,p,q, X = Xs[1] )
    
    for a in 2:size(Xs)[1]
        makeRM!(temp_sk, f,p,q, X = Xs[a] )
        product!(sk, temp_sk)
    end

end





function product!(sk1, sk2)

    if sk1.x != sk2.x
        println("Skyrmion grids are not equal. Aborting product.")
        return
    end

    x = sk1.x
    lp = sk1.lp
    ls = sk1.ls

    tempsk = Skyrmion(lp,ls)

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        tempsk.phi[i,j,k,4] = sk1.phi[i,j,k,4]*sk2.phi[i,j,k,4] 
        for a in 1:3
            tempsk.phi[i,j,k,a] = sk1.phi[i,j,k,a]*sk2.phi[i,j,k,4] + sk1.phi[i,j,k,4]*sk2.phi[i,j,k,a]
            tempsk.phi[i,j,k,4] -= sk1.phi[i,j,k,a]*sk2.phi[i,j,k,a] 
        end

    end

    normer(tempsk)

    sk1.phi .= tempsk.phi

end


function product(sk1, sk2)

    if sk1.x != sk2.x
        println("Skyrmion grids are not equal. Aborting product.")
        return
    end

    x = sk1.x
    lp = sk2.lp
    ls = sk2.ls

    tempsk = Skyrmion(lp,ls)

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        tempsk.phi[i,j,k,4] = sk1.phi[i,j,k,4]*sk2.phi[i,j,k,4] 
        for a in 1:3
            tempsk.phi[i,j,k,a] = sk1.phi[i,j,k,a]*sk2.phi[i,j,k,4] + sk1.phi[i,j,k,4]*sk2.phi[i,j,k,a]
            tempsk.phi[i,j,k,4] -= sk1.phi[i,j,k,a]*sk2.phi[i,j,k,a] 
        end

    end

    normer(tempsk)

    return tempsk

end







function shift(skyrmion,X)

    x = skyrmion.x
    lp = skyrmion.lp
    ls = skyrmion.ls

    sky_temp = Skyrmion(lp,ls)

    vac = [0.0,0.0,0.0,1.0]

    ϕinterp = [ linear_interpolation((x[1],x[2],x[3]), skyrmion.phi[:,:,:,a] )  for a in 1:4 ]

    for a in 1:4, i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
        if x[3][k] > x[3][1] + X[3] && x[3][k] < x[3][end] + X[3]
            sky_temp.phi[i,j,k,a] = ϕinterp[a](x[1][i], x[2][j], x[3][k] - X[3])
        else
            sky_temp.phi[i,j,k,a] = vac[a]
        end
    end

    normer(sky_temp)

    return sky_temp

end


function isorotate!(sky_temp,skyrmion,θ,n)

    rotation_matrix = [ cos(θ) sin(θ) 0
                        -sin(θ) cos(θ) 0
                        0 0 1 ]

    x = skyrmion.x
    lp = skyrmion.lp

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3], a in 1:3
        sky_temp.phi[i,j,k,a] = 0.0
    end

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        for a in 1:3, b in 1:3
                sky_temp.phi[i,j,k,a] += rotation_matrix[a,b]*skyrmion.phi[i,j,k,b]
        end

        sky_temp.phi[i,j,k,4] = skyrmion.phi[i,j,k,4]

    end

    #normer(sky_temp)
end