function make_RM_product!(sk, Xs)

    x = sk.x
    lp = sk.lp
    ls = sk.ls

    temp_sk = Skyrmion(lp,ls)

    a=1
    makeRM!(sk, Xs[a][1],Xs[a][2],Xs[a][3], X = Xs[a][4], iTH = Xs[a][5], i_n = Xs[a][6], jTH = Xs[a][7], j_n = Xs[a][8]  )
    
    for a in 2:size(Xs)[1]

        makeRM!(temp_sk, Xs[a][1],Xs[a][2],Xs[a][3], X = Xs[a][4], iTH = Xs[a][5], i_n = Xs[a][6], jTH = Xs[a][7], j_n = Xs[a][8]  )
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


function translate(skyrmion,X)

    x = skyrmion.x
    lp = skyrmion.lp
    ls = skyrmion.ls

    sky_temp = Skyrmion(lp,ls)

    vac = [0.0,0.0,0.0,1.0]

    ϕinterp = [ linear_interpolation((x[1],x[2],x[3]), skyrmion.phi[:,:,:,a] )  for a in 1:4 ]

    for a in 1:4, i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
        if x[1][i] > x[1][1] + X[1] && x[1][i] < x[1][end] + X[1] && x[2][j] > x[2][1] + X[2] && x[2][j] < x[2][end] + X[2] && x[3][k] > x[3][1] + X[3] && x[3][k] < x[3][end] + X[3]
            sky_temp.phi[i,j,k,a] = ϕinterp[a](x[1][i] - X[1], x[2][j] - X[2], x[3][k] - X[3])
        else
            sky_temp.phi[i,j,k,a] = vac[a]
        end
    end

    normer(sky_temp)

    return sky_temp

end

function translate!(skyrmion,X)

    x = skyrmion.x
    lp = skyrmion.lp
    ls = skyrmion.ls

    sky_temp = Skyrmion(lp,ls)

    vac = [0.0,0.0,0.0,1.0]

    ϕinterp = [ linear_interpolation((x[1],x[2],x[3]), skyrmion.phi[:,:,:,a] )  for a in 1:4 ]

    for a in 1:4, i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]
        if x[1][i] > x[1][1] + X[1] && x[1][i] < x[1][end] + X[1] && x[2][j] > x[2][1] + X[2] && x[2][j] < x[2][end] + X[2] && x[3][k] > x[3][1] + X[3] && x[3][k] < x[3][end] + X[3]
            sky_temp.phi[i,j,k,a] = ϕinterp[a](x[1][i] - X[1], x[2][j] - X[2], x[3][k] - X[3])
        else
            sky_temp.phi[i,j,k,a] = vac[a]
        end
    end

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3], a in 1:4
        skyrmion.phi[i,j,k,a] = sky_temp.phi[i,j,k,a]
    end

    normer(skyrmion)

end


function isorotate!(skyrmion,θ,n)

    rotation_matrix = R_from_axis_angle(θ,n)


    x = skyrmion.x
    lp = skyrmion.lp
    ls = skyrmion.ls

    tempsk = Skyrmion(lp,ls)

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3], a in 1:3
        tempsk.phi[i,j,k,a] = 0.0
    end

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        for a in 1:3, b in 1:3
            tempsk.phi[i,j,k,a] += rotation_matrix[a,b]*skyrmion.phi[i,j,k,b]
        end

        tempsk.phi[i,j,k,4] = skyrmion.phi[i,j,k,4]

    end

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3], a in 1:4
        skyrmion.phi[i,j,k,a] = tempsk.phi[i,j,k,a]
    end

    normer(skyrmion)

end


function isorotate(skyrmion,θ,n)

    rotation_matrix = R_from_axis_angle(θ,n)

    x = skyrmion.x
    lp = skyrmion.lp
    ls = skyrmion.ls

    tempsk = Skyrmion(lp,ls)

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3], a in 1:3
        tempsk.phi[i,j,k,a] = 0.0
    end

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        for a in 1:3, b in 1:3
            tempsk.phi[i,j,k,a] += rotation_matrix[a,b]*skyrmion.phi[i,j,k,b]
        end

        tempsk.phi[i,j,k,4] = skyrmion.phi[i,j,k,4]

    end

    normer(tempsk)

    return tempsk

end




function rotate!(skyrmion,θ,n)

    rotation_matrix = R_from_axis_angle(θ,n)

    x = skyrmion.x
    lp = skyrmion.lp
    ls = skyrmion.ls

    sky_temp = Skyrmion(lp,ls)

    vac = [0.0,0.0,0.0,1.0]

    ϕinterp = [ linear_interpolation((x[1],x[2],x[3]), skyrmion.phi[:,:,:,a] )  for a in 1:4 ]


    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        newx = rotation_matrix*[ x[1][i], x[2][j], x[3][k] ]

        if x[1][1] < newx[1] < x[1][end] && x[2][1] < newx[2] < x[2][end] && x[3][1] < newx[3] < x[3][end] 
            for a in 1:4    
                    sky_temp.phi[i,j,k,a] = ϕinterp[a](newx[1], newx[2], newx[3])
            end
        else 
            for a in 1:4    
                sky_temp.phi[i,j,k,a] = vac[a]
            end
        end
    end

    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3], a in 1:4
        skyrmion.phi[i,j,k,a] = sky_temp.phi[i,j,k,a]
    end

    normer(skyrmion)

end



function rotate(skyrmion,θ,n)

    rotation_matrix = R_from_axis_angle(θ,n)

    x = skyrmion.x
    lp = skyrmion.lp
    ls = skyrmion.ls

    sky_temp = Skyrmion(lp,ls)

    vac = [0.0,0.0,0.0,1.0]

    ϕinterp = [ linear_interpolation((x[1],x[2],x[3]), skyrmion.phi[:,:,:,a] )  for a in 1:4 ]


    for i in 1:lp[1], j in 1:lp[2], k in 1:lp[3]

        newx = rotation_matrix*[ x[1][i], x[2][j], x[3][k] ]

        if x[1][1] < newx[1] < x[1][end] && x[2][1] < newx[2] < x[2][end] && x[3][1] < newx[3] < x[3][end] 
            for a in 1:4    
                    sky_temp.phi[i,j,k,a] = ϕinterp[a](newx[1], newx[2], newx[3])
            end
        else 
            for a in 1:4    
                sky_temp.phi[i,j,k,a] = vac[a]
            end
        end
    end

    normer(sky_temp)

    return sky_temp

end