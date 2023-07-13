"""
    set_dirichlet!(skyrmion; vac = [0.0,0.0,0.0,1.0])

    Sets the boundary of `skyrmion` equal to `vac`, with default value `[0.0, 0.0, 0.0, 1.0]`

"""
function set_dirichlet!(sk; vac=[0.0,0.0,0.0,1.0])

    vac = [0.0,0.0,0.0,1.0]

    for i in 1:sk.lp[1], j in 1:sk.lp[2]
        sk.phi[i,j,1,:] .= vac
        sk.phi[i,j,2,:] .= vac
        sk.phi[i,j,sk.lp[3],:] .= vac
        sk.phi[i,j,sk.lp[3]-1,:] .= vac
    end
    
    for k in 1:sk.lp[3], j in 1:sk.lp[2]
        sk.phi[1,j,k,:] .= vac
        sk.phi[2,j,k,:] .= vac
        sk.phi[sk.lp[1],j,k,:] .= vac
        sk.phi[sk.lp[1]-1,j,k,:] .= vac
    end
    
    for k in 1:sk.lp[3], i in 1:sk.lp[1]
        sk.phi[i,1,k,:] .= vac
        sk.phi[i,2,k,:] .= vac
        sk.phi[i,sk.lp[2],k,:] .= vac
        sk.phi[i,sk.lp[2]-1,k,:] .= vac
    end
    
    

end

"""
    make_RM_product!(skyrmion, X_list) 

Makes a product approximation of many rational map skyrmions, determined through the  list `X_list`. The final field is written into `skyrmion`.

The formatting of the list is as follow:
`X_list = [ data_1, data_2, data_3, ... ]`
where
`data_1 = [ f(r), p(z), q(z), X, θiso, n_iso, θrot, n_rot ]`

See also [`product`]

# Example of list
```
p1(z) = z; q1(z) = 1; f1(r) = 4*atan(exp(-r));
p2(z) = z^2; q2(z) = 1; f2(r) = 4*atan(exp(-0.7*r));
X_list = [ [ f1, p1, q1, [0.0,0.0,1.5], 0.0, [0.0,0.0,1.0], 0.0, [0.0,0.0,1.0] ], [ f2, p2, q2, [0.0,0.0,-1.5], pi, [1.0,0.0,0.0], 0.0, [0.0,0.0,1.0] ] ]
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
    makeRM!(sk, Xs[a][1],Xs[a][2],Xs[a][3], X = Xs[a][4], iTH = Xs[a][5], i_n = Xs[a][6], jTH = Xs[a][7], j_n = Xs[a][8]  )
    
    for a in 2:size(Xs)[1]

        makeRM!(temp_sk, Xs[a][1],Xs[a][2],Xs[a][3], X = Xs[a][4], iTH = Xs[a][5], i_n = Xs[a][6], jTH = Xs[a][7], j_n = Xs[a][8]  )
        product_approx!(sk, temp_sk)
    end

end


"""
    product_approx!(skyrmion1,skyrmion2) 

Makes the symmetrised product approximation of `skyrmion1` and `skyrmion2`. The output is written in to `skyrmion1`. The returned field is normalised.

See also [`product_approx`]
"""
function product_approx!(sk1, sk2)

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

    normer!(tempsk)
    sk1.phi .= tempsk.phi

end

"""
    product_approx(skyrmion1,skyrmion2) -> product_skyrmion

Returns the symmetrised product approximation of `skyrmion1` and `skyrmion2`. The returned field is normalised.

See also [`product_approx!`]
"""
function product_approx(sk1, sk2)

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

    normer!(tempsk)

    return tempsk

end

"""
    translate_sk(skyrmion,X) -> translated_skyrmion

Returns `skyrmion` translated by 3-Vector `X`, e.g. `X = [1.0, 0.0, 0.0]`

See also [`translate_sk!`]
"""
function translate_sk(skyrmion,X)

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

    normer!(sky_temp)

    return sky_temp

end

"""
    translate_sk!(skyrmion,X)

Translates `skyrmion` by the 3-Vector `X`, e.g. `X = [1.0, 0.0, 0.0]`

See also [`translate_sk`]
"""
function translate_sk!(skyrmion,X)

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

    normer!(skyrmion)

end

"""
    isorotate_sk!(skyrmion,θ,n)

Isorotates `skyrmion` by `θ` around the vector `n`. The given vector is automatically normalised.

See also [`isorotate_sk!`]


"""
function isorotate_sk!(skyrmion,θ,n)

    if n == [0.0, 0.0, 0.0]
        println("ERROR: your vector is zero.")
        return
    end

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

    normer!(skyrmion)

end

"""
    isorotate_sk(skyrmion,θ,n) -> isorotated_skyrmion

Returns `skyrmion` isorotated by `θ` around the vector `n`. The given vector is automatically normalised.

See also [`isorotate_sk!`]
"""
function isorotate_sk(skyrmion,θ,n)

    if n == [0.0, 0.0, 0.0]
        println("ERROR: your vector is zero.")
        return
    end

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

    normer!(tempsk)

    return tempsk

end



"""
    rotate_sk!(skyrmion,θ,n)

Rotates `skyrmion` by `θ` around the vector `n`. The given vector is automatically normalised.

See also [`rotate_sk`]
"""
function rotate_sk!(skyrmion,θ,n)

    if n == [0.0, 0.0, 0.0]
        println("ERROR: your vector is zero.")
        return
    end

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

    normer!(skyrmion)

end


"""
    rotate_sk(skyrmion,θ,n) -> rotated_skyrmion

Returns `skyrmion` rotated by `θ` around the vector `n`. The given vector is automatically normalised.

See also [`rotate_sk!`]
"""
function rotate_sk(skyrmion,θ,n)

    if n == [0.0, 0.0, 0.0]
        println("ERROR: your vector is zero.")
        return
    end

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

    normer!(sky_temp)

    return sky_temp

end

"""
    center_skyrmion(my_skyrmion)

Translates `skyrmion' so that the center of mass is `(0,0,0)'.

"""
function center_skyrmion!(sk)


    for _ in 1:5
        current_CoM = center_of_mass(sk)
        translate_sk!(sk,-current_CoM)
    end
    
end




