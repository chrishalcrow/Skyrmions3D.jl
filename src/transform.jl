"""
    product_approx!(skyrmion1, skyrmion2) 

Makes the symmetrised product approximation of `skyrmion1` and `skyrmion2`. The output is written in to `skyrmion1`. The returned field is normalised.

See also [`product_approx`](@ref). 

"""
function product_approx!(sk1, sk2)
    # Assign via a temporary skyrmion made via product_approx. 
    tempsk = product_approx(sk1, sk2)
    # note here the .= is broadcast assignment 
    sk1.pion_field .= tempsk.pion_field

end

"""
    product_approx(skyrmion1, skyrmion2)

Returns the symmetrised product approximation of `skyrmion1` and `skyrmion2`. The returned field is normalised.

See also [`product_approx!`](@ref).

"""
function product_approx(sk1, sk2)

    check_grids(sk1, sk2)

    lp = sk2.grid.lp

    tempsk = deepcopy(sk1)

    for k = 1:lp[3], j = 1:lp[2], i = 1:lp[1]
        tempsk.pion_field[i, j, k, :] = product_approx_pt(sk1, sk2, i, j, k)
    end

    normer!(tempsk)

    return tempsk

end

function check_grids(sk1, sk2)
    if sk1.grid.x != sk2.grid.x
        error("skyrmion grids are not equal")
        return
    end
end

function product_approx_pt(sk1, sk2, i, j, k)

    temp_pt = MVector{4,Float64}(0.0, 0.0, 0.0, 0.0)

    temp_pt[4] = sk1.pion_field[i, j, k, 4]*sk2.pion_field[i, j, k, 4]
    for a = 1:3
        temp_pt[a] =
            sk1.pion_field[i, j, k, a]*sk2.pion_field[i, j, k, 4] +
            sk1.pion_field[i, j, k, 4]*sk2.pion_field[i, j, k, a]
        temp_pt[4] -= sk1.pion_field[i, j, k, a]*sk2.pion_field[i, j, k, a]
    end

    return temp_pt

end



"""
    translate_sk(skyrmion; X = [0.0, 0.0, 0.0])

Returns `skyrmion` translated by 3-Vector `X`, e.g. `X = [1.0, 0.0, 0.0]`

See also [`translate_sk!`](@ref). 

"""
function translate_sk(skyrmion; X = [0.0, 0.0, 0.0])

    x = skyrmion.grid.x
    lp = skyrmion.grid.lp

    sky_temp = deepcopy(skyrmion)

    vac = [0.0, 0.0, 0.0, 1.0]

    ϕinterp = quadratic_spline_interpolation(skyrmion.pion_field, x)

    for a = 1:4, k = 1:lp[3], j = 1:lp[2], i = 1:lp[1]
        if x[1][i] > x[1][1] + X[1] &&
           x[1][i] < x[1][end] + X[1] &&
           x[2][j] > x[2][1] + X[2] &&
           x[2][j] < x[2][end] + X[2] &&
           x[3][k] > x[3][1] + X[3] &&
           x[3][k] < x[3][end] + X[3]
            sky_temp.pion_field[i, j, k, a] =
                ϕinterp[a](x[1][i] - X[1], x[2][j] - X[2], x[3][k] - X[3])
        else
            sky_temp.pion_field[i, j, k, a] = vac[a]
        end
    end

    normer!(sky_temp)

    return sky_temp

end

"""
    translate_sk!(skyrmion; X = [0.0, 0.0, 0.0])

Translates `skyrmion` by the 3-Vector `X`, e.g. `X = [1.0, 0.0, 0.0]`

See also [`translate_sk`](@ref). 

"""
function translate_sk!(skyrmion; X = [0.0, 0.0, 0.0])

    tempsk = translate_sk(skyrmion, X=X)
    skyrmion.pion_field .= tempsk.pion_field

end

"""
    isorotate_sk!(skyrmion; theta = 0, n = [0, 0, 1])

Isorotates `skyrmion` by `theta` around the vector `n`. The given vector is automatically normalised.

See also [`isorotate_sk`](@ref). 

"""
function isorotate_sk!(skyrmion; theta = 0, n = [0, 0, 1])

    tempsk = isorotate_sk(skyrmion, theta=theta, n=n)
    skyrmion.pion_field .= tempsk.pion_field

end

"""
    isorotate_sk(skyrmion; theta = 0, n = [0, 0, 1])

Returns `skyrmion` isorotated by `theta` around the vector `n`. The given vector is automatically normalised.

See also [`isorotate_sk!`](@ref). 

"""
function isorotate_sk(skyrmion; theta = 0, n = [0, 0, 1])

    if n == [0.0, 0.0, 0.0]
        error("Your vector is zero.")
        return
    end

    rotation_matrix = R_from_axis_angle(theta, n)

    tempsk = deepcopy(skyrmion)

    lp = skyrmion.grid.lp
    for i = 1:lp[1], j = 1:lp[2], k = 1:lp[3], a = 1:3
        tempsk.pion_field[i, j, k, a] = 0.0
    end

    for i = 1:lp[1], j = 1:lp[2], k = 1:lp[3]

        for a = 1:3, b = 1:3
            tempsk.pion_field[i, j, k, a] +=
                rotation_matrix[a, b]*skyrmion.pion_field[i, j, k, b]
        end

        tempsk.pion_field[i, j, k, 4] = skyrmion.pion_field[i, j, k, 4]

    end

    normer!(tempsk)

    return tempsk

end



"""
    rotate_sk!(skyrmion; theta = 0, n = [0, 0, 1])

Rotates `skyrmion` by `theta` around the vector `n`. The given vector is automatically normalised.

See also [`rotate_sk`](@ref). 

"""
function rotate_sk!(skyrmion; theta = 0, n = [0, 0, 1])

    tempsk = rotate_sk(skyrmion, theta=theta, n=n)
    skyrmion.pion_field .= tempsk.pion_field

end


"""
    rotate_sk(skyrmion; theta = 0, n = [0, 0, 1])

Returns `skyrmion` rotated by `theta` around the vector `n`. The given vector is automatically normalised.

See also [`rotate_sk!`](@ref). 

"""
function rotate_sk(skyrmion; theta = 0, n = [0, 0, 1])

    if n == [0.0, 0.0, 0.0]
        error("Your vector is zero.")
        return
    end

    rotation_matrix = R_from_axis_angle(theta, n)

    lp = skyrmion.grid.lp
    x = skyrmion.grid.x

    sky_temp = deepcopy(skyrmion)

    vac = [0.0, 0.0, 0.0, 1.0]

    ϕinterp = quadratic_spline_interpolation(skyrmion.pion_field, x)

    for i = 1:lp[1], j = 1:lp[2], k = 1:lp[3]

        newx = rotation_matrix*[x[1][i], x[2][j], x[3][k]]

        if x[1][1] < newx[1] < x[1][end] &&
           x[2][1] < newx[2] < x[2][end] &&
           x[3][1] < newx[3] < x[3][end]
            for a = 1:4
                sky_temp.pion_field[i, j, k, a] = ϕinterp[a](newx[1], newx[2], newx[3])
            end
        else
            for a = 1:4
                sky_temp.pion_field[i, j, k, a] = vac[a]
            end
        end
    end

    normer!(sky_temp)

    return sky_temp

end

"""
    center_skyrmion!(skyrmion; steps = 10, tolerance = 1e-9)

Translates `skyrmion` so that the center of mass is `(0, 0, 0)`.

The method works by succesively finding the center of mass of `skyrmion` and translating by it, either until the L1 difference between the center of mass and `(0, 0, 0)` is less than `tolerance`, or until `steps` many translations have occurred. This process should converge provided the Skyrme field of `skyrmion` is small at the boundary of the grid. 

"""
function center_skyrmion!(sk; steps = 10, tolerance = 1e-9)
    # A single translation by the center of mass may not leave the new center 
    # sufficiently close to the origin. By using repeat translations this 
    # method is more robust, but succeptible to the same problems which haunt
    # center_of_mass. 
    tolerance < 0 && @warn "Tolerance is negative"

    current_CoM = center_of_mass(sk)
    counter = 0
    L1 = tolerance + 1 # dummy value so we enter the while loop
    while counter < steps && L1 >= tolerance 
        translate_sk!(sk, X = -current_CoM)
        current_CoM = center_of_mass(sk)
        L1 = maximum(abs, current_CoM)
        counter += 1
    end
    L1 >= tolerance && @warn "Centering failed to converge"

end

function set_dirichlet_boudary!(sk; vac = [0.0, 0.0, 0.0, 1.0])

    for i = 1:sk.grid.lp[1], j = 1:sk.grid.lp[2]
        sk.pion_field[i, j, 1, :] .= vac
        sk.pion_field[i, j, 2, :] .= vac
        sk.pion_field[i, j, sk.grid.lp[3], :] .= vac
        sk.pion_field[i, j, sk.grid.lp[3]-1, :] .= vac
    end

    for k = 1:sk.grid.lp[3], j = 1:sk.grid.lp[2]
        sk.pion_field[1, j, k, :] .= vac
        sk.pion_field[2, j, k, :] .= vac
        sk.pion_field[sk.grid.lp[1], j, k, :] .= vac
        sk.pion_field[sk.grid.lp[1]-1, j, k, :] .= vac
    end

    for k = 1:sk.grid.lp[3], i = 1:sk.grid.lp[1]
        sk.pion_field[i, 1, k, :] .= vac
        sk.pion_field[i, 2, k, :] .= vac
        sk.pion_field[i, sk.grid.lp[2], k, :] .= vac
        sk.pion_field[i, sk.grid.lp[2]-1, k, :] .= vac
    end

end

"""
    evaluate_sk(skyrmion, y)

Evaluates the Skyrme field at the spatial position `y`.

See also [`quadratic_spline_interpolation`](@ref). 

"""
function evaluate_sk(skyrmion, y)

    x = skyrmion.grid.x
    vac = [0.0, 0.0, 0.0, 1.0]
    phi=vac

    phiinterp = quadratic_spline_interpolation(skyrmion.pion_field, x)

    if x[1][1] < y[1] < x[1][end] &&
       x[2][1] < y[2] < x[2][end] &&
       x[3][1] < y[3] < x[3][end]
        for a = 1:4
            phi[a] = phiinterp[a](y[1], y[2], y[3])
        end
    end

    return phi/sqrt(phi'*phi)

end


"""
    quadratic_spline_interpolation(pion_field, x)

Provides a function which interpolates `pion_field` on the grid `x` 

In particular, the return object is an array of interpolation objects from the package `Interpolations`, one for each component of the pion field. These are constructed using a quadatic b-spline. 

"""
function quadratic_spline_interpolation(pion_field, x)

    phiinterp = [
        extrapolate(
            scale(
                interpolate(pion_field[:, :, :, a], BSpline(Quadratic())),
                (x[1], x[2], x[3]),
            ),
            Throw(),
        ) for a = 1:4
    ]

    return phiinterp

end
