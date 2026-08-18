"""
VelocityBoxKernels.jl

ParallelStencil kernels that write prescribed velocities from `VelBox2D`
boxes (see VelocityBoxes.jl) into the Stokes solver's velocity arrays, plus
`apply_vel_boxes!`, the driver that loops over all active boxes and calls
them.

ORDERING REQUIREMENT: this file uses `@parallel_indices` / `@parallel` /
`@idx`, which depend on `@init_parallel_stencil` having already run in your
main script. `include` this file AFTER `@init_parallel_stencil(...)` -- the
same place the equivalent inline code used to live in the main script. It
also needs the `VelBox2D` type from VelocityBoxes.jl, so include that file
first too.
"""

# Velocity boxes are applied on the same staggered coordinates as the Stokes
# solver. In the Geometry API these coordinates are stored in `grid.xi_vel`:
# - `grid.xi_vel[1]` are the coordinates for Vx (x-face, z)
# - `grid.xi_vel[2]` are the coordinates for Vy (x, z-face)
# so the box region is applied to the correct velocity DoFs.

@parallel_indices (i, j) function _apply_vel_box_Vx!(
        Vx,
        xvx,
        yvx,
        cenx,
        cenz,
        halfx,
        halfz,
        vx_val,
    )
    if i ≤ size(Vx, 1) && j ≤ size(Vx, 2)
        x = xvx[i]
        z = yvx[j]
        if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
            @inbounds Vx[i, j] = vx_val
        end
    end
    return nothing
end

@parallel_indices (i, j) function _apply_vel_box_Vy!(
        Vy,
        xvy,
        yvy,
        cenx,
        cenz,
        halfx,
        halfz,
        vy_val,
    )
    if i ≤ size(Vy, 1) && j ≤ size(Vy, 2)
        x = xvy[i]
        z = yvy[j]
        if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
            @inbounds Vy[i, j] = vy_val
        end
    end
    return nothing
end

@parallel_indices (i, j) function _mark_vbox_mask_Vx!(
        mask_vbox_x,
        xvx,
        yvx,
        cenx,
        cenz,
        halfx,
        halfz,
    )
    if i ≤ size(mask_vbox_x, 1) && j ≤ size(mask_vbox_x, 2)
        # mask indices (i,j) correspond to velocity DoFs at (i+1,j+1)
        ii = i + 1
        jj = j + 1
        if ii ≤ length(xvx) && jj ≤ length(yvx)
            x = xvx[ii]
            z = yvx[jj]
            if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
                @inbounds mask_vbox_x[i, j] = 1
            end
        end
    end
    return nothing
end

@parallel_indices (i, j) function _mark_vbox_mask_Vy!(
        mask_vbox_y,
        xvy,
        yvy,
        cenx,
        cenz,
        halfx,
        halfz,
    )
    if i ≤ size(mask_vbox_y, 1) && j ≤ size(mask_vbox_y, 2)
        # mask indices (i,j) correspond to velocity DoFs at (i+1,j+1)
        ii = i + 1
        jj = j + 1
        if ii ≤ length(xvy) && jj ≤ length(yvy)
            x = xvy[ii]
            z = yvy[jj]
            if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
                @inbounds mask_vbox_y[i, j] = 1
            end
        end
    end
    return nothing
end

@parallel_indices (i, j) function _mark_vbox_mask_center!(mask_vbox_c, xc, zc, cenx, cenz, halfx, halfz)
    if i ≤ size(mask_vbox_c, 1) && j ≤ size(mask_vbox_c, 2)
        x = xc[i]
        z = zc[j]
        if abs(x - cenx) ≤ halfx && abs(z - cenz) ≤ halfz
            @inbounds mask_vbox_c[i, j] = 1
        end
    end
    return nothing
end

"""
    apply_vel_boxes!(stokes, grid, boxes::Vector{VelBox2D}, mask_vbox_c)

Overwrite `stokes`'s velocity arrays inside each box in `boxes` with its
prescribed value(s), and set the corresponding entries of
`stokes.mask_vbox_x`, `stokes.mask_vbox_y`, and `mask_vbox_c` so the solver
treats those degrees of freedom as fixed. Resets all three masks to zero
first, so calling this with an empty `boxes` list clears any previously
applied constraints.
"""
function apply_vel_boxes!(
        stokes,
        grid,
        boxes::Vector{VelBox2D},
        mask_vbox_c,
    )
    isempty(boxes) && return nothing
    Vx, Vy = @velocity(stokes)
    grid_vx, grid_vy = grid.xi_vel
    xvx, yvx = grid_vx
    xvy, yvy = grid_vy
    xc, zc = grid.xci
    stokes.mask_vbox_x.mask .= 0
    stokes.mask_vbox_y.mask .= 0
    mask_vbox_c .= 0
    for box in boxes
        halfx = box.widthx / 2
        halfz = box.widthz / 2
        if box.has_vx
            nx = length(xvx)
            ny = length(yvx)
            @parallel (@idx (nx, ny)) _apply_vel_box_Vx!(Vx, xvx, yvx, box.cenx, box.cenz, halfx, halfz, box.vx)
            @parallel (@idx (nx, ny)) _mark_vbox_mask_Vx!(stokes.mask_vbox_x.mask, xvx, yvx, box.cenx, box.cenz, halfx, halfz)
        end
        if box.has_vy
            nx = length(xvy)
            ny = length(yvy)
            @parallel (@idx (nx, ny)) _apply_vel_box_Vy!(Vy, xvy, yvy, box.cenx, box.cenz, halfx, halfz, box.vy)
            @parallel (@idx (nx, ny)) _mark_vbox_mask_Vy!(stokes.mask_vbox_y.mask, xvy, yvy, box.cenx, box.cenz, halfx, halfz)
        end
        if box.has_vx || box.has_vy
            nxc = length(xc)
            nzc = length(zc)
            @parallel (@idx (nxc, nzc)) _mark_vbox_mask_center!(mask_vbox_c, xc, zc, box.cenx, box.cenz, halfx, halfz)
        end
    end
    return nothing
end
