"""
VelocityBoxKernels.jl

Builds a `VelocityBoundaryConditions` whose `dirichlet` field enforces the
currently-registered `VelBox2D` boxes (see VelocityBoxes.jl) as an internal
prescribed-velocity region.

This replaces the older `stokes.mask_vbox_x`/`mask_vbox_y` +
`apply_velocity_box` mechanism now that JustRelax's DYREL solver enforces
internal velocity Dirichlet regions itself via
`VelocityBoundaryConditions.dirichlet` (JustRelax.jl branch
`feature/velbox-dirichlet-bc`). The DYREL DR kernel re-asserts the prescribed
value at every pseudo-transient iteration on its own, so there is no separate
pre-solve velocity write anymore -- just build `flow_bcs` with
`velocity_box_flow_bcs` and hand it to `solve_DYREL!` as before.

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
#
# `mask`/`value` here are sized exactly like the velocity component array
# they belong to (Vx or Vy) -- unlike the old mask_vbox_x/y (sized like the
# interior-only residual array Rx/Ry), so there is no index offset to keep in
# sync with the DYREL kernels by hand.

@parallel_indices (i, j) function _mark_vbox!(mask, value, xv, yv, cenx, cenz, halfx, halfz, v_val)
    if i <= size(mask, 1) && j <= size(mask, 2)
        x = xv[i]; z = yv[j]
        if abs(x - cenx) <= halfx && abs(z - cenz) <= halfz
            @inbounds mask[i, j] = 1
            @inbounds value[i, j] = v_val
        end
    end
    return nothing
end

"""
    velocity_box_flow_bcs(stokes, grid, boxes::Vector{VelBox2D}; free_slip, free_surface = false)

Build a `VelocityBoundaryConditions` with the interior `dirichlet` region set
from `boxes`. Call this once per timestep before `solve_DYREL!` -- rebuilding
it is cheap (a couple of small array allocations + a `@parallel` mask pass),
and picks up any change to `boxes` (e.g. a ramped-up prescribed velocity)
immediately.

Uses `DirichletBoundaryCondition(value, Mask(mask))` directly (not the
`(; constant, mask)` shorthand) so that a box prescribing exactly zero
velocity is still correctly marked as constrained -- the shorthand's
array form infers its mask from non-zero value entries, which cannot
represent a prescribed value of exactly zero.
"""
function velocity_box_flow_bcs(
        stokes, grid, boxes::Vector{VelBox2D};
        free_slip = (left = true, right = true, top = true, bot = true),
        free_surface = false,
    )
    Vx, Vy = @velocity(stokes)
    xvx, yvx = grid.xi_vel[1]
    xvy, yvy = grid.xi_vel[2]

    mask_x, value_x = @zeros(size(Vx)...), @zeros(size(Vx)...)
    mask_y, value_y = @zeros(size(Vy)...), @zeros(size(Vy)...)

    for box in boxes
        halfx, halfz = box.widthx / 2, box.widthz / 2
        if box.has_vx
            @parallel (@idx size(mask_x)) _mark_vbox!(mask_x, value_x, xvx, yvx, box.cenx, box.cenz, halfx, halfz, box.vx)
        end
        if box.has_vy
            @parallel (@idx size(mask_y)) _mark_vbox!(mask_y, value_y, xvy, yvy, box.cenx, box.cenz, halfx, halfz, box.vy)
        end
    end

    dirichlet = (;
        Vx = JustRelax.DirichletBoundaryCondition(value_x, JustRelax.Mask(mask_x)),
        Vy = JustRelax.DirichletBoundaryCondition(value_y, JustRelax.Mask(mask_y)),
    )
    return VelocityBoundaryConditions(; free_slip = free_slip, free_surface = free_surface, dirichlet = dirichlet)
end
## END OF HELPER FUNCTION ------------------------------------------------------------
