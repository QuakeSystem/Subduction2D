"""
VelocityBoxes.jl

Defines a `VelBox2D`: a rectangular region of the domain where a velocity
component is prescribed directly (a Dirichlet-like constraint), plus the
list of boxes currently active in the model and a constructor to add one.

This file is pure Julia with no ParallelStencil/JustRelax dependency, so it
can be `include`d at any point in your script -- but note that
`GMG_subduction_2D_with_coords` (in Subduction2D_setup.jl) calls
`add_vel_box!`, so this file must be included BEFORE the setup file.

Applying these boxes to the actual solver arrays (the parallel kernels that
write into `stokes`) lives in VelocityBoxKernels.jl instead, since that code
needs ParallelStencil to already be initialized (`@init_parallel_stencil`)
before it can be defined
"""

struct VelBox2D
    cenx::Float64
    cenz::Float64
    widthx::Float64
    widthz::Float64
    vx::Float64
    vy::Float64
    has_vx::Bool
    has_vy::Bool
end

const vel_boxes_2D = VelBox2D[]

"""
    add_vel_box!(; cenx, cenz, widthx, widthz, vx=nothing, vy=nothing)

Register a velocity box centered at `(cenx, cenz)` with size `(widthx,
widthz)` (all in meters). Pass `vx` and/or `vy` (in m/s) to prescribe that
velocity component inside the box; omit one to leave it unconstrained there.
"""
function add_vel_box!(
        ; cenx,
        cenz,
        widthx,
        widthz,
        vx = nothing,
        vy = nothing,
    )
    vx_val = vx === nothing ? 0.0 : Float64(vx)
    vy_val = vy === nothing ? 0.0 : Float64(vy)
    has_vx = vx !== nothing
    has_vy = vy !== nothing
    push!(
        vel_boxes_2D,
        VelBox2D(
            Float64(cenx),
            Float64(cenz),
            Float64(widthx),
            Float64(widthz),
            vx_val,
            vy_val,
            has_vx,
            has_vy,
        ),
    )
    return nothing
end

clear_vel_boxes!() = empty!(vel_boxes_2D)
