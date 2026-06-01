using GeophysicalModelGenerator
using Statistics
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

# ------------------------------------------------------------
# Non-uniform grid helpers (logistic refinement)
# ------------------------------------------------------------

function logistic_vertices(
    n_cells::Int,
    L::Float64,
    x0::Float64;
    x_center::Float64,
    w_ref::Float64,
    refine_factor::Float64,
    k::Float64 = 4.0,
)
    # Compute cell widths with a logistic stretching weight, then integrate.
    widths_lin = fill(L / n_cells, n_cells)
    x_centers_u = [x0 + (i - 0.5) * (L / n_cells) for i in 1:n_cells]

    local_width_weight(x) = begin
        d = abs(x - x_center) / w_ref
        s = 1 / (1 + exp(-k * (d - 1.0)))
        1 + (refine_factor - 1) * s
    end

    widths = similar(widths_lin)
    @inbounds for i in 1:n_cells
        widths[i] = widths_lin[i] * local_width_weight(x_centers_u[i])
    end

    # Renormalize to preserve the total domain length exactly.
    widths .*= L / sum(widths)

    vertices = zeros(n_cells + 1)
    vertices[1] = x0
    @inbounds for i in 1:n_cells
        vertices[i + 1] = vertices[i] + widths[i]
    end

    return vertices, widths
end

function subduction_nonuniform_coords_1d(
    n_points::Int,
    x0::Float64,
    x1::Float64;
    ref_grid::Int,
    refine_factor::Float64,
    w_ref_ratio::Float64, # width of refined region as a fraction of total domain length
    x_center_frac::Float64,
    k::Float64 = 4.0, # logistic stretching parameter (higher k = sharper transition)
    verbose::Int = 0,
)
    if ref_grid == 0 || refine_factor == 1.0
        # Preserve the original type returned by `Geometry` (LinRange/StepRangeLen),
        # which downstream JustPIC advection currently dispatches on.
        return LinRange(x0, x1, n_points)
    end

    n_cells = n_points - 1
    L = x1 - x0
    x_center = x0 + x_center_frac * L
    w_ref = w_ref_ratio * abs(L)

    vertices, widths = logistic_vertices(
        n_cells,
        abs(L),
        x0;
        x_center = x_center,
        w_ref = w_ref,
        refine_factor = refine_factor,
        k = k,
    )

    if verbose == 1
        println("[grid verbose] non-uniform grid from $(x0) to $(x1):")
        println("  max cell size = $(maximum(widths))")
        println("  min cell size = $(minimum(widths))")
        println("  average cell size = $(mean(widths))")
        println("  refine_factor = $(refine_factor)")
        println("  w_ref_ratio = $(w_ref_ratio) -> w_ref = $(w_ref)")
        println("  center = $(x_center)")
    end

    return vertices
end

"""
Like `GMG_subduction_2D`, but also returns the non-uniform 1D coordinate vectors
needed by the DYREL grid setup:
`xvi` (staggered/vertex coordinates) and `xci` (cell-center coordinates) in meters.
"""
function GMG_subduction_2D_with_coords(
    nx_points::Int,
    ny_points::Int;
    ref_grid::Int = 0,
    refine_factor_x::Float64 = 30.0,
    refine_factor_y::Float64 = 10.0,
    w_ref_ratio_x::Float64 = 1/2, 
    w_ref_ratio_y::Float64 = 0.6, 
    k_x::Float64 = 8.0,
    k_y::Float64 = 12.0,
    x_center_frac::Float64 = 0.525,
    y_center_frac::Float64 = 0.9,
    verbose::Int = 1,
)
    model_depth = 260.0 # km
    Tsurface = 20
    Tbot = 1743 - 273.15 # K, 1445 C, from Katsura 2022 # increased by 25 since v0.271 for 260km depth
    x0_km, x1_km = 0.0, 1500.0
    air_thickness = 0.0
    z0_km, z1_km = -model_depth * 1.0, air_thickness

    # Our coordinate arrays are "points" for CartData: xvi has length nx_points.
    x = subduction_nonuniform_coords_1d(
        nx_points,
        x0_km,
        x1_km;
        ref_grid = ref_grid,
        refine_factor = refine_factor_x,
        w_ref_ratio = w_ref_ratio_x,
        x_center_frac = x_center_frac,
        k = k_x,
        verbose = verbose,
    )
    z = subduction_nonuniform_coords_1d(
        ny_points,
        z0_km,
        z1_km;
        ref_grid = ref_grid,
        refine_factor = refine_factor_y,
        w_ref_ratio = w_ref_ratio_y,
        x_center_frac = y_center_frac,
        k = k_y,
        verbose = verbose,
    )

    Grid2D = CartData(xyz_grid(x, 0, z))

    # Phases and temperature on the CartData grid points ------------------
    Phases = zeros(Int64, nx_points, 1, ny_points)
    Temp = fill(Tbot, nx_points, 1, ny_points)
    Tlab = 1300


    # phases
    # These were all unique up to v0.75 (28 jan 2026). Now made contiguous for simplicity (and prevent ABI overload).
    # 1: asthenosphere - mantle
    # 2: air
    # 3: mantle weakzone
    # 4: oceanic crust
    # 5: oceanic gabbro
    # 6: felsic crust
    # 7: sediments
    # 4: deflected oceanic crust
    # 5: deflected gabbro


    #### Temperature fields. Four total
    # Subducting plate temperature, halfspace cooling (80 Myr age)
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 940, 1035, 1100, x0_km, x0_km),
        zlim=(-12.5, -35, -80, -112.5, -112.5, -12.5),
        T=HalfspaceCoolingTemp(Tsurface=Tsurface, Tmantle=Tlab, Age=80, Adiabat=0.5)
    )

    # Overriding plate temperature, linear geotherm with T0=0 and Tbot=TLab=1300C
    # NOTE: slight error for the left top edge (-12.5km to -8 km gradient), as temperature here is 0-39C at the surface (gradient is vertical with -8 as surface, leaving left side exposed
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(954, 854, 940, 1035, 1100, x1_km, x1_km),
        zlim=(-8, -12.5, -35, -80, -112.5, -112.5, -8),
        T=LinearTemp(Ttop=Tsurface, Tbot=Tlab)
    )

    # Air temperature, constant zero C.
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km, x0_km, 854,  954, x1_km, x1_km),
        zlim=(z1_km,-12.5, -12.5,  -8, -8, z1_km),
        T=ConstantTemp(T=0)
    )

    # Mantle temperature, linear geotherm with TLab = 1300 and T = Tbot 1445 C, from Katsura 2022.
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km, x1_km),
        zlim=(-model_depth, -112.5),
        T=LinearTemp(Ttop=Tlab, Tbot=Tbot)
    )

    #### Material phases described with polygons similar to SZU2019.
    #  asthenosphere - mantle
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km, x1_km),
        zlim=(-model_depth, 0),
        phase=ConstantPhase(1),
    )
    # air
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km, x1_km),
        zlim=(-12.5, z1_km),
        phase=ConstantPhase(2),
    )

    # asthenosphere
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km, x0_km, x1_km, x1_km),
        zlim=(-19.5, -112.5, -112.5, -19.5),
        phase=ConstantPhase(1), # NOTE: Changed this visualisation field 23 January 2026 to relieve size of Rheology tuple
    )


    # mantle weakzone
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(920, 1035, 1050, 970),
        zlim=(-25, -80, -80, -38),
        phase=ConstantPhase(3) #4),#Making phases contiguous (28 jan, v0.76)
    )

    # oceanic crust/interface material
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km, x0_km, 906, 902),
        zlim=(-12.5, -14.5, -14.5, -12.5),
        phase=ConstantPhase(4),# 5)#Making phases contiguous (28 jan, v0.76)
    )

    # oceanic gabbro
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km, x0_km, 916, 906),
        zlim=(-14.5, -19.5, -19.5, -14.5),
        phase=ConstantPhase(5)# 6)#Making phases contiguous (28 jan, v0.76)
    )

    # felsic crust 1
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(954, 854, x1_km, x1_km),
        zlim=(-8, -19.5, -11.5, -8),
        phase=ConstantPhase(6) # 7),#Making phases contiguous (28 jan, v0.76) 
    )
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 890, x1_km, x1_km),
        zlim=(-19.5, -23, -23, -11.5),
        phase=ConstantPhase(6) # 7),#Making phases contiguous (28 jan, v0.76) 
    )

    # felsic crust 2
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(890, 960, x1_km, x1_km),
        zlim=(-23, -38, -38, -23),
        phase=ConstantPhase(6) # 7),#Making phases contiguous (28 jan, v0.76)  # NOTE: Changed this visualisation field 23 January 2026 to relieve size of Rheology tuple
    )

    # sediments
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(954, 854, 1007, 1030),
        zlim=(-8, -12.5, -19.5, -8),
        phase=ConstantPhase(7) # 9), #Making phases contiguous (28 jan, v0.76)
    )

    # sediments 2
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 960, 970, 1007),
        zlim=(-12.5, -38, -38, -19.5),
        phase=ConstantPhase(7) # 9), #Making phases contiguous (28 jan, v0.76) # NOTE: Changed this visualisation field 23 January 2026 to relieve size of Rheology tuple
    )


    #   /Deflected_Hydrated_fractured_top_oceanic_crust_Thrust_Interface
    # Same rheology as Oceanic Crust (4)
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 854, 967, 980),
        zlim=(-12.5, -18.5, -39.3, -33),
        phase=ConstantPhase(4)# 5),#Making phases contiguous (28 jan, v0.76)
    )

    #    /Deflected_Crust_(Gabbro)
    # Same rheology as Oceanic crust (Gabbro) (5)
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 854, 925, 900),
        zlim=(-14.5, -19.5, -32, -21),
        phase=ConstantPhase(5)# 6),#Making phases contiguous (28 jan, v0.76)
    )
    #    left boundary
    # low viscous left boundary
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(x0_km,x0_km+10),
        zlim=(-100, 0),
        phase=ConstantPhase(8),
        T=LinearTemp(Ttop=200, Tbot=Tlab+200)
        
    )
    add_vel_box!(
        cenx   = 180 * 1.0e3,  # m
        cenz   = -54.8 * 1.0e3,          # m
        widthx = 20 * 1.0e3,          # m
        widthz = 23.6 * 1.0e3,           # m
        vx     = 0 * 0.01 / (3600*24*365)             # m/s (optional)
    #    vy     = -4.0e-9,                    # m/s (optional)
    )

    Grid2D = addfield(Grid2D, (; Phases, Temp))
    # write_paraview(Grid2D, "Initial_Setup_Subduction_rank")

    # Convert coordinates/geometry to meters for DYREL
    li = (abs(last(x) - first(x)), abs(last(z) - first(z))) .* 1.0e3
    origin = (x[1], z[1]) .* 1.0e3

    ph = Phases[:, 1, :]
    T = Temp[:, 1, :] .+ 273

    # Staggered grid coordinate vectors in meters:
    # - xvi are vertices (length nx_points)
    # - xci are cell centers (length nx_points-1)
    xvi = (x .* 1.0e3, z .* 1.0e3)
    xci = (
        0.5 .* (x[1:end-1] .+ x[2:end]) .* 1.0e3,
        0.5 .* (z[1:end-1] .+ z[2:end]) .* 1.0e3,
    )

    return li, origin, ph, T, xvi, xci
end
