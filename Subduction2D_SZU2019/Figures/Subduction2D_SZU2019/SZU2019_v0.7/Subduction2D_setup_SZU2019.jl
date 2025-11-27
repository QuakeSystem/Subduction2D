using GeophysicalModelGenerator

function GMG_subduction_2D(nx, ny)
    model_depth = 300.0 # km
    # Our starting basis is the example above with ridge and overriding slab
    nx, nz = nx, ny
    Tbot = 1474.0
    xmin = 0
    xmax = 1500 # km

    x = range(0, xmax, nx)
    air_thickness = 12.5
    z = range(-model_depth, air_thickness, nz)
    Grid2D = CartData(xyz_grid(x, 0, z))
    Phases = zeros(Int64, nx, 1, nz)
    Temp = fill(Tbot, nx, 1, nz)
    Tlab = 1300
    # lith   = LithosphericPhases(Layers=[80], Phases=[1 0], Tlab=Tlab)

    # phases
    # 0: asthenosphere - mantle
    # 4: air
    # 0: asthenosphere - mantle
    # x: mantle weakzone
    # 3: oceanic crust
    # x: oceanic gabbro
    # x: felsic crust (two parts for visualization)
    # x: sediments
    # x: deflected oceanic crust
    # x: deflected gabbro

    # 0: asthenosphere - mantle
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-model_depth, 0),
        Origin=nothing, StrikeAngle=0, DipAngle=0,
        phase=LithosphericPhases(Layers=[], Phases=[0], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        T=HalfspaceCoolingTemp(Tsurface=20, Tmantle=Tbot, Age=50, Adiabat=0)
    )
    # 2: air
    # added explicitly at end?
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-12.5, air_thickness),
        Origin=nothing, StrikeAngle=0, DipAngle=0,
        phase=LithosphericPhases(Layers=[], Phases=[1], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        T=HalfspaceCoolingTemp(Tsurface=20, Tmantle=Tbot, Age=50, Adiabat=0)
    )

    # x: asthenosphere2 - for visualisation
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmin, xmax, xmax),
        zlim=(-19.5, -70, -70, -19.5),
        phase=LithosphericPhases(Layers=[], Phases=[2], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )


    # # x: mantle weakzone
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(920, 1035, 1050, 970),
        zlim=(-25, -80, -80, -38),
        phase=LithosphericPhases(Layers=[], Phases=[3], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )

    # # 3: oceanic crust/interface material
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(0, 0, 906, 902),
        zlim=(-12.5, -14.5, -14.5, -12.5),
        phase=LithosphericPhases(Layers=[], Phases=[4], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )

    # x: oceanic gabbro
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(0, 0, 916, 906),
        zlim=(-14.5, -19.5, -19.5, -14.5),
        phase=LithosphericPhases(Layers=[], Phases=[5], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )
    # x: felsic crust1
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(954, 854, xmax, xmax),
        zlim=(-8, -19.5, -11.5, -8),
        phase=LithosphericPhases(Layers=[], Phases=[6], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 890, xmax, xmax),
        zlim=(-19.5, -23, -23, -11.5),
        phase=LithosphericPhases(Layers=[], Phases=[6], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )

    # x: felsic crust 2
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(890, 960, xmax, xmax),
        zlim=(-23, -38, -38, -23),
        phase=LithosphericPhases(Layers=[], Phases=[7], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )

    # x: sediments
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(954, 854, 1007, 1030),
        zlim=(-8, -12.5, -19.5, -8),
        phase=LithosphericPhases(Layers=[], Phases=[8], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )

    # x: sediments 2
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 960, 970, 1007),
        zlim=(-12.5, -38, -38, -19.5),
        phase=LithosphericPhases(Layers=[], Phases=[9], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )


    #   /Deflected_Hydrated_fractured_top_oceanic_crust_Thrust_Interface
    # Same rheology as Oceanic Crust
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 854, 967, 980),
        zlim=(-12.5, -18.5, -39.3, -33),
        phase=LithosphericPhases(Layers=[], Phases=[4 + 1], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )

    #    /Deflected_Crust_(Gabbro)
    # Same rheology as Oceanic crust (Gabbro)
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 854, 925, 900),
        zlim=(-14.5, -19.5, -32, -21),
        phase=LithosphericPhases(Layers=[], Phases=[5 + 1], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )

    #       Fixed_Air
    # add_polygon!(
    #     Phases,
    #     Temp,
    #     Grid2D;
    #     xlim=(),
    #     zlim=(),
    #     phase=LithosphericPhases(Layers=[], Phases=[12], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
    #     # T=LinearTemp(Ttop=20, Tbot=Tbot)
    # )

    # #       Fixed_Asthenosphere
    # add_polygon!(
    #     Phases,
    #     Temp,
    #     Grid2D;
    #     xlim=(),
    #     zlim=(),
    #     phase=LithosphericPhases(Layers=[], Phases=[13], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
    #     # T=LinearTemp(Ttop=20, Tbot=Tbot)
    # )

    # surf = Grid2D.z.val .> 0.0
    # Temp[surf] .= 20.0
    # Phases[surf] .= 3

    Grid2D = addfield(Grid2D, (; Phases, Temp))
    # write_paraview(Grid2D, "Initial_Setup_Subduction_rank")
    li = (abs(last(x) - first(x)), abs(last(z) - first(z))) .* 1.0e3
    origin = (x[1], z[1]) .* 1.0e3

    ph = Phases[:, 1, :] .+ 1
    T = Temp[:, 1, :]

    return li, origin, ph, T
end
