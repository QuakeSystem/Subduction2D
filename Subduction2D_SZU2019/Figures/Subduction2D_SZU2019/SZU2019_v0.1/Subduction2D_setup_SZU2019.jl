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
        zlim=(-8, -19.5, -8, -11.5),
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
        xlim=(854, 890, xmax, xmax),
        zlim=(-19.5, -23, -23, -11.5),
        phase=LithosphericPhases(Layers=[], Phases=[8], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
        # T=LinearTemp(Ttop=20, Tbot=Tbot)
    )


    # # x: deflected oceanic crust
    # add_polygon!(
    #     Phases,
    #     Temp,
    #     Grid2D;
    #     xlim=(854, 854, 967, 980),
    #     zlim=(-12.5, -18.5, -39.3, -33),
    #     phase=LithosphericPhases(Layers=[], Phases=[9], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
    #     # T=LinearTemp(Ttop=20, Tbot=Tbot)
    # )

    # # x: deflected gabbro
    # add_polygon!(
    #     Phases,
    #     Temp,
    #     Grid2D;
    #     xlim=(854, 854, 925, 900),
    #     zlim=(-14.5, -19.5, -32, -21),
    #     phase=LithosphericPhases(Layers=[], Phases=[10], Tlab=Tlab), # √GeophysicalModelGenerator.jl/src/Setup_geometry.jl
    #     # T=LinearTemp(Ttop=20, Tbot=Tbot)
    # )


    # add_box!(
    #     Phases,
    #     Temp,
    #     Grid2D;
    #     xlim=(0,),
    #     zlim=(-model_depth, 0.0),
    #     Origin=nothing, StrikeAngle=0, DipAngle=0,
    #     phase=LithosphericPhases(Layers=[80], Phases=[1 0], Tlab=Tlab),
    #     T=HalfspaceCoolingTemp(Tsurface=20, Tmantle=Tbot, Age=50, Adiabat=0)
    # )

    # # Add right oceanic plate crust
    # add_box!(
    #     Phases,
    #     Temp,
    #     Grid2D;
    #     xlim=(1500 - 715, 1500 - 100),
    #     zlim=(-model_depth, 0.0),
    #     Origin=nothing, StrikeAngle=0, DipAngle=0,
    #     phase=LithosphericPhases(Layers=[8 72], Phases=[2 1 0], Tlab=Tlab),
    #     T=HalfspaceCoolingTemp(Tsurface=20, Tmantle=Tbot, Age=50, Adiabat=0)
    # )

    # # Add slab
    # add_box!(
    #     Phases,
    #     Temp,
    #     Grid2D;
    #     xlim=(800, 250),
    #     zlim=(-80, 0.0),
    #     Origin=nothing, StrikeAngle=0, DipAngle=-15,
    #     phase=LithosphericPhases(Layers=[8 80], Phases=[2 1 0], Tlab=Tlab),
    #     T=HalfspaceCoolingTemp(Tsurface=20, Tmantle=Tbot, Age=50, Adiabat=0)
    # )

    # surf = Grid2D.z.val .> 0.0
    # Temp[surf] .= 20.0
    # Phases[surf] .= 3

    Grid2D = addfield(Grid2D, (; Phases, Temp))
    write_paraview(Grid2D, "Initial_Setup_Subduction_rank")
    li = (abs(last(x) - first(x)), abs(last(z) - first(z))) .* 1.0e3
    origin = (x[1], z[1]) .* 1.0e3

    ph = Phases[:, 1, :] .+ 1
    T = Temp[:, 1, :]

    return li, origin, ph, T
end
