using GeophysicalModelGenerator

function GMG_subduction_2D(nx, ny)
    model_depth = 300.0 # km
    nx, nz = nx, ny
    Tbot = 1474.0
    xmin = 0
    xmax = 1500 # km
    extra_air_thickness = 0 #12.5
    x = range(0, xmax, nx)
    z = range(-model_depth, extra_air_thickness, nz)
    Grid2D = CartData(xyz_grid(x, 0, z))
    Phases = zeros(Int64, nx, 1, nz)
    Temp = fill(Tbot, nx, 1, nz)
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

    # Temperature field is simple layers + halfspace cooling. To be improved.
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-8, 0),
        phase=ConstantPhase(1),
        T=ConstantTemp(T=0)
    )

    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-112.5, -8),
        phase=ConstantPhase(1),
        T=HalfspaceCoolingTemp(Tsurface=0, Tmantle=Tlab, Age=40, Adiabat=0)
    )

    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-model_depth, -112.5),
        phase=ConstantPhase(1),
        T=LinearTemp(Ttop=Tlab, Tbot=Tbot)
    )

    # Material phases described with polygons similar to SZU2019.
    # 0: asthenosphere - mantle
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-model_depth, 0),
        Origin=nothing, StrikeAngle=0, DipAngle=0,
        phase=ConstantPhase(1),
    )
    # 2: air
    # added explicitly at end?
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-12.5, extra_air_thickness),
        Origin=nothing, StrikeAngle=0, DipAngle=0,
        phase=ConstantPhase(2),
    )

    # x: asthenosphere2 - for visualisation
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmin, xmax, xmax),
        zlim=(-19.5, -112.5, -112.5, -19.5),
        phase=ConstantPhase(1), # NOTE: Changed this visualisation field 23 January 2026 to relieve size of Rheology tuple
    )


    # # x: mantle weakzone
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(920, 1035, 1050, 970),
        zlim=(-25, -80, -80, -38),
        phase=ConstantPhase(3) #4),#Making phases contiguous (28 jan, v0.76)
    )

    # # 4: oceanic crust/interface material
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(0, 0, 916, 902),
        zlim=(-12.5, -19.5, -19.5, -12.5),
        phase=ConstantPhase(4),# 5)#Making phases contiguous (28 jan, v0.76)
    )

    # # 5: oceanic gabbro ## incorporated on top
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(0, 0, 916, 906),
        zlim=(-14.5, -19.5, -19.5, -14.5),
        phase=ConstantPhase(5)# 6)#Making phases contiguous (28 jan, v0.76)
    )

    ####
    # 6: felsic crust1
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(954, 854, xmax, xmax),
        zlim=(-8, -19.5, -11.5, -8),
        phase=ConstantPhase(6) # 7),#Making phases contiguous (28 jan, v0.76) 
    )
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 890, xmax, xmax),
        zlim=(-19.5, -23, -23, -11.5),
        phase=ConstantPhase(6) # 7),#Making phases contiguous (28 jan, v0.76) 
    )

    # 6: felsic crust 2
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(890, 960, xmax, xmax),
        zlim=(-23, -38, -38, -23),
        phase=ConstantPhase(6) # 7),#Making phases contiguous (28 jan, v0.76)  # NOTE: Changed this visualisation field 23 January 2026 to relieve size of Rheology tuple
    )

    # 7: sediments
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(954, 854, 1007, 1030),
        zlim=(-8, -12.5, -19.5, -8),
        phase=ConstantPhase(7) # 9), #Making phases contiguous (28 jan, v0.76)
    )

    # 7: sediments 2
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

    Grid2D = addfield(Grid2D, (; Phases, Temp))
    # write_paraview(Grid2D, "Initial_Setup_Subduction_rank")
    li = (abs(last(x) - first(x)), abs(last(z) - first(z))) .* 1.0e3
    origin = (x[1], z[1]) .* 1.0e3

    ph = Phases[:, 1, :]
    T = Temp[:, 1, :]

    return li, origin, ph, T
end
