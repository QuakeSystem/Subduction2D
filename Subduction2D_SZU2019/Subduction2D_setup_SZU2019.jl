using GeophysicalModelGenerator

function GMG_subduction_2D(nx, ny)
    model_depth = 175.0 # km
    nx, nz = nx, ny
    Tbot = 1718 - 273.15 # K, 1445 C, from Katsura 2022
    xmin = 0 # km, left edge of model
    xmax = 1500 # km, right edge of model
    extra_air_thickness = 0 # km, can be increased
    x = range(0, xmax, nx)
    z = range(-model_depth, extra_air_thickness, nz)
    Grid2D = CartData(xyz_grid(x, 0, z))
    Phases = zeros(Int64, nx, 1, nz)
    Temp = fill(0.0, nx, 1, nz)
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
        xlim=(854, 940, 1035, 1100, xmin, xmin),
        zlim=(-12.5, -35, -80, -112.5, -112.5, -12.5),
        T=HalfspaceCoolingTemp(Tsurface=0, Tmantle=Tlab, Age=80, Adiabat=0.5)
    )

    # Overriding plate temperature, linear geotherm with T0=0 and Tbot=TLab=1300C
    # NOTE: slight error for the left top edge (-12.5km to -8 km gradient), as temperature here is 0-39C at the surface (gradient is vertical with -8 as surface, leaving left side exposed
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(854, 940, 1035, 1100, xmax, xmax, 954),
        zlim=(-12.5, -35, -80, -112.5, -112.5, -8, -8),
        T=LinearTemp(Ttop=0, Tbot=Tlab)
    )

    # Air temperature, constant zero C.
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmin, 854,  954, xmax, xmax),
        zlim=(extra_air_thickness,-12.5, -12.5,  -8, -8, extra_air_thickness),
        T=ConstantTemp(T=0)
    )

    # Mantle temperature, linear geotherm with TLab = 1300 and T = Tbot 1445 C, from Katsura 2022.
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-model_depth, -112.5),
        T=LinearTemp(Ttop=Tlab, Tbot=Tbot)
    )

    #### Material phases described with polygons similar to SZU2019.
    #  asthenosphere - mantle
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-model_depth, 0),
        phase=ConstantPhase(1),
    )
    # air
    add_box!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmax),
        zlim=(-12.5, extra_air_thickness),
        phase=ConstantPhase(2),
    )

    # asthenosphere
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmin, xmax, xmax),
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
        xlim=(xmin, xmin, 916, 902),
        zlim=(-12.5, -19.5, -19.5, -12.5),
        phase=ConstantPhase(4),# 5)#Making phases contiguous (28 jan, v0.76)
    )

    # oceanic gabbro
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(xmin, xmin, 916, 906),
        zlim=(-14.5, -19.5, -19.5, -14.5),
        phase=ConstantPhase(5)# 6)#Making phases contiguous (28 jan, v0.76)
    )

    # felsic crust 1
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

    # felsic crust 2
    add_polygon!(
        Phases,
        Temp,
        Grid2D;
        xlim=(890, 960, xmax, xmax),
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

    Grid2D = addfield(Grid2D, (; Phases, Temp))
    # write_paraview(Grid2D, "Initial_Setup_Subduction_rank")
    li = (abs(last(x) - first(x)), abs(last(z) - first(z))) .* 1.0e3
    origin = (x[1], z[1]) .* 1.0e3

    ph = Phases[:, 1, :]
    T = Temp[:, 1, :] .+ 273

    return li, origin, ph, T
end
