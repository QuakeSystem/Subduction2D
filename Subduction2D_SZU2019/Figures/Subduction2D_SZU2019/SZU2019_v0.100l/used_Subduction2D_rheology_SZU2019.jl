using GeoParams.Dislocation
using GeoParams.Diffusion
using GeoParams.MaterialParameters
searchsorted

function init_rheology_nonNewtonian()
    #dislocation laws
    disl_wet_olivine = SetDislocationCreep(Dislocation.wet_olivine1_Hirth_2003)
    # diffusion laws
    diff_wet_olivine = SetDiffusionCreep(Diffusion.wet_olivine_Hirth_2003)

    el = ConstantElasticity(; G=40.0e9)

    lithosphere_rheology = CompositeRheology((el, disl_wet_olivine, diff_wet_olivine))
    return init_rheologies(lithosphere_rheology)
end

function init_rheology_nonNewtonian_plastic()
    #dislocation laws
    disl_wet_olivine = SetDislocationCreep(Dislocation.wet_olivine1_Hirth_2003)
    # diffusion laws
    diff_wet_olivine = SetDiffusionCreep(Diffusion.wet_olivine_Hirth_2003)
    # plasticity
    ϕ_wet_olivine = asind(0.1)
    C_wet_olivine = 1.0e6
    η_reg = 1.0e16
    el = ConstantElasticity(; G=40.0e9, ν=0.45)
    lithosphere_rheology = CompositeRheology(
        (
        el,
        disl_wet_olivine,
        diff_wet_olivine,
        DruckerPrager_regularised(; C=C_wet_olivine, ϕ=ϕ_wet_olivine, η_vp=η_reg, Ψ=0.0), # non-regularized plasticity
    )
    )
    return init_rheologies(lithosphere_rheology)
end

function init_rheology_linear()
    el = ConstantElasticity(; G=40.0e9, ν=0.45)
    # lithosphere_rheology = CompositeRheology( (LinearViscous(; η=1e23), ))
    lithosphere_rheology = CompositeRheology((LinearViscous(; η=1.0e23), el))
    return init_rheologies(lithosphere_rheology)
end


function init_rheologies()
    el = ConstantElasticity(; G=40.0e9, ν=0.45)
    # lithosphere_rheology = CompositeRheology( (LinearViscous(; η=1e23), ))
    lithosphere_rheology = CompositeRheology((LinearViscous(; η=1.0e23), el))
    # common physical properties
    α = 2.4e-5 # 1 / K
    Cp = 750    # J / kg K
    # Define rheolgy struct
    return rheology = (
        # Name = "Asthenoshpere",
        SetMaterialParams(;
            Phase=1,
            Density=ConstantDensity(; ρ=3.2e3),
            HeatCapacity=ConstantHeatCapacity(; Cp=Cp),
            Conductivity=ConstantConductivity(; k=2.5),
            CompositeRheology=CompositeRheology((LinearViscous(; η=1.0e20),)),
            Gravity=ConstantGravity(; g=9.81),
        ),
        # Name              = "Oceanic lithosphere",
        SetMaterialParams(;
            Phase=2,
            Density=PT_Density(; ρ0=3.2e3, α=α, β=0.0e0, T0=273 + 1474),
            HeatCapacity=ConstantHeatCapacity(; Cp=Cp),
            Conductivity=ConstantConductivity(; k=2.5),
            CompositeRheology=lithosphere_rheology,
        ),
        # Name              = "oceanic crust",
        SetMaterialParams(;
            Phase=3,
            Density=ConstantDensity(; ρ=3.2e3),
            HeatCapacity=ConstantHeatCapacity(; Cp=Cp),
            Conductivity=ConstantConductivity(; k=2.5),
            CompositeRheology=CompositeRheology((LinearViscous(; η=1.0e20),)),
        ),
        # Name              = "StickyAir",
        SetMaterialParams(;
            Phase=4,
            Density=ConstantDensity(; ρ=100), # water density
            HeatCapacity=ConstantHeatCapacity(; Cp=3.0e3),
            Conductivity=ConstantConductivity(; k=1.0),
            CompositeRheology=CompositeRheology((LinearViscous(; η=1.0e19),)),
        ),
    )
end

#       Mantle1_DRY_OL_Ranalli1995
#  /___NUM__NU(Pa^MM*s)DE(J)V(J/bar)SS(Pa)_MM(Power)_____LL(KOEF)_____RO(kg/M^3)_____bRo(1/K)_____aRo(1/kbar)CP(J/kg)Kt(Wt/(m*K))_Ht(Wt/kg)
#       num, minvis,maxvis,minsig,maxsig   1/prefac,     activ_energ.  act. vol       <diffcreep?>n(powerlaw)shrmod(Pa),fluidpres,coh1,coh2,  fric1,fric2,strweak?,strweak?,denscorr_markf0 unused, denscorr_markf1 unused.,density,      bb: tempdependence in adiab.          aa compr dep in adiab?    Cp????    base kt conduct. markkf temp dependence,markkp pressure dependence. ragiogenic heat prod. W/m3
#       10  1e+18 1e+26 0e+00 5e+29        3.98E+16      5.32E+05      0.80E+00      3.00E+04     3.50E+00  6.7E+10  1.00         1e+07 1e+07 0.600 0.600 0.5      1.5      0                       0                       3.30E+03      3.00E-05                              1.00E-03                  1.00E+03  0.73E+00         12.93e+02              4.00E-06                    2.20E-08

#=
Below function added to combine all rheologies used in SZU2019 model. No need for modularity currently.

    1. Function has been edited to make phase numbers contiguous (28 jan, v0.76)
    2. Conductivity parameter d is currently not sufficient. Some error in conversion between SZU2019 and GeoParams / Gerya's values.
    3. Same for Dislocation Creep. Dislocation creep parameters A, E, V, n are commented out for now (and are therefore some default). To be added later.
    4. Rheology names are kept the same as in SZU2019 for clarity.
    5. Density is currently constant, but should be expanded to PT_Density later. 

    Comments and commented out blocks have been cleaned up for clarity between v0.86 and v0.87.
=#
function init_rheologies()

    # common physical properties
    # α = 2.4e-5 # 1 / K
    # Cp = 750    # J / kg K
    # Define rheology struct
    return rheology = (
        SetMaterialParams(; Name="Mantle1_DRY_0",
            Phase=1,
            Density=ConstantDensity(ρ=3300), # Can and should be expanded to PT_Density
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=0.73e0,
                b=12.93e2,
                d=4e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=2.20E-08), # W/m3
            CompositeRheology=CompositeRheology((
                # elasticity
                ConstantElasticity(; G=6.7e10),
                # dislocation creep
                DislocationCreep(;
                # A=1 / 3.98e16,
                # A=2.5125e4,
                # E=5.32e5,
                # V=8.0e-6,
                # n=3.5,
                ),
                # Drucker–Prager plasticity
                DruckerPrager_regularised(;
                    C=1e7,
                    ϕ=asind(0.6),
                    η_vp=0
                )
            )),
        ),

        #       /Air
        #       0   1e+17 1e+17 0e+00 5e+04        1.00E+17      0.00E+05      0.00E+00      0.00E+04     1.00E+00  7.0E+11  0.00 0e+06 0e+06 0.000 0.000 0.0 1.0 0 0  1.00E+00      0.00E-05      0.00E-03      3.33E+06      2.00E+02   0.00E+00  0.00E+00  0.00E-10
        SetMaterialParams(;
            Name="Air",
            Phase=2,
            Density=ConstantDensity(; ρ=1.00E+02), # water density
            HeatCapacity=ConstantHeatCapacity(; Cp=3.0e3),
            Conductivity=ConstantConductivity(; k=1.0),
            CompositeRheology=CompositeRheology((LinearViscous(; η=1.00E+17),)),
        ),

        #       /Mantle_Weak_zone
        #       Shear_Zone_Mantle4_WET_OL_Ranalli1995
        #       12  1e+18 1e+26 0e+00 5e+29        5.01E+20      
        #       4.70E+05      0.80E+00      3.00E+04     4.00E+00  
        #       6.7E+10  1.00 1e+07 1e+07 0.100 0.100 0.5 1.5 0 0  
        #       3.30E+03      3.00E-05      1.00E-03      1.00E+03      
        #       0.73E+00  12.93e+02  4.00E-06  2.20E-08
        SetMaterialParams(; Name="Mantle_Weak_zone",
            Phase=3, #4, #Making phases contiguous (28 jan, v0.76)
            Density=ConstantDensity(ρ=3300), # Can and should be expanded to PT_Density
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=0.73e0,
                b=12.93e2,
                d=4.00e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=2.20e-08), # W/m3
            CompositeRheology=CompositeRheology((
                # elasticity
                ConstantElasticity(; G=6.7e10),

                # dislocation creep (wet olivine)
                DislocationCreep(;
                # A=1 / 5.01e20,
                # A=1.996e3,
                # E=4.70e5,
                # V=8.0e-6,
                # n=4.0,
                ),

                # Drucker–Prager plasticity
                DruckerPrager_regularised(;
                    C=1e7,
                    ϕ=asind(0.100),
                    η_vp=0)
            )),
        ),

        #       Oceanic crust-interface material
        #       /Hydrated_fractured_top_oceanic_crust_Thrust_Interface
        #       /Basalts_WET_QUARTZITE_RANALLI_1995
        #       7   1e+17 1e+26 1e+00 5e+29        1.97E+17      1.54E+05      0.80E+00      3.00E+04     2.30E+00  2.5E+10  1.00 6e+06 6e+06 10.00 10.00 0.5 1.5 0 0  3.00E+03      3.00E-05      1.00E-03      1.00E+03      1.18E+00   4.74E+02  4.00E-06  0.25E-06
        SetMaterialParams(; Name="Hydrated_fractured_top_oceanic_crust_Thrust_Interface",
            Phase=4, #5, #Making phases contiguous (28 jan, v0.76)
            Density=ConstantDensity(ρ=3000), # Can and should be expanded to PT_Density
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=1.18e0,
                b=4.74e2,
                d=4.00e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=0.25e-06), # W/m3
            CompositeRheology=CompositeRheology((
                ConstantElasticity(; G=2.5e10),
                DislocationCreep(;
                # A=3.2e-4,
                # E=1.54e5,
                # V=8.0e-6,
                # n=2.3,
                ),
                DruckerPrager_regularised(;
                    C=6e6,
                    ϕ=asind(0.5),
                    η_vp=0)
            )),
        ),


        #       /Oceanic_Crust_(Gabbro)
        #       /Basic_Crust2_An75_Ranalli1995
        #       8   1e+18 1e+26 0e-01 5e+29        4.80E+22      2.38E+05      0.80E+00      3.00E+04     3.20E+00  2.5E+10  1.00 1e+07 1e+07 0.850 0.850 0.5 1.5 0 0  3.00E+03      3.00E-05      1.00E-03      1.00E+03      1.18E+00   4.74E+02  4.00E-06  0.25E-06
        SetMaterialParams(; Name="Oceanic_Crust_(Gabbro)",
            Phase=5, # 6, #Making phases contiguous (28 jan, v0.76)
            Density=ConstantDensity(ρ=3000), # Can and should be expanded to PT_Density
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=1.18e0,
                b=4.74e2,
                d=4.00e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=0.25e-06), # W/m3
            CompositeRheology=CompositeRheology((
                ConstantElasticity(; G=2.5e10),
                DislocationCreep(;
                # A=1 / 4.80e22,
                # A=3.3e-4,
                # E=2.38e5,
                # V=8.0e-6,
                # n=3.2,
                ),
                DruckerPrager_regularised(;
                    C=1e7,
                    ϕ=asind(0.85),
                    η_vp=0)
            )),
        ),

        #       /Felsic_Crust1
        #       /Felsic_Crust1_WET_QUARTZITE_RANALLI_1995
        #       5   1e+18 1e+26 0e+00 5e+29        1.97E+17      1.54E+05      1.20E+00      3.00E+04     2.30E+00  2.5E+10  1.00 1e+07 1e+07 0.720 0.720 0.5 1.5 0 0  2.70E+03      3.00E-05      1.00E-03      1.00E+03      0.64E+00   8.07E+02  4.00E-06  1.00E-06

        SetMaterialParams(; Name="Felsic_Crust1",
            Phase=6, # 7, #Making phases contiguous (28 jan, v0.76)
            Density=ConstantDensity(ρ=2700), # Can and should be expanded to PT_Density
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=0.64e0,
                b=8.07e2,
                d=4.00e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=1.00e-06), # W/m3
            CompositeRheology=CompositeRheology((
                ConstantElasticity(; G=2.5e10),
                DislocationCreep(;
                # A=3.2e-4,
                # E=1.54e5,
                # V=12.0e-6,
                # n=2.3,
                ),
                DruckerPrager_regularised(;
                    C=1e7,
                    ϕ=asind(0.72),
                    η_vp=0)
            )),
        ),

        #       /Passive_margin_sediments
        #       /Sediments2_WET_QUARTZITE_RANALLI_1995
        #       3   1e+18 1e+26 1e+00 5e+29        1.97E+17      1.54E+05      0.80E+00      3.00E+04     2.30E+00  1.0E+10  1.00 1e+07 1e+07 0.350 0.350 0.5 1.5 0 0  2.60E+03      3.00E-05      1.00E-03      1.00E+03      0.64E+00   8.07E+02  4.00E-06  2.00E-06
        SetMaterialParams(; Name="Passive_margin_sediments2",
            Phase=7, # 9,#Making phases contiguous (28 jan, v0.76)
            Density=ConstantDensity(ρ=2600), # Can and should be expanded to PT_Density
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=0.64e0,
                b=8.07e2,
                d=4.00e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=2.00e-06), # W/m3
            CompositeRheology=CompositeRheology((
                # Elasticity
                ConstantElasticity(; G=1.0e10),

                # Dislocation creep
                DislocationCreep(;
                # A=1 / 1.97e17,
                # A=3.2e-4,
                # E=1.54e5,
                # V=8.0e-6,
                # n=2.3,
                ),

                # Drucker–Prager plasticity
                DruckerPrager_regularised(;
                    C=1e7,
                    ϕ=asind(0.35),
                    η_vp=0)
            )),
        ),
    )

end

function init_phases!(phases, phase_grid, particles, xvi)
    ni = size(phases)
    return @parallel (@idx ni) _init_phases!(phases, phase_grid, particles.coords, particles.index, xvi)
end

@parallel_indices (I...) function _init_phases!(phases, phase_grid, pcoords::NTuple{N,T}, index, xvi) where {N,T}

    ni = size(phases)

    for ip in cellaxes(phases)
        # quick escape
        @index(index[ip, I...]) == 0 && continue

        pᵢ = ntuple(Val(N)) do i
            @index pcoords[i][ip, I...]
        end

        d = Inf # distance to the nearest particle
        particle_phase = -1
        for offi in 0:1, offj in 0:1
            ii = I[1] + offi
            jj = I[2] + offj

            !(ii ≤ ni[1]) && continue
            !(jj ≤ ni[2]) && continue

            xvᵢ = (
                xvi[1][ii],
                xvi[2][jj],
            )
            d_ijk = √(sum((pᵢ[i] - xvᵢ[i])^2 for i in 1:N))
            if d_ijk < d
                d = d_ijk
                particle_phase = phase_grid[ii, jj]
            end
        end
        @index phases[ip, I...] = Float64(particle_phase)
    end

    return nothing
end
