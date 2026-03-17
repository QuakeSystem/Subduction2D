using GeoParams.Dislocation
using GeoParams.Diffusion
using GeoParams.MaterialParameters

#=
Below function added to combine all rheologies used in SZU2019 model. No need for modularity currently.

    1. Function has been edited to make phase numbers contiguous (28 jan, v0.76)
    2. Conductivity parameter d is currently not sufficient. Some error in conversion between SZU2019 and GeoParams / Gerya's values.
    3. Rheology names are kept the same as in SZU2019 for clarity.

=#
function init_rheologies()

    # common physical properties
    α = 2.4e-5 # 1 / K
    # Cp = 750    # J / kg K
    # Define rheology struct
    return rheology = (
        # SZU material 10
        #Mantle1_DRY_OL_Ranalli1995
        SetMaterialParams(; Name="Mantle1_DRY_0",
            Phase=1,
            Density = PT_Density(; ρ0 = 3.3e3, α = α, β = 0/kbar, T0 = 20C),
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=0.73e0,
                b=12.93e2,
                d=4e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=2.20E-08), # W/m3
            CompositeRheology=CompositeRheology((
                ConstantElasticity(; G=6.7e10),
                DislocationCreep(;
                A=2.513e-17,
                E=5.32e5,
                V=8.0e-6,
                n=3.5,
                ),
                DruckerPrager_regularised(;
                    C=1e7,
                    ϕ=asind(0.6),
                    η_vp=0
                )
            )),
        ),

        # SZU material 0
        #       /Air
        SetMaterialParams(;
            Name="Air",
            Phase=2,
            Density=ConstantDensity(; ρ=1.00E+02), # water density
            HeatCapacity=ConstantHeatCapacity(; Cp=3.0e3),
            Conductivity=ConstantConductivity(; k=1.0),
            CompositeRheology=CompositeRheology((LinearViscous(; η=1.00E+17),)),
        ),

        # SZU material 12
        #       /Mantle_Weak_zone
        #       Shear_Zone_Mantle4_WET_OL_Ranalli1995
        SetMaterialParams(; Name="Mantle_Weak_zone",
            Phase=3, #4, #Making phases contiguous (28 jan, v0.76)
            Density = PT_Density(; ρ0 = 3.3e3, α = α, β = 0/kbar, T0 = 20C),
            HeatCapacity=ConstantHeatCapacity(; Cp=1.00e3),
            Conductivity=TP_Conductivity(;
                a=0.73e0,
                b=12.93e2,
                d=4.00e-6 * 1e-5 * 1e-6,
            ),
            RadioactiveHeat=ConstantRadioactiveHeat(; H_r=2.20e-08), # W/m3
            CompositeRheology=CompositeRheology((
                ConstantElasticity(; G=6.7e10),
                DislocationCreep(;
                A=1.9960e-21,
                E=4.7e5,
                V=8.0e-6,
                n=4,
                ),
                DruckerPrager_regularised(;
                    C=2e6,#1e7, # edit v0.101
                    ϕ=asind(0.100),
                    η_vp=0)
            )),
        ),

        # SZU material 7
        #       Oceanic crust-interface material
        #       /Hydrated_fractured_top_oceanic_crust_Thrust_Interface
        #       /Basalts_WET_QUARTZITE_RANALLI_1995
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
                A=5.076142131979696e-18,
                E=1.54e5,
                V=8.0e-6,
                n=2.3,
                ),
                DruckerPrager_regularised(;
                    C=6e6,
                    ϕ=asind(0.025), # edit v0.101
                    η_vp=0)
            )),
        ),

        # SZU material 8
        #       /Oceanic_Crust_(Gabbro)
        #       /Basic_Crust2_An75_Ranalli1995
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
                A=2.0833333333333332e-23,
                E=2.38e5,
                V=8.0e-6,
                n=3.2,
                ),
                DruckerPrager_regularised(;
                    C=200e6,
                    ϕ=asind(0.85),
                    η_vp=0)
            )),
        ),

        # SZU material 5/6
        #       /Felsic_Crust1
        #       /Felsic_Crust1_WET_QUARTZITE_RANALLI_1995
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
                A=5.076142131979696e-18,
                E=1.54e5,
                V=12e-6,
                n=2.3,
                ),
                DruckerPrager_regularised(;
                    C=200e8,
                    ϕ=asind(0.72),
                    η_vp=0)
            )),
        ),

        # SZU material 3
        #       /Passive_margin_sediments
        #       /Sediments2_WET_QUARTZITE_RANALLI_1995
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
                ConstantElasticity(; G=1.0e10),
                DislocationCreep(;
                A=5.076142131979696e-18,
                E=1.54e5,
                V=8e-6,
                n=2.3,
                ),
                DruckerPrager_regularised(;
                    C=200e6,
                    ϕ=asind(0.35),
                    η_vp=0)
            )),),
        SetMaterialParams(; Name="left_boundary", # low viscosity boundary condition
            Phase=8,
            Density=ConstantDensity(; ρ=2600),
            HeatCapacity=ConstantHeatCapacity(; Cp=1.0e3),
            Conductivity=ConstantConductivity(; k=1.0),
            CompositeRheology=CompositeRheology((LinearViscous(; η=1.00E+20),))
        ))

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
