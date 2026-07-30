using GeoParams

dp = GeoParams.DruckerPrager_regularised()
τxx, τyy, τxy = 2.0e6, -1.3e6, 0.7e6
τij = (τxx, τyy, τxy)

CR = GeoParams.MaterialParameters.ConstitutiveRelationships
# @which CR.∂Q∂τxx(dp, τij)
τII = CR.second_invariant(τij)

∂xx = CR.∂Q∂τxx(dp, τij)
∂yy = CR.∂Q∂τyy(dp, τij)

implicit_dτzz_per_unit  = 2 * (∂xx + ∂yy)          # coefficient of η_ve*λ in implicit τzz update
explicit_dτzz_per_unit  = (τxx + τyy) / τII         # coefficient of η_ve*λ in correct 3D update

println("∂Q∂τxx = ", ∂xx)
println("∂Q∂τyy = ", ∂yy)
println("implicit factor = ", implicit_dτzz_per_unit)
println("explicit factor = ", explicit_dτzz_per_unit)
println("ratio = ", implicit_dτzz_per_unit / explicit_dτzz_per_unit, "  (should be 1.0)")