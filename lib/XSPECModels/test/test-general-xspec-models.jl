using Test
using SpectralFitting

# this model has been problematic so it gets its own explicit test case
model = XS_PhotoelectricAbsorption() * (XS_PowerLaw() + XS_Laor())

energy = collect(range(0.2, 8.0, 1000))
fluxes = zeros(Float64, (length(energy) - 1, objective_cache_count(model)))
flux = view(fluxes, :, 1)

# check model invokes okay
invokemodel!(fluxes, energy, model)
@test all(.!isnan.(flux))
@test all(.!isinf.(flux))
