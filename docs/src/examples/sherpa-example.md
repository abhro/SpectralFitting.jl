# A quick guide to modelling and fitting in SpectralFitting.jl

This is SpectralFitting.jl version of [A quick guide to modeling and fitting in Sherpa](https://sherpa.readthedocs.io/en/latest/quick.html).


```@example sherpa
using SpectralFitting, Plots
using Random
Random.seed!(0)

x = collect(range(-5, 5, 200))

A_true = 3.0
pos_true = 1.3
sigma_true = 0.8
err_true = 0.2

y = @. A_true * exp(-(x - pos_true)^2 / (2 * sigma_true^2))

y_noisy = y .+ (0.2 * randn(length(y)))

scatter(x, y_noisy)
```

To make this into a fittable dataset, we observe that our layout is injective (i.e. `length(x) == length(y)`). This is subtly different from how the majority of spectral models are implemented, which usually assume some kind of binning (`length(x) == length(y) + 1`). Fortunately, SpectralFitting.jl can track this for us, and do various conversion to make the models work correctly for the data. We need only tell the package what our [`AbstractDataLayout`](@ref) is:

```@example sherpa
data = InjectiveData(x, y_noisy; name = "example")
```

The data prints the _data card_, which provides us some high level information about our data at a glance. We can plot the data trivially using one of the Plots.jl recipes

```@example sherpa
plot(data, markersize = 3)
```

Next we want to specify a model to fit to this data. Models that are prefixed with `XS_` are models that are linked from the XSPEC model library, provided via [LibXSPEC_jll](https://github.com/astro-group-bristol/LibXSPEC_jll.jl). For a full list of the models, see [Models index](@ref).

!!! warning
    It is advised to **use the Julia implemented models**. This allows various calculations to benefit from automatic differentiation, efficient multi-threading, GPU offloading, and various other useful things, see [Why & how](@ref).


```@example sherpa
model = GaussianLine(μ = FitParam(0.0))
```

We can plot our model over the same domain range quite easily too:

```@example sherpa
plot(data.domain[1:end-1], invokemodel(data.domain, model))
```

Note that we've had to adjust the domain here. As stated before, most models are implemented for binned data, and therefore return one fewer bin than given.


SpectralFitting.jl adopts the SciML problem-solver abstraction, so to fit a model to data we specify a [`FittingProblem`](@ref):

```@example sherpa
prob = FittingProblem(model => data)
```

We fit problem then by calling [`fit`](@ref):

```@example sherpa
result = fit(prob, LevenbergMarquadt())
```

The result card tells us a little bit about how successful the fit was. We further inspect the fit by overplotting result on the data:

```@example sherpa
plot(data, markersize = 3)
plot!(result)
```

We can create a contour plot of the fit statistic by evaluating the result everywhere on the grid and measuring the statistic:

```@example sherpa
amps = range(50, 200, 50)
devs = range(0.5, 1.2, 50)

stats = [
    measure(ChiSquared(), result, [a, result.u[2], d])
    for d in devs, a in amps
]

# 1, 2, and 3 sigma contours
levels = [2.3, 4.61, 9.21]
contour(
    amps, 
    devs, 
    stats .- sum(result.stats), 
    levels = levels, 
    xlabel = "K", 
    ylabel = "σ"
)
scatter!([result.u[1]], [result.u[3]])
```