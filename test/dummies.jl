using SparseArrays

# dummy datasets
function make_dummy_dataset(
    shape_function;
    energy = collect(range(0.2, 10.0, 101)),
    kwargs...,
)
    # calculate flux from shape function
    flux = shape_function.(energy[1:end-1])
    make_dummy_dataset(flux, energy; kwargs...)
end

function make_dummy_dataset(
    flux,
    energy;
    exposure_time = 2.0,
    units = "counts",
    error_fraction = 0.1,
)
    bins_low = energy[1:end-1]
    bins_high = energy[2:end]
    N = length(flux)

    channels = collect(1:N)

    spectrum = Spectrum(
        channels,
        zeros(Int, N),
        ones(Int, N),
        flux,
        units,
        exposure_time,
        1.0, # background scale
        1.0, # area scale
        SpectralFitting.ErrorStatistics.Unknown,
        error_fraction .* flux,
        0.0,
        "test-telescope",
        "test-name",
    )

    matrix = SpectralFitting.sparse(Matrix(1.0 * SpectralFitting.I, N, N))
    response = ResponseMatrix(
        matrix,
        spectrum.channels,
        eltype(bins_low).(channels),
        eltype(bins_low).(push!(channels[2:end], N + 1)),
        bins_low,
        bins_high,
    )

    SpectralData(spectrum, response)
end

# standard julia models for testing
struct DummyAdditive{T,F} <: AbstractSpectralModel{T,Additive}
    K::T
    a::T
    b::T
end
function DummyAdditive(; K = FitParam(1.0), a = FitParam(1.0), b = FitParam(5.0))
    DummyAdditive{typeof(K),SpectralFitting.FreeParameters{(:K, :a)}}(K, a, b)
end
function SpectralFitting.invoke!(flux, energy, model::DummyAdditive)
    let a = model.a, b = model.b
        @. flux = a + b
    end
end

struct DummyMultiplicative{T,F} <: AbstractSpectralModel{T,Multiplicative}
    a::T
    b::T
end
function DummyMultiplicative(; a = FitParam(1.0), b = FitParam(5.0))
    DummyMultiplicative{typeof(a),SpectralFitting.FreeParameters{(:a,)}}(a, b)
end
function SpectralFitting.invoke!(flux, energy, model::DummyMultiplicative)
    let a = model.a, b = model.b
        @. flux = a * b
    end
end

# table models for testing

struct DummyAdditiveTableModel{D,T,F} <: AbstractTableModel{T,Additive}
    table::D
    K::T
    a::T
    b::T
end
function DummyAdditiveTableModel(K::T, a::T, b::T) where {T}
    # table is just an interpolation anywhere so lambda for tests
    table = (x) -> x^2
    DummyAdditiveTableModel{typeof(table),T,SpectralFitting.FreeParameters{(:K, :a)}}(
        table,
        K,
        a,
        b,
    )
end
function SpectralFitting.invoke!(flux, energy, model::DummyAdditiveTableModel)
    let a = model.a, b = model.b, table = model.table
        flux[:] .= table(a) + b
    end
end
function DummyAdditiveTableModel(; K = FitParam(1.0), a = FitParam(1.0), b = FitParam(2.0))
    DummyAdditiveTableModel(K, a, b)
end

struct DummyMultiplicativeTableModel{D,T,F} <: AbstractTableModel{T,Multiplicative}
    table::D
    a::T
    b::T
end
function DummyMultiplicativeTableModel(a::T, b::T) where {T}
    # table is just an interpolation anywhere so lambda for tests
    table = (x, k) -> k * x
    DummyMultiplicativeTableModel{typeof(table),T,SpectralFitting.FreeParameters{(:a,)}}(
        table,
        a,
        b,
    )
end
function SpectralFitting.invoke!(flux, energy, model::DummyMultiplicativeTableModel)
    let a = model.a, b = model.b, table = model.table
        @. flux = table(flux, a) + b
    end
end
function DummyMultiplicativeTableModel(; a = FitParam(1.0), b = FitParam(2.0))
    DummyMultiplicativeTableModel(a, b)
end

# standard julia models for testing
struct DummyAdditiveWithManyFrozen{T,F} <: AbstractSpectralModel{T,Additive}
    K::T
    a::T
    b::T
    c::T
    d::T
    e::T
    f::T
    g::T
    h::T
    i::T
    j::T
end
function DummyAdditiveWithManyFrozen(;
    K = FitParam(1.0),
    a = FitParam(1.0),
    b = FitParam(5.0),
    c = FitParam(2.0),
    d = FitParam(2.0),
    e = FitParam(2.0),
    f = FitParam(2.0),
    g = FitParam(2.0),
    h = FitParam(2.0),
    i = FitParam(2.0),
    j = FitParam(2.0),
)
    DummyAdditiveWithManyFrozen{typeof(K),SpectralFitting.FreeParameters{(:K, :h)}}(
        K,
        a,
        b,
        c,
        d,
        e,
        f,
        g,
        h,
        i,
        j,
    )
end
function SpectralFitting.invoke!(flux, energy, model::DummyAdditiveWithManyFrozen)
    for i in eachindex(flux)
        flux[i] =
            model.a * 2model.b +
            model.c +
            model.d +
            model.e +
            model.f +
            model.g +
            model.h +
            model.i +
            model.j
    end
end
