mutable struct Spectrum{T} <: AbstractDataset
    channels::Vector{Int}
    quality::Vector{Int}
    grouping::Vector{Int}

    # this could be counts or flux
    data::Vector{T}
    unit_string::String

    exposure_time::T
    background_scale::T
    area_scale::T

    error_statistics::SpectralFitting.ErrorStatistics.T
    errors::Union{Missing,Vector{T}}
    systematic_error::T

    telescope_name::String
    instrument::String
end

function normalize!(spectrum::Spectrum)
    if spectrum.unit_string == "counts"
        @. spectrum.data /= spectrum.exposure_time
        @. spectrum.errors /= spectrum.exposure_time
        spectrum.unit_string = "count / s"
    end
    spectrum
end

supports_contiguosly_binned(::Type{<:Spectrum}) = true

function make_objective(::ContiguouslyBinned, dataset::Spectrum)
    if dataset.unit_string == "counts"
        @warn "Spectrum is currently still in counts. Most models fit in rate (count / s). Use `normalize!(dataset)` to ensure the dataset is in a standard format."
    end
    dataset.data
end

function make_objective_variance(::ContiguouslyBinned, dataset::Spectrum)
    if dataset.unit_string == "counts"
        @warn "Spectrum is currently still in counts. Most models fit in rate (count / s). Use `normalize!(dataset)` to ensure the dataset is in a standard format."
    end
    @. dataset.errors^2
end

function make_domain(::ContiguouslyBinned, dataset::Spectrum)
    @warn "Spectrum doesn't know the energy values by default. Domain is channels. Proceed only if you know what you are doing."
    dataset.channels
end

isgrouped(spectrum::Spectrum) = all(==(1), spectrum.grouping)

regroup!(spectrum::Spectrum) = regroup!(spectrum, spectrum.grouping)
function regroup!(spectrum::Spectrum{T}, grouping) where {T}
    if all(==(1), grouping)
        # no op
        return spectrum
    end

    itt = GroupingIterator(grouping)
    for grp in itt
        spectrum.channels[grp[1]] = spectrum.channels[grp[2]]
        regroup_vector!(spectrum.data, grp)
        regroup_quality_vector!(spectrum.quality, grp)
        if !ismissing(spectrum.errors)
            vs = spectrum.data[grp[1]]
            if spectrum.unit_string == "counts"
                spectrum.errors[grp[1]] = count_error(vs, 1.0)
            elseif spectrum.unit_string == "count / s"
                vs = vs * spectrum.exposure_time
                es = count_error(vs, 1.0)
                spectrum.errors[grp[1]] = es / spectrum.exposure_time
            else
                error(
                    "No method for grouping errors with given spectral units ($(spectrum.unit_string)).",
                )
            end
        end
    end

    resize!(spectrum, length(itt))
    spectrum.grouping .= 1
    spectrum
end

function Base.resize!(spectrum::Spectrum, n::Int)
    resize!(spectrum.channels, n)
    resize!(spectrum.data, n)
    resize!(spectrum.grouping, n)
    resize!(spectrum.quality, n)
    if !isnothing(spectrum.errors)
        resize!(spectrum.errors, n)
    end
end

function drop_channels!(spectrum::Spectrum, indices)
    deleteat!(spectrum.channels, indices)
    deleteat!(spectrum.data, indices)
    deleteat!(spectrum.quality, indices)
    deleteat!(spectrum.grouping, indices)
    if !ismissing(spectrum.errors)
        deleteat!(spectrum.errors, indices)
    end
    length(indices)
end

_readable_boolean(b) = b ? "yes" : "no"

function _printinfo(io, spectrum::Spectrum)
    dmin, dmax = prettyfloat.(extrema(spectrum.data))
    is_grouped = isgrouped(spectrum) |> _readable_boolean
    num_bad = count(==(BAD_QUALITY), spectrum.quality)
    has_bad = num_bad > 0 ? "yes ($num_bad)" : "no"
    descr = """Spectrum: $(spectrum.telescope_name)[$(spectrum.instrument)]
      Units                 : $(spectrum.unit_string)
      . Exposure time       : $(spectrum.exposure_time)
      . Channels            : $(length(spectrum.channels))
      . Data (min/max)      : ($dmin, $dmax)
      . Grouped             : $is_grouped
      . Bad channels        : $has_bad
    """
    print(io, descr)
end

export Spectrum
