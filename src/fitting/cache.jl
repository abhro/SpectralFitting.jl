abstract type AbstractFittingCache end

_invoke_and_transform!(cache::AbstractFittingCache, domain, params) =
    error("Not implemented for $(typeof(cache))")

# one of these for each (mulit)model / data pair
struct SpectralCache{M,O,T,K,P,TransformerType} <: AbstractFittingCache
    model::M
    model_output::O
    calculated_objective::T
    output_cache::K
    parameter_cache::P
    transformer!!::TransformerType
    function SpectralCache(
        layout::AbstractLayout,
        model::M,
        domain,
        objective,
        transformer::XfmT;
        param_diff_cache_size = nothing,
    ) where {M,XfmT}
        model_output = DiffCache(construct_objective_cache(layout, model, domain))
        # fix for https://github.com/fjebaker/SpectralFitting.jl/issues/79
        # output must be a vector but can only give matrix to `mul!`, so we need to
        # unfortunately duplicate the array to ensure we have both types
        calc_obj = zeros(eltype(objective), (length(objective), 1))
        calc_obj_cache = DiffCache(calc_obj)
        # vector chache
        output = similar(objective)
        output .= 0
        output_cache = DiffCache(output)
        param_cache =
            make_diff_parameter_cache(model; param_diff_cache_size = param_diff_cache_size)
        new{
            M,
            typeof(model_output),
            typeof(calc_obj_cache),
            typeof(output_cache),
            typeof(param_cache),
            XfmT,
        }(
            model,
            model_output,
            calc_obj_cache,
            output_cache,
            param_cache,
            transformer,
        )
    end
end

function Base.show(io::IO, @nospecialize(config::SpectralCache{M})) where {M}
    descr = "SpectralCache{$(Base.typename(M).name)}"
    print(io, descr)
end

function Base.show(
    io::IO,
    ::MIME"text/plain",
    @nospecialize(config::SpectralCache{M})
) where {M}
    descr = "SpectralCache{$(Base.typename(M).name)}"
    print(io, descr)
end

function _invoke_and_transform!(cache::SpectralCache, domain, params)
    # read all caches
    model_output = get_tmp(cache.model_output, params)
    calc_obj = get_tmp(cache.calculated_objective, params)

    # update the free parameters, and then get all of them
    update_free_parameters!(cache.parameter_cache, params)
    parameters = _get_parameters(cache.parameter_cache, params)

    output = invokemodel!(model_output, domain, cache.model, parameters)
    cache.transformer!!(calc_obj, domain, output)

    output_vector = get_tmp(cache.output_cache, params)
    output_vector .= calc_obj
    output_vector
end

struct MultiModelCache{K,N,CacheTypes<:Tuple,ParameterMappingType} <: AbstractFittingCache
    caches::CacheTypes
    all_outputs::K
    domain_mapping::NTuple{N,Int}
    output_domain_mapping::NTuple{N,Int}
    objective_mapping::NTuple{N,Int}
    parameter_mapping::ParameterMappingType
end

function _get_range(mapping::NTuple, i)
    m_start = i == 1 ? 1 : mapping[i-1] + 1
    m_end = mapping[i]
    (m_start, m_end)
end

function _invoke_and_transform!(cache::MultiModelCache, domain, params)
    all_outputs = get_tmp(cache.all_outputs, params)

    for (i, ch) in enumerate(cache.caches)
        p = @views params[cache.parameter_mapping[i]]

        domain_start, domain_end = _get_range(cache.domain_mapping, i)
        objective_start, objective_end = _get_range(cache.objective_mapping, i)

        d = @views domain[domain_start:domain_end]
        all_outputs[objective_start:objective_end] .= _invoke_and_transform!(ch, d, p)
    end

    all_outputs
end

function _build_parameter_mapping(model::FittableMultiModel, bindings)
    parameters = map(m -> collect(filter(isfree, parameter_tuple(m))), model.m)
    parameters_counts = _accumulated_indices(map(length, parameters))

    all_parameters = reduce(vcat, parameters)

    parameter_mapping, remove = _construct_bound_mapping(bindings, parameters_counts)
    # remove duplicate parameters that are bound
    deleteat!(all_parameters, remove)

    all_parameters, parameter_mapping
end

function _build_mapping_length(f, itt::Tuple)
    values = map(f, itt)
    mapping = _accumulated_indices(map(length, values))
    values, mapping
end

_build_objective_mapping(layout::AbstractLayout, dataset::FittableMultiDataset) =
    _build_mapping_length(i -> make_objective(layout, i), dataset.d)
_build_domain_mapping(layout::AbstractLayout, dataset::FittableMultiDataset) =
    _build_mapping_length(i -> make_model_domain(layout, i), dataset.d)
_build_output_domain_mapping(layout::AbstractLayout, dataset::FittableMultiDataset) =
    _build_mapping_length(i -> make_output_domain(layout, i), dataset.d)
