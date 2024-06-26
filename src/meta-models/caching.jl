mutable struct CacheEntry{T}
    cache::Vector{T}
    params::Vector{T}
    domain_limits::Tuple{T,T}
    size_of_element::Int

    function CacheEntry(params::AbstractVector{<:Number})
        T = eltype(params)
        cache = zeros(T, 1)
        new{T}(cache, params, (zero(T), zero(T)), sizeof(T))
    end
end

struct AutoCache{M,T,K,C<:CacheEntry} <: AbstractModelWrapper{M,T,K}
    model::M
    cache::C
    abstol::Float64
    function AutoCache(
        model::AbstractSpectralModel{T,K},
        cache::CacheEntry,
        abstol,
    ) where {T,K}
        new{typeof(model),T,K,typeof(cache)}(model, cache, abstol)
    end
end

function AutoCache(model::AbstractSpectralModel{T,K}; abstol = 1e-3) where {T,K}
    params = [get_value.(parameter_tuple(model))...]
    cache = CacheEntry(params)
    AutoCache(model, cache, abstol)
end

function _reinterpret_dual(::Type, v::AbstractArray, n::Int)
    needs_resize = n > length(v)
    if needs_resize
        @warn "AutoCache: Growing dual buffer..."
        resize!(v, n)
    end
    view(v, 1:n), needs_resize
end
function _reinterpret_dual(
    DualType::Type{<:ForwardDiff.Dual},
    v::AbstractArray{T},
    n::Int,
) where {T}
    n_elems = div(sizeof(DualType), sizeof(T)) * n
    needs_resize = n_elems > length(v)
    if needs_resize
        @warn "AutoCache: Growing dual buffer..."
        resize!(v, n_elems)
    end
    reinterpret(DualType, view(v, 1:n_elems)), needs_resize
end

function invoke!(output, domain, model::AutoCache{M,T}) where {M,T}
    D = promote_type(eltype(domain), T)

    _new_params = parameter_tuple(model.model)
    _new_limits = (first(domain), last(domain))

    output_cache, out_resized = _reinterpret_dual(D, model.cache.cache, length(output))
    param_cache, _ = _reinterpret_dual(D, model.cache.params, length(_new_params))

    same_domain = model.cache.domain_limits == _new_limits

    # if the parameter size has changed, need to rerun the model
    if (!out_resized) && (model.cache.size_of_element == sizeof(D)) && (same_domain)
        # if all parameters within some tolerance, then just return the cache
        within_tolerance = all(zip(param_cache, _new_params)) do I
            p, pm = I
            abs((get_value(p) - get_value(pm)) / p) < model.abstol
        end

        if within_tolerance
            @. output = output_cache
            return output
        end
    end

    model.cache.size_of_element = sizeof(D)
    invoke!(output_cache, domain, model.model)
    # update the auto cache infos
    model.cache.domain_limits = _new_limits

    @. param_cache = get_value(_new_params)
    @. output = output_cache
end

export AutoCache
