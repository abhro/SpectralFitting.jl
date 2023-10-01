struct FittableMultiDataset{D}
    d::D
    FittableMultiDataset(data::Vararg{<:AbstractDataset}) = new{typeof(data)}(data)
end

function Base.getindex(multidata::FittableMultiDataset, i)
    multidata.d[i]
end

struct FittableMultiModel{M}
    m::M
    FittableMultiModel(model::Vararg{<:AbstractSpectralModel}) = new{typeof(model)}(model)
end

function Base.getindex(multimodel::FittableMultiModel, i)
    multimodel.m[i]
end

function _multi_constructor_wrapper(
    T::Union{<:Type{<:FittableMultiDataset},<:Type{<:FittableMultiModel}},
    args,
)
    if args isa T
        args
    else
        if args isa Tuple
            T(args...)
        else
            T(args)
        end
    end
end

"""
    _accumulated_indices(items)

`items` is a tuple or vector of lengths `n1, n2, ...`

Returns a tuple or array with same length as items, which gives the index boundaries of 
an array with size `n1 + n2 + ...`.
"""
function _accumulated_indices(items)
    total::Int = 0
    map(items) do item
        total = item + total
        total
    end
end

function _assemble_parameter_indices(bindings, n_params)
    remove = Int[]
    carry = Ref(0)
    parameter_indices = map((1:length(n_params)...,)) do i
        s = i == 1 ? 1 : n_params[i-1] + 1
        e = n_params[i]
        indices = collect(s:e)
        # assign bindings for all but the first models parameters
        if i > 1
            b = bindings[i]
            map(enumerate(indices)) do (j, index)
                # remove and replace index
                if j in b
                    push!(remove, index)
                    carry[] += 1
                    j
                else
                    # subtract the number of indices we've replaced / removed
                    index - carry[]
                end
            end
        else
            indices
        end
    end
    parameter_indices, remove
end


struct FittingProblem{M<:FittableMultiModel,D<:FittableMultiDataset,B}
    model::M
    data::D
    bindings::B
    function FittingProblem(m::FittableMultiModel, d::FittableMultiDataset)
        bindings = Vector{Pair{Int,Int}}[]
        new{typeof(m),typeof(d),typeof(bindings)}(m, d, bindings)
    end
end

FittingProblem(pairs::Vararg{<:Pair}) = FittingProblem(first.(pairs), last.(pairs))

function FittingProblem(m, d)
    _m = _multi_constructor_wrapper(FittableMultiModel, m)
    _d = _multi_constructor_wrapper(FittableMultiDataset, d)
    FittingProblem(_m, _d)
end

function model_count(prob::FittingProblem)
    return length(prob.model.m)
end

function data_count(prob::FittingProblem)
    return length(prob.data.d)
end

function _get_binding_indices(prob::FittingProblem, symbols::Vararg{Symbol})
    map(prob.model.m) do model
        param_symbols = all_parameter_symbols(model)
        map(symbols) do s
            i = findfirst(==(s), param_symbols)
            if isnothing(i)
                @warn "Model contains no symbol `$(Meta.quot(s))`"
                -1
            else
                i
            end
        end
    end
end

function old_bind!(prob::FittingProblem, symbols...)
    indices = _get_binding_indices(prob, symbols...)
    for (i, I) in enumerate(indices)
        if any(==(-1), I)
            j = findfirst(==(-1), I)
            throw("Could not bind symbol `$(Meta.quot(symbols[j]))`")
        end
        for index in I
            # avoid duplicating
            if index ∉ prob.bindings[i]
                push!(prob.bindings[i], index)
            end
        end
    end
    true
end


function Base.show(io::IO, ::MIME"text/plain", @nospecialize(prob::FittingProblem))
    buff = IOBuffer()
    println(buff, "FittingProblem:")
    println(buff, "  Models:")
    for model in prob.model.m
        println(buff, " "^4 * ". $(model)")
    end
    println(buff, "  Data:")
    for data in prob.data.d
        println(buff, " "^4 * ". $(data)")
    end
    print(io, encapsulate(String(take!(buff))))
end

export FittingProblem, bind!, FittableMultiModel, FittableMultiDataset
