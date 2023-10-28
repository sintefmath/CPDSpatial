using Roots
using SpecialFunctions
import SparsityTracing
import ForwardDiff


@inline function get_value(v::V) where V
    return value(v)
end

@inline function get_value(v::SparsityTracing.ADval{V}) where V
    return v.val
end

@inline function get_value(
    v::Union{Vector{SparsityTracing.ADval{V}},
             Vector{Vector{SparsityTracing.ADval{V}}}}) where V
    return [get_value(x) for x in v]
end

function find_zero_wrapper(fun,
                           param::Union{V,
                                        Vector{V},
                                        Vector{Vector{V}}},
                           interval) where V <: Real

    # the default is just a wrapper around Roots.find_zero
    paramval = get_value(param)

    try
        result = find_zero(x->fun(x, paramval...), interval)
        return attach_gradient(result, param, fun)
    catch ex
        println("error occured: ", ex.msg) 
        rethrow() 
    end
end

function attach_gradient(result, param, fun)
    # in the default (non-AD) case, there is no gradient
    return result
end

function attach_gradient(result, param::Union{V, Vector{V}}, fun) where V <: SparsityTracing.ADval
    # This is just to support the use of sparsity tracing.  Note that the actual
    # derivative value is not important (just the sparsity pattern), so we
    # "borrow" the derivative information from the input value.

    return SparsityTracing.ADval(result, sum(param).derivnode)
end

function attach_gradient(result, param::Vector{Vector{SparsityTracing.ADval}}, fun)
    # This is just to support the use of sparsity tracing.  Note that the actual
    # derivative value is not important (just the sparsity pattern), so we
    # "borrow" the derivative information from the input value.

    tmp = sum((sum(v) for v in param)) # combine sparsity of all input parameters
    
    return SparsityTracing.ADval(result, tmp.derivnode)
end

function attach_gradient(result,
                         param::Union{AD,
                                      Vector{AD},
                                      Vector{Vector{AD}}
                                      },
                         fun) where AD <: ForwardDiff.Dual{T, V, L} where {T, V, L}
    # collect all parameter values (not derivatives) in a single vector
    pvecAD = reduce(vcat, param)
    if isa(pvecAD, Real)
        pvecAD = [pvecAD]
    end
    pvec = value(pvecAD)

    # make a function wrapper that takes each vector of the tuple as a separate
    # argument, but using only the numerical values, not derivatives.  
    regroup = (isa(param[1], Real)) ?
        (v) -> v :
        (v) -> [x for x in eachcol(reshape(v, :, length(param)))]

    fun2 = (arg) -> fun(arg[1], regroup(arg[2:end])...)

    # computing gradient with respect to 'result' as well as all the parameter values
    grad = ForwardDiff.gradient(fun2, [result, pvec...])

    # grouping into columns, where each columns contains the partial derivatives associated
    # with a given parameter group
    ∂ = -grad[2:end] ./ grad[1] # compute dy/dx = -∂f/∂x / ∂f/∂y

    # multiply with existing derivatives and return new object
    partials = sum(∂ .* [x.partials for x in pvecAD])
    ForwardDiff.Dual{T, V, L}(result, partials)
end
