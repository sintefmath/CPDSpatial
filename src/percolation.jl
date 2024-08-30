# Copyright (c) 2023 Odd Andersen
# Copyright (c) 2023 SINTEF Digital

using Roots
using SpecialFunctions

export f_gas, f_tar, binned_molecular_weights

"""
   p⃰(p, σ)

Compute the `p⋆` value that is needed in the formula to compute the total
fraction of sites, F(p) contained in all the finite clusters.  

`p⃰` is defined as the root of the following equation in p:
```
        p⃰(1-p⃰)^(σ-1) = p(1-p)^(σ-1)
```
            
For `p <= 1/σ`, the solution is the trivial one of `p⃰ = p`.  Otherwise, a
root-finding function must be invoked.

# Arguments
- `p::Real`: Probability of a given bridge being intact. Should be in [0, 1].  
             Can be a vector with multiple values.
- `σ::Real`: (σ+1) is the coordination number, i.e. average number of bridges
             per site.  Should be a number >= 2.

# Returns
- `p⃰val:Real`: Return value, which will be in [0, 1/σ].  If `p` is a vector,
                the return value will also be.
"""
function p⃰(p::Real, σ)
    relation = function(p⃰, q)
        p⃰ * (1 - p⃰)^(σ-1) - q * (1-q)^(σ-1)
    end

    # if p is less than the threshold value, we already know the (trivial)
    # answer, and can avoid the call to the root-finding function
    (p <= 1/σ) ? p :
    (p >= 1.0 - eps()) ? 0.0 : # in reality, should never be larger than 1.0
                         find_zero_wrapper(relation, p, (eps()^2, 1/σ - eps()))
end

function p⃰(p::Vector, σ)
    return [p⃰(x, σ) for x in p]
end

# ----------------------------------------------------------------------------
"""
    F(p, σ)

Compute the total fraction of sites belonging to finite clusters.

# Arguments
- `p::Real`: Probability of a given bridge being intact. Should be in [0, 1].
- `σ::Real`: (σ+1) is the coordination number, i.e. average number of bridges
             per site.  Should be a number >= 2.
- `bins::Int64`: Number of bins to use for the finite clusters.  If 0, the
                 aggregate result is returned.  Otherwise, the result is
                 returned as a vector of length `bins`.
# Returns
- `Fval`: Return value, if `bins` equals 0, the aggregate result is returned 
          (i.e. the total fraction of sites belonging to finite clusters).  
          Otherwise, the result is returned as a vector of length `bins`, where 
          bin 'i' represents the fraction of sites in clusters of size i.

          Note that in case the total fraction is returned (bins==0), the 
          p⃰ value is also returned as a second return value.
"""
function F(p, σ; bins=0)
    @assert all(value(p) .<= 1.0)
    if bins == 0
        # compute the total fraction of sites belonging to finite clusters
        p⃰val = p⃰(p, σ) # store computed p⃰ value so that it can be returned
        # together with F.
        return ( p⃰val ./ p ).^( (σ + 1) / (σ - 1)), p⃰val
    end

    # If we got here, we want to return the value distributed over cluster sizes.

    # define generalized binomial function that can take non-integer arguments.
    # This implementation is good for η up to slighly above 40.
    bicoef = (η, μ) -> gamma.(η+1) ./ (gamma.(μ+1) .* gamma.(η-μ+1))

    nrange = 1:bins
    nbn = (σ+1) ./ (nrange .* σ .+ 1) .* bicoef.(nrange .* σ .+ 1, nrange .- 1)

    return nbn .* p .^ (nrange .- 1) .* (1 - p) .^ (nrange .* (σ - 1) .+ 2)
    
end

# ----------------------------------------------------------------------------
"""
    K(p, σ)

Referred to as the 'configuration generating function'.  It represents the 
total number of finite clusters, regardless of size.

# Arguments
- `p::Real`: Probability of a given bridge being intact. Should be in [0, 1].
- `σ::Real`: (σ+1) is the coordination number, i.e. average number of bridges
             per site.  Should be a number >= 2.

# Returns
- `Kval:Real`: Return value, the total number of finite clusters (expressed 
               as a fraction relative to the total number of sites available).
"""
function K(p, σ)

    Fval, p⃰val = F(p, σ)
    
    (1 .- (σ .+ 1) ./ 2 .* p⃰val) .* Fval
    
end
# ----------------------------------------------------------------------------
"""
    f_gas(r, g, σ, c0)

Compute the mass fraction of gas.

# Arguments
- `r::Float64`: The ratio of bridge mass `mb` to site mass `ma`.
- `g::Float64`: Number of light gas molecules (on a per-site basis).
- `σ::Float64`: (σ+1) is the coordination number.
- `c0::Float64`: Initial number of char bridges (per-site basis).

# Returns
- `result::Float64`: Mass fraction of gas.
"""
function f_gas(r, g, σ, c0)

    return  ( r * g * (σ + 1) ) /
            ( 4 + 2 * r * (1 - c0) * (σ + 1))
end

# ----------------------------------------------------------------------------
"""
    f_tar(r, p, σ, c0, δ, £)

Compute the mass fraction of tar.

# Arguments
- `r::Float64`: The ratio of bridge mass `mb` to site mass `ma`.
- `p::Float64`: Ratio of intact bridges (compared to total)
- `σ::Float64`: (σ+1) is the coordination number.
- `c0::Float64`: Initial number of char bridges (per-site basis).
- `δ::Float64`: Number of side chains, relative to the theoretical maximum number
                of labile bridges as defined by the coordination number (σ+1).
- `£::Float64`: Current number of labile bridges (relative to theoretical max).

# Keyword arguments
- `bins::Int64`: Number of bins to use for the tar distribution.  If 0, the
                 aggregate result is returned.  Otherwise, the result is
                 returned as a vector of length `bins`.

- `consistent_with_original_code::Bool`: if working with binned values, setting
                                         this value to `true` aims for
                                         consistent treatment with original
                                         cpdheat code, but the bins will not add
                                         up to the theoretical total.  If set to
                                         `false`, the bins will be enforced to
                                         add up to the total. (Default: true)
# Returns
- `result`: if `bins == 0`, the aggregate result is returned.  Otherwise,
            the result is returned as a vector of length `bins`, where 
            bin 'i' represents the fraction of tar in clusters of size i.
"""
function f_tar(r, p, σ, c0, δ, £; bins=0, consistent_with_original_code=true)
    @assert all(value(p) .<= 1.0)
    if p == 1.0
        # all bridges are intact, no tar possible
        return (bins <= 1) ? 0.0 : zeros(bins)
    elseif p == 0.0
        # no intact bridges, everything is turned into monomers
        return (bins <= 1) ? 1.0 : [1.0; zeros(bins-1)]
    end
    
    
    # local helper functions
    Φ = function(p, r, £, σ, δ)
        1 .+ r * (£./p .+ ((σ-1) * δ) ./ (4 * (1 .- p)))
    end

    Ω = function(p, £, δ)
        δ ./ (2 * (1 .- p)) .- £ ./ p
    end

    # intermediary factor
    fac = 2 / (2 + r * (1 - c0) * (σ + 1))

    total = fac * (Φ(p, r, £, σ, δ) .* F(p, σ)[1] .+ r * Ω(p, £, δ) .* K(p, σ))
    
    if bins <= 1
        # we only want the aggregate result, for which there is a closed form
        return bins == 0 ? total : [total]
    end

    # if we got here, we want the tar disaggregated into bins by cluster size

    if consistent_with_original_code
        # ensure treatment close to the originally published cpd code
        Fn = F(p, σ, bins = bins)
        Qn = Fn ./ collect(1:bins)

        res = fac * ( Φ(p, r, £, σ, δ) .* Fn .+ r * Ω(p, £, δ) .* Qn )
    else
        # ensure the total matches, but at the cost of diverging from the originally
        # published cpdheat code.
        Fn = F(p, σ, bins = bins-1)
        Qn = Fn ./ collect(1:bins-1)

        res = fac * ( Φ(p, r, £, σ, δ) .* Fn .+ r * Ω(p, £, δ) .* Qn )
        push!(res, total-sum(res))
    end

    


    #res[1:bins-1] = fac * ( Φ(p, r, £, σ, δ) .* Fn .+ r * Ω(p, £, δ) .* Qn )
    #res[bins] = total - sum(res[1:bins-1])
    
    return res
end

# ----------------------------------------------------------------------------

"""
    binned_molecular_weights(r, p, σ, c0, δ, £)
    
Compute the molecular weight of tar 1-mers, 2-mers up to n-mers.
Unit: same unit as used for `ma` in the input.

# Arguments
- `ma::Float64`: Mean weight of a single site.
- `r::Float64`: The ratio of bridge mass `mb` to mean site mass `ma`.
- `p::Float64`: Ratio of intact bridges (compared to total)
- `σ::Float64`: (σ+1) is the coordination number.
- `δ::Float64`: Number of side chains, relative to the theoretical maximum number
                of labile bridges as defined by the coordination number (σ+1). 
- `£::Float64`: Current number of labile bridges (relative to theoretical max).
- `bins::Int64`: Number of bins to use for the tar distribution. 
"""                 
function binned_molecular_weights(ma, r, p, σ, δ, £, bins)

    EPS = 1e-9

    # intermediary variables
    f1 = (p < EPS)   ? 0 : £ / p
    f2 = (1-p < EPS) ? 0 : δ / (1-p)
    
    nvec = collect(1:bins)

    tn = nvec .* (σ .- 1) .+ 2

    # return value: vector with molecular weights for each bin
    res = nvec .* ma +
        (nvec .- 1) .* r .* ma .* f1 +
        tn .* r .* ma ./ 4 .* f2

end

