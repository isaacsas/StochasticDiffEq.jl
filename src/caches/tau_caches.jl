struct TauLeapingConstantCache <: StochasticDiffEqConstantCache end

@cache struct TauLeapingCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    newrate::rateType
    EEstcache::rateType
end

function alg_cache(alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    TauLeapingConstantCache()
end

function alg_cache(alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    newrate = zero(jump_rate_prototype)
    EEstcache = zero(jump_rate_prototype)
    TauLeapingCache(u, uprev, tmp, newrate, EEstcache)
end

function alg_cache(alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    TauLeapingConstantCache()
end

function alg_cache(alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    TauLeapingCache(u, uprev, tmp, nothing, nothing)
end


@cache struct SplitTauLeapingCache{uType, tinvType, majType, ratefunType, 
        affectfunType, rngType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    majumps::majType
    ratefuns::ratefunType
    affects!::Any
    rng::rngType
end

function alg_cache(alg::SplitTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, jprob::JumpProblem) where {uEltypeNoUnits, 
        uBottomEltypeNoUnits, tTypeNoUnits}
    # TauLeapingConstantCache()
    error("Constant caches are not yet supported for SplitTauLeaping")
end

# note this should ultimately go in JumpProcesses as an API function, 
# this is temporary here while testing
function get_unpacked_jumps(jprob)
    maj = jprob.massaction_jump

    # ConstantRateJumps
    rates = jprob.discrete_jump_aggregation.rates
    affects! = jprob.discrete_jump_aggregation.affects!
    
    # VariableRateJumps -- these are not yet FunctionWrapped
    vrjs = jprob.variables_jumps
    if !isempty(vrjs)
        RateWrapper = eltype(rates)
        vrjrates = RateWrapper[RateWrapper(vrj.rate) for vrj in vrjs]
        vrjaffects = Any[(x -> (vrj.affect!(x); nothing)) for vrj in vrjs]
        append!(rates, vrjrates)
        append!(affects!, vrjaffects)
    end

    maj, rates, affects!
end

function alg_cache(alg::SplitTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, jprob::JumpProblem) where {uEltypeNoUnits, 
        uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)    
    maj, rates, affects! = get_unpacked_jumps(jprob)
    SplitTauLeapingCache(u, uprev, tmp, maj, rates, affects!, jprob.rng)
end