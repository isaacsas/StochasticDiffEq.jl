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


@cache struct SplitTauLeapingCache{uType, majType, ratefunCollectType, 
        affectfunCollectType, rngType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    majumps::majType
    ratefuns::ratefunCollectType
    affects!::affectfunCollectType
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
function get_unpacked_jumps(jprob, u::utype, p::ptype, t::ttype) where {utype, ptype, ttype}
    maj = jprob.massaction_jump

    RateType = FunctionWrappers.FunctionWrapper{ttype, Tuple{utype, ptype, ttype}}
    AffectType = FunctionWrappers.FunctionWrapper{Nothing, Tuple{utype, ptype, ttype, Int}}
   
    # ConstantRateJumps
    crjs = jprob.constant_jumps
    if !isempty(crjs)
        rates = RateType[RateType(crj.rate) for crj in crjs]
        affects! = AffectType[AffectType((u,p,t,c) -> (crj.affect!(u,p,t,c); nothing)) for crj in crjs]
    else
        rates = RateType[]
        affects! = AffectType[]
    end
    
    # VariableRateJumps -- these are not yet FunctionWrapped
    vrjs = jprob.variable_jumps
    if !isempty(vrjs)
        vrjrates = RateType[RateType(vrj.rate) for vrj in vrjs]
        append!(rates, vrjrates)
        vrjaffects! = AffectType[AffectType((u,p,t,c) -> (vrj.affect!(u,p,t,c); nothing)) for vrj in vrjs]
        append!(affects!, vrjaffects!)
    end

    maj, rates, affects!
end

function alg_cache(alg::SplitTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, jprob::JumpProblem) where {uEltypeNoUnits, 
        uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)    
    maj, rates, affects! = get_unpacked_jumps(jprob, u, p, t)
    SplitTauLeapingCache{typeof(u),typeof(maj),typeof(rates),typeof(affects!),
        typeof(jprob.rng)}(u, uprev, tmp, maj, rates, affects!, jprob.rng)
end