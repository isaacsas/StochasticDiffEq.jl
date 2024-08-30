@muladd function perform_step!(integrator,cache::TauLeapingConstantCache)
  @unpack t,dt,uprev,u,W,p,P,c = integrator
  tmp = c(uprev, p, t, P.dW, nothing)
  integrator.u = uprev .+ tmp

  if integrator.opts.adaptive
    if integrator.alg isa TauLeaping
      oldrate = P.cache.currate
      newrate = P.cache.rate(integrator.u,p,t+dt)
      EEstcache = @. abs(newrate - oldrate) / max(50integrator.opts.reltol*oldrate,integrator.rate_constants/integrator.dt)
      integrator.EEst = maximum(EEstcache)
      if integrator.EEst <= 1
        P.cache.currate = newrate
      end
    elseif integrator.alg isa CaoTauLeaping
      # Calculate τ as EEst
    end
  end
end

@muladd function perform_step!(integrator,cache::TauLeapingCache)
  @unpack t,dt,uprev,u,W,p,P,c = integrator
  @unpack tmp, newrate, EEstcache = cache
  c(tmp, uprev, p, t, P.dW, nothing)
  @.. u = uprev + tmp

  if integrator.opts.adaptive
    if integrator.alg isa TauLeaping
      oldrate = P.cache.currate
      P.cache.rate(newrate,u,p,t+dt)
      @.. EEstcache =  abs(newrate - oldrate) / max(50integrator.opts.reltol*oldrate,integrator.rate_constants/integrator.dt)
      integrator.EEst = maximum(EEstcache)
      if integrator.EEst <= 1
        P.cache.currate .= newrate
      end
    elseif integrator.alg isa CaoTauLeaping
      # Calculate τ as EEst
    end
  end
end

############# split tau leaping ###########

@inline function concretize_affects!(c::SplitTauLeapingCache, 
        ::I) where {I <: SDEIntegrator}
    if (c.affects! isa Vector) && 
            !(c.affects! isa Vector{FunctionWrappers.FunctionWrapper{Nothing, Tuple{I}}})
        AffectWrapper = FunctionWrappers.FunctionWrapper{Nothing, Tuple{I}}
        c.affects! = AffectWrapper[JumpProcesses.makewrapper(AffectWrapper, aff) for aff in c.affects!]
    end
    nothing
end

@inline function initialize!(integrator, cache::SplitTauLeapingCache) 
    concretize_affects!(cache, integrator)
    nothing
end

# TEMPORARY, should move to JumpProcesses 
@inline function executerx_n_times!(speciesvec::AbstractVector{T}, rxidx::S,
        majump::M, n) where {T, S, M <: AbstractMassActionJump}
    @inbounds net_stoch = majump.net_stoch[rxidx]
    @inbounds for specstoch in net_stoch
        speciesvec[specstoch[1]] += n * specstoch[2]
    end
    nothing
end

@muladd function _perform_step!(integrator::SDEIntegrator{SplitTauLeaping}, 
        cache::SplitTauLeapingCache, affects!)
    @unpack t, dt, uprev, u, p = integrator
    @unpack tmp, newrate, majumps, rates, rng = cache
    @.. u = uprev

    @inbounds for i in 1:get_num_majumps(majumps)
        integrated_intensity = dt * evalrxrate(u, i, majumps)
        num_jumps = PoissonRandom.pois_rand(rng, integrated_intensity)
        executerx_n_times!(u, i, majumps, num_jumps)
    end

    @inbounds for (i, rate) in enumerate(rates)
        integrated_intensity = dt * rate(u, p, t)
        num_jumps = PoissonRandom.pois_rand(rng, integrated_intensity)
        
    end
end

@muladd function perform_step!(integrator::I, cache::SplitTauLeapingCache) where {
        I <: SDEIntegrator{SplitTauLeaping}}
    # for type stability we need to pull out the affects!
    affects! = cache.affects!
    if affects! isa Vector{FunctionWrappers.FunctionWrapper{Nothing, Tuple{I}}}
        _perform_step!(integrator, cache, affects!)
    else
        error("Error, invalid affects! type. Expected a vector of function wrappers and got $(typeof(affects!))")
    end
end
  