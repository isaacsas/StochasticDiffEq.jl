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

# @inline function initialize!(integrator, cache::SplitTauLeapingCache) 
#     nothing
# end

@muladd function perform_step!(integrator::I, cache::SplitTauLeapingCache) where {
        I <: SDEIntegrator{SplitTauLeaping}}

    @unpack t, dt, uprev, u, p = integrator
    @unpack majumps, ratefuns, affects!, rng = cache
    @.. u = uprev
    
    @inbounds for i in 1:get_num_majumps(majumps)
        integrated_intensity = dt * evalrxrate(u, i, majumps)
        num_jumps = PoissonRandom.pois_rand(rng, integrated_intensity)
        executerx!(u, i, majumps, num_jumps)
    end

    @inbounds for (i, rate) in enumerate(ratefuns)
        integrated_intensity = dt * rate(u, p, t)
        num_jumps = PoissonRandom.pois_rand(rng, integrated_intensity)
        affects![i](u, p, t, num_jumps)
    end

    nothing
end
  