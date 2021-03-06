################################################################################
# Roessler SRK for second order weak approx

struct DRI1ConstantCache{T,T2} <: StochasticDiffEqConstantCache
  # hard-coded version
  a021::T
  a031::T
  a032::T

  a121::T
  a131::T
  #a132::T

  #a221::T
  #a231::T
  #a232::T

  b021::T
  b031::T
  #b032::T

  b121::T
  b131::T
  #b132::T

  b221::T
  b222::T
  b223::T
  b231::T
  b232::T
  b233::T

  α1::T
  α2::T
  α3::T

  c02::T2
  c03::T2

  #c11::T2
  c12::T2
  c13::T2

  #c21::T2
  #c22::T2
  #c23::T2

  beta11::T
  beta12::T
  beta13::T

  #beta21::T
  beta22::T
  beta23::T

  beta31::T
  beta32::T
  beta33::T

  #beta41::T
  beta42::T
  beta43::T

  #quantile(Normal(),1/6)
  NORMAL_ONESIX_QUANTILE::T
end

function DRI1ConstantCache(::Type{T}, ::Type{T2}) where {T,T2}

  a021 = convert(T, 1//2)
  a031 = convert(T, -1)
  a032 = convert(T, 2)

  a121 = convert(T, 342//491)
  a131 = convert(T, 342//491)

  b021 = convert(T, (6-sqrt(6))/10)
  b031 = convert(T, (3+2*sqrt(6))/5)

  b121 = convert(T, 3*sqrt(38//491))
  b131 = convert(T, -3*sqrt(38/491))

  b221 = convert(T, -214//513*sqrt(1105//991))
  b222 = convert(T, -491//513*sqrt(221//4955))
  b223 = convert(T, -491//513*sqrt(221//4955))
  b231 = convert(T, 214//513*sqrt(1105//991))
  b232 = convert(T, 491//513*sqrt(221//4955))
  b233 = convert(T, 491//513*sqrt(221//4955))

  α1 = convert(T, 1//6)
  α2 = convert(T, 2//3)
  α3 = convert(T, 1//6)

  c02 = convert(T2, 1//2)
  c03 = convert(T2, 1)

  c12 = convert(T2, 342//491)
  c13 = convert(T2, 342//491)

  beta11 = convert(T, 193//684)
  beta12 = convert(T, 491//1368)
  beta13 = convert(T, 491//1368)

  beta22 = convert(T, 1//6*sqrt(491//38))
  beta23 = convert(T, -1//6*sqrt(491//38))

  beta31 = convert(T, -4955//7072)
  beta32 = convert(T, 4955//14144)
  beta33 = convert(T, 4955//14144)

  beta42 = convert(T, -1//8*sqrt(4955//221))
  beta43 = convert(T, 1//8*sqrt(4955//221))

  NORMAL_ONESIX_QUANTILE = convert(T, -0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::DRI1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  DRI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end



function RI1ConstantCache(T::Type, T2::Type)
  a021 = convert(T, 2//3) # convert(T, 2//3)
  a031 = convert(T, -1//3)
  a032 = convert(T, 1)

  a121 = convert(T, 1)
  a131 = convert(T, 1)

  b021 = convert(T, 1)
  b031 = convert(T, 0)

  b121 = convert(T, 1)
  b131 = convert(T, -1)

  b221 = convert(T, 1)
  b222 = convert(T, 0)
  b223 = convert(T, 0)
  b231 = convert(T, -1)
  b232 = convert(T, 0)
  b233 = convert(T, 0)

  α1 = convert(T, 1//4)
  α2 = convert(T, 1//2)
  α3 = convert(T, 1//4)

  c02 = convert(T2, 2//3)
  c03 = convert(T2, 2//3)

  c12 = convert(T2, 1)
  c13 = convert(T2, 1)

  beta11 = convert(T, 1//2)
  beta12 = convert(T, 1//4)
  beta13 = convert(T, 1//4)

  beta22 = convert(T, 1//2)
  beta23 = convert(T, -1//2)

  beta31 = convert(T, -1//2)
  beta32 = convert(T, 1//4)
  beta33 = convert(T, 1//4)

  beta42 = convert(T, 1//2)
  beta43 = convert(T, -1//2)

  NORMAL_ONESIX_QUANTILE = convert(T,-0.9674215661017014)

  DRI1ConstantCache(a021,a031,a032,a121,a131,b021,b031,b121,b131,b221,b222,b223,b231,b232,b233,α1,α2,α3,c02,c03,c12,c13,beta11,beta12,beta13,beta22,beta23,beta31,beta32,beta33,beta42,beta43,NORMAL_ONESIX_QUANTILE)
end


function alg_cache(alg::RI1,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{false}})
  RI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
end


@cache struct DRI1Cache{uType,randType,MType1,tabType,rateNoiseType,rateType,possibleRateType} <: StochasticDiffEqMutableCache
  u::uType
  uprev::uType

  _dW::randType
  _dZ::randType
  chi1::randType
  Ihat2::MType1

  tab::tabType

  g1::rateNoiseType
  g2::Vector{rateNoiseType}
  g3::Vector{rateNoiseType}

  k1::rateType
  k2::rateType
  k3::rateType

  H02::uType
  H03::uType
  H12::Vector{uType}
  H13::Vector{uType}
  H22::Vector{uType}
  H23::Vector{uType}

  tmp1::possibleRateType
  tmpg::rateNoiseType

end


function alg_cache(alg::DRI1,prob,u,ΔW,ΔZ,p,rate_prototype,
                   noise_rate_prototype,jump_rate_prototype,uEltypeNoUnits,
                   uBottomEltypeNoUnits,tTypeNoUnits,uprev,f,t,dt,::Type{Val{true}})
  if typeof(ΔW) <: Union{SArray,Number}
    _dW = copy(ΔW)
    _dZ = copy(ΔW)
    chi1 = copy(ΔW)
  else
    _dW = zero(ΔW)
    _dZ = zero(ΔW)
    chi1 = zero(ΔW)
  end
  m = length(ΔW)
  Ihat2 = zeros(eltype(ΔW), m, m)
  tab = DRI1ConstantCache(real(uBottomEltypeNoUnits), real(tTypeNoUnits))
  g1 = zero(noise_rate_prototype)
  g2 = [zero(noise_rate_prototype) for k=1:m]
  g3 = [zero(noise_rate_prototype) for k=1:m]
  k1 = zero(rate_prototype); k2 = zero(rate_prototype); k3 = zero(rate_prototype)

  H02 = zero(u)
  H03 = zero(u)
  H12 = Vector{typeof(u)}()
  H13 = Vector{typeof(u)}()
  H22 = Vector{typeof(u)}()
  H23 = Vector{typeof(u)}()

  for k=1:m
    push!(H12,zero(u))
    push!(H13,zero(u))
    push!(H22,zero(u))
    push!(H23,zero(u))
  end

  tmp1 = zero(rate_prototype)
  tmpg = zero(noise_rate_prototype)

  DRI1Cache(u,uprev,_dW,_dZ,chi1,Ihat2,tab,g1,g2,g3,k1,k2,k3,H02,H03,H12,H13,H22,H23,tmp1,tmpg)
end
