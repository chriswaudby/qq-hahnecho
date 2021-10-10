## set up functions for multiplet calculations

"Return a matrix of simulated intensities of size (δH x δC x τ)"
function multiplet(param, δobsH, δobsC, τ, bfH, bfC, windowH, windowC)
  δH = param[1]
  δC = param[2]
  R2H = param[3]
  R2C = param[4]
  S2tc = param[5]
  csa = param[6]*1e-6
  J = param[7]

  # convert chemical shift axes to angular frequencies
  ωobsH = reshape(2π * δobsH * bfH, :, 1, 1)
  ωobsC = reshape(2π * δobsC * bfC, 1, :, 1)
  τobs = reshape(τ, 1, 1, :)

  # convert peak positions to angular frequencies
  ωH = 2π * δH * bfH
  ωC = 2π * δC * bfC

  # physical parameters
  μ₀ = 4π * 1e-7
  ħ = 1.055e-34
  γH = 2.675e8
  γC = 6.726e7
  rCH = 1.117e-10  # Tugarinov 2004
  rHH = sqrt(3) * rCH * sin(110.4*π/180)  # Tugarinov 2004
  ω0C = 2π * bfC * 1e6
  P2cosβ = -1/3.  # 109 degree tetrahedral geometry
  TAU = 1.0 / (2*125)  # coherence transfer period

  # relaxation rates (as function of S2tc)
  σ = 2/45 * (μ₀*ħ*γC*γH / (4π*rCH^3))^2 * S2tc * 1e-9
  η = 2/5 * (μ₀/(4π)) * P2cosβ/rCH^3 * ħ * γH * γC * ω0C * csa * S2tc * 1e-9
  ηHHHC = 1/5 * (μ₀/(4π))^2 * rCH^-3 * rHH^-3 * ħ^2 * γH^3 * γC * S2tc * 1e-9
  Δ = exp(-TAU * (9/20 * (μ₀*ħ*γH^2 / (4π*rHH^3))^2 * S2tc * 1e-9)) * cosh(-TAU * ηHHHC)
  ηHHHH = 9/40 * (μ₀/(4π))^2 * rHH^-6 * ħ^2 * γH^4 * S2tc * 1e-9

  # initial intensities (observing slow relaxing 1H coherences only)
  Iouter = 3 + 3Δ
  Iinner = 3 - Δ

  # relaxation rates for 13C transitions
  Rααα = R2C + 3σ + 2η
  Rααβ = R2C - σ + 2/3 * η
  Rαββ = R2C - σ - 2/3 * η
  Rβββ = R2C + 3σ - 2η

  return @. (
        Iouter * exp(-Rααα * τobs) * lineshape(ωC + 3π*J, Rααα, ωobsC, windowC) +
        Iinner * exp(-Rααβ * τobs) * lineshape(ωC + π*J, Rααβ, ωobsC, windowC) +
        Iinner * exp(-Rαββ * τobs) * lineshape(ωC - π*J, Rαββ, ωobsC, windowC) +
        Iouter * exp(-Rβββ * τobs) * lineshape(ωC - 3π*J, Rβββ, ωobsC, windowC)
        ) * lineshape(ωH, R2H, ωobsH, windowH)
end

## window functions
function lineshape(ω, R, ωobs, w::ExponentialWindow)
  x = @. -R + 1im*(ωobs - ω) - π*w.lb
  T = w.tmax
  return @. real((1 - exp(T*x)) / x)
end

function lineshape(ω, R, ωobs, w::CosWindow)
  x = @. -R + 1im*(ωobs - ω)
  T = w.tmax
  Tx = T * x
  return @. real(T * (π*exp(Tx) + 2*Tx) / (π^2 + 4*Tx^2))
end

function lineshape(ω, R, ωobs, w::Cos²Window)
  x = @. -R + 1im*(ωobs - ω)
  Tx = w.tmax * x
  return @. real((π^2*(1-exp(Tx)) + 2*Tx^2) / ((π^2 + Tx^2) * x))
end

Base.Broadcast.broadcastable(w::WindowFunction) = Ref(w)
