abstract type ExptData end

struct R2Data <: ExptData
    ZQ
    DQ
    QQ⁻
    QQ⁺
    B0
    R2Data(ZQ, DQ, QQ⁻, QQ⁺, B0=nothing) = new(ZQ, DQ, QQ⁻, QQ⁺, B0)
end

abstract type CPMGData <: ExptData end

struct CPMGData_MQ <: CPMGData
    νCPMG
    R2eff
    B0
    Trelax
end

struct CPMGData_SQ1H <: CPMGData
    νCPMG
    R2eff
    B0
    Trelax
end


struct modelpars_2state
    kex
    pB
end