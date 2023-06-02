module BoltzmannSolverSymbolicsExt
    using BoltzmannSolvers

    using Symbolics: Num, unwrap, SymbolicUtils

    (i::BoltzmannSolvers.NamedInterpolation)(x::Num) = SymbolicUtils.term(i, unwrap(x))
    SymbolicUtils.promote_symtype(::BoltzmannSolvers.NamedInterpolation, _...) = Real
end
