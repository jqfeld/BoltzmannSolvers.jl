module BoltzmannSolversSymbolicsExt

using BoltzmannSolvers
using Symbolics

# Proper registration (not a raw SymbolicUtils.term): gives the calls the
# correct scalar symtype *and shape* — raw terms carry Unknown shape, which
# breaks additions inside downstream equation assembly (e.g. Catalyst).
Symbolics.@register_symbolic (i::BoltzmannSolvers.NamedInterpolation)(x)
Symbolics.@register_symbolic (d::BoltzmannSolvers.NamedInterpolationDerivative)(x)

# Symbolic derivative of an interpolation. Without it, differentiating
# through a table yields *zero* (silently, in Symbolics 7 — unknown
# operations fall back to a zero rule), so an implicit solver that has to
# invert a tabulated quantity — e.g. E/N from a prescribed current via the
# drift velocity — sees no gradient and never leaves its initial guess.
#
# Symbolics 7 registers derivative rules with `@register_derivative`;
# version 6 used `Symbolics.derivative(f, args, ::Val{i})`. Support both,
# since the package allows either.
@static if isdefined(Symbolics, Symbol("@register_derivative"))
    Symbolics.@register_derivative (i::BoltzmannSolvers.NamedInterpolation)(x) 1 begin
        Symbolics.SConst(Symbolics.unwrap(BoltzmannSolvers.derivative_of(i)(Symbolics.wrap(x))))
    end
else
    Symbolics.derivative(i::BoltzmannSolvers.NamedInterpolation, args::NTuple{1,Any}, ::Val{1}) =
        BoltzmannSolvers.derivative_of(i)(args[1])
end

end
