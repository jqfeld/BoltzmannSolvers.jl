using DataInterpolations
import DataInterpolations: LinearInterpolation, AbstractInterpolation
struct NamedInterpolation{I} 
    name::Symbol
    itp::I
end
NamedInterpolation(name, y, x) = NamedInterpolation(name,LinearInterpolation(y,x))

Base.nameof(i::NamedInterpolation) = i.name
Base.show(io::IO, i::NamedInterpolation) = Base.print(io::IO, "[$(i.name)]")

using Symbolics: Num, unwrap, SymbolicUtils
(i::NamedInterpolation)(x::Num) = SymbolicUtils.term(i, unwrap(x))
SymbolicUtils.promote_symtype(::NamedInterpolation, _...) = Real

(i::NamedInterpolation)(x) = i.itp(x)

# for later when I implement the reduced interpolation scheme
findclosest(a, A) = argmin(x -> abs(x - a), A)
exprange(start, stop, length) = [start; exp.(range(0, log(stop), length))]

function create_interpolation(df::DataFrame, y_name, x_name; name=y_name)
    return NamedInterpolation(Symbol(name), df[!,y_name], df[!,x_name])
end

function create_interpolation(df::DataFrame, y_names::T, x_name; 
    name=reduce((x, y) -> String(x) * "," * String(y), y_names)) where T<:AbstractVector
    ys = zero(df[!,y_names[1]])
    for y_name in y_names
        ys += df[!,y_name] 
    end
    return NamedInterpolation(Symbol(name), ys, df[!,x_name])
end