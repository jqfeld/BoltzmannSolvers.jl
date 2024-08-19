using DataInterpolations
import DataInterpolations: LinearInterpolation, AbstractInterpolation
struct NamedInterpolation{I} 
    name::Symbol
    itp::I
end
NamedInterpolation(name, y, x) = NamedInterpolation(name,LinearInterpolation(y,x; extrapolate=true))

Base.nameof(i::NamedInterpolation) = i.name
Base.show(io::IO, i::NamedInterpolation) = Base.print(io::IO, "[$(i.name)]")


(i::NamedInterpolation)(x) = i.itp(x)

# for later when I implement the reduced interpolation scheme
# TODO: implement reduced interpolation scheme (see LoKIHelpers.jl)
findclosest(a, A) = argmin(x -> abs(x - a), A)
exprange(start, stop, length) = [start; exp.(range(0, log(stop), length))]


"""
    `create_interpolation(df::DataFrame, y_name, x_name; name=y_name)`

Create interpolation object `y(x)` from `DataFrame` `df`. `y_name` is the column of
the dependent variable, `x_name` is the column name of the independent variable.

Returns a `NamedInterpolation` which works with `Symbolics.jl` and can be used in 
`ModelingToolkit.jl` models.

"""
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
