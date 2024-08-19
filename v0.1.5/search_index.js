var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = BoltzmannSolvers","category":"page"},{"location":"#BoltzmannSolvers","page":"Home","title":"BoltzmannSolvers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for BoltzmannSolvers.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [BoltzmannSolvers]","category":"page"},{"location":"#BoltzmannSolvers.create_interpolation-Tuple{DataFrames.DataFrame, Any, Any}","page":"Home","title":"BoltzmannSolvers.create_interpolation","text":"`create_interpolation(df::DataFrame, y_name, x_name; name=y_name)`\n\nCreate interpolation object y(x) from DataFrame df. y_name is the column of the dependent variable, x_name is the column name of the independent variable.\n\nReturns a NamedInterpolation which works with Symbolics.jl and can be used in  ModelingToolkit.jl models.\n\n\n\n\n\n","category":"method"},{"location":"#BoltzmannSolvers.normalize_reaction_names-Tuple{Any}","page":"Home","title":"BoltzmannSolvers.normalize_reaction_names","text":"normalize_reaction_name!(x)\n\nNormalize the reaction name using PlasmaSpecies.jl.  If PlasmaSpecies.jl is not loaded, this is a NOP. PlasmaSpecies is NOT loaded!\n\n\n\n\n\n","category":"method"},{"location":"#BoltzmannSolvers.parse_reaction_names-Tuple{MultiBolt, Any}","page":"Home","title":"BoltzmannSolvers.parse_reaction_names","text":"parse_reaction_names(::Multibolt, source)\n\nDo nothing. For MultiBolt it is more convenient to parse the reaction process already in load_raw_dataframe.\n\n\n\n\n\n","category":"method"}]
}
