using CSV, DataFrames
using Dates

function yyyymmdd()
    return Dates.format(Dates.today(), "yyyymmdd")
end

function dated_base(base::AbstractString)
    return string(yyyymmdd(),"_",base)
end

"""
    save_timeseries_csv(path; comps::Dict, time=nothing)

Save multiple time-series (vector-of-vectors) into a single CSV file.

`comps` should be a Dict mapping a short series name (e.g. "X", "Y", "Z" or "force")
to a vector of length `nsteps`, where each element is a vector of length `nnodes` (or a scalar for node-insensitive series).

The produced CSV has `nsteps` rows (one per time step). If `time` is provided it is saved as the first column `time`.
For each series `name` with `nnodes > 1` the CSV will contain columns `name_1, name_2, ..., name_nnodes` (per-node values).
For scalar series (no underscore in column naming) `name` will be written as a single column.

Returns `true` on success.
"""
function save_timeseries_csv(path; metadata, comps::Dict=Dict(),  time=nothing)
    isempty(comps) && error("`comps` must be a Dict mapping series names to vector-of-vectors, e.g. Dict(\"X\"=>x_, ...)")

    # Basic validation: all series must have same number of time steps
    nsteps = nothing
    for (name, vecs) in comps
        if nsteps === nothing
            nsteps = length(vecs)
        else
            @assert length(vecs) == nsteps "All series must have same number of time steps; mismatch for $name"
        end
    end

    # Convert each series into a (nnodes, nsteps) matrix where possible
    comp_mats = Dict{String, Matrix{Float64}}()
    comp_scalars = Dict{String, Vector{Float64}}()
    for (name, vecs) in comps
        # Determine if inner elements are vectors (per-node) or scalars
        first_inner = first(vecs)
        if isa(first_inner, AbstractVector)
            nnodes = length(first_inner)
            for (i, vv) in enumerate(vecs)
                @assert length(vv) == nnodes "Inconsistent node length in $name at step $i"
            end
            M = hcat([vec(v) for v in vecs]...)   # (nnodes, nsteps)
            comp_mats[name] = M
        else
            # treat as scalar-per-time series
            comp_scalars[name] = [float(v) for v in vecs]
        end
    end

    base = replace(path, r"\.csv$" => "")
    df = DataFrame()
    if time !== nothing
        @assert length(time) == nsteps "time vector must match number of steps"
        df.time = time
    end

    # Add scalar series first
    for (name, vec) in comp_scalars
        df[!, Symbol(name)] = vec
    end

    # Add per-node columns for each matrix series
    for (name, M) in comp_mats
        nnodes = size(M, 1)
        for node in 1:nnodes
            colname = Symbol(string(name, "_", node))
            df[!, colname] = vec(M[node, :])
        end
    end

    save_metadata(base, metadata)

    CSV.write(base * ".csv", df ; writeheader = true, append=true)
    return true
end


"""
    load_timeseries_csv(path)

Load CSV file previously written by `save_timeseries_csv`.
- Accepts `path` which can be the base prefix (e.g. "results_combined").
- Returns a NamedTuple with fields:
  - `time`: time vector (if present in CSV), or nothing
  - `series`: Dict mapping series names (e.g., "X", "Y", "Z") to (nnodes, nsteps) matrices
  - `scalar_series`: Dict mapping scalar series names to vectors of length nsteps
  - `shape`: (nnodes, nsteps) tuple for per-node series (or nothing if no per-node series)
"""
function load_timeseries_csv(path)
    base = replace(path, r"\.csv$" => "")
    csvfile = base * ".csv"
    
    !isfile(csvfile) && error("CSV file not found: $csvfile")
    
    # Extract metadata
    first_lines = readlines(csvfile)[1:4]

    # Extract metadata, headers, and values
    headers = split(first_lines[2], ",")
    values = split(first_lines[3], ",")

    # Build dictionary
    metadata_dict = Dict(h => to_value(v) for (h, v) in zip(headers, values))

    # Extract data
    df = CSV.File(csvfile; skipto=6, header = 5) |> DataFrame
    nsteps = size(df, 1)
    
    # Extract time column if present
    time = nothing
    colnames = String.(names(df))
    if "time" in colnames
        time = Vector(df[:, :time])
        colnames = filter(x -> x != "time", colnames)
    end
    
    # Parse column names to identify per-node and scalar series
    # Per-node columns: "X_1", "X_2", ...; Scalar columns: "force", etc.
    series_dict = Dict{String, Matrix{Float64}}()
    scalar_dict = Dict{String, Vector{Float64}}()
    
    per_node_cols = filter(c -> occursin("_", c), colnames)
    scalar_cols = filter(c -> !occursin("_", c), colnames)
    
    # Group per-node columns by series name
    for col in per_node_cols
        parts = split(col, "_")
        if length(parts) == 2
            series_name = parts[1]
            node_idx = parse(Int, parts[2])
            
            if !haskey(series_dict, series_name)
                series_dict[series_name] = zeros(0, nsteps)
            end
        end
    end
    
    # Reconstruct matrices for each per-node series
    for (series_name, _) in series_dict
        cols_for_series = filter(c -> startswith(c, series_name * "_"), per_node_cols)
        nnodes = length(cols_for_series)
        # Sort by numeric index, not lexicographic (to get X_1, X_2, ..., X_10, X_11 not X_1, X_10, X_11, ...)
        cols_sorted = sort(cols_for_series; by=c -> parse(Int, split(c, "_")[2]))
        M = Matrix{Float64}(undef, nnodes, nsteps)
        for (node_idx, col) in enumerate(cols_sorted)
            M[node_idx, :] = vec(Vector(df[:, Symbol(col)]))
        end
        series_dict[series_name] = M
    end
    

    # Load scalar series
    for col in scalar_cols
        scalar_dict[col] = Vector(df[:, Symbol(col)])
    end
    
    return (time=time, series=series_dict, scalar_series=scalar_dict, metadata=metadata_dict)
end

function add_struct_to_dict(dict, structure)
    for field in fieldnames(typeof(structure))
        dict[string(field)] = getfield(structure, field)
    end
    return dict
end

function save_metadata(base, meta)
    # Collect field names and values
    fields = []
    values = []
    for k in meta 
        push!(fields, k[1])
        push!(values, k[2])
    end

    # Convert everything to strings safely
    title_row = ["METADATA"]
    header_row = string.(fields)
    metadata_row = string.(values)
    data_row = ["DATA"]

    # Open file and write properly quoted CSV
    open(base * ".csv", "w") do io
        println(io, join(title_row, ","))   # First row
        println(io, join(header_row, ","))     # Second row
        println(io, join(metadata_row, ","))       # Third row
        println(io, join(data_row, ","))   # Fourth row
    end

end

# Tests
#------------------------------------------
