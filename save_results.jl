using CSV, DataFrames
using Dates

function yyyymmdd()
    return Dates.format(Dates.today(), "yyyymmdd")
end

function dated_base(base::AbstractString)
    return string(yyyymmdd(),"_",base)
end

"""
    parse_metadata_csv(path::AbstractString) -> Dict{String,Any}

Read the first three CSV records (METADATA, header row, metadata row) using CSV.jl
and convert each metadata value to Float64, Vector{Float64}, or Matrix{Float64}
when possible. If conversion fails the raw string is returned.
"""
function parse_metadata_csv(csvfile::AbstractString)
    !isfile(csvfile) && error("CSV file not found: $csvfile")

    rows = collect(CSV.File(csvfile; header=false, limit=3))
    length(rows) < 3 && error("CSV must contain at least 3 rows (METADATA, header, values)")

    hdr_row = rows[2]
    val_row = rows[3]

    headers = [string(x) for x in values(hdr_row)]
    rawvals = [string(x) for x in values(val_row)]

    # (reconstruction of headers/raw values moved below after helper defs)

    function try_parse_num(tok::AbstractString)
        tok = strip(tok)
        isempty(tok) && return nothing
        try
            return parse(Float64, tok)
        catch
            return nothing
        end
    end

    # Split string at top-level separators (ignore separators inside nested brackets)
    function split_top_level(s::AbstractString, sep::Char)
        res = String[]
        buf = IOBuffer()
        depth = 0
        for c in collect(s)
            if c == '['
                depth += 1
                print(buf, c)
            elseif c == ']'
                depth -= 1
                print(buf, c)
            elseif c == sep && depth == 0
                push!(res, String(take!(buf)))
            else
                print(buf, c)
            end
        end
        leftover = String(take!(buf))
        if !isempty(strip(leftover))
            push!(res, leftover)
        end
        return res
    end

        # Reconstruct header/value lists if CSV.jl fragmented the metadata row
        # into multiple columns (e.g. produced many "missing" placeholders).
        # Detect a combined header token anywhere (not just at headers[1]).
        if any(h -> occursin(',', h), headers) && (length(headers) == 1 || any(h -> lowercase(strip(h)) == "missing", headers))
            # find the element that contains the combined header string
            idx = findfirst(h -> occursin(',', h), headers)
            combined = idx === nothing ? headers[1] : headers[idx]
            headers = [strip(h) for h in split_top_level(combined, ',')]
            # reconstruct rawvals by joining fragments (treat literal "missing" as empty)
            raw_join = join([v == "missing" ? "" : v for v in rawvals], ' ')
            rawvals = [strip(v) for v in split_top_level(raw_join, ',')]
        end

        # Merge any bracketed values that were split across multiple fields
        # (e.g. a long "[...; ...; ...]" that CSV widened). This concatenates
        # subsequent rawvals until the closing ']' is found.
        i = 1
        merged = String[]
        while i <= length(rawvals)
            s = rawvals[i]
            s_strip = strip(s)
            if startswith(s_strip, "[") && !endswith(s_strip, "]")
                buf = s
                j = i + 1
                found = false
                while j <= length(rawvals)
                    buf *= ";" * rawvals[j]
                    if occursin("]", rawvals[j])
                        found = true
                        break
                    end
                    j += 1
                end
                push!(merged, strip(buf))
                i = found ? j + 1 : j + 1
            else
                push!(merged, s_strip)
                i += 1
            end
        end
        rawvals = merged

        # (reconstructed header/value lists)

    function parse_bracketed(s::AbstractString)
        inner = strip(s[2:end-1])
        # Nested bracketed arrays like [[1,2],[3,4]]
        if startswith(inner, "[")
            parts = split_top_level(inner, ',')
            parsed_parts = [begin
                p = strip(x)
                if startswith(p, "[") && endswith(p, "]")
                    parse_bracketed(p)
                else
                    # try to parse as simple vector
                    toks = filter(x->x!="", split(p, r"[,\s]+"))
                    try
                        [parse(Float64,t) for t in toks]
                    catch
                        toks
                    end
                end
            end for x in parts]

            # If all parts are numeric vectors of equal length -> matrix
            if all(x->isa(x, AbstractVector{Float64}), parsed_parts)
                lengths = unique(length.(parsed_parts))
                if length(lengths) == 1
                    nrows = length(parsed_parts)
                    ncols = lengths[1]
                    M = Array{Float64}(undef, nrows, ncols)
                    for i in 1:nrows
                        M[i, :] = parsed_parts[i]
                    end
                    return M
                end
            end
            return parsed_parts
        end

        # matrix using semicolons to separate rows
        if occursin(";", inner)
            rows_s = split(inner, ';')
            parsed = [filter(x->x!="", split(strip(r), r"[,\s]+")) for r in rows_s]
            parsed_nums = Vector{Vector{Float64}}()
            for r in parsed
                row_nums = Float64[]
                for t in r
                    try
                        push!(row_nums, parse(Float64, t))
                    catch
                        # fallback to strings if any element fails to parse
                        return parsed
                    end
                end
                push!(parsed_nums, row_nums)
            end
            nrows = length(parsed_nums); ncols = length(parsed_nums[1])
            M = Array{Float64}(undef, nrows, ncols)
            for i in 1:nrows
                length(parsed_nums[i]) == ncols || error("inconsistent columns in matrix metadata")
                M[i, :] = parsed_nums[i]
            end
            return M
        end

        # simple bracketed vector or flattened matrix
        toks = filter(x->x!="", split(inner, r"[,\s]+"))
        try
            nums = [parse(Float64,t) for t in toks]
            # If flattened node coordinates (triplets), reshape to n x 3 matrix
            if length(nums) >= 3 && length(nums) % 3 == 0
                nrows = length(nums) ÷ 3
                M = Array{Float64}(undef, nrows, 3)
                for i in 1:nrows
                    M[i, :] = nums[(3*(i-1)+1):(3*i)]
                end
                return M
            end
            return nums
        catch
            return toks
        end
    end

    function convert_value(s::AbstractString)
        s = strip(s)
        isempty(s) && return s
        # bracketed list or matrix
        if startswith(s, "[") && endswith(s, "]")
            return parse_bracketed(s)
        end
        # whitespace- or comma-separated vector without brackets
        if occursin(',', s) || occursin(' ', s)
            toks = filter(x->x!="", split(s, r"[,\s]+"))
            nums = map(t -> try_parse_num(t), toks)
            if all(x->x !== nothing, nums)
                return [Float64(x) for x in nums]
            else
                return toks
            end
        end
        # scalar number
        n = try_parse_num(s)
        n !== nothing && return n
        # fallback: return original string
        return s
    end

    meta = Dict{String,Any}()
    for (h,v) in zip(headers, rawvals)
        meta[h] = convert_value(v)
    end
    return meta
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
    raw_values = split(first_lines[3], ",")

    # Recombine raw_values in case some metadata fields contained commas
    function unquote_and_unescape(s)
        s = strip(s)
        if startswith(s, '"') && endswith(s, '"')
            inner = s[2:end-1]
                return replace(inner, "\"\"" => "\"")
        else
            return s
        end
    end

    values = String[]
    i = 1
    while i <= length(raw_values)
        v = raw_values[i]
        vstr = strip(v)
        if startswith(vstr, "[") && !endswith(vstr, "]")
            # accumulate until closing bracket found
            acc = v
            j = i + 1
            while j <= length(raw_values) && !occursin("]", raw_values[j])
                acc *= "," * raw_values[j]
                j += 1
            end
            if j <= length(raw_values)
                acc *= "," * raw_values[j]
                push!(values, unquote_and_unescape(acc))
                i = j + 1
            else
                push!(values, unquote_and_unescape(acc))
                break
            end
        else
            push!(values, unquote_and_unescape(v))
            i += 1
        end
    end

    # Build dictionary
    metadata_dict = parse_metadata_csv(path)
    # metadata_dict = Dict(h => convert_to_value(v) for (h, v) in zip(headers, values))


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

    # Convert everything to strings
    title_row = ["METADATA"]
    header_row = string.(fields)
    metadata_row = string.(values)
    # Quote fields that contain commas or quotes so CSV splitting is safe
    function quote_csv_field(s)
        s = string(s)
        if occursin('"', s)
            s = replace(s, '"' => "\"\"")
        end
        if occursin(',', s) || occursin('"', s) || occursin('\n', s)
            return '"' * s * '"'
        else
            return s
        end
    end
    metadata_row_quoted = quote_csv_field.(metadata_row)
    data_row = ["DATA"]

    # Open file and write properly quoted CSV
    open(base * ".csv", "w") do io
        println(io, join(title_row, ","))   # First row
        println(io, join(header_row, ","))     # Second row
        println(io, join(metadata_row_quoted, ","))       # Third row (quoted as needed)
        println(io, join(data_row, ","))   # Fourth row
    end
end



# Tests
#------------------------------------------
