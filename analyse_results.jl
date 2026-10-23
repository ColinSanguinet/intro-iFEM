using CSV, DataFrames, GLMakie, Statistics, Interpolations

# Load your save_results.jl module
include("save_results.jl")

function load_analysis(file_dir, file_inv)
    # Load direct and inverse results
    direct = load_timeseries_csv(file_dir)
    inverse = load_timeseries_csv(file_inv)

    println("======================\nSuccessfully loaded data for analysis\n=====================\n")

    return(direct, inverse)
end


# Function to calculate RMSE
function rmse(a, b)
    return sqrt(mean((a - b).^2))
end

function calculate_RMSE(direct, inverse, nSteps)
    # Calculate RMSE
    rmses = Dict{String, Vector{Float64}}()
    comp = "Y"
    rmse_y = [rmse(direct.series[comp][:, j], inverse.series[comp][:, j]) for j in 1:nSteps]
    max_val, max_idx = findmax(rmse_y)
    println("RMSE maximum $max_val for time index $max_idx .")

    return max_val, max_idx
end

function plot_eachComponent(direct, inverse, node_number, component_dict, output_dir, output_name)
    # Create a figure for each component on a chosen node
    for (comp,lab) in component_dict
        local fig = Figure(size = (1000, 600))
        local ax = Axis(
            fig[1, 1], 
            # title = "$comp @ node $node_number", 
            # titlesize = 25,
            xlabel = "Time (s)", 
            ylabel = comp*lab,
            xlabelsize = 25,       # Increases X label font size
            ylabelsize = 25,       # Increases Y label font size
            xticklabelsize = 20,   # Increases X tick label font size
            yticklabelsize = 20    # Increases Y tick label font size
            )
        lines!(ax, direct.time, direct.series[comp][node_number, :], label = "Direct (Node $node_number)", linewidth = 2)
        lines!(ax, inverse.time, inverse.series[comp][node_number, :], label = "Inverse (Node $node_number)", linestyle = :dash, linewidth = 2)
        axislegend(ax, labelsize = 25)
        mkpath(output_dir)
        save(joinpath(output_dir,"$comp--" * output_name * ".png"), fig)
    end
end

function plot_boxplot_disp(direct, inverse, component_dict, output_dir, output_name)
    # Plot boxplot of the error
    fig_rmse = Figure(size = (1000, 600))
    ax_rmse = Axis(
        fig_rmse[1, 1],
        # title = "Boxplot", 
        # titlesize = 25,
        xlabel = "Position [m]", 
        ylabel = "absolute error",
        xlabelsize = 25,       # Increases X label font size
        ylabelsize = 25,       # Increases Y label font size
        xticklabelsize = 20,   # Increases X tick label font size
        yticklabelsize = 20    # Increases Y tick label font size
    )
    nnodes = size(direct.series["X"], 1)
    nodes_pos_dir = [direct.metadata["nodeCoord"][node,1] for node in 1:nnodes]
    nnodes_inv = size(inverse.series["X"], 1)
    nodes_pos_inv = [inverse.metadata["nodeCoord"][node,1] for node in 1:nnodes_inv]
    nSteps = length(direct.time)
    abs_error = zeros((nSteps,nnodes))
    for (comp,lab) in component_dict
        for t in 1:nSteps
            comp_serie_inv =  linear_interpolation(nodes_pos_inv, inverse.series[comp][:,t])
            # flatten the per-node vectors into one long vector
            local abs_error_t = vcat([abs.(direct.series[comp][i, t] - comp_serie_inv(nodes_pos_dir[i])) for i in 1:nnodes]...)
            abs_error[t,:] = abs_error_t
        end
        # do the same for the x positions (repeat node index nSteps times)
        local nnodes_list = vcat([fill(nodes_pos_dir[i], nSteps) for i in 1:nnodes]...)
        if length(abs_error[:,1]) == 1
            scatter!(ax_rmse, nnodes_list, abs_error, label=comp*lab, markersize = 15)
        else
            w = nodes_pos_dir[2]-nodes_pos_dir[1]
            boxplot!(ax_rmse, nnodes_list, abs_error, label=comp*lab, width = w)
        end

    end
    # axislegend(ax_rmse, labelsize = 25)
    fig_rmse[1,2] = Legend(fig_rmse, ax_rmse, labelsize = 25)
    mkpath(output_dir)
    save(joinpath(output_dir, output_name * ".png"), fig_rmse)
end

function plot_deflection(direct, inverse, rng, output_dir, output_name)
    nnodes = size(direct.series["X"], 1)
    nodes_pos_dir = [direct.metadata["nodeCoord"][node,1] for node in 1:nnodes]
    nnodes_inv = size(inverse.series["X"], 1)
    nodes_pos_inv = [inverse.metadata["nodeCoord"][node,1] for node in 1:nnodes_inv]
    # Create a figure for beam y displacements
    fig = Figure(size = (1000, 600))
    ax_disp = Axis(
        fig[1, 1], 
        # title = "y displacement", 
        # titlesize = 25,
        xlabel = "Position [m]", 
        ylabel = "y disp [m]",
        xlabelsize = 25,       # Increases X label font size
        ylabelsize = 25,       # Increases Y label font size
        xticklabelsize = 20,   # Increases X tick label font size
        yticklabelsize = 20    # Increases Y tick label font size
        )
    for state_id in rng
        global elem1 = GLMakie.scatter!(ax_disp, nodes_pos_dir, direct.series["Y"][:, state_id], label = "Forward : loadstep = "*string(state_id), color = (:blue,(state_id/nSteps*9+1)*0.1), markersize = 10)
        # GLMakie.scatter!(ax, 1:nnodes, inverse.metadata["Y_input"][:, state_id], label = "Noised input", color = (:purple, (state_id/nSteps*9+1)*0.1), markersize = 10)
        global elem2 = GLMakie.scatter!(ax_disp, nodes_pos_inv, inverse.series["Y"][:, state_id], label = "Inverse : loadstep = "*string(state_id), color = (:red, (state_id/nSteps*9+1)*0.1),markersize = 10)
    end
    # axislegend(ax_disp, labelsize = 25, position = :rb)
    Legend(fig[1,2], [elem1, elem2], ["forward", "inverse"], labelsize = 25)
    mkpath(output_dir)
    save(joinpath(output_dir,output_name * ".png"), fig)
end

function plot_forceComparison(direct, inverse, rng, output_dir, output_name)
    nnodes = size(direct.series["X"], 1)
    nodes_pos_dir = [direct.metadata["nodeCoord"][node,1] for node in 1:nnodes]
    nnodes_inv = size(inverse.series["X"], 1)
    nodes_pos_inv = [inverse.metadata["nodeCoord"][node,1] for node in 1:nnodes_inv]
    # Create a figure for comparing unknown forces reconstruction
    fig = Figure(size = (1000, 600))
    ax_force = Axis(
        fig[1, 1], 
        # title = "Comparison: Fy", 
        # titlesize = 25,
        xlabel = "Position [m]", 
        ylabel = "Fy [N]",
        xlabelsize = 25,       # Increases X label font size
        ylabelsize = 25,       # Increases Y label font size
        xticklabelsize = 20,   # Increases X tick label font size
        yticklabelsize = 20    # Increases Y tick label font size
        )
    for state_id in rng
        local ps_forward = Point2f.(nodes_pos_dir, 0)
        local ps_inverse = Point2f.([pos + 0.02 for pos in nodes_pos_inv], 0)
        local vs_forward = Vec2f.(0,direct.series["Fdir"][:, end])
        local vs_inverse = Vec2f.(0,inverse.series["Fy"][:, end])
        global elem1 = arrows2d!(ax_force, ps_forward, vs_forward, minshaftlength = 1, color = (:blue, (state_id/nSteps*9+1)*0.1))
        global elem2 = arrows2d!(ax_force, ps_inverse, vs_inverse, minshaftlength = 1, color = (:red, (state_id/nSteps*9+1)*0.1) )
    end
    elem1, elem2 = LineElement(color = :blue),  LineElement(color = :red)
    Legend(fig[1,2], [elem1, elem2], ["forward", "inverse"], labelsize = 25)
    GLMakie.scatter!(ax_force, nodes_pos_dir, zeros(nnodes), color = :black, markersize = 10)
    GLMakie.scatter!(ax_force, nodes_pos_inv, zeros(nnodes_inv), color = :purple, markersize = 10)
    mkpath(output_dir)
    save(joinpath(output_dir,output_name * ".png"), fig)
end


function plot_boxplot_force(direct, inverse, output_dir, output_name)
    nnodes = size(direct.series["X"], 1)
    nodes_pos_dir = [direct.metadata["nodeCoord"][node,1] for node in 1:nnodes]
    nnodes_inv = size(inverse.series["X"], 1)
    nodes_pos_inv = [inverse.metadata["nodeCoord"][node,1] for node in 1:nnodes_inv]
    # Plot boxplot of the error
    fig_ = Figure(size = (1000, 600))
    ax_ = Axis(
        fig_[1, 1],
        # title = "Boxplot", 
        # titlesize = 25,
        xlabel = "Position [m]", 
        ylabel = "absolute error",
        xlabelsize = 25,       # Increases X label font size
        ylabelsize = 25,       # Increases Y label font size
        xticklabelsize = 20,   # Increases X tick label font size
        yticklabelsize = 20    # Increases Y tick label font size
        )
    nSteps = length(direct.time)
    abs_error = zeros((nSteps,nnodes))
    lab = " [N]"
    for t in 1:nSteps
        comp_serie_inv =  linear_interpolation(nodes_pos_inv, inverse.series["Fy"][:,t])
        # flatten the per-node vectors into one long vector
        local abs_error_t = vcat([abs.(direct.series["Fdir"][i, t] - comp_serie_inv(nodes_pos_dir[i])) for i in 1:nnodes]...)
        abs_error[t,:] = abs_error_t
    end
    # do the same for the x positions (repeat node index nSteps times)
    local nnodes_list = vcat([fill(nodes_pos_dir[i], nSteps) for i in 1:nnodes]...)
    if length(abs_error[:,1]) == 1
        scatter!(ax_, nnodes_list, abs_error, label="Fy"*lab, markersize = 15)
    else
        w = nodes_pos_dir[2]-nodes_pos_dir[1]
        boxplot!(ax_, nnodes_list, abs_error, label="Fy"*lab, width = w)
    end
    axislegend(ax_, labelsize = 25)
    mkpath(output_dir)
    save(joinpath(output_dir,output_name * ".png"), fig_)
end

function plot_forces(direct, inverse, node_number, output_dir, output_name)
    fig = Figure(size = (1000, 600))
    ax = Axis(
        fig[1, 1], 
        # title = "$comp @ node $node_number", 
        # titlesize = 25,
        xlabel = "Time (s)", 
        ylabel = "force y [N]",
        xlabelsize = 25,       # Increases X label font size
        ylabelsize = 25,       # Increases Y label font size
        xticklabelsize = 20,   # Increases X tick label font size
        yticklabelsize = 20    # Increases Y tick label font size
        )
    lines!(ax, direct.time, direct.series["Fdir"][node_number, :], label = "Direct (Node $node_number)", linewidth = 2)
    lines!(ax, inverse.time, inverse.series["Fy"][node_number, :], label = "Inverse (Node $node_number)", linestyle = :dash, linewidth = 2)
    axislegend(ax, labelsize = 25)
    mkpath(output_dir)
    save(joinpath(output_dir,"Fy--" * output_name * ".png"), fig)
end

function calculate_forceReconstructionAccuracy(direct, inverse, output_dir, output_name)
    Wdir = sum(direct.series["Fdir"][:,:] .* direct.series["Y"][:,:], dims=1)
    Winv = sum(inverse.series["Fy"][:,:] .* inverse.series["Y"][:,:], dims=1)
    epsilon = vec(abs.((Wdir .- Winv) ./ Wdir) * 100)
    epsilon = replace!(epsilon, Inf=>NaN)
    fig = Figure(size = (1000, 600))
    ax = Axis(
        fig[1, 1], 
        # title = "$comp @ node $node_number", 
        # titlesize = 25,
        xlabel = "Time (s)", 
        ylabel = "force_error [%]",
        xlabelsize = 25,       # Increases X label font size
        ylabelsize = 25,       # Increases Y label font size
        xticklabelsize = 20,   # Increases X tick label font size
        yticklabelsize = 20    # Increases Y tick label font size
        )
    scatter!(ax, direct.time, epsilon, markersize = 15)
    mkpath(output_dir)
    save(joinpath(output_dir,output_name * ".png"), fig)
end