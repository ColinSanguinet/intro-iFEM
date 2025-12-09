using CSV, DataFrames, GLMakie, Statistics

# Load your save_results.jl module
include("save_results.jl")

# Load direct and inverse results
direct = load_timeseries_csv("20251128_direct_static_100elts.csv")
inverse = load_timeseries_csv("20251209_inverse_static_100eltsDirect100eltsInverse.csv")

# Get the number of nodes from the series (assuming all components have the same number of nodes)
nnodes = size(direct.series["X"], 1)
@assert nnodes == size(inverse.series["X"], 1) "Direct and inverse results must have the same number of nodes"

# Function to calculate RMSE
function rmse(a, b)
    return sqrt(mean((a - b).^2))
end

# Calculate RMSE for each node and component
components = ["X", "Y", "Z", "R3"]
rmses = Dict{String, Vector{Float64}}()
for comp in components
    rmses[comp] = [rmse(direct.series[comp][i, :], inverse.series[comp][i, :]) for i in 1:nnodes]
end

# Create a figure for each component
for comp in components
    fig = Figure(size = (1000, 600))
    ax = Axis(fig[1, 1], title = "Comparison: $comp", xlabel = "Time", ylabel = comp)
    for node in 1:nnodes
        lines!(ax, direct.time, direct.series[comp][node, :], label = "Direct (Node $node)", linewidth = 2)
        lines!(ax, inverse.time, inverse.series[comp][node, :], label = "Inverse (Node $node)", linestyle = :dash, linewidth = 2)
    end
    axislegend(ax)
    save("comparison_$comp.png", fig)
end

# Create a figure for each component
fig = Figure(size = (1000, 600))
ax = Axis(fig[1, 1], title = "Comparison: Fy", xlabel = "nodenumber", ylabel = "Fy")
GLMakie.scatter!(ax, 1:nnodes, inverse.series["Fy"][:, end], label = "Inverse Fy", color = :blue, markersize = 10)
axislegend(ax)
save("comparison_Fy.png", fig)

# Plot RMSE for each component
fig_rmse = Figure(size = (1000, 600))
ax_rmse = Axis(fig_rmse[1, 1], title = "RMSE per Node", xlabel = "Node", ylabel = "RMSE")
for comp in components
    GLMakie.scatterlines!(ax_rmse, 1:nnodes, rmses[comp], label = comp, marker = :circle, markersize = 10, linewidth = 2)
end
axislegend(ax_rmse)
save("rmse_per_node.png", fig_rmse)

println("Plots saved to current directory.")

# Select 3 timesteps to plot (e.g., first, middle, last)
ntimesteps = length(direct.time)
timesteps_to_plot = [1, div(ntimesteps, 2), ntimesteps]

# Create a 3D figure for each selected timestep
for (i, t) in enumerate(timesteps_to_plot)
    fig = Figure(size = (1200, 800))
    ax = Axis3(fig[1, 1],
        title = "Beam at timestep $t (time = $(round(direct.time[t], digits=3)) s)",
        xlabel = "X", ylabel = "Y", zlabel = "Z",
        perspectiveness = 0.5, azimuth = 1.5π, elevation = 0.2π
    )

    # Plot direct results
    direct_positions = Point3f[
        (direct.series["X"][n, t], direct.series["Y"][n, t], direct.series["Z"][n, t])
        for n in 1:nnodes
    ]
    lines!(ax, direct_positions, color = :blue, linewidth = 4, label = "Direct")

    # Plot inverse results
    inverse_positions = Point3f[
        (inverse.series["X"][n, t], inverse.series["Y"][n, t], inverse.series["Z"][n, t])
        for n in 1:nnodes
    ]
    lines!(ax, inverse_positions, color = :red, linestyle = :dash, linewidth = 4, label = "Inverse")

    # Add nodes as scatter points for both
    GLMakie.scatter!(ax, direct_positions, color = :blue, markersize = 10, label = "Direct Nodes")
    GLMakie.scatter!(ax, inverse_positions, color = :red, markersize = 10, label = "Inverse Nodes")

    # Add legend
    axislegend(ax, position = :rb)

    save("beam_3d_timestep_$t.png", fig)
end

println("3D beam plots saved to current directory.")