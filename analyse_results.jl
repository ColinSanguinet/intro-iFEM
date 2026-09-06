using CSV, DataFrames, GLMakie, Statistics

# Load your save_results.jl module
include("save_results.jl")

path_direct = "20260906_direct_static_middleForce_30elts.csv"
path_inverse = "20260906_inverse_static_noised_30elts.csv"

base_direct = split(path_direct, '/')[end][10:end-4]
base_inverse = split(path_inverse, '/')[end][10:end-4]

# Load direct and inverse results
direct = load_timeseries_csv(path_direct)
inverse = load_timeseries_csv(path_inverse)

# Create results folder path
output_dir = joinpath(@__DIR__, "results")
mkpath(output_dir)

println("\n===============")

# Get the number of nodes from the series (assuming all components have the same number of nodes)
nnodes = size(direct.series["X"], 1)
@assert nnodes == size(inverse.series["X"], 1) "Direct and inverse results must have the same number of nodes"

nLoadSteps = length(direct.time)
rng = nLoadSteps:1:nLoadSteps 

components = ["X", "Y", "Z", "R3"]

# Function to calculate RMSE
function rmse(a, b)
    return sqrt(mean((a - b).^2))
end

# Calculate RMSE
rmses = Dict{String, Vector{Float64}}()
comp = "Y"
rmse_y = [rmse(direct.series[comp][:, j], inverse.series[comp][:, j]) for j in 1:nLoadSteps]
max_val, max_idx = findmax(rmse_y)
println("RMSE maximum $max_val for time index $max_idx .")

# # Create a figure for each component on a chosen node
# node_number = 8
# for comp in components
#     local lab = (comp == "R3") ? " [rad]" : " [m]"
#     local fig = Figure(size = (1000, 600))
#     local ax = Axis(
#         fig[1, 1], 
#         title = "$comp @ node $node_number", 
#         xlabel = "Time (s)", 
#         ylabel = comp*lab,
#         titlesize = 25,
#         xlabelsize = 25,       # Increases X label font size
#         ylabelsize = 25,       # Increases Y label font size
#         xticklabelsize = 20,   # Increases X tick label font size
#         yticklabelsize = 20    # Increases Y tick label font size
#         )
#     lines!(ax, direct.time, direct.series[comp][node_number, :], label = "Direct (Node $node_number)", linewidth = 2)
#     lines!(ax, inverse.time, inverse.series[comp][node_number, :], label = "Inverse (Node $node_number)", linestyle = :dash, linewidth = 2)
#     axislegend(ax, labelsize = 25)
#     save(joinpath(output_dir,"$comp@node$node_number-" * base_inverse * ".png"), fig)
# end

# Plot boxplot of the error
fig_rmse = Figure(size = (1000, 600))
ax_rmse = Axis(
    fig_rmse[1, 1],
    title = "Boxplot", 
    titlesize = 25,
    xlabel = "node number", 
    ylabel = "absolute error",
    xlabelsize = 25,       # Increases X label font size
    ylabelsize = 25,       # Increases Y label font size
    xticklabelsize = 20,   # Increases X tick label font size
    yticklabelsize = 20    # Increases Y tick label font size
    )

for comp in components
# for comp in ["Y"]
    # flatten the per-node vectors into one long vector
    local abs_error = vcat([abs.(direct.series[comp][i, :] - inverse.series[comp][i, :]) for i in 1:nnodes]...)
    # do the same for the x positions (repeat node index nLoadSteps times)
    local nnodes_list = vcat([fill(i, nLoadSteps) for i in 1:nnodes]...)
    local lab = (comp == "R3") ? " [rad]" : " [m]"
    if length(abs_error[1]) == 1
        scatter!(ax_rmse, nnodes_list, abs_error, label=comp*lab, markersize = 15)
    else
        boxplot!(ax_rmse, nnodes_list, abs_error, label=comp*lab)
    end
end
# axislegend(ax_rmse, labelsize = 25)
fig_rmse[1,2] = Legend(fig_rmse, ax_rmse, labelsize = 25)
save(joinpath(output_dir,"error_boxplot-" * base_inverse * ".png"), fig_rmse)

# Create a figure for beam y displacements
fig = Figure(size = (1000, 600))
ax_disp = Axis(
    fig[1, 1], 
    title = "y displacement", 
    titlesize = 25,
    xlabel = "node number", 
    ylabel = "y disp [m]",
    xlabelsize = 25,       # Increases X label font size
    ylabelsize = 25,       # Increases Y label font size
    xticklabelsize = 20,   # Increases X tick label font size
    yticklabelsize = 20    # Increases Y tick label font size
    )
for state_id in rng
    global elem1 = GLMakie.scatter!(ax_disp, 1:nnodes, direct.series["Y"][:, state_id], label = "Forward : loadstep = "*string(state_id), color = (:blue,(state_id/nLoadSteps*9+1)*0.1), markersize = 10)
    # GLMakie.scatter!(ax, 1:nnodes, inverse.metadata["Y_input"][:, state_id], label = "Noised input", color = (:purple, (state_id/nLoadSteps*9+1)*0.1), markersize = 10)
    global elem2 = GLMakie.scatter!(ax_disp, 1:nnodes, inverse.series["Y"][:, state_id], label = "Inverse : loadstep = "*string(state_id), color = (:red, (state_id/nLoadSteps*9+1)*0.1),markersize = 10)
end
# axislegend(ax_disp, labelsize = 25, position = :rb)
Legend(fig[1,2], [elem1, elem2], ["forward", "inverse"], labelsize = 25)
save(joinpath(output_dir,"comparison_ydisp-" * base_inverse * ".png"), fig)


# Print the sum of reconstructed forces
println("Sum of all the forces components : sum(f) = ", sum(inverse.series["Fy"][:, end]))

# # Create a figure for unknown forces
# fig = Figure(size = (1000, 600))
# ax_force = Axis(
#     fig[1, 1], 
#     title = "Unknown force reconstruction", 
#     titlesize = 25,
#     xlabel = "node number", 
#     ylabel = "Fy [N]",
#     xlabelsize = 25,       # Increases X label font size
#     ylabelsize = 25,       # Increases Y label font size
#     xticklabelsize = 20,   # Increases X tick label font size
#     yticklabelsize = 20    # Increases Y tick label font size
#     )
# for state_id in rng
#     ps = Point2f.(1:nnodes, 0)
#     vs = Vec2f.(0,inverse.series["Fy"][:, end])
#     arrows2d!(ax_force, ps, vs, minshaftlength = 1, color = (:blue, (state_id/nLoadSteps*9+1)*0.1), label = "loadstep = "*string(state_id))
# end
# GLMakie.scatter!(ax_force, 1:nnodes, zeros(nnodes), color = :black, markersize = 10)
# # axislegend(ax_force, labelsize = 25, position = :rb)
# save(joinpath(output_dir,"Fy-" * base_inverse * ".png"), fig)

# Create a figure for comparing unknown forces reconstruction
fig = Figure(size = (1000, 600))
ax_force = Axis(
    fig[1, 1], 
    title = "Comparison: Fy", 
    titlesize = 25,
    xlabel = "node number", 
    ylabel = "Fy [N]",
    xlabelsize = 25,       # Increases X label font size
    ylabelsize = 25,       # Increases Y label font size
    xticklabelsize = 20,   # Increases X tick label font size
    yticklabelsize = 20    # Increases Y tick label font size
    )
for state_id in rng
    local ps_forward = Point2f.(1:nnodes, 0)
    local ps_inverse = Point2f.([i + 0.25 for i in 1:nnodes], 0)
    local vs_forward = Vec2f.(0,direct.series["Fdir"][:, end])
    local vs_inverse = Vec2f.(0,inverse.series["Fy"][:, end])
    global elem1 = arrows2d!(ax_force, ps_forward, vs_forward, minshaftlength = 1, color = (:blue, (state_id/nLoadSteps*9+1)*0.1))
    global elem2 = arrows2d!(ax_force, ps_inverse, vs_inverse, minshaftlength = 1, color = (:red, (state_id/nLoadSteps*9+1)*0.1) )
end
elem1, elem2 = LineElement(color = :blue),  LineElement(color = :red)
Legend(fig[1,2], [elem1, elem2], ["forward", "inverse"], labelsize = 25)
GLMakie.scatter!(ax_force, 1:nnodes, zeros(nnodes), color = :black, markersize = 10)
save(joinpath(output_dir,"comparison_Fy-" * base_inverse * ".png"), fig)

# Plot boxplot of the error
fig_ = Figure(size = (1000, 600))
ax_ = Axis(
    fig_[1, 1],
    title = "Boxplot", 
    titlesize = 25,
    xlabel = "node number", 
    ylabel = "absolute error",
    xlabelsize = 25,       # Increases X label font size
    ylabelsize = 25,       # Increases Y label font size
    xticklabelsize = 20,   # Increases X tick label font size
    yticklabelsize = 20    # Increases Y tick label font size
    )


# flatten the per-node vectors into one long vector
abs_error = vcat([abs.(direct.series["Fdir"][i, :] - inverse.series["Fy"][i, :]) for i in 1:nnodes]...)
# do the same for the x positions (repeat node index nLoadSteps times)
nnodes_list = vcat([fill(i, nLoadSteps) for i in 1:nnodes]...)
lab = " [N]"
if length(abs_error[1]) == 1
    scatter!(ax_, nnodes_list, abs_error, label="Fy"*lab, markersize = 15)
else
    boxplot!(ax_, nnodes_list, abs_error, label="Fy"*lab)
end
axislegend(ax_, labelsize = 25)
save(joinpath(output_dir,"error_boxplot_force-" * base_inverse * ".png"), fig_)


println("Plots saved to current directory.")


# # Select 3 timesteps to plot (e.g., first, middle, last)
# ntimesteps = length(direct.time)
# timesteps_to_plot = [1, div(ntimesteps, 2), ntimesteps]

# # Compute power spectral density for Y displacement of direct results
# using DSP

# # Get Y displacement data for all nodes
# y_displacements = direct.series["Y"]
# nnodes = size(y_displacements, 1)

# # Compute PSD for each node
# fig = Figure(size = (1200, 700))
# ax = Axis(fig[1, 1], title = "Power Spectral Density - Y Displacement", 
#           xlabel = "Frequency (Hz)", ylabel = "PSD ()", yscale = log10)

# for node in 1:nnodes
#     # Compute power spectral density using periodogram
#     psd = DSP.periodogram(y_displacements[node, :])
#     # Convert to dB scale with floor to prevent -Inf from log of zero/negative values
#     # psd_db = 10 .* log10.(max.(psd.power, 1e-12))  # Floor at 1e-12 to avoid log of zero
#     # # Filter out any remaining invalid values (NaN, Inf)
#     # valid_idx = isfinite.(psd_db)

#     lines!(ax, psd.freq, psd.power, label = "Node $node", linewidth = 2, color = (:red, (2*node-nnodes)^2/nnodes^2))
# end

# axislegend(ax, position = :rb)
# display(fig)
# save(joinpath(output_dir,"psd_y_displacement.png"), fig)
# println("PSD plot saved to current directory.")