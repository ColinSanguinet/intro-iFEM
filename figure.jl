# Juste un fichier pour homogénéiser le format des figures:


### Template 2D
fig_disp = Figure(size = (1000,1000)) # Ratio 1:1 préférable
ax = Axis(fig_disp[1, 1], xlabel = "Time [s]", xlabelsize = 22.0f0, xticklabelsize = 18.0f0, ylabel = "Displacement, y [m]", ylabelsize = 22.0f0, yticklabelsize = 18.0f0)
[lines!(ax, time_dir, dir_timeseries; color = :blue, label = "Direct Result")] # Couleurs froides pour direct
[lines!(ax, time_inv, inv_timeseries; color = :tomato, label = "Inverse Result")] # Couleurs chaudes pour sol inverse
[scatter!(ax, time_inv, inv_timeseries; color = :tomato, label = "Inverse Result")] 
axislegend(ax, merge = false, unique = false)


### Template 3D
# A venir 