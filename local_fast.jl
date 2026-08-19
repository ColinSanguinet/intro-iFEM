### Beam dynamic analysis with sinusoidal load and inverse analysis

# The dev branch of Muscade should be used. 2 solutions are available:
# 1) ] add Muscade#dev
# 2) ] dev Muscade ; open the Muscade repo (.julia/dev/Muscade) and checkout the dev branch.
#      If changes are made to the Muscade repo locally the use of Revise is advised.

using Muscade, StaticArrays, GLMakie, CSV, DataFrames, Interpolations, Statistics, FFTW
using Muscade.Toolbox, CairoMakie

##########################################
# Inputs
##########################################


# Material Properties and Beam Geometry
R   = 0.0;          # Radius of the bend [m]
EI₂ = 1e5;          # Bending stiffness [Nm²]
EI₃ = 1e5;          # Bending stiffness [Nm²]
EA  = 1e6;          # Axial stiffness [N]
GJ  = 1e6;          # Torsional stiffness [Nm²]
L   = 10.;          # Length of the beam [m]
μ   = 1.;           # Linear mass along main axis [kg/m]
ι₁  = 1.;           # Mass moment of inertia around main axis [kg·m³]
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁,w=μ*9.81);

# General properties
t₁ = 0.0
t₂ = 2.0
maxiter     = 10
maxiter_inv = 20
maxΔx       = 1e-6
maxΔu       = 1e-6
maxΔa       = 1e-6
maxΔx_inv       = 1e-4
maxΔu_inv       = 9e-4
maxΔa_inv       = 1e-4
maxΔλ       = Inf

# Mesh Definition Direct
nel         = 50                    # Number of elements
nnodes      = nel+1                 # Number of nodes
ϵ = 1e-3 # Small offset to avoid zero coordinates for interpolation
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,nnodes,2)); # Node coordinates

# Mesh Definition Inverse
nel_inv         = 50                # Number of elements
nnodes_inv      = nel_inv+1         # Number of nodes
nodeCoord_inv   = hcat((0:L/nel_inv:L),zeros(Float64,nnodes_inv,2)); # Node coordinates

# Direct Dynamic analysis properties
Δt₀                 = 0.005         # Initial time step [s]
dyn_time = t₁+Δt₀:Δt₀:t₂
γ = 0.505                             # Newmark parameter
β = 1/3                             # Newmark parameter


ID ="DynLC3"
bNodalForceImpulse  = false
bNodalForceSin      = false
bNodalStaticForce   = false
bDistributedForceSinStatic = false
bDistributedForceSinDyn    = true
node_number        = Int(floor(nnodes/2))+1 # Node where the static load is applied

bEigenAnalysis      = false         # Enable/disable eigenvalue analysis
bPlanar             = true          # Planar motion constraint for eigenvalue analysis
stat_or_dyn         = 2             # 0 for static analysis, 2 for dynamic analysis
DirectSolver        = SweepX{stat_or_dyn}   

if stat_or_dyn == 2
    time = dyn_time
elseif stat_or_dyn == 1
    time = dyn_time
else
    time = t₁+Δt₀:Δt₀:t₂
end
nLoadSteps = length(time)

# Load properties
F           = -10.0      # Amplitude of the load [N] -10 for distributed load, -1000 for nodal static load
q           = 0.0          # Uniform lateral load [N/m]
t_impulse   = 0.1        # Duration of the impulse load [s]

# Inverse analysis properties
Δtᵢₙᵥ             = 0.005               # Time step for the inverse analysis [s]
time_inv = dyn_time
nLoadSteps_inv = length(time)
bInverseAnalysis = true
InvSolver        = DirectXUA{stat_or_dyn,0,0}   # Dynamic solver for the inverse analysis

# Noise properties
bNoise = false                         # Enable/disable noise
ϵ = 0.00                             # Noise level (fraction of the maximum displacement)

# Post-processing of the direct analysis
bSaveFigures = false
bShow3DBeam = false
titlesize = 55.0f0
labelsize = 50.0f0
legendsize = 34.0f0
ticksize = 36.0f0


##########################################
# Direct model
##########################################


# Direct model
name        = :SinusoidalLoadEigen
model       = Model(name)
nodid       = addnode!(model, nodeCoord)
mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
eleid       = addelement!(model, EulerBeam3D, mesh;mat=mat, orient2=SVector(0.,1.,0.))

# Boundary conditions
[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]];           # Support at one end
[addelement!(model,Hold,[nodid[nnodes]]  ;field) for field∈[:t1,:t2,:t3,:r1]];     # Support at the other end
if bPlanar 
    [[addelement!(model,Hold,[nodid[i]] ;field) for field∈[:t3]] for i in 2:nnodes-1] # Planar motion constraint for eigenvalue analysis
end

# Loading conditions
if bNodalForceSin
    @functor with(F, t₂) NodalSin(t) = F * sin(0.5*π*t/t₂)
    @functor with(node_number) costU(u, t, node_num, factor_node = 1.0, factor_other = 1.0) = node_num == node_number ? factor_node*(F * sin(2*π*t/t₂) - u)^2 : factor_other*u^2
    requestable_load = [addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalSin )]
end

if bDistributedForceSinStatic
    @functor with(F, L, nodeCoord) DistSinStat(t, x_pos) = F * t * sin(2*π*x_pos.second/L)
    for node in 1:nnodes
        addelement!(model,DofLoad,[nodid[node]];field=:t2, value= DistSinStat, valueargs=(x_pos = nodeCoord[node,1]) )
    end
end

if bDistributedForceSinDyn
    @functor with(F, t_impulse, L) NodalSin(t, x_pos) = t <= t_impulse ? - (t / t_impulse) * F * sin(2*π*x_pos.second/L) : 0.0
    @functor with(nodeCoord) costU(u, t, node_num, factor_node = 1.0) = node_num != -1 ? factor_node*( NodalSin(t, "zz" => nodeCoord[node_num, 1])-u)^2 : factor_node*u^2
    requestable_load = [addelement!(model,DofLoad,[nodid[node]];field=:t2,value= NodalSin, valueargs=(x_pos = nodeCoord[node,1]) ) for node in 1:nnodes]
end

if bNodalForceImpulse
    @functor with(F, t_impulse) NodalImpulse(t) = t <= t_impulse ? F * (1 - cos(2*π*t/t_impulse))/2 : 0.0
    @functor with(node_number) costU(u, t, node_num, factor_node = 1.0, factor_other = 1.0) = node_num == node_number ? factor_node*( NodalImpulse(t)-u)^2 : factor_other*u^2
    addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalImpulse )
end

if bNodalStaticForce
    @functor with(F) NodalStatic(t) = t*F
    @functor with(node_number) costU(u, t, node_num, factor_node = 1.0, factor_other = 1.0, F_est = -1000, F_other = 0.0) = node_num == node_number ? factor_node*(NodalStatic(t)-u)^2 : factor_other*(F_other-u)^2
    requestable_load = [addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalStatic )]
end

# Initializing the model (sets all DoF to 0 at t=t₁)
initialstate                = initialize!(model; time=t₁);

if bEigenAnalysis
    fₙ(k) = √(EI₂/μ)*(k^2*π)/(2*L^2)
    Φₙ(k,x) = sin.(k*π/L.*x)
    state_stat = solve(SweepX{0}; initialstate, time=[t₁])

    nmod = 5  # request a few extra modes
    res  = solve(EigX{ℝ}; state=state_stat[1], nmod)

    fig  = Figure(size = (1600, 900))
    ax   = Axis(fig[1, 1],
                xlabel = "Position, x (m)",
                xlabelsize = labelsize,
                xticklabelsize = ticksize,
                ylabel = "Displacement, y (m)",
                ylabelsize = labelsize,
                yticklabelsize = ticksize,
                title = "Solver: EigX{ℝ}, Modes: 1-3",
                titlefont = :italic,
                titlesize = titlesize)
    mode_prec = Vector{Float64}([0.0 for i in 1:nnodes])
    for idxMod in 1:3
        eigres = increment(state_stat[1], res, [idxMod], [1.0])
        t2_eig = getdof(eigres; field=:t2, nodID=nodid[1:nnodes])
        δ = maximum(abs.(t2_eig[2:end-1]))
        if δ ≈ 0
            @warn "Mode $idxMod has near-zero amplitude; skipping"
            continue
        end
        sgn     = sign((Φₙ(idxMod, 0:L/nel:L)' * t2_eig))
        t2_norm = sgn .* (t2_eig/δ)
        labelStr = "Mode $idxMod, Muscade: $(round(res.ω[idxMod]/(2π), digits=3)) Hz, Analyt.: $(round(fₙ(idxMod), digits=3)) Hz"
        scatter!(ax, 0:L/nel:L, t2_norm; label=labelStr)
        lines!(ax, 0:L/nel:L, Φₙ(idxMod, 0:L/nel:L))
    end
    
    xlims!(ax, 0, L); ylims!(ax, -1.2, 1.2)
    axislegend(ax, merge = true, unique = false, position = :lb, labelsize = legendsize)
    display(fig)
    save("output/Eig_Beam.png", fig)
end

# Solving in direct mode
state = solve(DirectSolver;initialstate,time,verbose=false,maxΔx, maxiter, γ, β);

x_dir = [getdof(state[1:nLoadSteps];field=:t1,nodID=[nodid[node]]) for node ∈ 1:nnodes]
y_dir = [getdof(state[1:nLoadSteps];field=:t2,nodID=[nodid[node]]) for node ∈ 1:nnodes]
z_dir = [getdof(state[1:nLoadSteps];field=:t3,nodID=[nodid[node]]) for node ∈ 1:nnodes]
r1_dir = [getdof(state[1:nLoadSteps];field=:r1,nodID=[nodid[node]]) for node ∈ 1:nnodes] 
r2_dir = [getdof(state[1:nLoadSteps];field=:r2,nodID=[nodid[node]]) for node ∈ 1:nnodes] 
r3_dir = [getdof(state[1:nLoadSteps];field=:r3,nodID=[nodid[node]]) for node ∈ 1:nnodes] 

req = @request F
out = getresult(state[1:nLoadSteps], req, requestable_load)



##########################################
# Decay Plots
##########################################

plot_figures_dir = true

if plot_figures_dir
    node_2_plot = 13 # 26 is center node for 50 elements
    time_2_plot = 0.4
    timestep_2_plot = findfirst(t -> t >= time_2_plot, time)
    act_time_2_plot = time[timestep_2_plot]

    # Displacement of one node over time
    fig_decay = Figure(size = (2000,1000))
    ax1 = Axis(fig_decay[1, 1],
                xlabel = "Time [s]",
                xlabelsize = labelsize,
                xticklabelsize = ticksize,
                ylabel = "Displacement, y [m]",
                ylabelsize = labelsize,
                yticklabelsize = ticksize,
                title = "Solver: $(DirectSolver), β = $β, γ = $γ, Δt = $Δt₀ s, Node: $node_2_plot",
                titlefont = :italic,
                titlesize = titlesize)
    lines!(ax1, time, vcat(y_dir[node_2_plot]...); color = :blue, label = "Direct Analysis")
    axislegend(ax1, merge = true, unique = false)
    display(fig_decay)
    save("output/Dir_Beam_Decay_$(ID).pdf", fig_decay; backend = CairoMakie)

    # Beam visualization at a given timestep
    fig_beam = Figure(size = (2000,1000))
    ax1 = Axis(fig_beam[1, 1],
                xlabel = "Displacement, x [m]",
                xlabelsize = labelsize,
                xticklabelsize = ticksize,
                ylabel = "Displacement, y [m]",
                ylabelsize = labelsize,
                yticklabelsize = ticksize,
                title = "Solver: $(DirectSolver), β = $β, γ = $γ, Δt = $Δt₀ s, Time: $(round(act_time_2_plot, digits=3)) s",
                titlefont = :italic,
                titlesize = titlesize)
    draw!(ax1,state[timestep_2_plot];EulerBeam3D=(;nseg=20,  line_color= RGBf(0.0, 0.0, 1.0 )))
    display(fig_beam)
    save("output/Dir_Beam_T_$(timestep_2_plot)_$(ID).png", fig_beam)

    # FFT of the displacement of a given node after release
    t_fft = 0.2:Δt₀:t₂
    n_impulse = length(time) - length(t_fft) + 1
    fs = 1/Δt₀
    F = fftshift(fft(y_dir[node_2_plot][n_impulse:end]))
    freqs =  fftshift(fftfreq(length(t_fft), fs))
    fig_ot = Figure(size = (2000,1000))
    ax1 = Axis(fig_ot[1, 1], xlabel = "Time [s]", xlabelsize = labelsize, xticklabelsize = ticksize, ylabel = "Displacement, y [m]", ylabelsize = labelsize, yticklabelsize = ticksize  )
    time_domain = lines!(ax1, t_fft, y_dir[node_2_plot][n_impulse:end])
    ax2 = Axis(fig_ot[1, 2], xlabel = "Frequency [Hz]", xlabelsize = labelsize, xticklabelsize = ticksize, ylabel = "Amplitude, [m.s]", ylabelsize = labelsize, yticklabelsize = ticksize)
    xlims!(ax2, 0, min(50,1/Δt₀/2))
    freq_domain = lines!(ax2, freqs, abs.(F))
    ax0 = Label(fig_ot[0, :], "Solver: $(DirectSolver), β = $β, γ = $γ, Δt = $Δt₀ s, Node: $node_2_plot", font = :italic, fontsize = titlesize)
    display(fig_ot)
    save("output/Dir_FFT_N_$(node_2_plot)_$(ID).png", fig_ot)
end

######################################
#           Inverse model            #
######################################

if bInverseAnalysis
    nnodes = nnodes_inv
    nodeCoord = nodeCoord_inv

    name        = :BeamDynSinusoidalLoads
    inv_model   = Model(name)
    nodid       = addnode!(inv_model, nodeCoord)
    mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
    eleid       = addelement!(inv_model, EulerBeam3D, mesh;mat=mat, orient2=SVector(0.,1.,0.) )

    ####################################
    #     Set boundary conditions      #
    ####################################

    [addelement!(inv_model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]];           # Support at one end
    [addelement!(inv_model,Hold,[nodid[nnodes]]  ;field) for field∈[:t1, :t2,:t3,:r1]];     # Support at the other end
    if bPlanar
        [[addelement!(inv_model,Hold,[nodid[i]] ;field) for field∈[:t3]] for i in 2:nnodes-1] # Planar motion constraint
    end
    

    ####################################
    #          Input Data              #
    ####################################
    
    if bNoise
        using Random
        Random.seed!(1234)  # Set seed for reproducibility (optional)
        x_noisy = [vcat(x_dir[node]...) .+ ϵ .* randn(nLoadSteps) for node in 2:nnodes-1]
        y_noisy = [vcat(y_dir[node]...) .+ ϵ .* randn(nLoadSteps) for node in 2:nnodes-1]
        z_noisy = [vcat(z_dir[node]...) for node in 1:nnodes]
        r3_noisy = [vcat(r3_dir[node]...) for node in 1:nnodes]
        
        x_int = [linear_interpolation(time, vcat(x_dir[1]...)), [linear_interpolation(time, vcat(x_noisy[node])) for node in 1:nnodes-2]..., linear_interpolation(time, vcat(x_dir[end]...))]
        y_int = [linear_interpolation(time, vcat(y_dir[1]...)), [linear_interpolation(time, vcat(y_noisy[node])) for node in 1:nnodes-2]..., linear_interpolation(time, vcat(y_dir[end]...))]
        z_int = [linear_interpolation(time, z_noisy[node]) for node in 1:nnodes]
        r3_int = [linear_interpolation(time, r3_noisy[node]) for node in 1:nnodes]
    else
        x_int = [linear_interpolation(time, vcat(x_dir[node]...)) for node in 1:nnodes]
        y_int = [linear_interpolation(time, vcat(y_dir[node]...)) for node in 1:nnodes]
        z_int = [linear_interpolation(time, vcat(z_dir[node]...)) for node in 1:nnodes]
        r3_int = [linear_interpolation(time, vcat(r3_dir[node]...)) for node in 1:nnodes]
    end

    ####################################
    #             Costs                #
    ####################################

    @functor with() costX(x, t, meas) = 1000000000 * (meas(t)-x)^2
    @functor with() costXother(x, t, meas) = 1000 * (meas(t)-x)^2
    e5             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t1,    cost= costXother, costargs= (meas = x_int[node],) ) for node in 1:nnodes]
    e6             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t2,    cost= costX, costargs= (meas = y_int[node],) ) for node in 1:nnodes]
    e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t3,    cost= costXother, costargs= (meas = z_int[node],) ) for node in 1:nnodes];
    e8             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:r3,    cost= costXother, costargs= (meas = r3_int[node],) ) for node in 1:nnodes];
    e2             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t2,Ufield=:t2           ,    cost=costU, costargs=(node_num=node, factor_node=1000.0) )  for node in 2:nnodes-1];
    e3             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t3,Ufield=:t3           ,    cost=costU, costargs=(node_num=-1, factor_node=1000.0) )  for node in 2:nnodes-1];
    e4             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t1,Ufield=:t1           ,    cost=costU, costargs=(node_num=-1, factor_node=1000.0) )  for node in 2:nnodes-1];

    ####################################
    #             Solving              #
    ####################################

    initialstate    = initialize!(inv_model;time=t₁)
    stateXUA         = solve( InvSolver; initialstate=[initialstate], time= [time],maxiter= maxiter_inv, maxΔx= maxΔx_inv, maxΔu= maxΔu_inv, maxΔa= maxΔa_inv, maxΔλ= maxΔλ, verbose=true);

    ####################################
    # Data Processing for plotting     #
    ####################################

    x_inv =  [getdof(stateXUA[1][1:nLoadSteps_inv];field=:t1,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    y_inv =  [getdof(stateXUA[1][1:nLoadSteps_inv];field=:t2,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    z_inv =  [getdof(stateXUA[1][1:nLoadSteps_inv];field=:t3,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    r1_inv = [getdof(stateXUA[1][1:nLoadSteps_inv];field=:r1,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    r2_inv = [getdof(stateXUA[1][1:nLoadSteps_inv];field=:r2,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    r3_inv = [getdof(stateXUA[1][1:nLoadSteps_inv];field=:r3,nodID=[nodid[node]]) for node ∈ 1:nnodes]

    U_t_inv = [getdof(stateXUA[1][idxLoad];class = :U, field=:t2,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps_inv]
    
    
    ###################################
    #              Draw               #
    ###################################

    time_2_plot = 0.09
    time_plot = findfirst(t -> t >= time_2_plot, time_inv)
    act_time_plotted = time_inv[time_plot]
    node_2_plot = 13 # 26 is center node for 50 elements
    U_inv_nodes = [[0.0, U_t_inv[i]..., 0.0] for i in 1:nLoadSteps_inv]


    fig_2D = Figure(size = (2000,1000))
    cmap = cgrad([RGBf(1.0, 0.9, 0.0), RGBf(1.0, 0.2, 0.0)])
    crange = (0.0, 1.0)
    title = Label(fig_2D[0, 1:3], L"Solver: %$InvSolver, ϵ = %$ϵ, Time: %$(round(act_time_plotted, digits=3)) s", fontsize = 50.0f0)
    ax = Axis(  fig_2D[1, 1],
                xlabel = "Position, x [m]",
                xlabelsize = labelsize,
                xticklabelsize = ticksize,
                ylabel = "Displacement, y [m]",
                ylabelsize = labelsize,
                yticklabelsize = ticksize)
    ax2 = Axis( fig_2D[1, 1],
                ylabelsize = labelsize,
                yticklabelsize = ticksize,
                yaxisposition = :right,
                ylabel = L"\text{Relative Error, } \Delta y_{\mathrm{rel}}\ [-]",
                yticklabelcolor = :green)
    ax3 = Axis(  fig_2D[1, 2],
                xlabel = "Position, x (m)",
                xlabelsize = labelsize,
                xticklabelsize = ticksize,
                ylabel = "Force amplitude, Uₜ₂ [N]",
                ylabelsize = labelsize,
                yticklabelsize = ticksize)
    lines!(ax, nodeCoord[1:nnodes, 1],  [y_inv[i][time_plot] for i in 1:nnodes]; color = :tomato, label = "Inverse Analysis")
    scatter!(ax, nodeCoord[1:nnodes, 1],  [y_inv[i][time_plot] for i in 1:nnodes]; color = :tomato, label = "Inverse Analysis")
    scatter!(ax, nodeCoord[1:nnodes, 1],  [y_int[i](time_inv[time_plot]) for i in 1:nnodes]; color = :blue, label = "Inverse Input")
    lines!(ax, nodeCoord[1:nnodes, 1],  [y_dir[i][time_plot] for i in 1:nnodes]; color = :blue, label = "Direct Analysis")
    lines!(ax2, nodeCoord[1:nnodes, 1],  [(y_inv[i][time_plot]-y_dir[i][time_plot])/(y_dir[i][time_plot]+1e-12) for i in 1:nnodes]; color = :green, label = "Difference")
    
    hideydecorations!(ax2, grid = true, ticks = false, ticklabels = false, label = false)
    hidexdecorations!(ax2, grid = true)
    axislegend(ax, merge = true, unique = false, position = :lb, labelsize = legendsize)
    axislegend(ax2, merge = true, unique = false, position = :rb, labelsize = legendsize)
    ylims!(ax2, -1, 1)
    sc_for_cb = nothing
    for i in 1:nLoadSteps_inv
            τ = nLoadSteps_inv == 1 ? 0.0 : (i - 1) / (nLoadSteps_inv - 1)
            sc_for_cb = scatter!(ax3, 0:L/nel_inv:L, U_inv_nodes[i];
                color = fill(τ, length(U_inv_nodes[i])), colormap = cmap, colorrange = crange)
        end
        Colorbar(fig_2D[1, 3], sc_for_cb,
            flipaxis = false, label = "Time, t [s]", labelsize = legendsize, ticklabelsize = ticksize)
    display(fig_2D; backend = GLMakie)
    save("output/Inv_Beam_Disp_$(ID).pdf", fig_2D)
    
    if false
        fig_disp = Figure(size = (2000,1000))
        ax = Axis(fig_disp[1, 1], xlabel = "Time [s]", xlabelsize = 22.0f0, xticklabelsize = 18.0f0, ylabel = "Displacement, y [m]", ylabelsize = 22.0f0, yticklabelsize = 18.0f0)
        [lines!(ax, time_inv, vcat(y_inv[25]...); color = :tomato, label = "Inverse Analysis")]
        [lines!(ax, time_inv,[y_int[25](i) for i in time_inv]; color = :blue, label = "Direct Analysis")]
        axislegend(ax, merge = true, unique = false)
        display(fig_disp; backend = GLMakie)
    end


    if false
        figure     = Figure(size = (2000,1000))
        ax      = Axis3(figure[1,1],xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",aspect=:equal)
        for to_draw in 1:10:nLoadSteps
        draw!(ax,stateXUA[1][to_draw];EulerBeam3D=(;nseg=20,  line_color= RGBf(1.0, to_draw/nLoadSteps, 0.)))
        draw!(ax,state[to_draw];EulerBeam3D=(;nseg=20,  line_color= RGBf(0.0, 0.0, 1.0 )))
        end
        display(figure; backend = GLMakie)
        figure
    end

    if true
        U_t_inv_nodes = [U_t_inv[i][node_2_plot-1] for i in 1:nLoadSteps_inv]
        U_t_dir = [out[node_2_plot,i].F for i in 1:nLoadSteps_inv]
        fig_disp = Figure(size = (2000,1000))
        ax = Axis(  fig_disp[1, 1],
                    xlabel = "Time, t [s]",
                    xlabelsize = labelsize,
                    xticklabelsize = ticksize,
                    ylabel = "Force amplitude, Uₜ₂ [N]",
                    ylabelsize = labelsize,
                    yticklabelsize = ticksize,
                    title = "Solver: $InvSolver, ϵ = $ϵ, Node: $node_2_plot",
                    titlefont = :italic,
                    titlesize = titlesize)
        [scatter!(ax, time, U_t_inv_nodes; color = :tomato, label = "Inverse Analysis")]
        [lines!(ax, time, U_t_dir; color = :blue, label = "Direct Analysis")]
        axislegend(ax, merge = true, unique = false, position = :rt, labelsize = legendsize)
        display(fig_disp)
        save("output/Inv_Beam_Load_Evol_$(ID).png", fig_disp)
    end

    if true
        U_inv_nodes = [[0.0, U_t_inv[i]..., 0.0] for i in 1:nLoadSteps_inv]
        fig_disp = Figure(size = (2000,1000))
        cmap = cgrad([RGBf(1.0, 0.9, 0.0), RGBf(1.0, 0.2, 0.0)])
        crange = (0.0, 1.0)
        ax = Axis(  fig_disp[1, 1],
                    xlabel = "Position, x (m)",
                    xlabelsize = labelsize,
                    xticklabelsize = ticksize,
                    ylabel = "Force amplitude, Uₜ₂ [N]",
                    ylabelsize = labelsize,
                    yticklabelsize = ticksize,
                    title = "Solver: $InvSolver, ϵ = $ϵ",
                    titlefont = :italic,
                    titlesize = titlesize)
        sc_for_cb = nothing
        for i in 1:nLoadSteps_inv
            τ = nLoadSteps_inv == 1 ? 0.0 : (i - 1) / (nLoadSteps_inv - 1)
            sc_for_cb = scatter!(ax, 0:L/nel_inv:L, U_inv_nodes[i];
                color = fill(τ, length(U_inv_nodes[i])), colormap = cmap, colorrange = crange)
        end
        Colorbar(fig_disp[1, 2], sc_for_cb,
            flipaxis = false, label = "Time, t [s]", labelsize = legendsize, ticklabelsize = ticksize)
        display(fig_disp)
        save("output/Inv_Beam_Loads_$(ID).png", fig_disp)
    end
end

req = @request gp(resultants(mᵢ))
o = getresult(stateXUA[1][400],req,eleid)
Fgp1_ = [ o[idxEl].gp[1][:resultants][:mᵢ] for idxEl ∈ 1:nel]
Fgp2_ = [ o[idxEl].gp[2][:resultants][:mᵢ] for idxEl ∈ 1:nel]
Fgp3_ = [ o[idxEl].gp[3][:resultants][:mᵢ] for idxEl ∈ 1:nel]
Fgp4_ = [ o[idxEl].gp[4][:resultants][:mᵢ] for idxEl ∈ 1:nel];

fig      = Figure(size = (1000,1000))
xgp1 = (1. /nel)*( (0.5-1/2*sqrt(3/7+2/7*sqrt(6/5))) :1:nel)
xgp2 = (1. /nel)*( (0.5-1/2*sqrt(3/7-2/7*sqrt(6/5))) :1:nel)
xgp3 = (1. /nel)*( (0.5+1/2*sqrt(3/7-2/7*sqrt(6/5))) :1:nel)
xgp4 = (1. /nel)*( (0.5+1/2*sqrt(3/7+2/7*sqrt(6/5))) :1:nel)
xgps = [xgp1;xgp2;xgp3;xgp4];
ax=Axis(fig[1,1], ylabel="Forces F [N]",       yminorgridvisible = true,xminorgridvisible = true,xticks = (0:1. /nel:1))
scatter!(ax,xgps,  [[Fgp1_[iel][2] for iel=1:nel] ; [Fgp2_[iel][2] for iel=1:nel] ; [Fgp3_[iel][2] for iel=1:nel] ; [Fgp4_[iel][2] for iel=1:nel]],          label="F₁");
display(fig)