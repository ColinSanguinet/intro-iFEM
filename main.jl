### Beam dynamic analysis with sinusoidal load and inverse analysis

# The dev branch of Muscade should be used. 2 solutions are available:
# 1) ] add Muscade#dev
# 2) ] dev Muscade ; open the Muscade repo (.julia/dev/Muscade) and checkout the dev branch.
#      If changes are made to the Muscade repo locally the use of Revise is advised.

using Muscade, StaticArrays, GLMakie, CSV, DataFrames, Interpolations, Statistics
using Muscade.Toolbox
include("save_results.jl")


##########################################
# Inputs    
##########################################


# Material Properties and Beam Geometry
#R   = 0.0;          # Radius of the bend [m]
EI₂ = 1e5;          # Bending stiffness [Nm²]
EI₃ = 1e5;          # Bending stiffness [Nm²]
EA  = 1e6;          # Axial stiffness [N]
GJ  = 1e6;          # Torsional stiffness [Nm²]
L   = 10.;          # Length of the beam [m]
μ   = 1.;           # Linear mass along main axis [kg/m]
ι₁  = 1.;           # Mass moment of inertia around main axis [kg·m³]
@functor with() g̃_null(t) = SVector(0.,0.,0.); # Inexistant gravity
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁ , g̃=g̃_null);

# General properties
t₁ = 0.0
t₂ = 1.0
mode = 2 #mode
# ω = (mode*π)^2 * (EI₂/(μ*L^4))^(.5) #(EI₂*mode^4*π^4/(μ*L^4))^(.5) #1 * 2*π/(t₂-t₁)
ω = 2*π/(t₂-t₁) * 2
σ = 0.001   # [m] ; noise level OR displacement sensor measurements precision 
β = 0.005     # [m] ; bias on displacement measurements

noise(x) = (rand()-0.5)*σ
maxiter     = 50
maxΔx       = 1e-5
maxΔu       = 1e-4
maxΔa       = 1e-6
maxΔλ       = Inf

# Mesh Definition Direct
nel         = 30    # Number of elements
nnodes      = nel+1 # Number of nodes
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,nnodes,2)); # Node coordinates

# Mesh Definition Inverse
nel_inv         = 40    # Number of elements
nnodes_inv      = nel_inv+1 # Number of nodes
nodeCoord_inv   = hcat((0:L/nel_inv:L),zeros(Float64,nnodes_inv,2)); # Node coordinates
# @show sensornodes     = [Int(floor((nnodes_inv+1)/4)), Int(floor((nnodes_inv+1)/4)) +2,  Int(floor((nnodes_inv+1)*3/4))]
@show sensornodes     = 1:nnodes_inv

# Costs properties
σₓ      = 1e-3 # [m]
σₓ₂     = 1e-2 # [m]
σᵤ      = 1e0 # [N]
σᵤ₂     = 1e-3 # [N]
Fᵤ      = -0.   # [N]
Fᵤ₂     = -0.   # [N]
Cx      = 1/(2*σₓ^2)
Cxother = 1/(2*σₓ₂^2)
Cu      = 1/(2*σᵤ^2)
Cuother = 1/(2*σᵤ₂^2)

# Direct Dynamic analysis properties
bDynamic_else_static = true
bDirectAnalysis = true
Δt₀                 = 0.01          # Initial time step [s]
time = t₁+Δt₀:Δt₀:t₂
nLoadSteps = length(time)
bNodalStaticForce   = true
bDistributedForceSinStatic = false
bNodalForceSin      = false
bNodalForceImpulse  = false
bDistributedForceSinDyn = false
bEigenAnalysis      = false
bPlanar             = bEigenAnalysis # Planar motion constraint for eigenvalue analysis
DirectSolver        = bDynamic_else_static ? SweepX{2} : SweepX{0}    # 2 for Dynamic solver, 1 for Static solver

# Load properties
point_load_position = L/2
_,node_number         = findmin(abs.(nodeCoord[:,1].-point_load_position)) # Node where the static load is applied
_,node_number_inv        = findmin(abs.(nodeCoord_inv[:,1].-point_load_position)) # Node where the static load is applied
F           = -1000.0      # Amplitude of the load [N]
q           = 0.0          # Uniform lateral load [N/m]
t_impulse   = 0.01        # Duration of the impulse load [s]

# Inverse analysis properties
Δtᵢₙᵥ             = 0.01               # Time step for the inverse analysis [s]
time_inv = bDynamic_else_static ? (t₁+Δt₀:Δtᵢₙᵥ:t₂) : (t₂-Δtᵢₙᵥ:Δtᵢₙᵥ:t₂)
nLoadSteps_inv = bDynamic_else_static ? length(time_inv) : 2
bInverseAnalysis = true
bNoise = false
InvSolver = bDynamic_else_static ? DirectXUA{2,0,0} : DirectXUA{0,0,0}

# Saving config
analysis_type = bDynamic_else_static ? "dynamic_" : "static_"
saveDirect = true
basename_direct = "direct_" * analysis_type * "fullSensing_" * string(nel) * "elts"
path_to_data = "20261023_" * basename_direct * ".csv"
# metadata_direct = add_struct_to_dict(
    # Dict{}(
    #     "nel" => nel,
    #     "nodeCoord" => nodeCoord,
    #     "typeOfLoad" => "pointload",
    #     "F" => F,
    #     # "nodenumber_of_force" => node_number
    # ), 
    # mat)
metadata_direct = Dict{}(
    "nel" => nel,
    "nodeCoord" => nodeCoord,
    "F" => F,
    "pointLoadPosition" => point_load_position,
)

saveInverse = true
basename_inverse = "inverse_" * analysis_type * "fullSensing_" * string(nel_inv) * "elts"
# metadata_inverse = add_struct_to_dict(
    # Dict{}(
    #     "nel" => nel_inv,
    #     "nodeCoord" => nodeCoord_inv,
    #     "typeOfLoad" => "pointload",
    #     "pointLoadPosition" => point_load_position,
    #     "estimated_force" => 0/2,
    #     "CostXfactor" =>        Cx,
    #     "CostXotherfactor" =>   Cxother,
    #     "CostUfactor" =>        Cu,
    #     "CostUotherfactor" =>   Cuother,
    #     "nodenumber_of_sensors" => sensornodes,

    # ), 
    # mat)
metadata_inverse = Dict{}(
    "nel" => nel_inv,
    "nodeCoord" => nodeCoord_inv,
    "pointLoadPosition" => point_load_position,
    "SigmaX" =>        σₓ ,
    "SigmaXother" =>   σₓ₂,
    "SigmaU" =>        σᵤ ,
    "SigmaUother" =>   σᵤ₂,
    "FU" =>        Fᵤ ,
    "FUother" =>   Fᵤ₂,
    "nodenumber_of_sensors" => sensornodes,
)


# Scaling
@show scale = (
    X=(
        t1=1e0, 
        t2=1e0, 
        t3=1e0, 
        r1=1e0, 
        r2=1e0, 
        r3=1e0, 
    ),
    U=(
        t1=1e0,
        t2=1e0,
        t3=1e0,
    )

)
@show Λscale = 1e0 # The bigger the more Importance given to Residual over Costs, so the more importance given to Physical Model on Data Model 

##########################################
# Solving
#########################################

# Direct model
#------------------------------------------

if bDirectAnalysis
    # Direct model
    name        = :directmodel
    model       = Model(name)
    nodid       = addnode!(model, nodeCoord)
    mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
    eleid       = addelement!(model, EulerBeam3D, mesh;mat=mat, orient2=SVector(0.,1.,0.))

    # Boundary conditions
    [addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]];           # Support at one end
    [addelement!(model,Hold,[nodid[nnodes]]  ;field) for field∈[:t1, :t2,:t3,:r1]];     # Support at the other end
    if bPlanar 
        [[addelement!(model,Hold,[nodid[i]] ;field) for field∈[:t3]] for i in 2:nnodes-1] # Planar motion constraint for eigenvalue analysis
    end

    # Loading conditions
    #------------------------------------------

    if bNodalForceSin
        bDynamic_else_static ? println("") : println("\n================================================\nWARNING : the chosen loading condition is appropriate for a dynamic case and you chose a static analysis.\n================================================")
        @functor with(F, t₂) NodalSin(t) = F * sin(2*π*t/t₂)
        @functor with() null(t) = 0.
        @functor with(node_number) costU(u, t, node_num, factor_node = 1.0, factor_other = 1.0) = node_num == node_number ? factor_node*(F * sin(2*π*t/t₂) - u)^2 : factor_other*u^2
        # addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalSin )
        requestable_load = [addelement!(model,DofLoad,[nodid[node]];field=:t2, value= node==node_number ? NodalSin : null) for node in 1:nnodes]
    end

    if bDistributedForceSinStatic
        bDynamic_else_static ? println("\n================================================\nWARNING : the chosen loading condition is appropriate for a static case and you chose a dynamic analysis.\n================================================") : println("")
        @functor with(F, L) DistSinStat(t, x_pos) = F * t * sin(2*π*x_pos/L)
        global applied_load = zeros(nnodes, length(time))
        requestable_load = [addelement!(model,DofLoad,[nodid[node]];field=:t2, value= DistSinStat, args=(x_pos = nodeCoord[node,1]) ) for node in 1:nnodes]
    end

    if bDistributedForceSinDyn
        bDynamic_else_static ? println("") : println("\n================================================\nWARNING : the chosen loading condition is appropriate for a dynamic case and you chose a static analysis.\n================================================")
        @functor with(F, t_impulse, L) NodalSin(t, x_pos) = t <= t_impulse ? (t / t_impulse) * F * sin(2*π*x_pos/L) : 0.0
        # @functor with(nodeCoord) costU(u, t, node_num, factor_node = 1.0) = node_num == -1 ? factor_node*( NodalSin(t, "zz" => nodeCoord[node_num, 1])-u)^2 : factor_node*u^2
        requestable_load = [addelement!(model,DofLoad,[nodid[node]];field=:t2,value= NodalSin, args=(x_pos = nodeCoord[node,1]) ) for node in 1:nnodes]
    end
    
    if bNodalForceImpulse
        bDynamic_else_static ? println("") : println("\n================================================\nWARNING : the chosen loading condition is appropriate for a dynamic case and you chose a static analysis.\n================================================")
        @functor with(F, t_impulse) NodalImpulse(t) = t <= t_impulse ? F * (1 - cos(π*t/t_impulse))/2 : 0.0
        # @functor with(F, t_impulse) NodalImpulse(t) = t <= t_impulse ? F * (1 - cos(2*π*t/t_impulse))/2 : 0.0
        @functor with() null(t) = 0.
        @functor with(node_number) costU(u, t, node_num, factor_node = 1.0, factor_other = 1.0) = node_num == node_number ? factor_node*( NodalImpulse(t)-u)^2 : factor_other*u^2
        # requestable_load = [addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalImpulse )]
        requestable_load = [addelement!(model,DofLoad,[nodid[node]];field=:t2,value= node==node_number ? NodalImpulse : null) for node in 1:nnodes]
    end
    
    if bNodalStaticForce
        bDynamic_else_static ? println("\n================================================\nWARNING : the chosen loading condition is appropriate for a static case and you chose a dynamic analysis.\n================================================") : println("")
        @functor with(F) NodalStatic(t) = F
        @functor with() null(t) = 0.
        @functor with(node_number) costU(u, t, node_num, factor_node = 1.0, factor_other = 1.0) = node_num == node_number ? factor_node*(F-u)^2 : factor_other*u^2
        # addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalStatic)
        requestable_load = [addelement!(model,DofLoad,[nodid[node]];field=:t2,value= node==node_number ? NodalStatic : null) for node in 1:nnodes]
    end

    # initialize
    #------------------------------------------
    
    # Initializing the model (sets all DoF to 0 at t=t₁)
    initialstate                = initialize!(model; time=t₁);

    if bEigenAnalysis
        state_stat   = solve(SweepX{0};initialstate,time=[t₁]);
        # Solve eigenvalue problem
        nmod            = 10
        res             = solve(EigX{ℝ}; state=state_stat[1],nmod);
    end

    # Solving in direct mode
    state                       = solve(DirectSolver;initialstate,time,verbose=true,maxΔx, maxiter);

    x_dir = [getdof(state[idxLoad];field=:t1,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps]
    y_dir = [getdof(state[idxLoad];field=:t2,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps]
    z_dir = [getdof(state[idxLoad];field=:t3,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps]
    r1_dir = [getdof(state[idxLoad];field=:r1,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps] 
    r2_dir = [getdof(state[idxLoad];field=:r2,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps] 
    r3_dir = [getdof(state[idxLoad];field=:r3,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps] 
    
    req = @request(F)
    F_dir = [[getresult(state[idxLoad], req, [requestable_load[node]])[1][1] for node in 1:nnodes] for idxLoad ∈ 1:nLoadSteps]

    time_saved = time
    
    if !bDynamic_else_static
        x_dir =  [x_dir[end] ]
        y_dir =  [y_dir[end] ]
        z_dir =  [z_dir[end] ]
        r1_dir = [r1_dir[end]]
        r2_dir = [r2_dir[end]]
        r3_dir = [r3_dir[end]]
        
        F_dir = [F_dir[end]]

        time_saved = [time_saved[end]]
    end

    # Save
    if saveDirect

        timeseries = Dict(
            "X" => x_dir,
            "Y" => y_dir,
            "Z" => z_dir,
            "R1" => r1_dir,
            "R2" => r2_dir,
            "R3" => r3_dir,
            "Fdir" => F_dir,
        )

        save_timeseries_csv(dated_base(basename_direct); metadata = metadata_direct, comps=timeseries, time=time_saved)
    end

    # # Draw
    # figure     = Figure(size = (1000,1000))
    # ax      = Axis3(figure[1,1],xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",aspect=:equal)
    # for to_draw in 1:1:nLoadSteps
    #     draw!(ax,state[to_draw];EulerBeam3D=(;nseg=20,  line_color= RGBf(1.0, to_draw/nLoadSteps, 0.)))
    # end
    # display(figure)
    # figure
end


# Inverse model
#------------------------------------------

if bInverseAnalysis
    name        = :inversemodel
    inv_model       = Model(name)
    data = load_timeseries_csv(path_to_data)
    x_dir = data.series["X"]
    y_dir = data.series["Y"]
    z_dir = data.series["Z"]
    r3_dir = data.series["R3"]
    time = data.time;
    nodeCoord = data.metadata["nodeCoord"]

    # Add noise to raw data before interpolation
    x_dir_noised =  [x_dir[node,:]  .+ (bNoise ? noise.(x_dir[node,:] ) .+ β : 0.0) for node in 1:nnodes]
    y_dir_noised =  [y_dir[node,:]  .+ (bNoise ? noise.(y_dir[node,:] ) .+ β : 0.0) for node in 1:nnodes]
    z_dir_noised =  [z_dir[node,:]  .+ (bNoise ? noise.(z_dir[node,:] ) .+ β : 0.0) for node in 1:nnodes]
    r3_dir_noised = [r3_dir[node,:] .+ (bNoise ? noise.(r3_dir[node,:]) .+ β : 0.0) for node in 1:nnodes]
    
    # Each entry in *_dir_noised is a vector (time series) for one node.
    # For quick plotting we extract the value at the final time for each node.
    x_dir_noised_last = [x_dir_noised[node][end] for node in 1:nnodes]
    y_dir_noised_last = [y_dir_noised[node][end] for node in 1:nnodes]
    z_dir_noised_last = [z_dir_noised[node][end] for node in 1:nnodes]
    r3_dir_noised_last = [r3_dir_noised[node][end] for node in 1:nnodes]

    # fig = Figure(size = (1000, 600))
    # ax = Axis(fig[1, 1], title = "y displacement", xlabel = "nodenumber", ylabel = "y disp")
    # GLMakie.scatter!(ax, 1:nnodes, y_dir_noised_last, label = "noised input (final time)", color = :purple, markersize = 10)
    # GLMakie.scatter!(ax, 1:nnodes, y_dir[:,end], label = "original input (final time)", color = :blue, markersize = 10)
    # axislegend(ax)
    # display(fig)

    # # Interpolate the noised dat
    # if bDynamic_else_static
    #     x_int  = [linear_interpolation(time_inv, x_dir_noised[node]) for node in 1:nnodes]
    #     y_int  = [linear_interpolation(time_inv, y_dir_noised[node]) for node in 1:nnodes]
    #     z_int  = [linear_interpolation(time_inv, z_dir_noised[node]) for node in 1:nnodes]
    #     r3_int = [linear_interpolation(time_inv, r3_dir_noised[node]) for node in 1:nnodes]
    # else
    #     x_int  = [linear_interpolation(time_inv, [x_dir_noised[node][1]  for i in 1:length(time_inv)]) for node in 1:nnodes]
    #     y_int  = [linear_interpolation(time_inv, [y_dir_noised[node][1]  for i in 1:length(time_inv)]) for node in 1:nnodes]
    #     z_int  = [linear_interpolation(time_inv, [z_dir_noised[node][1]  for i in 1:length(time_inv)]) for node in 1:nnodes]
    #     r3_int = [linear_interpolation(time_inv, [r3_dir_noised[node][1] for i in 1:length(time_inv)]) for node in 1:nnodes]
    # end

    # Interpolate the noised dat
    if bDynamic_else_static
        x_int   = interpolate((time_inv,nodeCoord[:,1]),  hcat([x_dir_noised[i] for i in 1:nnodes ]...),Gridded(Linear())) 
        y_int   = interpolate((time_inv,nodeCoord[:,1]),  hcat([y_dir_noised[i] for i in 1:nnodes ]...),Gridded(Linear())) 
        z_int   = interpolate((time_inv,nodeCoord[:,1]),  hcat([z_dir_noised[i] for i in 1:nnodes ]...),Gridded(Linear())) 
        r3_int  = interpolate((time_inv,nodeCoord[:,1]),  hcat([r3_dir_noised[i] for i in 1:nnodes ]...),Gridded(Linear()))
    else
        x_int   = interpolate((time_inv,nodeCoord[:,1]), vcat(hcat([x_dir_noised[i] for i in 1:nnodes ]...),hcat([x_dir_noised[i] for i in 1:nnodes ]...)),Gridded(Linear())) 
        y_int   = interpolate((time_inv,nodeCoord[:,1]), vcat(hcat([y_dir_noised[i]  for i in 1:nnodes ]...),hcat([y_dir_noised[i]  for i in 1:nnodes ]...)),Gridded(Linear())) 
        z_int   = interpolate((time_inv,nodeCoord[:,1]), vcat(hcat([z_dir_noised[i]  for i in 1:nnodes ]...),hcat([z_dir_noised[i]  for i in 1:nnodes ]...)),Gridded(Linear())) 
        r3_int  = interpolate((time_inv,nodeCoord[:,1]), vcat(hcat([r3_dir_noised[i] for i in 1:nnodes ]...),hcat([r3_dir_noised[i] for i in 1:nnodes ]...)),Gridded(Linear()))
    end

    # Mesh
    nodid       = addnode!(inv_model, nodeCoord_inv)
    mesh        = hcat(nodid[1:nnodes_inv-1],nodid[2:nnodes_inv])
    eleid       = addelement!(inv_model, EulerBeam3D, mesh;mat=mat, orient2=SVector(0.,1.,0.))

    # Boundary conditions
    [addelement!(inv_model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]];           # Support at one end
    [addelement!(inv_model,Hold,[nodid[nnodes_inv]]  ;field) for field∈[:t1, :t2,:t3,:r1]];     # Support at the other end
    if bPlanar
        [[addelement!(inv_model,Hold,[nodid[i]] ;field) for field∈[:t3]] for i in 2:nnodes_inv-1] # Planar motion constraint for eigenvalue analysis
    end

    # Costs
    @functor with() costX(x, t, meas) = Cx * (meas(t)-x)^2
    @functor with() costXother(x, t, meas) = Cxother * (meas(t)-x)^2
    @functor with() costU(u, t) = Cu *(Fᵤ-u)^2
    @functor with() costUother(u, t) = Cuother *(Fᵤ₂-u)^2
    e5             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t1,    cost= costXother, costargs= (meas = t -> x_int(t,nodeCoord_inv[node,1]) ,) ) for node in sensornodes]
    e6             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t2,    cost= costX, costargs= (meas = t -> y_int(t,nodeCoord_inv[node,1]),) ) for node in sensornodes]
    e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t3,    cost= costXother, costargs= (meas = t -> z_int(t,nodeCoord_inv[node,1]),) ) for node in sensornodes];
    e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:r3,    cost= costXother, costargs= (meas = t -> r3_int(t,nodeCoord_inv[node,1]),) ) for node in sensornodes];
    # e5             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t1,    cost= costXother, costargs= (meas = x_int[node] ,) ) for node in sensornodes]
    # e6             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t2,    cost= costX, costargs= (meas = y_int[node],) ) for node in sensornodes]
    # e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t3,    cost= costXother, costargs= (meas = z_int[node],) ) for node in sensornodes];
    # e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:r3,    cost= costXother, costargs= (meas = r3_int[node],) ) for node in sensornodes];
    e2             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t3,Ufield=:t3           ,    cost=costUother )  for node in 1:nnodes_inv];
    e3             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t2,Ufield=:t2           ,    cost=costUother )  for node in 1:node_number_inv-1];
    e3             = [addelement!(inv_model,SingleUdof,[nodid[node_number_inv]]; Xfield=:t2,Ufield=:t2    ,    cost=costU )];
    e3             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t2,Ufield=:t2           ,    cost=costUother )  for node in node_number_inv+1:nnodes_inv];
    e4             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t1,Ufield=:t1           ,    cost=costUother )  for node in 1:nnodes_inv];


    [[addelement!(inv_model,Hold,[nodid[i]] ;field) for field∈[:t3, :r2, :r1]] for i in 1:nnodes_inv]

    setscale!(inv_model; scale = scale, Λscale = Λscale)

    initialstate    = initialize!(inv_model;time=t₁)
    stateXUA         = solve(InvSolver;primerstate=[initialstate], time=[time_inv],verbose=true,maxiter,maxΔx,maxΔλ,maxΔu,maxΔa);

    x_inv =  [getdof(stateXUA[1][idxLoad];field=:t1,nodID=nodid[1:nnodes_inv]) for idxLoad ∈ 1:nLoadSteps_inv]
    y_inv =  [getdof(stateXUA[1][idxLoad];field=:t2,nodID=nodid[1:nnodes_inv]) for idxLoad ∈ 1:nLoadSteps_inv]
    z_inv =  [getdof(stateXUA[1][idxLoad];field=:t3,nodID=nodid[1:nnodes_inv]) for idxLoad ∈ 1:nLoadSteps_inv]
    r1_inv = [getdof(stateXUA[1][idxLoad];field=:r1,nodID=nodid[1:nnodes_inv]) for idxLoad ∈ 1:nLoadSteps_inv]
    r2_inv = [getdof(stateXUA[1][idxLoad];field=:r2,nodID=nodid[1:nnodes_inv]) for idxLoad ∈ 1:nLoadSteps_inv]
    r3_inv = [getdof(stateXUA[1][idxLoad];field=:r3,nodID=nodid[1:nnodes_inv]) for idxLoad ∈ 1:nLoadSteps_inv]

    U_t_inv = [getdof(stateXUA[1][idxLoad];class = :U, field=:t2,nodID=nodid[1:nnodes_inv]) for idxLoad ∈ 1:nLoadSteps_inv]

    time_inv_saved = time_inv

    if !bDynamic_else_static
        x_inv   = [  x_inv[end]] 
        y_inv   = [  y_inv[end]] 
        z_inv   = [  z_inv[end]] 
        r1_inv  = [ r1_inv[end]]
        r2_inv  = [ r2_inv[end]]
        r3_inv  = [ r3_inv[end]]
        U_t_inv = [U_t_inv[end]] 

        time_inv_saved = [time_inv_saved[end]]
    end

    # Save
    if saveInverse
        metadata_inverse["X_input"] = x_dir_noised
        metadata_inverse["Y_input"] =  y_dir_noised
        metadata_inverse["Z_input"] = z_dir_noised
        metadata_inverse["R3_input"] = r3_dir_noised     
        

        timeseries = Dict(
            "X" => x_inv,
            "Y" => y_inv,
            "Z" => z_inv,
            "R1" => r1_inv,
            "R2" => r2_inv,
            "R3" => r3_inv,
            "Fy" => U_t_inv,
        )

        save_timeseries_csv(dated_base(basename_inverse); metadata = metadata_inverse, comps=timeseries, time=time_inv_saved)
    end
end



### Plotting
include("analyse_results.jl")
date = "20261023"
@show filename_direct = date * "_" * basename_direct * ".csv"
@show filename_inverse = date * "_" * basename_inverse * ".csv"
# Create results folder path
base_direct = split(filename_direct, '/')[end][10:end-4]
base_inverse = split(filename_inverse, '/')[end][10:end-4]
output_dir = joinpath(@__DIR__, "results", base_inverse)

# Components to plot for the boxplot
component_dict = Dict(
    "X" => " [m]",
    "Y" => " [m]",
    "Z" => " [m]",
    "R3" => " [rad]",
)

direct, inverse = load_analysis(filename_direct, filename_inverse)

nSteps = length(direct.time)
rng = nSteps:1:nSteps
 
nodeToLookAt = 16

# calculate_RMSE(direct, inverse, nSteps)
# plot_eachComponent(direct, inverse, nodeToLookAt, component_dict, output_dir, "node" * "$nodeToLookAt" * "_withTime_" * base_inverse )
plot_boxplot_disp(direct, inverse, component_dict, output_dir, "boxplotDisp_" * base_inverse)
plot_deflection(direct, inverse, rng, output_dir, "deflection_" * base_inverse)
plot_forceComparison(direct, inverse, rng, output_dir, "forces_" * base_inverse)
plot_boxplot_force(direct, inverse, output_dir, "boxplotForce_" * base_inverse)
plot_forces(direct, inverse, nodeToLookAt, output_dir, "forcewithTime_" * base_inverse)
calculate_forceReconstructionAccuracy(direct, inverse, output_dir, "errorForce_" * base_inverse)

println("Plots saved to current directory.")
