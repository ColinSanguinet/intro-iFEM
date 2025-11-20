### Beam dynamic analysis with sinusoidal load and inverse analysis

# The dev branch of Muscade should be used. 2 solutions are available:
# 1) ] add Muscade#dev
# 2) ] dev Muscade ; open the Muscade repo (.julia/dev/Muscade) and checkout the dev branch.
#      If changes are made to the Muscade repo locally the use of Revise is advised.

using Muscade, StaticArrays, GLMakie, CSV, DataFrames, Interpolations
using Muscade.Toolbox

# Material Properties and Beam Geometry
#R   = 0.0;          # Radius of the bend [m]
EI₂ = 1e5;          # Bending stiffness [Nm²]
EI₃ = 1e5;          # Bending stiffness [Nm²]
EA  = 1e6;          # Axial stiffness [N]
GJ  = 1e6;          # Torsional stiffness [Nm²]
L   = 10.;          # Length of the beam [m]
μ   = 1.;           # Linear mass along main axis [kg/m]
ι₁  = 1.;           # Mass moment of inertia around main axis [kg·m³]
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=μ,ι₁=ι₁);

# Mesh Definition
nel         = 50    # Number of elements
nnodes      = nel+1 # Number of nodes
nodeCoord   = hcat((0:L/nel:L),zeros(Float64,nnodes,2)); # Node coordinates

# Direct Dynamic analysis properties
bNodalForceImpulse  = false
bNodalForceSin      = false
bNodalStaticForce   = true
node_number        = Int(floor(nnodes/2)) # Node where the static load is applied
bEigenAnalysis      = false
bPlanar             = bEigenAnalysis # Planar motion constraint for eigenvalue analysis
DirectSolver        = SweepX{0}     # 2 for Dynamic solver, 1 for Static solver
Δt₀                 = 0.01          # Initial time step [s]

# Load properties
F           = -1000.0      # Amplitude of the load [N]
q           = 0.0          # Uniform lateral load [N/m]
t_impulse   = 0.1        # Duration of the impulse load [s]

# Inverse analysis properties
bInverseAnalysis = true
InvSolver        = DirectXUA{0,0,0}   # Dynamic solver for the inverse analysis
Δtᵢₙᵥ             = 0.01               # Time step for the inverse analysis [s]

# General properties
t₁ = 0.0
t₂ = 1.0
time = t₁+Δt₀:Δt₀:t₂
nLoadSteps = length(time)
time_inv = t₁+Δt₀:Δtᵢₙᵥ:t₂
nLoadSteps_inv = length(time_inv)
maxiter     = 10
maxΔx       = 1e-6
maxΔu       = 1e-6
maxΔa       = 1e-6
maxΔλ       = Inf

# Direct model
name        = :BeamDynSinusoidalLoad
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
if bNodalForceSin
    @functor with(F, t₂) NodalSin(t) = F * sin(2*π*t/t₂)
    addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalSin )
end

if bNodalForceImpulse
    @functor with(F, t_impulse) NodalImpulse(t) = t <= t_impulse ? F * (1 - cos(2*π*t/t_impulse))/2 : 0.0
    addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalImpulse )
end

if bNodalStaticForce
    @functor with(F) NodalStatic(t) = F
    addelement!(model,DofLoad,[nodid[node_number]];field=:t2,value= NodalStatic )
end

# Initializing the model (sets all DoF to 0 at t=t₁)
initialstate                = initialize!(model; time=0.);

if bEigenAnalysis
    state_stat   = solve(SweepX{0};initialstate,time=[t₁]);
    # Solve eigenvalue problem
    nmod            = 10
    res             = solve(EigX{ℝ}; state=state_stat[1],nmod);
end

# Solving in direct mode
state                       = solve(DirectSolver;initialstate,time,verbose=true,maxΔx, maxiter);




# Post-processing of the direct analysis
bSaveFigures = false
bShow3DBeam = true

x_sin = [[getdof(state[idxLoad];field=:t1,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
y_sin = [[getdof(state[idxLoad];field=:t2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
z_sin = [[getdof(state[idxLoad];field=:t3,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
r1_sin = [[getdof(state[idxLoad];field=:r1,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
r2_sin = [[getdof(state[idxLoad];field=:r2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
r3_sin = [[getdof(state[idxLoad];field=:r3,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]

figure     = Figure(size = (1000,1000))
ax      = Axis3(figure[1,1],xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",aspect=:equal)
for to_draw in 1:10:nLoadSteps
    draw!(ax,state[to_draw];EulerBeam3D=(;nseg=20,  line_color= RGBf(1.0, to_draw/nLoadSteps, 0.)))
end
display(figure)
figure


if bInverseAnalysis
    name        = :BeamDynSinusoidalLoads
    inv_model       = Model(name)
    nodid       = addnode!(inv_model, nodeCoord)
    mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
    eleid       = addelement!(inv_model, EulerBeam3D, mesh;mat=mat, orient2=SVector(0.,1.,0.))

    # Boundary conditions
    [addelement!(inv_model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]];           # Support at one end
    [addelement!(inv_model,Hold,[nodid[nnodes]]  ;field) for field∈[:t1, :t2,:t3,:r1]];     # Support at the other end
    if bPlanar
        [[addelement!(inv_model,Hold,[nodid[i]] ;field) for field∈[:t3]] for i in 2:nnodes-1] # Planar motion constraint for eigenvalue analysis
    end

    x_int = [linear_interpolation(time_inv, vcat(x_sin[node]...)) for node in 1:nnodes]
    y_int = [linear_interpolation(time_inv, vcat(y_sin[node]...)) for node in 1:nnodes]
    z_int = [linear_interpolation(time_inv, vcat(z_sin[node]...)) for node in 1:nnodes]
    r3_int = [linear_interpolation(time_inv, vcat(r3_sin[node]...)) for node in 1:nnodes]


    @functor with() costX(x, t, meas) = 1 * (meas(t)-x)^2
    @functor with() costXother(x, t, meas) = 1 * (meas(t)-x)^2
    @functor with() costU(u, t) = 1*(F-u)^2
    @functor with() costUother(u, t) = 1*u^2
    e5             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t1,    cost= costXother, costargs= (meas = x_int[node],) ) for node in 1:nnodes]
    e6             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t2,    cost= costX, costargs= (meas = y_int[node],) ) for node in 1:nnodes]
    e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t3,    cost= costXother, costargs= (meas = z_int[node],) ) for node in 1:nnodes];
    e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:r3,    cost= costXother, costargs= (meas = r3_int[node],) ) for node in 1:nnodes];
    e2             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t3,Ufield=:t3           ,    cost=costUother )  for node in 1:nnodes-1];
    e3             = [addelement!(inv_model,SingleUdof,[nodid[node_number]]; Xfield=:t2,Ufield=:t2    ,    cost=costU )];
    e3             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t2,Ufield=:t2           ,    cost=costUother )  for node in 1:nnodes-1 if node != node_number];
    e4             = [addelement!(inv_model,SingleUdof,[nodid[node]]; Xfield=:t1,Ufield=:t1           ,    cost=costUother )  for node in 1:nnodes-1];


    [[addelement!(inv_model,Hold,[nodid[i]] ;field) for field∈[:t3, :r2, :r1]] for i in 1:nnodes]

    initialstate    = initialize!(inv_model;time=t₁)
    stateXUA         = solve(InvSolver;initialstate=[initialstate], time=[time_inv],verbose=true,maxiter,maxΔx,maxΔλ,maxΔu,maxΔa);


    x_inv = [[getdof(stateXUA[1][idxLoad];field=:t1,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
    y_inv = [[getdof(stateXUA[1][idxLoad];field=:t2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
    z_inv = [[getdof(stateXUA[1][idxLoad];field=:t3,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
    r1_inv = [[getdof(stateXUA[1][idxLoad];field=:r1,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
    r2_inv = [[getdof(stateXUA[1][idxLoad];field=:r2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
    r3_inv = [[getdof(stateXUA[1][idxLoad];field=:r3,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes] 


    U_t = [[getdof(stateXUA[1][idxLoad];class = :U, field=:t2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]

    fig_l = Figure(size = (1000,1000))
    axl = Axis(fig_l[1, 1])
    lines!(axl, 1:100, vcat(U_t[node_number+3]...); color = :tomato)
    lines!(axl, 1:100, vcat([F for t in 1:nLoadSteps_inv]...); color = :blue)
    display(fig_l)

    fig_disp = Figure(size = (1000,1000))
    ax = Axis(fig_disp[1, 1])
    [lines!(ax, time_inv, vcat(y_inv[25]...); color = :blue)]
    [lines!(ax, time_inv,[y_int[25](i) for i in time_inv]; color = :tomato)]
    display(fig_disp)

    figure     = Figure(size = (1000,1000))
    ax      = Axis3(figure[1,1],xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",aspect=:equal)
    for to_draw in 1:10:nLoadSteps
        draw!(ax,stateXUA[1][to_draw];EulerBeam3D=(;nseg=20,  line_color= RGBf(1.0, to_draw/nLoadSteps, 0.)))
        draw!(ax,state[to_draw];EulerBeam3D=(;nseg=20,  line_color= RGBf(0.0, 0.0, 1.0 )))
    end
    display(figure)
    figure
end

