using Muscade, StaticArrays, GLMakie, CSV, DataFrames, Interpolations
include("BeamElement.jl");
R   = 0.0;          # Radius of the bend [m]
EI₂ = 833.33e3;     # Bending stiffness [Nm²]
EI₃ = 833.33e3;     # Bending stiffness [Nm²]
EA  = 1e8;          # Axial stiffness [N]
GJ  = 705e3;        # Torsional stiffness [Nm²]
L   = 10.;           # Length of the beam [m]

nel         = 50
nnodes      = nel+1
nodeCoord   = hcat( -5. .+ ((1:nnodes).-1)/(nnodes-1)*L,
                     0  .+ zeros(Float64, nnodes, 1),
                     0  .+ zeros(Float64, nnodes, 1))
mat         = DampedBeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=1.,ι₁=1., β=-0.005)

function createSimplySupportedBeam(name::Symbol; bPlanar=false)
    model       = Model(name)
    nodid       = addnode!(model, nodeCoord)
    mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
    eleid       = addelement!(model, EulerBeam3D, mesh;mat=mat, orient2=SVector(0.,1.,0.))
    Muscade.describe(model, :eletyp)
    [addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]];   # Support at one end
    [addelement!(model,Hold,[nodid[nnodes]]  ;field) for field∈[:t1, :t2,:t3,:r1]];      # Support at the other end
    if bPlanar 
        [[addelement!(model,Hold,[nodid[i]] ;field) for field∈[:t3]] for i in 2:nnodes-1] 
    end
    return model, nodid, nnodes, eleid
end

function sinus_load(A, t, T, ϕ)
    sinus_load = A * sin(2*π*t/T + ϕ)
end

function distributed_sinus_load(q, t, T_max, L, x)
    distributed_sinus_load = q*(t/T_max)*sin(2*π*x/L)
end
model_dyn, nodid, nnodes, eleid   = createSimplySupportedBeam(:DynAnalysis_sin)

@functor (;) load(t) = sinus_load(100., t, 20., 0.)
[addelement!(model_dyn,DofLoad,[nodid[node]];field=:t2,value= load ) for node in 1:nnodes]

initialstate                = initialize!(model_dyn; time = 0.);
Tsin                        = 1.:0.05:20.
nLoadSteps                  = length(Tsin)
state                       = solve(SweepX{2};initialstate,time=Tsin,verbose=true,maxΔx=1e-8, maxiter = 80);

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

req = @request R
loads = getresult(state[250],req,[eleid[10]])
inertiaLoads = loads[1]



inv_model, nodid, nnodes, _ = createSimplySupportedBeam(:InverseModel_bis)

T = Tsin
x_int = [linear_interpolation(T, vcat(x_sin[node]...)) for node in 1:nnodes]
y_int = [linear_interpolation(T, vcat(y_sin[node]...)) for node in 1:nnodes]
z_int = [linear_interpolation(T, vcat(z_sin[node]...)) for node in 1:nnodes]
r3_int = [linear_interpolation(T, vcat(r3_sin[node]...)) for node in 1:nnodes]


@functor (;) costX(x, t, meas) = 100000 * (meas(t)-x)^2
@functor (;) costXother(x, t, meas) = 100 * (meas(t)-x)^2
@functor (;) costU(u, t) = 0.0005*(sinus_load(100., t, 20., 0.)-u)^2
@functor (;) costUother(u, t) = 10*u^2
e5             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t1,    cost= costXother, costargs= (meas = x_int[node],) ) for node in 1:nnodes]
e6             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t2,    cost= costX, costargs= (meas = y_int[node],) ) for node in 1:nnodes]
e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t3,    cost= costXother, costargs= (meas = z_int[node],) ) for node in 1:nnodes];
e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:r3,    cost= costXother, costargs= (meas = r3_int[node],) ) for node in 1:nnodes];
e2             = [addelement!(inv_model,SingleDofCost,[nodid[node]]; class=:U,field=:t3           ,    cost=costUother )  for node in 1:nnodes-1];
e3             = [addelement!(inv_model,SingleDofCost,[nodid[node]]; class=:U,field=:t2           ,    cost=costU )  for node in 1:nnodes-1];
e4             = [addelement!(inv_model,SingleDofCost,[nodid[node]]; class=:U,field=:t1           ,    cost=costUother )  for node in 1:nnodes-1];


[[addelement!(inv_model,Hold,[nodid[i]] ;field) for field∈[:t3, :r2, :r1]] for i in 1:nnodes]

initialstateXUA    = initialize!(inv_model;time=0.)
stateXUA         = solve(DirectXUA{2,0,0};initialstate=[initialstateXUA],time=[T],
                        maxiter=20,
                        maxΔx=1e-5,maxΔλ=Inf,maxΔu=1e-5,maxΔa=1e-5);


x_inv = [[getdof(stateXUA[1][idxLoad];field=:t1,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
y_inv = [[getdof(stateXUA[1][idxLoad];field=:t2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
z_inv = [[getdof(stateXUA[1][idxLoad];field=:t3,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
r1_inv = [[getdof(stateXUA[1][idxLoad];field=:r1,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
r2_inv = [[getdof(stateXUA[1][idxLoad];field=:r2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]
r3_inv = [[getdof(stateXUA[1][idxLoad];field=:r3,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]

figure     = Figure(size = (1000,1000))
ax      = Axis3(figure[1,1],xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",aspect=:equal)
for to_draw in 1:10:nLoadSteps
    draw!(ax,stateXUA[1][to_draw];EulerBeam3D=(;nseg=20,  line_color= RGBf(1.0, to_draw/nLoadSteps, 0.)))
end
display(figure)
figure

U_t = [[getdof(stateXUA[1][idxLoad];class = :U, field=:t2,nodID=[nodid[node]]) for idxLoad ∈ 1:nLoadSteps] for node in 1:nnodes]

fig_disp = Figure(size = (1000,1000))
ax = Axis(fig_disp[1, 1])
[lines!(ax, 1:381, vcat(y_inv[i]...); color = :blue) for i in 1:50]
display(fig_disp)

fig_l = Figure(size = (1000,1000))
axl = Axis(fig_l[1, 1])
lines!(axl, 1:381, vcat(U_t[25]...); color = :tomato)
lines!(axl, 1:381, vcat([sinus_load(100., t, 20., 0.) for t in 1.:0.05:20.]...); color = :blue)
display(fig_l)