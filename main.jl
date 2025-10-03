using Muscade, StaticArrays, GLMakie, CSV, DataFrames, Interpolations
include("BeamElement.jl");
R   = 0.0;          # Radius of the bend [m]
EI₂ = 833.33e3;     # Bending stiffness [Nm²]
EI₃ = 833.33e3;     # Bending stiffness [Nm²]
EA  = 1e8;          # Axial stiffness [N]
GJ  = 705e3;        # Torsional stiffness [Nm²]
L   = 10.;           # Length of the beam [m]

nel         = 2
nnodes      = nel+1
nodeCoord   = hcat( -5. .+ ((1:nnodes).-1)/(nnodes-1)*L,
                     0  .+ zeros(Float64, nnodes, 1),
                     0  .+ zeros(Float64, nnodes, 1))
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=1.,ι₁=1.)

function createSimplySupportedBeam(name::Symbol; bPlanar=false)
    model       = Model(name)
    nodid       = addnode!(model, nodeCoord)
    mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes], nodid[1:nnodes-1])
    eleid       = addelement!(model, EulerBeam3D{true}, mesh;mat=mat, orient2=SVector(0.,1.,0.))
    [addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3,:r1]];   # Support at one end
    [addelement!(model,Hold,[nodid[nnodes]]  ;field) for field∈[:t1, :t2,:t3,:r1]];      # Support at the other end
    if bPlanar 
        [[addelement!(model,Hold,[nodid[i]] ;field) for field∈[:t3]] for i in 2:nnodes-1] 
    end
    return model, nodid, nnodes
end

function sinus_load(A, t, T, ϕ)
    sinus_load = A * sin(2*π*t/T + ϕ)
end

function distributed_sinus_load(q, t, T_max, L, x)
    distributed_sinus_load = q*(t/T_max)*sin(2*π*x/L)
end
model_dyn, nodid, nnodes   = createSimplySupportedBeam(:DynAnalysis_sin)


addelement!(model_dyn,DofLoad,[nodid[floor(Int,nnodes/2)]];field=:t2,value= t -> sinus_load(1e4, t, 10., 0.) )

initialstate                = initialize!(model_dyn; time = 0.);
Tsin                        = 1.:0.05:10.
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







inv_model, nodid, nnodes = createSimplySupportedBeam(:InverseModel_bis)

T = Tsin
x_int = [linear_interpolation(T, vcat(x_sin[node]...)) for node in 1:nnodes]
y_int = [linear_interpolation(T, vcat(y_sin[node]...)) for node in 1:nnodes]
z_int = [linear_interpolation(T, vcat(z_sin[node]...)) for node in 1:nnodes]
r3_int = [linear_interpolation(T, vcat(r3_sin[node]...)) for node in 1:nnodes]

e5             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t1,    cost= (x,t) -> 10 * (x_int[node](t)-x)^2 ) for node in 1:nnodes];
e6             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t2,    cost= (x,t) -> 10 * (y_int[node](t)-x)^2 ) for node in 1:nnodes];
e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:t3,    cost= (x,t) -> 10 * (z_int[node](t)-x)^2 ) for node in 1:nnodes];
e7             = [addelement!(inv_model,SingleDofCost,[nodid[node]];class=:X,field=:r3,    cost= (x,t) -> 10 * (r3_int[node](t)-x)^2 ) for node in 1:nnodes];
e2             = [addelement!(inv_model,SingleDofCost,[nodid[node]]; class=:U,field=:t3           ,    cost=(u,t) -> node == 1 ? 10*(sinus_load(1e4, t, 10., 0.)-u)^2 : 10*u^2 )  for node in 1:nnodes-1];
e3             = [addelement!(inv_model,SingleDofCost,[nodid[node]]; class=:U,field=:t2           ,    cost=(u,t) -> 10*u^2 )  for node in 1:nnodes-1];
e4             = [addelement!(inv_model,SingleDofCost,[nodid[node]]; class=:U,field=:t1           ,    cost=(u,t) -> 10*u^2 )  for node in 1:nnodes-1];


[[addelement!(inv_model,Hold,[nodid[i]] ;field) for field∈[:t3, :r2, :r1]] for i in 1:nnodes]

initialstateXUA    = initialize!(inv_model;time=0.)
stateXUA         = solve(DirectXUA{2,0,0};initialstate=initialstateXUA,time=T,
                        maxiter=100,saveiter=true,
                        maxΔx=1e-5,maxΔλ=Inf,maxΔu=1e-5,maxΔa=1e-5);

Muscade.describe(inv_model, nodid[9])