using Muscade, StaticArrays, GLMakie
include("BeamElement.jl");


R = 0.0;  # Radius of the bend [m]
EI₂ = 833.33e3;  # Bending stiffness [Nm²]
EI₃ = 833.33e3;  # Bending stiffness [Nm²]
EA = 1e9;  # Axial stiffness [N]
GJ = 705e3;  # Torsional stiffness [Nm²]
L = 10;  # Length of the beam [m]

nel         = 9
nnodes      = nel+1
nodeCoord   = hcat( -5. .+ ((1:nnodes).-1)/(nnodes-1)*L,
                    0 .+ zeros(Float64,nnodes,1),
                    0 .+ zeros(Float64,nnodes,1))
mat         = BeamCrossSection(EA=EA,EI₂=EI₂,EI₃=EI₃,GJ=GJ,μ=1.,ι₁=1.)
model       = Model(:TestModel)
nodid       = addnode!(model,nodeCoord)
mesh        = hcat(nodid[1:nnodes-1],nodid[2:nnodes])
eleid       = addelement!(model,EulerBeam3D,mesh;mat=mat,orient2=SVector(0.,1.,0.))
[addelement!(model,Hold,[nodid[1]]  ;field) for field∈[:t1,:t2,:t3]]; # Support at one end
[addelement!(model,Hold,[nodid[10]]  ;field) for field∈[:t3]]; # Support at the other end

function load(t)
    a = -20
    t<=1. ? load=a*t*300. :
    t>1. && t<=2. ? load=a*(300. +(t-1)*150.) :
    load=a*(4500. +(t-2)*1500.)
end
addelement!(model,DofLoad,[nodid[5]];field=:t3,value=t->load(t));


initialstate    = initialize!(model);
loadSteps = [0.,1.,2.,3.];
nLoadSteps = length(loadSteps)
state           = solve(SweepX{0};initialstate,time=loadSteps,verbose=true,maxΔx=1e-9);

x_ = [getdof(state[idxLoad];field=:t1,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps]
y_ = [getdof(state[idxLoad];field=:t2,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps]
z_ = [getdof(state[idxLoad];field=:t3,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps]

fig     = Figure(size = (1000,1000))
ax      = Axis3(fig[1,1],xlabel="x [m]", ylabel="y [m]", zlabel="z [m]",aspect=:equal)
clr = [:black,:blue,:green,:red]
for idxLoad ∈ 1:nLoadSteps
    draw!(ax,state[idxLoad];EulerBeam3D=(;nseg=10,line_color=clr[idxLoad]))
end
xlims!(ax, -5,5); ylims!(ax, -5,5); zlims!(ax, -5,5);
save("first_plot.png", fig)