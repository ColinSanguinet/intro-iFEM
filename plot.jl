using Muscade, Serialization

loaded_data = open("output/data.bin", "r") do io
    deserialize(io)
end

stateXUA = loaded_data
nnodes = 51
x_inv =  [getdof(stateXUA[1][1:end];field=:t1,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    y_inv =  [getdof(stateXUA[1][1:end];field=:t2,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    z_inv =  [getdof(stateXUA[1][1:end];field=:t3,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    r1_inv = [getdof(stateXUA[1][1:end];field=:r1,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    r2_inv = [getdof(stateXUA[1][1:end];field=:r2,nodID=[nodid[node]]) for node ∈ 1:nnodes]
    r3_inv = [getdof(stateXUA[1][1:end];field=:r3,nodID=[nodid[node]]) for node ∈ 1:nnodes]

    U_t_inv = [getdof(stateXUA[1][idxLoad];class = :U, field=:t2,nodID=nodid[1:nnodes]) for idxLoad ∈ 1:nLoadSteps_inv]