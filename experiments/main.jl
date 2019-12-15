include("count.jl")
using LightGraphs



# input graph
G = LightGraphs.SimpleGraphs.erdos_renyi(100, 0.1);

n = nv(G); # number of nodes
m = ne(G); # number of edges

# adding updates
n_updates = 10;
updates = Array{BPDTriangleCount.Update, 1}(undef, n_updates);

for i in 1:n_updates
    u = rand(1:n)
    v = rand(1:n)
    if (u>v) swap(u,v) end
    udt_i = BPDTriangleCount.Update(BPDTriangleCount.LGEdge(u,v),true)
    updates[i] = udt_i
end



# count updated triangles
state = BPDTriangleCount.init_BPDTCState(G)
last_update_ind, inserts, deletes = BPDTriangleCount.process_updates(updates)
C = BPDTriangleCount.process(G, n, m, n_updates, updates, C, state, inserts, deletes )

for i in last_update_ind
    if updates[i].add
        add_edge!(G, updates[i].e.src, updates[i].e.dst)
    else
        rem_edge!(G, updates[i].e.src, updates[i].e.dst)
    end
end

#12) major rebalance if needed
if (m + len(inserts)-len(removes)) > M ||  (m + len(inserts)-len(removes)) < M/4
    state = BPDTriangleCount.init_BPDTCState(G)
end
