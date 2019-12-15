module BPDTriangleCount
using LightGraphs
using StaticArrays
include("util.jl")
include("bench.jl")
# PERSISTENT dictionary
const NO_KEY = -3 #TYPE?

# const LGEdge{T} = LightGraphs.SimpleGraphs.SimpleEdge{T}???

const LGEdge = LightGraphs.SimpleGraphs.SimpleEdge

# tools for users to define updates

# [] parametrize
# [] SVector

"""
data structure representing the updates
`e` edge to udpate
`add` true if adding, false if deleting
"""
struct Update
    e::LGEdge
    add::Bool
end



"""
designed to keep static global variables for the algorithm

in T, we only have T[u][v] u < v
"""
# preallocate?
# parametrize integer
# funvtions better to have non strict typed
struct BPDTCState{T, R}
    n::T#Integer
    m::T#Integer
    M::T#Integer
    t1::R #Real
    t2::R #Real
    T:: Dict{T, Dict{T, SizedArray{Tuple{5}, T}}} #static array
    HH::Dict{T, Dict{T, T}}
    HL::Dict{T, Dict{T, T}}
    LH::Dict{T, Dict{T, T}}
    LL::Dict{T, Dict{T, T}}
    D::Array{T, 1}
end


"""
helper functions
"""

# @inline function get_from(
#     table::Dict{Integer, Dict{Integer, Integer}},
#     u::Integer,
#     v::Integer
#     )
#     if (!haskey(table, u)) return NO_KEY end
#
#     return  get(table[u], v, NO_KEY)
# end

@inline function get_from(
    table::Dict{T, Dict{T, T}},
    u::T,
    v::T
    ) where {T}
    if (!haskey(table, u)) return T(NO_KEY) end

    return  get(table[u], v, NO_KEY)
end


@inline function add_to(
    table,
    u,
    v,
    val
    ) # shoudl I add type?
    table[u] = get(table, u, Dict{Integer, Integer}())
    table[u][v] = val
end

"""
return the value deleted and delete the pair from the table

return NO_KEY if not in table
"""
# @inline function delete_from(
#     table::Dict{Integer, Dict{Integer, Integer}},
#     u::Integer,
#     v::Integer
#     )
#     val = get_from(table, u, v)
#     if (val == NO_KEY) return NO_KEY end
#     delete!(table[u], v)
#     if (isempty(table[u])) delete!(table, u) end
#     return val
# end

@inline function delete_from(
    table,
    u::T,
    v
    ) where {T}
    val = get_from(table, u, v)
    if (val == NO_KEY) return T(NO_KEY) end
    delete!(table[u], v)
    if (isempty(table[u])) delete!(table, u) end
    return val
end

@inline in_high(u::Integer, state) = state.D[u] > state.t2 #::BPDTCState???
@inline in_low(u::Integer, state) = !in_high(u, state)

@inline function sameEdge(e1, e2)
    return e1 == e2
    # return e1.u==e2.u && e1.v==e2.v
end


# steps in algorithms

struct CleanFilter{A, G, T}
    updates::A
    G::G
    n::T
end

function (p::CleanFilter)(i)
    last_update = (i == p.n || !sameEdge(p.updates[i].e, p.updates[i+1].e))
    if (last_update)
        already_in = p.updates[i].add && !has_edge(p.G, p.updates[i].e.src, p.updates[i].e.dst)
        not_in = !p.updates[i].add && has_edge(p.G, p.updates[i].e.src, p.updates[i].e.dst)
        return already_in || not_in
    else
        return false
    end
end


function process_updates(G, updates, n_updates)
    #only leave last updates, need a stable sort!
    t1 = @elapsed sort!(updates, by=update->(update.e.src, update.e.dst), alg = MergeSort)#faster sort?
    # last_update_ind = filter(i->(i == n_updates || !sameEdge(updates[i].e, updates[i+1].e)), 1:n_updates)

    #remove edges already in graph/ not in graph
    f = CleanFilter(updates, G,  n_updates)
    # filter!(f, last_update_ind)
    t2 = @elapsed last_update_ind = filter(f, 1:n_updates)

    #inserts
    t3 = @elapsed inserts = filter(i->(updates[i].add), last_update_ind) # better memory usage?
    #removes
    t4 = @elapsed removes = filter(i->(!updates[i].add), last_update_ind)
    println(t1)
    println(t2)
    println(t3)
    println(t4)
    return last_update_ind, inserts, removes
end



function process_updates_par(G, updates, n_updates)
    #only leave last updates, need a stable sort!
    t = similar(updates, n_updates)
    t1 = @elapsed UTIL.SortParallel!(updates, t, update->(update.e.src, update.e.dst))

    #remove edges already in graph/ not in graph
    f = CleanFilter(updates, G,  n_updates)
    t2 = @elapsed last_update_ind = UTIL.FilterParallel(f, 1:n_updates)

    #inserts
    # t3 = @elapsed inserts, removes = UTIL.PartitionParallel(i->(updates[i].add), last_update_ind) # better memory usage?
    # t3 = @elapsed inserts, removes = BENCH.PartitionParallel(i->(updates[i].add), last_update_ind) # better memory usage?
    t3 = @elapsed inserts = filter(i->(updates[i].add), last_update_ind) # better memory usage?
    #removes
    t4 = @elapsed removes = filter(i->(!updates[i].add), last_update_ind)


    println(t1)
    println(t2)
    println(t3)
    println(t4)

    return last_update_ind, inserts, removes
end




"""
initialize the data structures and constants
including HH, HL, LH, LL, D, T

"""
function init_BPDTCState(G::SimpleGraph{TY}) where {TY}
    n = nv(G);
    m = ne(G);
    M = 2 * m + 1;
    t1 = sqrt(M)/2;
    t2 = 3 * t1;

    D = zeros(TY,floor(TY, 1.5 * n));
    # D = Array{TY}(undef, floor(Integer, 1.5 * n))
    # fill!(D, 0);
    HH = Dict{TY, Dict{TY, TY}}();
    HL = Dict{TY, Dict{TY, TY}}();
    LH = Dict{TY, Dict{TY, TY}}();
    LL = Dict{TY, Dict{TY, TY}}();
    T = Dict{TY, Dict{TY, SizedArray{Tuple{5}, TY}}}();

    state = BPDTCState(n,m,M,t1,t2,T, HH, HL, LH, LL, D)

    # init degrees and tables
    for e in edges(G)
        u,v = e.dst, e.src;
        if (v<u) u,v = v,u end
        state.D[u] += 1;
        state.D[v] += 1;


        # thread safe dictionary?
        if (in_high(u, state) && in_high(v, state))
            add_to(state.HH, u, v, 0)
            add_to(state.HH, v, u, 0)
        elseif (in_high(u, state))
            add_to(state.HL, u, v, 0)
            add_to(state.LH, v, u, 0)
        elseif (in_low(u, state) && in_low(v, state))
            add_to(state.LL, u, v, 0)
            add_to(state.LL, v, u, 0)
        else
            add_to(state.LH, u, v, 0)
            add_to(state.HL, v, u, 0)
        end
    end

    # init T
    for u in keys(state.HL)
        for w in keys(state.HL[u])
            for v in keys(state.LH[w])
                if (u < v)
                    increment_T(state.T, u, v, 1)
                else
                    increment_T(state.T, v, u, 1)
                end # same wedge or different: should be different?
            end
        end
    end


    return state
end

"""
#2,5) mark inserted/deleted edges
# update HH,HL,LH,LL
# treat in between vertices as low nodes
"""
function mark_edges(updates, inds, state, mark)
    for i in inds
        u,v = updates[i].e.src, updates[i].e.dst
        if (in_high(u, state) && in_high(v, state)) #HH
            add_to(state.HH, u, v, mark)
            add_to(state.HH, v, u, mark)
        elseif (in_high(u, state)) #HL, LH
            add_to(state.HL, u, v, mark)
            add_to(state.LH, v, u, mark)
        elseif (in_low(u, state) && in_high(v, state)) #LH, HL
            add_to(state.LH, u, v, mark)
            add_to(state.HL, v, u, mark)
        else #LL
            add_to(state.LL, u, v, mark)
            add_to(state.LL, v, u, mark)
        end
    end
end

"""
4) minor rebalance due to inserts
Perform all minor rebalancing due to insertions for all vertices v that exceed t2
in degree. This type of minor rebalance updates the
structure by making a low-degree vertex high-degree
 Move all edges in LH and LL to HH and HL.

 # also move from HL to HH ??
 # better method?
"""
function minor_rebalancing_inserts(state)
    LH, HL, HH, LL, D= state.LH, state.HL, state.HH, state.LL, state.D
    for u in keys(LH)
        if in_high(u, state)
            tmp = LH[u]
            delete!(LH, u)# move to HH
            HH[u] = tmp
            for v in keys(HH[u])
                val2 = delete_from(HL,v,u)
                add_to(HH,v,u,val2)
            end
        end
    end

    for u in keys(LL)
        if in_high(u, state)
            tmp = LL[u]
            delete!(LL, u)# move to HL
            HL[u] = tmp
            for v in keys(HL[u])
                val1 = delete_from(LL,v,u)
                add_to(LH,v,u,val1)
            end
        end
    end
end


"""
11) minor rebalance due to deletions
"""
# TODO: no_key val?
function minor_rebalancing_deletes(state)
    LH, HL, HH, LL, D= state.LH, state.HL, state.HH, state.LL, state.D
    for u in keys(HL)
        if in_low(u, state)
            tmp = HL[u]
            delete!(HL, u)# move to LL
            LL[u] = tmp
            for v in keys(LL[u])
                val2 = delete_from(LH,v,u)
                add_to(LL,v,u,val2)
            end
        end
    end

    for u in keys(HH)
        if in_low(u, state)
            tmp = HH[u]
            delete!(HH, u)# move to LH
            LH[u] = tmp
            for v in keys(LH[u])
                val1 = delete_from(HH,v,u)
                add_to(HL,v,u,val1)
            end
        end
    end
end

"""
#6) update table insertions
"""
function increment_T(T, u, v, i)
    print("incrementing")
    if (!haskey(T, u))
        temp = SizedArray{Tuple{5}, Integer}(zeros(5));
        temp[i] = 1;
        T[u] =  Dict{Integer, SizedArray{Tuple{5}, Integer}}(v=>temp)
        return;
    end
    if (!haskey(T[u],v))
        temp = SizedArray{Tuple{5}, Integer}(zeros(5));
        temp[i] = 1;
        T[u][v] = temp
        return;
    end
    T[u][v][i] += 1
end


function update_T_insertions(state, updates, inserts)
    LH, HL, HH, LL, T, D = state.LH, state.HL, state.HH, state.LL, state.T, state.D
    for i in inserts
        u,v = updates[i].e.src, updates[i].e.dst
        if (D[u] < D[v]) u,v = v,u end

        if (in_high(u, state) && in_low(v, state)) #HL
            for (w, val) in pairs(IndexStyle(LH[v]), LH[v])
                if (val == 0)
                    if (u < w)
                        increment_T(T, u, w, 2)
                    elseif (u > w)
                        increment_T(T, w, u, 2)
                    end
                elseif (val == 1)
                    if (u < w)
                        increment_T(T, u, w, 3)
                    elseif (u > w)
                        increment_T(T, w, u, 3)
                    end
                end
            end
        end

    end
end

"""
#7) update table deletions
"""
function update_T_deletions(state, updates, deletes)
    LH, HL, HH, LL, T, D = state.LH, state.HL, state.HH, state.LL, state.T, state.D
    for i in deletes
        u,v = updates[i].e.src, updates[i].e.dst
        if (D[u] < D[v]) u,v = v,u end
        if in_high(u, state) && in_low(v, state) #HL
            for (w, val) in pairs(IndexStyle(LH[v]), LH[v])
                if val == 0
                    if (u < w)
                        increment_T(T, u, w, 4)
                    elseif (u > w)
                        increment_T(T, w, u, 4)
                    end
                elseif val == 2
                    if (u < w)
                        increment_T(T, u, w, 5)
                    elseif (u > w)
                        increment_T(T, w, u, 5)
                    end
                end
            end
        end
    end
end


"""
### step 8 helpers
"""

@inline function count_triangles_helper_add(t_uw, t_wv, u, w, v, c)
    if (!haskey(t_wv[w], v)) return end # no edge wv
    if (t_uw[u][w] == 2 || t_wv[w][v] == 2)
        return
    end
    ind = t_uw[u][w] + t_wv[w][v] + 1
    c[ind] += 1
end

@inline function count_triangles_helper_remove(t_uw, t_wv, u, w, v, c)
    if (!haskey(t_wv[w], v)) return end # no edge wv
    if (t_uw[u][w] == 1 || t_wv[w][v] == 1)
        return
    end
    ind = (t_uw[u][w] + t_wv[w][v])/2 + 1
    c[Int(3+ind)] += 1
end

# c = zeros(6)
function count_triangles_add(e, c, state)
    LH, HL, HH, LL, T, D= state.LH, state.HL, state.HH, state.LL, state.T, state.D
    u,v = e.src, e.dst
    if (D[u] > D[v]) u,v = v,u end
    if (in_low(u, state)) # at least one in low
        if (haskey(LH,u))
            for w in keys(LH[u])
                if (in_high(v, state))
                    count_triangles_helper_add(LH, HH, u, w, v, c)
                else
                    count_triangles_helper_add(LH, HL, u, w, v, c)
                end
            end
        end
        if (haskey(LL,u))
            for w in keys(LL[u])
                if (in_high(v, state))
                    count_triangles_helper_add(LL, LH, u, w, v, c)
                else
                    count_triangles_helper_add(LL, LL, u, w, v, c)
                end
            end
        end
    else # both in high
        if (u > v) u,v = v,u  end
        for w in keys(HH[u])
            count_triangles_helper_add(HH, HH, u, w, v, c)
        end
        if (haskey(T, u) && haskey(T[u],v))
            c[1:3] .+= T[u][v][1:3]
        end
          # uwv ans vwu should be 2 different wedges TODO: check
    end
end

function count_triangles_remove(e, c, state)
    LH, HL, HH, LL, T, D = state.LH, state.HL, state.HH, state.LL, state.T, state.D#just pointers?
    u,v = e.src, e.dst
    if (D[u] > D[v]) u,v = v,u end
    if (in_low(u, state)) # at least one in low
        if (haskey(LH,u))
            for w in keys(LH[u])
                if (in_high(v, state))
                    count_triangles_helper_remove(LH, HH, u, w, v, c)
                else
                    count_triangles_helper_remove(LH, HL, u, w, v, c)
                end
            end
        end
        if (haskey(LL,u))
            for w in keys(LL[u])
                if (in_high(v, state))
                    count_triangles_helper_remove(LL, LH, u, w, v, c)
                else
                    count_triangles_helper_remove(LL, LL, u, w, v, c)
                end
            end
        end
    else # both in high
        if (u > v) u,v = v,u  end
        for w in keys(HH[u])
            count_triangles_helper_remove(HH, HH, u, w, v, c)
        end
        if (haskey(T, u) && haskey(T[u],v) ) c[4:6] .+= T[u][v][4:6] end
    end
end


"""
9) reset tables
"""
# u is in high degree, w is in low degree
function update_T_helper(u::Integer, w::Integer, T::Dict{Integer, Dict{Integer, SizedArray{Tuple{5}, Integer}}}, LH)
    for v in keys(LH[w])
        T[u][v][1] = sum(T[u][v][1:3]) - sum(T[u][v][4:5])
        T[u][v][2:5]  .= 0
    end
end

function reset_tables(state, updates, inserts, deletes)
    LH, HL, HH, LL, T = state.LH, state.HL, state.HH, state.LL, state.T
    for i in inserts
        u,v = updates[i].e.src, updates[i].e.dst
        if (in_low(u, state) && in_low(v, state))
            LL[u][v] = 0
            LL[v][u] = 0
        elseif (in_low(u, state) && in_high(v, state))
            LH[u][v] = 0
            HL[v][u] = 0
            update_T_helper(v, u, LH)
        elseif (in_high(u, state) && in_low(v, state))
            LH[v][u] = 0
            HL[u][v] = 0
            update_T_helper(u, v, LH)
        else
            HH[u][v] = 0
            HH[v][u] = 0
        end
    end
    for i in deletes
        u,v = updates[i].e.src, updates[i].e.dst
        if (in_low(u, state) && in_low(v, state))
            delete_from(LL,u,v)
            delete_from(LL,v,u)
        elseif (in_low(u, state) && in_high(v, state))
            delete_from(LH,u,v)
            delete_from(HL,v,u)
            update_T_helper(v, u, LH)
        elseif (in_high(u, state) && in_low(v, state))
            delete_from(LH,v,u)
            delete_from(HL,u,v)
            update_T_helper(u, v, LH)
        else
            delete_from(HH,u,v)
            delete_from(HH,v,u)
        end
    end
end


### main algorithms

function process(G, n, m, n_updates, updates, C, state, inserts, deletes)

    ##### inserts
    #2) mark inserted edges
    mark_edges(updates, inserts, state, 1)
    #3) update degree
    for i in inserts
        u,v = updates[i].e.src, updates[i].e.dst
        state.D[u] += 1
        state.D[v] += 1
    end
    #4) minor rebalancing on tables
    minor_rebalancing_inserts(state)

    ##### deletes
    #5) mark deleted edges
    mark_edges(updates, deletes, state, 2)

    #6) update table insertions
    update_T_insertions(state, updates, inserts)

    #7) update table deletions
    update_T_deletions(state, updates, deletes)

    #8) count triangles
    c = zeros(6)
    for i in inserts
        count_triangles_add(updates[i].e, c, state)
    end
    for i in deletes
        count_triangles_remove(updates[i].e, c, state)
    end

    C = C + c[1] + c[2]/2 + c[3]/3 - c[4] - c[5]/2 - c[6]/3

    #9) restore table
    reset_tables(state, updates, inserts, deletes)

    #10) update D due to deletions
    for i in deletes
        u,v = updates[i].e.src, updates[i].e.dst
        state.D[u] -= 1
        state.D[v] -= 1
    end

    #11） minor rebalance due to deletions
    minor_rebalancing_deletes(state)

    # update degree why updating degree of deletions last?

    print("C = ", C)
    return C
end

end
