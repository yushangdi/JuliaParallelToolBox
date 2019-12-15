include("util.jl")
include("bench.jl")
using TimerOutputs

P = Threads.nthreads() # number of processors

function partitionserial(N, f, g, itr = 20)
    println("===Partition Serial, N = ", N, " itr = ", itr)
    # compile
    arr = g(N)
    UTIL.PartitionSerial(f, arr)
    UTIL.PartitionSerial!(f, arr)
    BENCH.PartitionSerial(f, arr)

    t1 = 100000
    for _ = 1:itr
        arr = g(N)
        GC.gc()
        t1 = min(t1, @elapsed UTIL.PartitionSerial(f, arr))
    end

    t3 = 100000
    for _ = 1:itr
        arr = g(N)
        GC.gc()
        t3 = min(t3, @elapsed UTIL.PartitionSerial!(f, arr))
    end

    t2 = 100000
    for _ = 1:itr
        arr = g(N)
        GC.gc()
        t2 = min(t2, @elapsed BENCH.PartitionSerial(f, arr))
    end

    println("our time: ", t1)
    println("our time inplace: ", t3)
    println("bench time: ", t2)
end


function partitionparallel(N, f, g, itr = 20)
    println("===Partition Parallel, N = ", N, " itr = ", itr)
    # compile
    arr = g(N)
    flags = zeros(Int64, 2*N)
    UTIL.PartitionParallel(f, arr)
    BENCH.PartitionParallel(f, arr)

    t1 = 100000
    for _ = 1:itr
        arr = g(N)
        flags = zeros(Int64, 2*N)
        GC.gc()
        t1 = min(t1, @elapsed UTIL.PartitionParallel(f, arr, flags))
    end

    t2 = 100000
    for _ = 1:itr
        arr = g(N)
        flags = zeros(Int64, 2*N)
        GC.gc()
        t2 = min(t2, @elapsed BENCH.PartitionParallel(f, arr))
    end

    println("our time: ", t1)
    println("bench time: ", t2)
end

println("============ EXPERIMENT WITH ", P, " WORKERS ===================")
println("==Partition==")
function f(x)
     ct = 0
     for _ in 1:10
       ct += rand()
     end
     ct > 5
end
g(N) = collect(1:N)

for i = 6:7
    for j = 1:4:9
        partitionserial(j * 10^i, f, g, 3)
        partitionparallel(j * 10^i, f, g, 3)
    end
end
