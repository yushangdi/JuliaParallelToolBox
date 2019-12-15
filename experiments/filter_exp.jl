include("util.jl")
include("bench.jl")
using TimerOutputs

P = Threads.nthreads() # number of processors

function filterserial(N, f, g, itr = 20)
    println("===Filter Serial, N = ", N, " itr = ", itr)
    # compile
    arr = g(N)
    out = zeros(N)
    UTIL.FilterSerial(f, arr, out)
    BENCH.FilterSerial(f, arr)
    filter(f, arr)
    # test start
    t1 = 10000
    for _ = 1:itr
        arr = g(N)
        GC.gc()
        out = zeros(N)
        t1= min(t1, @elapsed UTIL.FilterSerial(f, arr, out))
    end

    #if (N >= 10^7)
    #    println("our time: ", t1)
    #    return t1, -1
    #end

    t2 = 10000
    for _ = 1:itr
        arr = g(N)
        GC.gc()
        t2 = min(t2, @elapsed BENCH.FilterSerial(f, arr))
    end

    t3 = 10000
    for _ = 1:itr
        arr = g(N)
        GC.gc()
        t3 = min(t3, @elapsed filter(f, arr))
    end

    println("our time: ", t1)
    println("bench time 1: ", t2)
    println("bench time 2: ", t3)
    return t1, t2, t3
end


function filterparallel(N, f, g, itr = 20)
    println("===Filter Parallel, N = ", N, " itr = ", itr)
    # compile
    arr = g(N)
    flags = zeros(Int64, 2*N)
    UTIL.FilterParallel(f, arr, flags)
    BENCH.FilterParallelD(f, arr)
    # test start


    t1 = 10000
    for _ = 1:itr
        arr = g(N)
        flags = zeros(Int64, 2*N)
        GC.gc()
        t1 = min(t1, @elapsed UTIL.FilterParallel(f, arr, flags))
    end

    if (N >= 10^7)
        println("our time: ", t1)
        return t1, -1
    end

    t2 = 10000
    for _ = 1:itr
        arr = g(N)
        GC.gc()
        t2 = min(t2, @elapsed BENCH.FilterParallelD(f, arr))
    end
    t2 /= itr

    println("our time: ", t1)
    println("bench time 1: ", t2)
    return t1, t2
end

println("============ EXPERIMENT WITH ", P, " WORKERS ===================")
println("==Filter==")

function f(x)
     ct = 0
     for _ in 1:10
       ct += rand()
     end
     ct > 5
end

g(N) = 1:N 


for i = 6:7
    for j = 1:4:9
        filterserial(j * 10^i, f, g, 3)
        filterparallel(j * 10^i, f, g, 3)
    end
end
