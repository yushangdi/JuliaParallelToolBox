include("util.jl")
include("bench.jl")
using TimerOutputs


function prefixserial(N, itr = 20)
    println("===PrefixSum Serial, N = ", N, " itr = ", itr)
    # compile
    arr = rand(10)
    UTIL.PrefixSumSerial!(arr)
    BENCH.PrefixSumSerial!(arr)

    t1 = 10000
    for _ = 1:itr
        arr = rand(N)
        GC.gc()
        t1 = min(t1, @elapsed UTIL.PrefixSumSerial!(arr))
    end

    # t2 = 0
    # for _ = 1:itr
    #     arr = rand(N)
    #     GC.gc()
    #     t2 += @elapsed BENCH.PrefixSumSerial!(arr)
    # end
    # t2 /= itr
    println("our time: ", t1)
    # println("bench time: ", t2)
end

function prefixparallel(N, itr = 20)
    println("===PrefixSum Parallel, N = ", N, " itr = ", itr)
    # compile
    arr = rand(N)
    Sums = zeros(N)
    UTIL.PrefixSumParallel_unsafe!(arr, Sums)
    BENCH.PrefixSumParallel!(arr)
    # test start
    t1 = 10000
    for _ = 1:itr
        arr = rand(N)
        Sums = zeros(N)
        GC.gc()
        t1 = min(t1, @elapsed UTIL.PrefixSumParallel_unsafe!(arr, Sums))
    end
    t1 /= itr

    t2 = 0
    for _ = 1:itr
        arr = rand(N)
        GC.gc()
        t2 += @elapsed BENCH.PrefixSumParallel!(arr)
    end
    t2 /= itr

    println("our time: ", t1)
    println("bench time: ", t2)
    return t1, t2
end


for i = 6:7
    for j = 1:2:9
        prefixserial(j * 10^i, 3)
        prefixparallel(j * 10^i, 3)
    end
end
