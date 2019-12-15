using UnsafeArrays

P = Threads.nthreads()

@inline function nblocks(_n,_bsize)
    return ceil(Int64,(_n/_bsize))
end

@inline function PrefixSumSerial(arr, out, s, e,  init = 0)
    @inbounds out[s] = arr[s] + init
    @simd ivdep for i = (s+1):e
        @inbounds out[i] = out[i-1] + arr[i]
    end
    return ;
end

@inline  function PrefixSumSerial(arr, out, init = 0)
    PrefixSumSerial(arr, out,1, length(arr),  init)
end

@inline  function PrefixSumSerial!(arr, init = 0)
    PrefixSumSerial( arr, arr, 1, length(arr), init)
end


function PrefixSumParallel_unsafe!(arr, Sums, _SCAN_SIZE::Integer = 1024)
    n = length(arr)
    l = nblocks(n, _SCAN_SIZE)

    if (l <= 2)
        PrefixSumSerial!(arr)
        return
    end

    Threads.@threads for i = 1:l
        # @show @which Base.mapfoldl_impl(identity,Base.add_sum, Base.Generator(identity,uview(arr, (i-1)*_SCAN_SIZE+1:min(i * _SCAN_SIZE, n))))
        let i = i, arr = arr, _SCAN_SIZE = _SCAN_SIZE, n = n, Sums = Sums
            @inbounds Sums[i+1] = sum(uview(arr, (i-1)*_SCAN_SIZE+1:min(i * _SCAN_SIZE, n)))
        end
    end
    PrefixSumParallel_unsafe!(uview(Sums, 1:(l+1)), uview(Sums, (l+2):length(Sums)), _SCAN_SIZE)
    Threads.@threads for i = 1:l #
        let i = i, arr = arr, _SCAN_SIZE = _SCAN_SIZE, n = n, Sums = Sums
            @inbounds PrefixSumSerial!(uview(arr, ((i-1)*_SCAN_SIZE+1):(min(i * _SCAN_SIZE, n))), Sums[i])
        end
    end

end

# check if that closes the gap with c
# the boxing bug -  no influence
# use indexing instead - does not help



# used Coverage.jl, the only extra allocations are in the threads


L = 6
N = 10^L

B = 1024

println(L)

arr = collect(1:300)
GC.gc()
@time PrefixSumSerial!(arr) # always  (4 allocations: 160 bytes)

arr = collect(1:N)
Sums = zeros(N)
GC.gc()
@time PrefixSumParallel_unsafe!(arr, Sums, B)

println("=====")
println("=====",P)

#
arr = collect(1:N)
GC.gc()
@time PrefixSumSerial!(arr) # always  (4 allocations: 160 bytes)

arr = collect(1:N)
Sums = zeros(N)
GC.gc()
@time PrefixSumParallel_unsafe!(arr, Sums, B)

arr = collect(1:N)
Sums = zeros(N)
GC.gc()
@time PrefixSumParallel_unsafe!(arr, Sums, B)

arr = collect(1:N)
Sums = zeros(N)
GC.gc()
@time PrefixSumParallel_unsafe!(arr, Sums, B)

arr = collect(1:N)
Sums = zeros(N)
GC.gc()
@time PrefixSumParallel_unsafe!(arr, Sums, B)
#
# # should run it on a machine that do not have other works
# repl is not stable
# but sometimes run the script is also not stable
