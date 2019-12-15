__precompile__()
module UTIL

using UnsafeArrays

Threads.nthreads()

@inline function nblocks(_n,_bsize)
    return ceil(Int64,(_n/_bsize))
end

@inline function get_par_loop_range(i, _SCAN_SIZE, n)
    return (i-1)*_SCAN_SIZE+1:min(i * _SCAN_SIZE, n)
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


function PrefixSumParallel!(arr::A, Sums::A, _SCAN_SIZE = 1024) where {A}
    n = length(arr)
    l = nblocks(n, _SCAN_SIZE)

    if (l <= 2)
        PrefixSumSerial!(arr)
        return
    end

    Threads.@threads  for i = 1:l
        @inbounds Sums[i+1] = sum(uview(arr, (i-1)*_SCAN_SIZE+1:min(i * _SCAN_SIZE, n)))
    end
    PrefixSumParallel!(@view(Sums[1:l+1]), @view(Sums[l+2:end]))
    Threads.@threads for i = 1:l
        @inbounds PrefixSumSerial!(@view(arr[(i-1)*_SCAN_SIZE+1:min(i * _SCAN_SIZE, n)]), Sums[i])
    end

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

function PrefixSumParallel!(arr, _SCAN_SIZE = 1024)
    Sums = zeros(length(arr))
    PrefixSumParallel_unsafe!(arr, Sums, _SCAN_SIZE)
end


function FilterSerial(f, a, out)
    j = 1
    for ai in a
        @inbounds out[j] = ai
        j = ifelse(f(ai), j+1, j)
    end
    resize!(out, j-1)
    sizehint!(out, length(out))
    length(out)
end

function FilterSerial!(f, arr)
    return FilterSerial(f, arr, arr)
end

function FilterParallel(f, arr, out::AbstractVector, flags::AbstractVector, _SCAN_SIZE = 1024)
    n = length(arr)
    l = nblocks(n, _SCAN_SIZE)

    if (l <= 2)
        return filter(f, arr, out)
    end

    #prefixsum to get index
    @inbounds Threads.@threads for i = 1:l
        r = get_par_loop_range(i, _SCAN_SIZE, n)
        @inbounds map!(f, uview(flags, r), uview(arr, r))
    end
    PrefixSumParallel_unsafe!(uview(flags,1:n), uview(flags,(n+1):(2*n)), _SCAN_SIZE)
    # PrefixSumSerial!(flags)
    if (flags[1] != 0)
        out[flags[1]] = arr[1]
    end
    @inbounds Threads.@threads for i = 1:l
        @simd ivdep for j = get_par_loop_range(i, _SCAN_SIZE, n)
            if (j<n && flags[j] != flags[j+1])
                out[flags[j+1]] = arr[j+1]
            end
        end
    end
    return flags[n]
end

"""
Return a vector of elements in arr that satisfies the predicate `f`.
`flags` is the auxiliary memory that the function will use. It should be a 0 vector of length 2* length(arr).
"""
function FilterParallel(f, arr, flags::AbstractVector, _SCAN_SIZE = 1024)
    n = length(arr)
    l = nblocks(n, _SCAN_SIZE)

    if (l <= 2)
        return filter(f, arr)
    end

    #prefixsum to get index
    Threads.@threads for i = 1:l
        r = get_par_loop_range(i, _SCAN_SIZE, n)
        @inbounds map!(f, uview(flags, r), uview(arr, r))
    end
    PrefixSumParallel_unsafe!(uview(flags,1:n), uview(flags,(n+1):(2*n)), _SCAN_SIZE)
    # PrefixSumSerial!(flags)
    out = zeros(typeof(arr[1]), flags[n])
    if (flags[1] != 0)
        out[flags[1]] = arr[1]
    end
    Threads.@threads for i = 1:l
        @simd ivdep for j = get_par_loop_range(i, _SCAN_SIZE, n)
            if (j<n && flags[j] != flags[j+1])
                @inbounds out[flags[j+1]] = arr[j+1]
            end
        end
    end
    return out
end


function FilterParallel(f, arr)
    flags = zeros(Int64, 2*length(arr))
    # flags = zeros(Int64, length(arr))
    # n_out = FilterParallel(f, arr, flags)
    # resize!(out, n_out)
    # sizehint!(out, n_out)
    return FilterParallel(f, arr, flags)
end


function PartitionSerial(f, a::A) where {A}
    j = 1
    k = 1
    b = A(undef, length(a))
    c = A(undef, length(a))
    @simd ivdep for ai in a
        if (f(ai))
            @inbounds b[j] = ai
            j += 1
        else
            @inbounds c[k] = ai
            k += 1
        end
    end
    resize!(b, j-1)
    sizehint!(b, length(b))
    resize!(c, k-1)
    sizehint!(c, length(c))
    b, c
end

"""
in place partition, but not stable
"""
function PartitionSerial!(f, a)
  i = 1
  n = length(a)
  j = n

  while ( true )

    while ( j>=1 && !f(a[j]) )
        j -=1
    end

    while ( i<= n && f(a[i]))
        i += 1
    end

    if ( i < j )
      a[[j i]] = a[[i j]]
    else
      return j
    end

  end
end

"""
Return two vectors of elements in arr that satisfies/fails the predicate `f`.
`flags` is the auxiliary memory that the function will use. It should be a 0 vector of length 2* length(arr).
"""
function PartitionParallel(f, arr, flags::AbstractVector, _SCAN_SIZE = 1024)
    n = length(arr)
    l = nblocks(n, _SCAN_SIZE)

    if (l <= 2)
        return PartitionSerial(f, arr)
    end

    #prefixsum to get index
    Threads.@threads for i = 1:l
        r = get_par_loop_range(i, _SCAN_SIZE, n)
        @inbounds map!(f, uview(flags, r), uview(arr, r))
    end
    PrefixSumParallel_unsafe!(uview(flags,1:n), uview(flags,(n+1):(2*n)), _SCAN_SIZE)
    # PrefixSumSerial!(flags)
    out = zeros(typeof(arr[1]), flags[n])
    outfalse = zeros(typeof(arr[1]), n-flags[n])
    if (flags[1] != 0)
        out[flags[1]] = arr[1]
    else
        outfalse[1] = arr[1]
    end
    Threads.@threads for i = 1:l
        @simd ivdep for j = get_par_loop_range(i, _SCAN_SIZE, n-1)
            if (flags[j] != flags[j+1])
                @inbounds out[flags[j+1]] = arr[j+1]
            else
                @inbounds outfalse[j+1-flags[j+1]] = arr[j+1]
            end
        end
    end
    return out, outfalse
end

function PartitionParallel(f, arr::AbstractVector, _SCAN_SIZE::Int = 1024)
    flags = zeros(Int64, length(arr) * 2)
    # flags = zeros(Int64, length(arr))
    return PartitionParallel(f, arr, flags, _SCAN_SIZE)
end


const SMALL_THRESHOLD  = 2000

function SortParallel!(v::AbstractVector, lo::Int, hi::Int, f = identity, t=similar(v,0))
    @inbounds if lo < hi
        hi-lo <= SMALL_THRESHOLD && return sort!(uview(v, lo: hi),  by=f, alg = MergeSort)

        m = (lo+hi)>>>1
        # (length(t) < m-lo+1) && resize!(t, m-lo+1)

        half = Threads.@spawn SortParallel!(v, m+1, hi, f, uview(t, 1:hi-m))
        SortParallel!(v, lo,  m, f,  uview(t, 1:m-lo+1))
        wait(half)

        j = lo
        while j <= m
            t[j] = v[j]
            j += 1
        end

        i = lo
        k = lo
        while k < j <= hi
            if  f(v[j]) <  f(t[i])
                v[k] = v[j]
                j += 1
            else
                v[k] = t[i]
                i += 1
            end
            k += 1
        end
        while k < j
            v[k] = t[i]
            k += 1
            i += 1
        end
    end

    return v
end
function SortParallel!(v::AbstractVector, t=similar(v,0), f = identity)
    return SortParallel!(v, 1, length(v), f , t)
end



end
