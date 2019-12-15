__precompile__()
module BENCH
using Distributed
using UnsafeArrays

function PrefixSumSerial!(arr::A) where {A}
    l = length(arr)
    k = ceil(Int, log2(l))
    @inbounds for j=1:k
         for i=2^j:2^j:min(l, 2^k)
            arr[i] = arr[i-2^(j-1)] + arr[i]
        end
    end
    @inbounds for j=(k-1):-1:1
        for i=3*2^(j-1):2^j:min(l, 2^k)
            arr[i] = arr[i-2^(j-1)] + arr[i]
        end
    end
end

# https://simondanisch.github.io/ReferenceImages/gallery/parallel_prefix_sum/index.html
function PrefixSumParallel!(arr::A) where {A}
    l = length(arr)
    k = ceil(Int, log2(l))
    for j=1:k
         Threads.@threads for i=2^j:2^j:min(l, 2^k)
             let
                @inbounds  arr[i] += arr[i-2^(j-1)]
            end
        end
    end
    for j=(k-1):-1:1
        Threads.@threads for i=3*2^(j-1):2^j:min(l, 2^k)
            let
                @inbounds  arr[i] += arr[i-2^(j-1)]
            end
        end
    end
end


function mapfilter(pred, f, itr, res)
    for x in itr
        pred(x) && f(res, x)
    end
    res
end

function FilterSerial(f, a)
    mapfilter(f, push!, a, empty(a))
end
# 
# function FilterParallel(f, arr, flags, _SCAN_SIZE = 1024)
#     n = length(arr)
#     l = nblocks(n, _SCAN_SIZE)
#
#     if (l <= 2)
#         return FilterSerial(f, arr)
#     end
#
#     #prefixsum to get index
#     Threads.@threads for i = 1:l
#         # let
#         for j = get_par_loop_range(i, _SCAN_SIZE, n)
#             @inbounds flags[j] = f(arr[j])
#         end
#         # end
#     end
#     PrefixSumParallel!(flags)
#     # PrefixSumSerial!(flags)
#     out = zeros(flags[n])
#     if (flags[1] != 0)
#         out[flags[1]] = arr[1]
#     end
#     Threads.@threads for i = 1:l
#         # let
#         for j = get_par_loop_range(i, _SCAN_SIZE, n)
#             if (j<n && flags[j] != flags[j+1])
#                 @inbounds out[flags[j+1]] = arr[j+1]
#             end
#         end
#         # end
#     end
#     return out
# end

function FilterParallelD(f, a)
    return a[pmap(f, a)]
end

function PartitionSerial(f, a)
    filter(f,a), filter(x->!f(x),a)
end

function PartitionParallel(f, a)
    v = Threads.@spawn filter(x->!f(x),a)
    filter(f,a), fetch(v)
end

# http://www.csd.uwo.ca/~moreno/cs2101a_moreno/Parallel_computing_with_Julia.pdf
function merge(a, b, istart, mid, iend)
    n = iend - istart + 1
    nb = iend - mid
    na = mid - istart + 1
    c = zeros(n)
    s = 1
    m = 1
    for tem = 1:n
    if s <= na && (m > nb || a[s] <= b[m])
    c[tem] = a[s]
    s=s+1
    else
    c[tem] = b[m]
    m=m+1
    end
    end
    c
end

function mergesort_parallel(data, istart, iend)
    if (iend - istart <= 2500000)
        sort(uview(data, istart: iend),  alg = MergeSort)
    else
        mid = ifloor((istart + iend)/2)
        a = Threads.@spawn mergesort_parallel(data, istart, mid)
        b = mergesort_parallel(data,mid+1, iend)
        c = merge(fetch(a), b, istart, mid, iend)
    end
end

function mergesort_parallel(data)
    mergesort_parallel(data, 1, length(data))
end

# https://stackoverflow.com/questions/47235390/how-can-i-parallelize-sorting/47235391#47235391
function blockranges(nblocks, total_len)
    rem = total_len % nblocks
    main_len = div(total_len, nblocks)

    starts=Int[1]
    ends=Int[]
    for ii in 1:nblocks
        len = main_len
        if rem>0
            len+=1
            rem-=1
        end
        push!(ends, starts[end]+len-1)
        push!(starts, ends[end] + 1)
    end
    @assert ends[end] == total_len
    starts[1:end-1], ends
end


function threadedsort!(data::Vector, f = identity)
    starts, ends = blockranges(Threads.nthreads(), length(data))

    # Sort each block
    Threads.@threads for (ss, ee) in collect(zip(starts, ends))
        @inbounds sort!(@views(data[ss:ee]),  by=f, alg = MergeSort)
    end


    # Go through each sorted block taking out the smallest item and putting it in the new array
    # This code could maybe be optimised. see https://stackoverflow.com/a/22057372/179081
    ret = similar(data) # main bit of allocation right here. avoiding it seems expensive.
    # Need to not overwrite data we haven't read yet
    @inbounds for ii in eachindex(ret)
        minblock_id = 1
        ret[ii]=data[starts[1]]
        @inbounds for blockid in 2:length(starts) # findmin allocates a lot for some reason, so do the find by hand. (maybe use findmin! ?)
            ele = data[starts[blockid]]
            if f(ret[ii]) > f(ele)
                ret[ii] = ele
                minblock_id = blockid
            end
        end
        starts[minblock_id]+=1 # move the start point forward
        if starts[minblock_id] > ends[minblock_id]
            deleteat!(starts, minblock_id)
            deleteat!(ends, minblock_id)
        end
    end
    data.=ret  # copy back into orignal as we said we would do it inplace
    return data
end

end
