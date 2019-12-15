# JuliaParallelToolBox

To run experiments for prefixsum with P threads:

```bash
./run_julia P prefixsum_exp.jl
```


To run experiments for filter with P threads:

```bash
./run_julia P filter_exp.jl
```


To run experiments for partition with P threads:

```bash
./run_julia P partition_exp.jl
```

To run experiments comparing Julia and C++:

```bash
./run_julia P prefixsum.jl
./run_cpp P
```
