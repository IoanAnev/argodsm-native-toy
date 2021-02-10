# Himeno

The Himeno benchmark [1] version 3.0 is a test program measuring the cpu performance in MFLOPS. Point-Jacobi method is employed in this Pressure Poisson equation solver as this method can be easily vectorized and be parallelized.

The size of this kernel benchmark can be choosen from the following sets of

`[mimax][mjmax][mkmax]:`\
`ext. small`&nbsp;&nbsp;&nbsp;&nbsp;: 33,33,65\
`small`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: 65,65,129\
`medium`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: 129,129,257\
`large`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: 257,257,513\
`ext.large`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;: 513,513,1025

## Compilation Flags:

`-O3 -mavx2`

Note: Include the CFLAG -mavx2 only if the machine provides AVX2 support (To check AVX2 support use the command : grep avx2 /proc/cpuinfo).

## References:

[1] http://accc.riken.jp/en/supercom/himenobmt, download at
http://accc.riken.jp/en/supercom/himenobmt/download/98-source/
