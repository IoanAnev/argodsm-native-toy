# Toy ArgoDSM benchmarks

Contains toy benchmarks that are used to evaluate the performance of ArgoDSM.

Current benchmarks include:
1. [daxpy](./daxpy/)
2. [himeno](./himeno/)
3. [matmul](./matmul/)
4. [matvec](./matvec/)
5. [stream](./stream/)

## Building

Build individually each application using `make`.

This assumes you have already installed [argodsm](https://github.com/etascale/argodsm) and set its path in the Makefile.

## Benchmarks

### **daxpy**

Performs the computation: `y += a * x` where `x`, `y`, are two vectors of `double`, of size `N` and `a` is a scalar.

#### **[Usage]**

```sh
mpirun $OMPIFLAGS ./daxpy N ITER [CHECK]

where:

N       the size of the vectors x and y
ITER    number of iterations for each to execute the computation
CHECK   an optional parameter that enables checks to make sure the computation is correct
```

### **himeno**

It measures the speed of major loops for solving Poisson’s equation using the Jacobi iteration method. Choose between different problem sizes by setting `SIZE`.

#### **[Usage]**

```sh
mpirun $OMPIFLAGS ./himenobmtxpa SIZE

where:

SIZE    XS, S, M, L, XL
```

> **Note:** For verification, output to a file and diff with the serial version. Code is included.

### **matmul**

Calculates the product of two matrices. Choose matrix dimensions by setting `N`, and select between the different implementations based on loop ordering by setting `I`.

#### **[Usage]**

```sh
mpirun $OMPIFLAGS ./matmul -s N -i I [-v]

where:

N       NxN matrix dimensions
I       0 (jik), 1 (ijk), or 2 (ikj)
```

> **Note:** For verification, simply pass `-v` to the executable.

### **matvec**

Calculates the matrix-vector product: `y = A * x`, where `A` is a matrix of `M` rows and `N` columns.

#### **[Usage]**

```sh
mpirun $OMPIFLAGS ./matvec M N ITER [CHECK]

where:

M       the rows of matrix A
N       the columns of matrix A
ITER    number of iterations for which to execute the computations
CHECK   an optional parameter that enables checks to make sure the computation is correct
```

### **stream**

It is a set of multiple kernel operations,
that is, sequential accesses over array data with simple arithmetic. Kernel operations include:

+ `Copy`&nbsp;&nbsp;&nbsp;&nbsp;-> a(i) = b(i)
+ `Scale`&nbsp;&nbsp;-> a(i) = q*b(i)
+ `Add`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-> a(i) = b(i) + c(i)
+ `Triad`&nbsp;&nbsp;-> a(i) = b(i) = q*c(i)

The benchmark options are chosen in compilation and thus, are set in the Makefile.

#### **[Usage]**

```sh
mpirun $OMPIFLAGS ./stream

compilation options:

STREAM_ARRAY_SIZE the size of the working arrays
RANDACC           0 or 1 to use random indexes to access the arrays
TEAMINIT          0 or 1 to use team process initialization on the arrays
```
