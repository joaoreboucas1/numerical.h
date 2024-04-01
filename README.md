# Numerical Tools in C
`numerical.h` is a header-only library written in C containing some of the simplest numerical algorithms.

**WARNING**: This is a personal project just for fun, so the public API is completely unstable and for now I don't really want to optimize the routines. Use this library at your own risk.

## Implemented routines
  - Linear algebra utilities: a `Matrix` struct with matrix multiplication, LU decomposition, QR decomposition, determinant, inverse matrix, solving linear systems of equations, eigenvalue finding
  - Numerical integration: trapezoid rule for functions of one variable (that may accept additional parameters)
  - ODEs: ODE system solver via Euler's algorithm
  - Interpolation: Linear interpolator and Cubic Spline interpolator
  - Root finding: secant method for functions of a single variable
  - Optimization: minimization of functions of a single variable using gradient descent with finite difference
  - FFT

## Examples
We provide C examples for each numerical routine in the `demos/` folder. We also provide a `Makefile` that builds executables for demonstration. To build the demos, simply go to the `demos/` folder and use `make`:

```console
    $ cd ./demos
    $ make
```

Run the demos to view the different features.

## Documentation
I follow the philosophy of Raylib: the `demos/` folder *is* the documentation. Read the demos source code to study the library API. 
