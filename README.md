# Numerical Tools in C
`numerical.h` is a header-only library written in C containing some of the simplest numerical algorithms.

**WARNING**: This is a personal project just for fun, so the public API is completely unstable. Use this library at your own risk.

## Implemented routines
  - Solving systems of linear algebraic equations
  - Trapezoid rule for numerical integration
  - 1D ODE solver via Euler's algorithm
  - Linear interpolator and Cubic Spline interpolator

## Examples
We provide C examples for each numerical routine in the `demos/` folder. We also provide a `Makefile` that builds executables for demonstration. To build the demos, simply go to the `demos/` folder and use `make`:

```console
    $ cd ./demos
    $ make
```

Run the demos to view the different features.

## Documentation
I follow the philosophy of Raylib: the `demos/` folder *is* the documentation. Read the demos source code to study the library API. 
