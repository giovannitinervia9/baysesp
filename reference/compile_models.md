# Compile the Stan programs for the cmdstanr backend

The default backend, `rstan`, uses programs compiled when the package is
installed. With the optional `cmdstanr` backend each Stan program is
compiled the first time it is used and the executable is saved in the
user cache directory. Call this function once after installing or
updating the package to compile all programs in advance.

## Usage

``` r
compile_models(programs = stan_programs())
```

## Arguments

- programs:

  Names of the Stan programs to compile. By default, all.

## Value

The paths of the compiled executables, invisibly.
