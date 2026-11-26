# Compile the Stan programs for the cmdstanr backend

The default backend, `rstan` (Stan Development Team, 2025), uses
programs written in Stan (Carpenter et al., 2017) and compiled when the
package is installed. With the optional `cmdstanr` backend each Stan
program is compiled the first time it is used and the executable is
saved in the user cache directory. Call this function once after
installing or updating the package to compile all programs in advance.

## Usage

``` r
compile_models(programs = stan_programs())
```

## Arguments

- programs:

  Names of the Stan programs to compile. By default, all.

## Value

The paths of the compiled executables, invisibly.

## References

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, 76(1), 1–32.
[doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)

Stan Development Team (2025). RStan: The R interface to Stan. R package.
<https://CRAN.R-project.org/package=rstan>
