# Bayesian Accuracy of Multiple Diagnostic Tests under Partial Verification

`{baysesp}` estimates the sensitivity, specificity and predictive values
of one or more binary diagnostic tests performed on the same subjects,
when the disease status is ascertained by the reference standard only
for a subset of them. Verification that depends on the test results
produces verification bias if it is ignored (Begg and Greenes, 1983).
The package takes a Bayesian approach, which places priors on the
parameters that are not identified (Buzoianu and Kadane, 2008). The
models are computed in Stan (Carpenter et al., 2017).

## References

Begg, C. B. and Greenes, R. A. (1983). Assessment of diagnostic tests
when disease verification is subject to selection bias. *Biometrics*,
39(1), 207–215. [doi:10.2307/2530820](https://doi.org/10.2307/2530820)

Buzoianu, M. and Kadane, J. B. (2008). Adjusting for verification bias
in diagnostic test evaluation: A Bayesian approach. *Statistics in
Medicine*, 27(13), 2453–2473.
[doi:10.1002/sim.3099](https://doi.org/10.1002/sim.3099)

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P. and Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, 76(1), 1–32.
[doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)

## See also

Useful links:

- <https://github.com/giovannitinervia9/baysesp>

- <https://giovannitinervia9.github.io/baysesp/>

- Report bugs at <https://github.com/giovannitinervia9/baysesp/issues>

## Author

**Maintainer**: Giovanni Tinervia <giovannitinervia9@gmail.com>
