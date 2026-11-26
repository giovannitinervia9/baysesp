# Generics of base R with S7 methods in this package --------------------------
#
# S7::method(f, class) <- fn assigns to the variable `f`. Using print() and
# plot() directly would create objects called `print` and `plot` in the
# namespace, which hide the generics of base R: the S3 methods of the package
# (such as print.prior_posterior) would then be registered on them and never
# be found. The S7 methods are therefore assigned to these copies instead.
print_generic <- base::print
plot_generic <- base::plot
update_generic <- stats::update
summary_generic <- base::summary
