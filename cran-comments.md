## Submission

This is the first submission of baysesp.

## Test environments

* local Windows 11, R 4.6.0
* GitHub Actions: windows-latest, macos-latest (R release); ubuntu-latest
  (R devel, release and oldrel-1)
* R-hub: linux, macos-arm64 and windows (R devel)
* win-builder: R devel and release

## R CMD check results

0 errors | 0 warnings | 1 note on Windows; 3 notes expected on Linux and macOS

* New submission.

* Suggests or Enhances not in mainstream repositories: cmdstanr.
  'cmdstanr' is an optional backend, available from the repository declared in
  `Additional_repositories` (https://stan-dev.r-universe.dev). It is used only
  conditionally (`requireNamespace()`); examples, tests and the vignette run
  without it.

* installed size is large; sub-directory 'libs' is large. (Linux and macOS)
  The package contains three Stan models precompiled with 'rstantools', as
  usual for packages built on 'rstan'; the size of the compiled code depends
  on the platform.

* GNU make is a SystemRequirements. (Linux and macOS)
  Required by 'rstantools' to compile the Stan models.

## Run time

Examples that fit Stan models are wrapped in `\donttest{}`. The tests that fit
Stan models are skipped on CRAN with `testthat::skip_on_cran()`. Examples, tests
and the vignette run the chains sequentially (no more than one core). On the
local machine the examples (with `--run-donttest`) take about 16 s, the tests
about 5 s and the vignette about 1 minute.
