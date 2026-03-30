# Contributing to serieshaz

Thank you for your interest in contributing to serieshaz!

## Reporting Bugs

Open an issue at <https://github.com/queelius/serieshaz/issues> with:

- A minimal reproducible example
- Your R version (`sessionInfo()`)
- The error message or unexpected behavior

## Suggesting Features

Open an issue describing the feature and its use case. For series system
extensions (new component types, new methods), consider whether they belong
in serieshaz or in the upstream [flexhaz](https://github.com/queelius/flexhaz)
package.

## Pull Requests

1. Fork the repository and create a feature branch
2. Follow the existing code style (S3 classes, closure-returning pattern)
3. Add tests for new functionality (`tests/testthat/`)
4. Run `devtools::check()` and ensure 0 errors, 0 warnings
5. Update documentation if needed (`roxygen2::roxygenise()`)
6. Submit a pull request against `main`

## Code of Conduct

Please note that this project follows the
[Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md).
