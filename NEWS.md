
# spaero 0.1.1

- Correct autocorrelation calculation. The previous version divided
  the autocovariance by the variance at the most recent time
  point. The current version divides by the geometric mean of the
  variance at each of the two time points, matching standard
  practice. The formula in the vignette for the autocorrelation has
  been changed accordingly.

- Clean up sloppy usage of the term statistic in the documentation.

# spaero 0.1.0

- Initial release
