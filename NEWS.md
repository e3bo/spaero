# spaero 0.2.0

- Add transmission argument to create_simulator to allow for
  frequency-dependent transmission. Density-dependent transmission
  remains the default model.

- Add vector of first difference of the variance vector produced by
  get_stats. This change makes it easier to use the convexity of the
  variance time series as an early warning signal. The name of the
  vector in the stats list is variance_first_diff. Note that this
  change makes the abbreviation stats$var ambiguous. Code using that
  abbreviation to obtain the vector of variance estimates should
  substitute in stats$variance.

- To the output of get_stats(), add list taus containing Kendall's
  correlation coefficient of the elements of each time series in the
  stats list in the output with time.

- Ensure variance and kurtosis esimates are non-negative. When using
  local linear for estimating statistics, it was possible in previous
  versions for negative values to occur.

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
