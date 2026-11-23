calculate_r2 <- function(beta, se, eaf, n) {
  numerator <- 2 * beta^2 * eaf * (1 - eaf)
  denominator <- numerator + 2 * n * eaf * (1 - eaf) * se^2
  r2 <- numerator / denominator
  return(r2)
}

calculate_f_from_r2 <- function(r2, n, k = 1) {
  if (k == 1) {
    f_stat <- (r2 * (n - 2)) / (1 - r2)
  } else {
    f_stat <- (r2 * (n - k - 1)) / (k * (1 - r2))
  }
  return(f_stat)
}
