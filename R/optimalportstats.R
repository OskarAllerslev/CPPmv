# optimalportstats.R

# Læs optimale vægte, forventede afkast og kovariansmatrix (ingen headers)

# 1) Optimal weights
opt_w <- as.matrix(read.csv("optimal_weights.csv", header = FALSE))
opt_w <- as.numeric(opt_w[,1])   # konverter til ren vektor

# 2) Expected returns
exp_ret <- as.matrix(read.csv("expected_returns.csv", header = FALSE))
exp_ret <- as.numeric(exp_ret[,1])  # vektor

# 3) Covariance matrix
cov_mat <- as.matrix(read.csv("cov_matrix.csv", header = FALSE))

# Tjek dimensioner
cat("Dimension af opt_w: ", length(opt_w), "\n")
cat("Dimension af exp_ret: ", length(exp_ret), "\n")
cat("Dimension af cov_mat: ", dim(cov_mat), "\n")

# Forvent, at length(opt_w) == length(exp_ret) == ncol(cov_mat) == nrow(cov_mat)
if (length(opt_w) != length(exp_ret)) {
  stop("Fejl: opt_w og exp_ret har forskellige længder!")
}
if (nrow(cov_mat) != ncol(cov_mat)) {
  stop("Fejl: Kovariansmatrix ikke kvadratisk!")
}
if (nrow(cov_mat) != length(opt_w)) {
  stop("Fejl: Dimensioner i kovariansmatrix og vægte/afkast matcher ikke!")
}

# Beregn porteføljens forventede afkast
port_return <- sum(opt_w * exp_ret)

# Beregn porteføljens varians
port_variance <- t(opt_w) %*% cov_mat %*% opt_w
port_sd <- sqrt(port_variance)

cat("Porteføljens forventede afkast:", port_return, "\n")
cat("Porteføljens varians:", port_variance, "\n")
cat("Porteføljens standardafvigelse:", port_sd, "\n")
