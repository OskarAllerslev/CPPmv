# datafetch.R
library(quantmod)

# Funktion til at hente data og beregne input til mean-variance optimering
get_mv_data <- function(tickers, start_date, end_date) {
  prices <- list()
  failed_tickers <- c()
  
  for (ticker in tickers) {
    stock_data <- tryCatch(
      getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(stock_data)) {
      # Lukkepris
      prices[[ticker]] <- Cl(stock_data)
    } else {
      failed_tickers <- c(failed_tickers, ticker)
    }
  }
  
  # Hvis ingen gyldige aktier blev hentet
  if (length(prices) == 0) {
    stop("Ingen gyldige aktier blev hentet fra Yahoo Finance.")
  }

  # Sammensæt data i én xts
  prices_df <- do.call(merge, prices)
  colnames(prices_df) <- names(prices)
  
  # Fjern eventuelle NA-rækker
  prices_df <- na.omit(prices_df)

  # Beregn log-afkast
  returns <- diff(log(prices_df))

  # Forventede afkast (middelværdien af hver kolonne)
  expected_returns <- colMeans(returns, na.rm = TRUE)
  
  # Kovariansmatrix
  cov_matrix <- cov(returns, use = "pairwise.complete.obs")
  
  # GEM TIL CSV - uden header og row.names
  write.csv(expected_returns, "expected_returns.csv", row.names = FALSE, col.names = FALSE)
  write.csv(cov_matrix, "cov_matrix.csv", row.names = FALSE, col.names = FALSE)

  if (length(failed_tickers) > 0) {
    cat("Følgende tickers blev ikke fundet på Yahoo Finance:\n",
        paste(failed_tickers, collapse = ", "), "\n")
  }
  
  return(list(expected_returns = expected_returns, cov_matrix = cov_matrix))
}

# Eksempel tickers
tickers <- c("AAPL","MSFT","GOOGL","AMZN")

start_date <- "2020-01-01"
end_date   <- "2025-01-01"

data <- get_mv_data(tickers, start_date, end_date)
print("CSV-filer: expected_returns.csv og cov_matrix.csv er nu oprettet uden headers.")
