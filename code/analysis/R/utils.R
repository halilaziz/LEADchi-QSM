# LEAD-chi analysis: small utilities (CLI parsing, logging, hashing, Monte Carlo precision, formatting).

parse_args <- function(args) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a)
    key <- gsub("-", "_", sub("^--", "", a))
    if (i < length(args) && !startsWith(args[i + 1], "--")) { out[[key]] <- args[i + 1]; i <- i + 2 }
    else { out[[key]] <- TRUE; i <- i + 1 }
  }
  out
}

.LOG <- new.env()
.LOG$lines <- character(0)
log_msg <- function(...) {
  s <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), "  ", paste0(..., collapse = ""))
  .LOG$lines <- c(.LOG$lines, s)
  message(s)
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

# Monte Carlo precision of a resampling p-value (v1.2 sec. 7.3): p = (k + 1) / (N + 1).
mc_info <- function(p, N, band = REG$mc_band) {
  k <- round(p * (N + 1) - 1)
  ci <- if (k >= 0 && k <= N) stats::binom.test(max(k, 0), N)$conf.int else c(NA, NA)
  list(p = p, N = N, mc_ci95 = as.numeric(ci), in_band = isTRUE(p >= band[1] && p <= band[2]))
}

# Holm step-down over a named vector of two-sided p values.
holm <- function(p) {
  adj <- stats::p.adjust(p, method = "holm")
  list(p = p, p_holm = adj, reject = adj < REG$alpha)
}

f3 <- function(x, d = 3) {
  if (is.null(x) || length(x) == 0) return("NA")
  ifelse(is.na(x), "NA", formatC(as.numeric(x), format = "f", digits = d))
}
fp <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
  if (p < 1e-4) return("<.0001")
  formatC(p, format = "f", digits = 4)
}

# Evaluate an expression; on error return NULL and log it (fits must never abort the run silently).
try_fit <- function(expr, label) {
  res <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    log_msg("warning in ", label, ": ", conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { log_msg("ERROR in ", label, ": ", conditionMessage(e)); NULL })
  res
}

cell_counts <- function(d) {
  tab <- table(factor(d$cell, levels = CELLS))
  c(setNames(as.integer(tab), CELLS), total = nrow(d))
}
cell_str <- function(d) {
  cc <- cell_counts(d)
  paste0(paste(cc[CELLS], collapse = " / "), " (total ", cc["total"], ")")
}
