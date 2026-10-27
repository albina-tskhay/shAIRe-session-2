# =============================================================================
# 00_install_packages.R
# Installs (only if missing) and checks all packages for the workshop.
# Usage: open in RStudio and click "Source", or run source("00_install_packages.R")
# =============================================================================

install_workshop_packages <- function(include_optional = FALSE,
                                      repos = "https://cloud.r-project.org") {
  required <- c(
    tidyverse = "data manipulation and visualization",
    tableone  = "descriptive statistics (Table 1)",
    cobalt    = "covariate balance assessment",
    sandwich  = "robust (sandwich) standard errors",
    lmtest    = "coefficient testing with robust SEs",
    boot      = "bootstrap confidence intervals",
    survey    = "weighted regression (for IPTW)"
  )
  optional <- c(MatchIt = "matching methods (optional)")
  pkgs <- if (include_optional) c(required, optional) else required
  
  # 1. Install only the packages that are not already installed
  missing <- setdiff(names(pkgs), rownames(installed.packages()))
  if (length(missing) > 0) {
    message("Installing ", length(missing), " package(s): ",
            paste(missing, collapse = ", "), "\nThis may take a few minutes...")
    tryCatch(
      install.packages(missing, repos = repos),
      error = function(e) message("Installation error: ", conditionMessage(e))
    )
  } else {
    message("All packages are already installed.")
  }
  
  # 2. Check that every package can be loaded
  ok <- vapply(names(pkgs), requireNamespace, logical(1), quietly = TRUE)
  status <- data.frame(
    package = names(pkgs),
    purpose = unname(pkgs),
    status  = ifelse(ok, "OK", "FAILED")
  )
  print(status, row.names = FALSE)
  
  if (all(ok)) {
    message("\nAll set for the workshop!")
  } else {
    message("\nThese packages failed: ", paste(names(pkgs)[!ok], collapse = ", "),
            "\nRun the function again, or see the Troubleshooting section of the guide.")
  }
  invisible(status)
}

install_workshop_packages()



