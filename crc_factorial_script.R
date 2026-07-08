# =============================================================================
# CRC factorial GLS analysis
# =============================================================================
# Author:  Kamil Myszczynski
# Date:    2026-07-07
# Version: 1.0

# Purpose: Run pooled and per-cell-line factorial models for CRC assay outcomes.
# Input:   crc_long.csv with columns assay, value, cell_line, inhibitor, ASA, FU
# Models:  Reduced 2-way models for 7AAD/CFSE/Necrosis/ATP/Caspase/Lactate/Size.
# Notes:   Pooled models use GLS with cell-line-specific residual variances
#          for consistency; variance LRTs are diagnostic only.
# =============================================================================

suppressPackageStartupMessages({
  library(nlme)      
  library(emmeans)   
  library(car)       
  library(dplyr)
  library(stringr)
})

# -----------------------------------------------------------------------------
# 0. READ & TRANSFORM
# -----------------------------------------------------------------------------

sample_tab <- read.csv("crc_long_simulated.csv")

options(contrasts = c("contr.sum", "contr.poly"))
sample_tab$cell_line <- factor(sample_tab$cell_line, levels = c("HCT116", "HT29"))
sample_tab$inhibitor <- factor(sample_tab$inhibitor, levels = c("None", "3-MA", "CQ"))
sample_tab$ASA <- ordered(sample_tab$ASA, levels = c(0, 2.5, 5))
sample_tab$FU  <- ordered(sample_tab$FU,  levels = c(0, 50, 75))

# Transform responses per assay; Caspase and Size are log-transformed.
TRANSFORM <- c("7AAD" = "identity", "CFSE" = "identity", "Necrosis" = "identity",
               "ATP" = "identity", "Lactate" = "identity",
               "Caspase" = "log", "Size" = "log")
resp <- function(x, how) if (how == "log") log(x) else x

# -----------------------------------------------------------------------------
# 1. MODEL FITTERS
# -----------------------------------------------------------------------------

# Diagnostic LRT for equal vs cell-line-specific variance; final pooled models retain varIdent for consistency.
check_variance <- function(df, reduced = FALSE) {
  f <- if (reduced) {
    y ~ (cell_line + ASA + FU + inhibitor)^2
  } else {
    y ~ cell_line * ASA * FU * inhibitor
  }
  
  m_eq  <- gls(f, data = df, method = "ML")
  m_var <- gls(f, data = df, weights = varIdent(form = ~1 | cell_line),
               method = "ML")
  
  cat("  variance check (LRT equal vs per-line):\n")
  print(anova(m_eq, m_var))
}

# Fit pooled GLS model; reduced models test main effects plus all 2-way interactions.
fit_pooled <- function(df, reduced = FALSE) {
  f <- if (reduced) y ~ (cell_line + ASA + FU + inhibitor)^2
  else         y ~ cell_line * ASA * FU * inhibitor
  m <- gls(f, data = df, weights = varIdent(form = ~1 | cell_line),
           method = "REML")
  cat("\n  Type III (marginal) tests -- pooled model:\n")
  print(anova(m, type = "marginal"))
  m
}

# Fit separate per-line LM models for cell-line-specific Type III tests.
fit_per_line <- function(df, reduced = FALSE) {
  out <- list()
  for (cl in levels(df$cell_line)) {
    sub <- droplevels(subset(df, cell_line == cl))
    f <- if (reduced) y ~ (ASA + FU + inhibitor)^2 else y ~ ASA * FU * inhibitor
    m <- lm(f, data = sub)
    cat(sprintf("\n  Type III -- %s only:\n", cl))
    print(Anova(m, type = 3))
    out[[cl]] <- m
  }
  out
}

# Build emmeans reference grids directly to avoid recover_data issues with varIdent-weighted GLS models.
emm_gls <- function(model, df, specs) {
  rg <- qdrg(formula = formula(model)[-2L],          # RHS only (drop response)
             data    = df,
             coef    = coef(model),
             vcov    = vcov(model),
             df      = length(residuals(model)) - length(coef(model)))
  emmeans(rg, specs)
}

# Inhibitor contrasts are by cell line; ASA x inhibitor contrasts are pooled for reduced models lacking the 3-way term.
contrasts_inhibitor <- function(model, df, by_cell = TRUE) {
  cat("\n  inhibitor vs None, by cell line (averaged over ASA and FU):\n")
  emm <- emm_gls(model, df, ~ inhibitor | cell_line)
  print(contrast(emm, method = "trt.vs.ctrl", ref = 1))   # ref = None (level 1)
  
  cat(sprintf("\n  inhibitor x ASA interaction contrasts (%s):\n",
              if (by_cell) "by cell line" else "pooled -- reduced model, no 3-way"))
  spec <- if (by_cell) ~ inhibitor * ASA | cell_line else ~ inhibitor * ASA
  emm2 <- emm_gls(model, df, spec)
  print(contrast(emm2, interaction = c("trt.vs.ctrl", "poly")))
}

# Run one assay: transform response, check variance, fit pooled and per-line models, then compute contrasts.
run_inferential <- function(assay_name, reduced = FALSE) {
  cat("\n==================================================================\n")
  cat(sprintf("ASSAY: %s  (%s)\n", assay_name,
              if (reduced) "REDUCED: mains + 2-way" else "FULL 4-way"))
  cat("==================================================================\n")
  df <- droplevels(subset(sample_tab, assay == assay_name))
  df$y <- resp(df$value, TRANSFORM[[assay_name]])   
  check_variance(df, reduced = reduced)
  mp <- fit_pooled(df, reduced = reduced)
  contrasts_inhibitor(mp, df, by_cell = !reduced)   
  ml <- fit_per_line(df, reduced = reduced)
  invisible(list(pooled = mp, per_line = ml))
}

# Descriptive summaries by full treatment combination.
describe <- function(assay_name) {
  df <- subset(sample_tab, assay == assay_name)
  aggregate(value ~ cell_line + inhibitor + ASA + FU, data = df,
            FUN = function(v) c(n = length(v), median = median(v),
                                q25 = quantile(v, .25), q75 = quantile(v, .75)))
}

# -----------------------------------------------------------------------------
# 2. RUN
# -----------------------------------------------------------------------------

# Reduced models
fit_7AAD     <- run_inferential("7AAD", reduced = TRUE)
fit_CFSE     <- run_inferential("CFSE", reduced = TRUE)
fit_Necrosis <- run_inferential("Necrosis", reduced = TRUE)
fit_ATP      <- run_inferential("ATP",     reduced = TRUE)
fit_Caspase  <- run_inferential("Caspase", reduced = TRUE)
fit_Lactate  <- run_inferential("Lactate", reduced = TRUE)
fit_Size     <- run_inferential("Size",    reduced = TRUE)
