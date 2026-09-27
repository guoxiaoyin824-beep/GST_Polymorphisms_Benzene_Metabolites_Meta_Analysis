# ============================================================
# Benzene urinary biomarkers and GST genotypes meta-analysis
# Input: biomarkers Excel workbook (one gene-metabolite per sheet)
# Effect measure: standardized mean difference (Hedges' g)
# Subgroup variable: Region
#
# Effect direction:
#   GSTM1 / GSTT1: Positive minus Null
#   GSTP1: Ile/Ile-related group minus variant genotype/carrier
# Positive SMD means a higher biomarker concentration in the
# first genotype named in the comparison.
# ============================================================

# Install once if needed:
# install.packages(c("meta", "readxl"))

required_packages <- c("meta", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Please install the missing package(s) first: ",
    paste(missing_packages, collapse = ", "),
    "\nRun: install.packages(c(\"meta\", \"readxl\"))"
  )
}

library(meta)
library(readxl)

# ---------- 1. Paths and analysis settings ----------
# Change this path to your edited workbook. If it is not found,
# R will open a file-selection window.
data_file <- "D:/科研/meta/benzene/biomarkers.xlsx"

if (!file.exists(data_file)) {
  message("Excel file not found. Please choose the biomarkers workbook manually.")
  data_file <- file.choose()
}

base_dir <- dirname(data_file)
output_root <- file.path(base_dir, "SMD_meta_results")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# If both All subjects and smoker/non-smoker strata are present for the same
# article, use the two non-overlapping strata and remove the overlapping total.
# Set FALSE if you prefer the overall row instead (then remove strata manually).
prefer_smoking_strata <- TRUE

# Same decision rule as the reference script:
# heterogeneity P < 0.10 -> random-effects result; otherwise common-effect result.
heterogeneity_cutoff <- 0.10

# Random-effects settings. REML estimates tau^2; Hartung-Knapp is used for
# the random-effects confidence interval.
tau_method <- "REML"
random_ci_method <- "HK"

# Width used for every forest plot PDF (overall, Region subgroup,
# and leave-one-out sensitivity plots).
forest_pdf_width <- 20

expected_sheets <- c(
  "GSTM1_SPMA", "GSTM1_ttMA",
  "GSTT1_SPMA", "GSTT1_ttMA",
  "GSTP1_SPMA", "GSTP1_ttMA"
)

# ---------- 2. General helper functions ----------
fill_down <- function(x) {
  x <- as.character(x)
  for (i in seq_along(x)) {
    if (i > 1 && (is.na(x[i]) || trimws(x[i]) == "")) x[i] <- x[i - 1]
  }
  x
}

num_clean <- function(x) {
  suppressWarnings(as.numeric(trimws(as.character(x))))
}

clean_column_names <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

format_year <- function(x) {
  sub("\\.0$", "", trimws(as.character(x)))
}

fmt_num <- function(x, digits = 2) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("")
  sprintf(paste0("%.", digits, "f"), x)
}

fmt_p <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_p4 <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("")
  if (p < 0.0001) return("<0.0001")
  sprintf("%.4f", p)
}

fmt_i2 <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("")
  if (x <= 1) x <- x * 100
  sprintf("%.1f%%", x)
}

fmt_ci <- function(te, lower, upper) {
  if (any(is.na(c(te, lower, upper)))) return("")
  sprintf("%.2f (%.2f, %.2f)", te, lower, upper)
}

get_first <- function(obj, candidates) {
  for (nm in candidates) {
    if (!is.null(obj[[nm]]) && length(obj[[nm]]) > 0) return(obj[[nm]][1])
  }
  NA_real_
}

write_csv_safe <- function(x, file) {
  tryCatch(
    write.csv(x, file, row.names = FALSE, na = "", fileEncoding = "UTF-8"),
    error = function(e) write.csv(x, file, row.names = FALSE, na = "")
  )
}

rbind_fill <- function(x, y) {
  all_names <- union(names(x), names(y))
  for (nm in setdiff(all_names, names(x))) x[[nm]] <- rep(NA, nrow(x))
  for (nm in setdiff(all_names, names(y))) y[[nm]] <- rep(NA, nrow(y))
  rbind(x[, all_names, drop = FALSE], y[, all_names, drop = FALSE])
}

safe_pdf <- function(file, expr, height = 7, width = 13) {
  pdf(file = file, height = height, width = width, onefile = TRUE)
  on.exit(dev.off(), add = TRUE)
  tryCatch(
    force(expr),
    error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Plot error:", e$message), cex = 0.8)
    }
  )
}

make_study_label <- function(dat) {
  base <- paste(dat$Author, dat$Year)
  stratum <- trimws(as.character(dat$Population_Stratum))
  show_stratum <- !is.na(stratum) & stratum != "" &
    !tolower(stratum) %in% c("all subjects", "overall", "total")
  ifelse(show_stratum, paste0(base, "(", stratum, ")"), base)
}

# ---------- 3. Read and clean each worksheet ----------
read_sheet_clean <- function(sheet_name) {
  dat <- as.data.frame(read_excel(data_file, sheet = sheet_name))
  names(dat) <- clean_column_names(names(dat))

  required_metadata <- c("Author", "Year", "Country", "Region")
  missing_metadata <- setdiff(required_metadata, names(dat))
  if (length(missing_metadata) > 0) {
    stop(sheet_name, " is missing: ", paste(missing_metadata, collapse = ", "))
  }

  if (!"Population_Stratum" %in% names(dat)) dat$Population_Stratum <- "All subjects"

  fill_cols <- intersect(
    c(
      "Author", "Year", "Country", "Region", "Ethnicity",
      "Detecting_method", "Source", "Population_Stratum", "Unit"
    ),
    names(dat)
  )
  for (cc in fill_cols) dat[[cc]] <- fill_down(dat[[cc]])

  dat$Author <- trimws(as.character(dat$Author))
  dat$Year <- format_year(dat$Year)
  dat$Region <- trimws(as.character(dat$Region))
  dat$Region[is.na(dat$Region) | dat$Region == ""] <- "Not reported"
  dat$Population_Stratum <- trimws(as.character(dat$Population_Stratum))

  dat
}

split_overlapping_totals <- function(dat) {
  empty_excluded <- dat[0, , drop = FALSE]
  empty_excluded$Exclusion_reason <- character(0)
  if (!prefer_smoking_strata || nrow(dat) == 0) {
    return(list(included = dat, excluded = empty_excluded))
  }

  key <- paste(dat$Author, dat$Year, sep = "__")
  pop <- tolower(dat$Population_Stratum)
  is_smoker <- grepl("smoker", pop) & !grepl("non[- ]?smoker", pop)
  is_nonsmoker <- grepl("non[- ]?smoker", pop)
  is_total <- pop %in% c("all subjects", "overall", "total")

  drop <- rep(FALSE, nrow(dat))
  for (k in unique(key)) {
    idx <- which(key == k)
    if (any(is_smoker[idx]) && any(is_nonsmoker[idx])) {
      drop[idx[is_total[idx]]] <- TRUE
    }
  }
  included <- dat[!drop, , drop = FALSE]
  excluded <- dat[drop, , drop = FALSE]
  excluded$Exclusion_reason <- "Overlapping overall row removed because smoking strata were used"
  list(included = included, excluded = excluded)
}

# ---------- 4. Convert worksheets into two-group comparisons ----------
prepare_two_group <- function(dat, n_e, mean_e, sd_e, n_c, mean_c, sd_c,
                              experimental_label, control_label) {
  required <- c(n_e, mean_e, sd_e, n_c, mean_c, sd_c)
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing analysis columns: ", paste(missing_cols, collapse = ", "))
  }

  for (cc in required) dat[[cc]] <- num_clean(dat[[cc]])

  dat$n_e <- dat[[n_e]]
  dat$mean_e <- dat[[mean_e]]
  dat$sd_e <- dat[[sd_e]]
  dat$n_c <- dat[[n_c]]
  dat$mean_c <- dat[[mean_c]]
  dat$sd_c <- dat[[sd_c]]
  dat$Experimental <- experimental_label
  dat$Control <- control_label

  valid <- complete.cases(dat[, c("n_e", "mean_e", "sd_e", "n_c", "mean_c", "sd_c")]) &
    dat$n_e >= 2 & dat$n_c >= 2 & dat$sd_e > 0 & dat$sd_c > 0

  excluded <- dat[!valid, , drop = FALSE]
  excluded$Exclusion_reason <- rep(
    "Incomplete n/mean/SD, n < 2, or SD <= 0",
    nrow(excluded)
  )

  included <- dat[valid, , drop = FALSE]
  overlap_split <- split_overlapping_totals(included)
  included <- overlap_split$included
  excluded <- rbind_fill(excluded, overlap_split$excluded)
  included$Study <- make_study_label(included)

  list(included = included, excluded = excluded)
}

combine_two_groups <- function(n1, mean1, sd1, n2, mean2, sd2) {
  n <- n1 + n2
  mean <- (n1 * mean1 + n2 * mean2) / n
  sd <- sqrt(
    (
      (n1 - 1) * sd1^2 + (n2 - 1) * sd2^2 +
        n1 * (mean1 - mean)^2 + n2 * (mean2 - mean)^2
    ) / (n - 1)
  )
  list(n = n, mean = mean, sd = sd)
}

prepare_gstp1_carrier <- function(dat) {
  cols <- c(
    "n_Ile_Ile", "Mean_Ile_Ile", "SD_Ile_Ile",
    "n_Ile_Val", "Mean_Ile_Val", "SD_Ile_Val",
    "n_Val_Val", "Mean_Val_Val", "SD_Val_Val"
  )
  missing_cols <- setdiff(cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing GSTP1 carrier columns: ", paste(missing_cols, collapse = ", "))
  }
  for (cc in cols) dat[[cc]] <- num_clean(dat[[cc]])

  valid <- complete.cases(dat[, cols]) &
    dat$n_Ile_Ile >= 2 & dat$n_Ile_Val >= 2 & dat$n_Val_Val >= 2 &
    dat$SD_Ile_Ile > 0 & dat$SD_Ile_Val > 0 & dat$SD_Val_Val > 0

  excluded <- dat[!valid, , drop = FALSE]
  excluded$Exclusion_reason <- rep(
    "Incomplete Ile/Val or Val/Val data for carrier pooling",
    nrow(excluded)
  )

  dat <- dat[valid, , drop = FALSE]
  if (nrow(dat) > 0) {
    combined <- mapply(
      combine_two_groups,
      dat$n_Ile_Val, dat$Mean_Ile_Val, dat$SD_Ile_Val,
      dat$n_Val_Val, dat$Mean_Val_Val, dat$SD_Val_Val,
      SIMPLIFY = FALSE
    )
    dat$n_e <- dat$n_Ile_Ile
    dat$mean_e <- dat$Mean_Ile_Ile
    dat$sd_e <- dat$SD_Ile_Ile
    dat$n_c <- vapply(combined, function(x) x$n, numeric(1))
    dat$mean_c <- vapply(combined, function(x) x$mean, numeric(1))
    dat$sd_c <- vapply(combined, function(x) x$sd, numeric(1))
    dat$Experimental <- "Ile/Ile"
    dat$Control <- "Ile/Val + Val/Val"
    overlap_split <- split_overlapping_totals(dat)
    dat <- overlap_split$included
    excluded <- rbind_fill(excluded, overlap_split$excluded)
    dat$Study <- make_study_label(dat)
  }

  list(included = dat, excluded = excluded)
}

prepare_gstp1_recessive <- function(dat) {
  cols <- c(
    "n_Ile_Ile", "Mean_Ile_Ile", "SD_Ile_Ile",
    "n_Ile_Val", "Mean_Ile_Val", "SD_Ile_Val",
    "n_Val_Val", "Mean_Val_Val", "SD_Val_Val"
  )
  missing_cols <- setdiff(cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing GSTP1 recessive-model columns: ", paste(missing_cols, collapse = ", "))
  }
  for (cc in cols) dat[[cc]] <- num_clean(dat[[cc]])

  valid <- complete.cases(dat[, cols]) &
    dat$n_Ile_Ile >= 2 & dat$n_Ile_Val >= 2 & dat$n_Val_Val >= 2 &
    dat$SD_Ile_Ile > 0 & dat$SD_Ile_Val > 0 & dat$SD_Val_Val > 0

  excluded <- dat[!valid, , drop = FALSE]
  excluded$Exclusion_reason <- rep(
    "Incomplete genotype data for recessive-model pooling",
    nrow(excluded)
  )

  dat <- dat[valid, , drop = FALSE]
  if (nrow(dat) > 0) {
    combined <- mapply(
      combine_two_groups,
      dat$n_Ile_Ile, dat$Mean_Ile_Ile, dat$SD_Ile_Ile,
      dat$n_Ile_Val, dat$Mean_Ile_Val, dat$SD_Ile_Val,
      SIMPLIFY = FALSE
    )
    dat$n_e <- vapply(combined, function(x) x$n, numeric(1))
    dat$mean_e <- vapply(combined, function(x) x$mean, numeric(1))
    dat$sd_e <- vapply(combined, function(x) x$sd, numeric(1))
    dat$n_c <- dat$n_Val_Val
    dat$mean_c <- dat$Mean_Val_Val
    dat$sd_c <- dat$SD_Val_Val
    dat$Experimental <- "Ile/Ile + Ile/Val"
    dat$Control <- "Val/Val"
    overlap_split <- split_overlapping_totals(dat)
    dat <- overlap_split$included
    excluded <- rbind_fill(excluded, overlap_split$excluded)
    dat$Study <- make_study_label(dat)
  }

  list(included = dat, excluded = excluded)
}

build_tasks <- function(sheet_name, raw_dat) {
  if (grepl("^(GSTM1|GSTT1)_", sheet_name)) {
    return(list(
      Positive_vs_Null = prepare_two_group(
        raw_dat,
        "n_Positive", "Mean_Positive", "SD_Positive",
        "n_Null", "Mean_Null", "SD_Null",
        "Positive", "Null"
      )
    ))
  }

  if (grepl("^GSTP1_", sheet_name)) {
    return(list(
      Heterozygote_IleIle_vs_IleVal = prepare_two_group(
        raw_dat,
        "n_Ile_Ile", "Mean_Ile_Ile", "SD_Ile_Ile",
        "n_Ile_Val", "Mean_Ile_Val", "SD_Ile_Val",
        "Ile/Ile", "Ile/Val"
      ),
      Homozygote_IleIle_vs_ValVal = prepare_two_group(
        raw_dat,
        "n_Ile_Ile", "Mean_Ile_Ile", "SD_Ile_Ile",
        "n_Val_Val", "Mean_Val_Val", "SD_Val_Val",
        "Ile/Ile", "Val/Val"
      ),
      Dominant_IleIle_vs_ValCarrier = prepare_gstp1_carrier(raw_dat),
      Recessive_IleCarrier_vs_ValVal = prepare_gstp1_recessive(raw_dat)
    ))
  }

  list()
}

# ---------- 5. Meta-analysis and result extraction ----------
make_meta <- function(dat, with_region = FALSE) {
  allow_random <- nrow(dat) >= 2
  args <- list(
    n.e = dat$n_e,
    mean.e = dat$mean_e,
    sd.e = dat$sd_e,
    n.c = dat$n_c,
    mean.c = dat$mean_c,
    sd.c = dat$sd_c,
    studlab = dat$Study,
    data = dat,
    sm = "SMD",
    method.smd = "Hedges",
    method.tau = tau_method,
    method.random.ci = if (allow_random) random_ci_method else "classic",
    common = TRUE,
    random = allow_random,
    prediction = FALSE
  )
  if (with_region) {
    args$subgroup <- dat$Region
    args$subgroup.name <- "Region"
    args$test.subgroup.common <- TRUE
    args$test.subgroup.random <- TRUE
  }
  do.call(metacont, args)
}

extract_result <- function(m) {
  k <- get_first(m, c("k"))
  ph <- get_first(m, c("pval.Q"))
  i2 <- get_first(m, c("I2"))
  tau2 <- get_first(m, c("tau2"))

  te_common <- get_first(m, c("TE.common", "TE.fixed"))
  lo_common <- get_first(m, c("lower.common", "lower.fixed"))
  up_common <- get_first(m, c("upper.common", "upper.fixed"))
  p_common <- get_first(m, c("pval.common", "pval.fixed"))

  te_random <- get_first(m, c("TE.random"))
  lo_random <- get_first(m, c("lower.random"))
  up_random <- get_first(m, c("upper.random"))
  p_random <- get_first(m, c("pval.random"))

  use_random <- !is.na(ph) && ph < heterogeneity_cutoff && k >= 2
  if (use_random) {
    selected_ci <- fmt_ci(te_random, lo_random, up_random)
    selected_p <- fmt_p(p_random)
    model_used <- paste0("Random (", tau_method, "+", random_ci_method, ")")
  } else {
    selected_ci <- fmt_ci(te_common, lo_common, up_common)
    selected_p <- fmt_p(p_common)
    model_used <- "Common-effect"
  }

  data.frame(
    K = k,
    SMD_common_95CI = fmt_ci(te_common, lo_common, up_common),
    P_common = fmt_p(p_common),
    SMD_random_95CI = fmt_ci(te_random, lo_random, up_random),
    P_random = fmt_p(p_random),
    Ph = fmt_p(ph),
    I2 = fmt_i2(i2),
    Tau2 = fmt_num(tau2, 4),
    Model_used = model_used,
    SMD_95CI = selected_ci,
    P = selected_p,
    stringsAsFactors = FALSE
  )
}

extract_region_results <- function(dat, m_subgroup) {
  regions <- unique(dat$Region)
  out <- do.call(
    rbind,
    lapply(regions, function(region_name) {
      region_dat <- dat[dat$Region == region_name, , drop = FALSE]
      region_meta <- make_meta(region_dat, with_region = FALSE)
      data.frame(
        Region = region_name,
        extract_result(region_meta),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })
  )

  p_between_common <- get_first(
    m_subgroup,
    c("pval.Q.b.common", "pval.Q.b.fixed", "pval.Q.b")
  )
  p_between_random <- get_first(m_subgroup, c("pval.Q.b.random"))
  out$P_between_regions_common <- fmt_p(p_between_common)
  out$P_between_regions_random <- fmt_p(p_between_random)
  out
}

get_bias_p <- function(m, method_bias) {
  b <- tryCatch(
    metabias(m, method.bias = method_bias, correct = TRUE, k.min = 3),
    error = function(e) NULL
  )
  if (is.null(b)) return(NA_real_)
  get_first(b, c("p.value", "pval", "pval.bias"))
}

extract_publication_bias <- function(m) {
  k <- get_first(m, c("k"))
  data.frame(
    K = k,
    Begg_P = if (k >= 3) fmt_p4(get_bias_p(m, "Begg")) else "",
    Egger_P = if (k >= 3) fmt_p4(get_bias_p(m, "Egger")) else "",
    Note = if (k < 10) "Interpret cautiously because fewer than 10 datasets" else "",
    stringsAsFactors = FALSE
  )
}

extract_sensitivity <- function(m, pooled_model) {
  inf <- tryCatch(metainf(m, pooled = pooled_model), error = function(e) NULL)
  if (is.null(inf)) return(data.frame())
  data.frame(
    Removed_dataset = inf$studlab,
    SMD_95CI = sprintf("%.2f (%.2f, %.2f)", inf$TE, inf$lower, inf$upper),
    P = vapply(inf$pval, fmt_p, character(1)),
    Model = pooled_model,
    stringsAsFactors = FALSE
  )
}

# ---------- 6. Plotting ----------
forest_style <- function(m, experimental_label, control_label,
                         year_sort, subgroup_plot = FALSE) {
  forest(
    m,
    common = TRUE,
    random = TRUE,
    prediction = FALSE,
    overall = TRUE,
    overall.hetstat = TRUE,
    sortvar = year_sort,
    leftcols = c("studlab", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c"),
    leftlabs = c("Dataset", experimental_label, "Mean", "SD", control_label, "Mean", "SD"),
    rightcols = c("effect", "ci", "w.common", "w.random"),
    rightlabs = c("SMD", "95% CI", "Weight (common)", "Weight (random)"),
    smlab = "Hedges' g (SMD)",
    xlab = paste0(
      "Hedges' g: ", experimental_label, " minus ", control_label,
      "  (positive = higher in ", experimental_label, ")"
    ),
    digits = 2,
    digits.I2 = 1,
    col.square = "#2C7FB8",
    col.square.lines = "#2C7FB8",
    col.diamond = "#D95F0E",
    col.diamond.lines = "#D95F0E",
    col.predict = "#636363",
    print.subgroup.name = subgroup_plot,
    test.subgroup = subgroup_plot
  )
}

plot_analysis <- function(dat, m, m_region, fig_dir, prefix) {
  experimental_label <- unique(dat$Experimental)[1]
  control_label <- unique(dat$Control)[1]
  forest_height <- max(6, 3 + 0.55 * nrow(dat))
  region_height <- max(8, 5 + 0.75 * nrow(dat))
  year_sort <- suppressWarnings(as.numeric(dat$Year))
  year_sort[is.na(year_sort)] <- Inf

  safe_pdf(file.path(fig_dir, paste0(prefix, "_overall_forest.PDF")), {
    forest_style(
      m, experimental_label, control_label,
      year_sort = year_sort, subgroup_plot = FALSE
    )
  }, height = forest_height, width = forest_pdf_width)

  safe_pdf(file.path(fig_dir, paste0(prefix, "_region_subgroup_forest.PDF")), {
    forest_style(
      m_region, experimental_label, control_label,
      year_sort = year_sort, subgroup_plot = TRUE
    )
  }, height = region_height, width = forest_pdf_width)

  if (nrow(dat) >= 3) {
    safe_pdf(file.path(fig_dir, paste0(prefix, "_sensitivity_common.PDF")), {
      forest(metainf(m, pooled = "common"), smlab = "Leave-one-out: common effect")
    }, height = forest_height, width = forest_pdf_width)

    safe_pdf(file.path(fig_dir, paste0(prefix, "_sensitivity_random.PDF")), {
      forest(metainf(m, pooled = "random"), smlab = "Leave-one-out: random effects")
    }, height = forest_height, width = forest_pdf_width)
  }

  if (nrow(dat) >= 2) {
    safe_pdf(file.path(fig_dir, paste0(prefix, "_funnel.PDF")), {
      funnel(m, pch = 21, xlab = "Hedges' g", studlab = TRUE)
    }, height = 6, width = 6)

    safe_pdf(file.path(fig_dir, paste0(prefix, "_galbraith_name.PDF")), {
      radial(m, text = dat$Study, level = 0.95)
    }, height = 7, width = 7)

    safe_pdf(file.path(fig_dir, paste0(prefix, "_galbraith_point.PDF")), {
      radial(m, pch = 19, level = 0.95)
    }, height = 7, width = 7)
  }
}

# ---------- 7. Run all worksheets and comparisons ----------
available_sheets <- excel_sheets(data_file)
sheets_to_run <- intersect(expected_sheets, available_sheets)
missing_sheets <- setdiff(expected_sheets, available_sheets)

if (length(sheets_to_run) == 0) {
  stop("None of the expected worksheets were found: ", paste(expected_sheets, collapse = ", "))
}
if (length(missing_sheets) > 0) {
  warning("These worksheets were not found and will be skipped: ", paste(missing_sheets, collapse = ", "))
}

master_overall <- list()
master_region <- list()

for (sheet_name in sheets_to_run) {
  message("Reading sheet: ", sheet_name)
  raw_dat <- read_sheet_clean(sheet_name)
  tasks <- build_tasks(sheet_name, raw_dat)

  for (comparison_name in names(tasks)) {
    task <- tasks[[comparison_name]]
    dat <- task$included

    analysis_dir <- file.path(output_root, sheet_name, comparison_name)
    fig_dir <- file.path(analysis_dir, "figures")
    res_dir <- file.path(analysis_dir, "results")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

    write_csv_safe(dat, file.path(res_dir, "included_datasets.csv"))
    write_csv_safe(task$excluded, file.path(res_dir, "excluded_rows.csv"))

    if (nrow(dat) < 2) {
      warning(sheet_name, " / ", comparison_name, ": fewer than 2 valid datasets; skipped.")
      next
    }

    m <- make_meta(dat, with_region = FALSE)
    m_region <- make_meta(dat, with_region = TRUE)

    overall <- data.frame(
      Sheet = sheet_name,
      Comparison = comparison_name,
      Direction = paste0(unique(dat$Experimental)[1], " minus ", unique(dat$Control)[1]),
      extract_result(m),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    region <- data.frame(
      Sheet = sheet_name,
      Comparison = comparison_name,
      extract_region_results(dat, m_region),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    bias <- data.frame(
      Sheet = sheet_name,
      Comparison = comparison_name,
      extract_publication_bias(m),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    sens_common <- extract_sensitivity(m, "common")
    sens_random <- extract_sensitivity(m, "random")
    sensitivity <- rbind(sens_common, sens_random)

    write_csv_safe(overall, file.path(res_dir, "Table1_overall_results.csv"))
    write_csv_safe(region, file.path(res_dir, "Table2_region_subgroup_results.csv"))
    write_csv_safe(bias, file.path(res_dir, "Table3_publication_bias.csv"))
    write_csv_safe(sensitivity, file.path(res_dir, "Table4_leave_one_out.csv"))

    prefix <- paste(sheet_name, comparison_name, sep = "_")
    plot_analysis(dat, m, m_region, fig_dir, prefix)

    master_overall[[prefix]] <- overall
    master_region[[prefix]] <- region
    message("Finished: ", sheet_name, " / ", comparison_name, "; datasets = ", nrow(dat))
  }
}

if (length(master_overall) > 0) {
  write_csv_safe(
    do.call(rbind, master_overall),
    file.path(output_root, "ALL_overall_results.csv")
  )
}
if (length(master_region) > 0) {
  write_csv_safe(
    do.call(rbind, master_region),
    file.path(output_root, "ALL_region_subgroup_results.csv")
  )
}

message("All analyses completed. Results saved in: ", output_root)
message("Effect measure: Hedges' g (SMD).")
message("Region subgroup results are in ALL_region_subgroup_results.csv and each analysis folder.")
