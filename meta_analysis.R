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
 
required_packages <- c("meta", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Please install: ", paste(missing_packages, collapse = ", "))
}

library(meta)
library(readxl)

# ---------- 1. Paths and settings ----------
data_file <- paste0(
  "/Users/roseguo/Desktop/科研/课题6_Benzene/step7_datatable/",
  "meta_analysis_data_2.xlsx"
)

if (!file.exists(data_file)) {
  if (!interactive()) stop("data_file not found: ", data_file)
  message("Excel file not found. Please choose manually.")
  data_file <- file.choose()
}

base_dir <- dirname(dirname(data_file))
output_root <- file.path(base_dir, "gxy_new result")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

prefer_smoking_strata <- TRUE
heterogeneity_cutoff <- 0.10
tau_method <- "REML"
random_ci_method <- "HK"
hk_min_k <- 2
sample_size_cutoff <- 100
forest_pdf_width <- 20

expected_sheets <- c(
  "GSTM1_SPMA", "GSTM1_ttMA",
  "GSTT1_SPMA", "GSTT1_ttMA"
)

# ---------- 2. Helpers ----------
fill_down <- function(x) {
  x <- as.character(x)
  for (i in seq_along(x)) {
    if (i > 1 && (is.na(x[i]) || trimws(x[i]) == "")) x[i] <- x[i - 1]
  }
  x
}

num_clean <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))

clean_column_names <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

format_year <- function(x) sub("\\.0$", "", trimws(as.character(x)))

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
  force(expr)
}

make_study_label <- function(dat) {
  author <- trimws(as.character(dat$Author))
  has_suffix <- grepl("^.+\\s*\\(.+\\)$", author)
  author_base <- sub("\\s*\\((.*)\\)$", "", author)
  author_suffix <- sub("^.+\\s*\\((.*)\\)$", "\\1", author)
  base <- ifelse(
    has_suffix,
    paste0(author_base, " ", dat$Year, "(", author_suffix, ")"),
    paste(author, dat$Year)
  )
  stratum <- trimws(as.character(dat$Population_Stratum))
  show_stratum <- !is.na(stratum) & stratum != "" &
    !tolower(stratum) %in% c("all subjects", "overall", "total")
  ifelse(show_stratum & !has_suffix, paste0(base, "(", stratum, ")"), base)
}

# ---------- 3. Read / clean ----------
read_sheet_clean <- function(sheet_name) {
  dat <- as.data.frame(read_excel(data_file, sheet = sheet_name))
  raw_names <- names(dat)
  unit_index <- match(c("单位", "Unit", "unit"), raw_names, nomatch = 0L)
  unit_index <- unit_index[unit_index > 0L][1]
  unit_values <- if (length(unit_index) == 0L || is.na(unit_index)) {
    rep(NA_character_, nrow(dat))
  } else {
    dat[[unit_index]]
  }
  names(dat) <- clean_column_names(names(dat))
  blank_names <- which(is.na(names(dat)) | names(dat) == "")
  if (length(blank_names) > 0) {
    names(dat)[blank_names] <- paste0("unnamed_", blank_names)
  }
  dat$Unit <- unit_values

  required_metadata <- c("Author", "Year", "Country", "Region")
  missing_metadata <- setdiff(required_metadata, names(dat))
  if (length(missing_metadata) > 0) {
    stop(sheet_name, " is missing: ", paste(missing_metadata, collapse = ", "))
  }

  if (!"Population_Stratum" %in% names(dat)) dat$Population_Stratum <- "All subjects"

  fill_cols <- intersect(
    c("Author", "Year", "Country", "Region", "Ethnicity",
      "Detecting_method", "Source", "Population_Stratum", "Unit"),
    names(dat)
  )
  for (cc in fill_cols) dat[[cc]] <- fill_down(dat[[cc]])

  dat$Author <- trimws(as.character(dat$Author))
  dat$Year <- format_year(dat$Year)
  dat$Region <- trimws(as.character(dat$Region))
  dat$Region[is.na(dat$Region) | dat$Region == ""] <- "Not reported"
  dat$Population_Stratum <- trimws(as.character(dat$Population_Stratum))
  dat$Unit <- trimws(as.character(dat$Unit))
  dat$Unit[is.na(dat$Unit) | dat$Unit == ""] <- "Not reported"

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
  # 修复：用 rep(..., nrow(excluded)) 保持行数一致
  excluded$Exclusion_reason <- rep(
    "Overlapping overall row removed because smoking strata were used",
    nrow(excluded)
  )
  list(included = included, excluded = excluded)
}

# ---------- 4. Prepare comparisons ----------
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
  # 修复：0 行时也能安全赋值
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
  list()
}

# ---------- 5. Meta-analysis ----------
make_meta <- function(dat, with_subgroup = NULL, random_ci_override = NULL) {
  allow_random <- nrow(dat) >= 2
  ci_method <- if (!is.null(random_ci_override)) {
    random_ci_override
  } else if (allow_random && nrow(dat) >= hk_min_k) {
    random_ci_method
  } else {
    "classic"
  }

  args <- list(
    n.e = dat$n_e, mean.e = dat$mean_e, sd.e = dat$sd_e,
    n.c = dat$n_c, mean.c = dat$mean_c, sd.c = dat$sd_c,
    studlab = dat$Study, data = dat,
    sm = "SMD", method.smd = "Hedges",
    method.tau = tau_method,
    method.random.ci = ci_method,
    common = TRUE,
    random = allow_random,
    prediction = FALSE
  )
  if (!is.null(with_subgroup)) {
    args$subgroup <- dat[[with_subgroup]]
    args$subgroup.name <- with_subgroup
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

  if (k >= 2 && is.na(ph)) {
    warning("Heterogeneity p is NA; defaulting to common-effect.")
  }

  use_random <- !is.na(ph) && ph <= heterogeneity_cutoff && k >= 2
  if (use_random) {
    selected_ci <- fmt_ci(te_random, lo_random, up_random)
    selected_p <- fmt_p(p_random)
    model_used <- paste0("Random (", tau_method, "+",
                         ifelse(k >= hk_min_k, random_ci_method, "classic"), ")")
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
    Ph = fmt_p(ph), I2 = fmt_i2(i2), Tau2 = fmt_num(tau2, 4),
    Model_used = model_used,
    SMD_95CI = selected_ci, P = selected_p,
    stringsAsFactors = FALSE
  )
}

extract_subgroup_results <- function(dat, subgroup_var, m_subgroup) {
  groups <- unique(dat[[subgroup_var]])
  out <- do.call(rbind, lapply(groups, function(g) {
    sub_dat <- dat[dat[[subgroup_var]] == g, , drop = FALSE]
    if (nrow(sub_dat) < 2) {
      # 修复：单研究分支也带上 P_between_* 两列，避免 rbind 列错位
      return(data.frame(
        Group = as.character(g), K = nrow(sub_dat),
        SMD_common_95CI = "", P_common = "",
        SMD_random_95CI = "", P_random = "",
        Ph = "", I2 = "", Tau2 = "",
        Model_used = "Single study (not pooled)",
        SMD_95CI = "", P = "",
        P_between_common = "", P_between_random = "",
        check.names = FALSE, stringsAsFactors = FALSE
      ))
    }
            ci_method <- if (nrow(sub_dat) >= hk_min_k) random_ci_method else "classic"
    m_g <- make_meta(
      sub_dat,
      with_subgroup = NULL,
      random_ci_override = ci_method
    )
    data.frame(
      Group = as.character(g),
      extract_result(m_g),
      P_between_common = "", P_between_random = "",
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }))
  names(out)[names(out) == "Group"] <- subgroup_var

  p_between_common <- get_first(m_subgroup,
    c("pval.Q.b.common", "pval.Q.b.fixed", "pval.Q.b"))
  p_between_random <- get_first(m_subgroup, c("pval.Q.b.random"))
  out$P_between_common <- fmt_p(p_between_common)
  out$P_between_random <- fmt_p(p_between_random)
  out
}

extract_region_results <- function(dat, m_subgroup) {
  extract_subgroup_results(dat, "Region", m_subgroup)
}

get_bias_p <- function(m, method_bias) {
  b <- tryCatch(
    metabias(m, method.bias = method_bias,
             correct = (method_bias == "Begg"), k.min = 3),
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

# 修复：失败时返回带列名的 0 行 data.frame，避免 rbind/write 出问题
empty_sensitivity <- function() {
  data.frame(
    Removed_dataset = character(0),
    SMD_95CI = character(0),
    P = character(0),
    Model = character(0),
    stringsAsFactors = FALSE
  )
}

extract_sensitivity <- function(m, pooled_model) {
  inf <- tryCatch(metainf(m, pooled = pooled_model), error = function(e) NULL)
  if (is.null(inf) || is.null(inf$studlab) || length(inf$studlab) == 0) {
    return(empty_sensitivity())
  }
  data.frame(
    Removed_dataset = inf$studlab,
    SMD_95CI = sprintf("%.2f (%.2f, %.2f)", inf$TE, inf$lower, inf$upper),
    P = vapply(inf$pval, fmt_p, character(1)),
    Model = pooled_model,
    stringsAsFactors = FALSE
  )
}

# ---------- 6. Plot ----------
forest_style <- function(m, experimental_label, control_label,
                         year_sort, subgroup_plot = FALSE) {
  forest(
    m,
    common = TRUE, random = TRUE, prediction = FALSE,
    overall = TRUE, overall.hetstat = TRUE,
    sortvar = year_sort,
    leftcols = c("studlab", "n.e", "mean.e", "sd.e", "n.c", "mean.c", "sd.c"),
    leftlabs = c("Dataset", experimental_label, "Mean", "SD",
                  control_label, "Mean", "SD"),
    rightcols = c("effect", "ci", "w.common", "w.random"),
    rightlabs = c("SMD", "95% CI", "Weight (common)", "Weight (random)"),
    smlab = "Hedges' g (SMD)",
    xlab = paste0("Hedges' g: ", experimental_label, " minus ", control_label,
                  "  (positive = higher in ", experimental_label, ")"),
    digits = 2, digits.I2 = 1,
    col.square = "#2C7FB8", col.square.lines = "#2C7FB8",
    col.diamond = "#D95F0E", col.diamond.lines = "#D95F0E",
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
    forest_style(m, experimental_label, control_label, year_sort, FALSE)
  }, height = forest_height, width = forest_pdf_width)

  safe_pdf(file.path(fig_dir, paste0(prefix, "_region_subgroup_forest.PDF")), {
    forest_style(m_region, experimental_label, control_label, year_sort, TRUE)
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

plot_subgroup <- function(dat, m_subgroup, subgroup_var, fig_dir, prefix,
                          experimental_label, control_label, year_sort,
                          height_mult = 0.75) {
  height <- max(8, 5 + height_mult * nrow(dat))
  safe_pdf(file.path(fig_dir, paste0(prefix, "_", subgroup_var, "_subgroup_forest.PDF")), {
    forest_style(m_subgroup, experimental_label, control_label,
                  year_sort, subgroup_plot = TRUE)
  }, height = height, width = forest_pdf_width)
}

# ---------- 7. Main loop ----------
available_sheets <- excel_sheets(data_file)
sheets_to_run <- intersect(expected_sheets, available_sheets)
missing_sheets <- setdiff(expected_sheets, available_sheets)

if (length(sheets_to_run) == 0) {
  stop("None of the expected worksheets were found.")
}
if (length(missing_sheets) > 0) {
  warning("Missing worksheets skipped: ", paste(missing_sheets, collapse = ", "))
}

master_overall <- list()
master_region <- list()
master_unit <- list()
master_size <- list()

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
      warning(sheet_name, " / ", comparison_name, ": <2 valid datasets; skipped.")
      next
    }

    dat$Total_N <- dat$n_e + dat$n_c
    dat$SampleSizeGroup <- ifelse(
      dat$Total_N < sample_size_cutoff,
      paste0("<", sample_size_cutoff),
      paste0(">=", sample_size_cutoff)
    )
    dat$SampleSizeGroup <- factor(
      dat$SampleSizeGroup,
      levels = c(paste0("<", sample_size_cutoff), paste0(">=", sample_size_cutoff))
    )

    experimental_label <- unique(dat$Experimental)[1]
    control_label <- unique(dat$Control)[1]
    year_sort <- suppressWarnings(as.numeric(dat$Year))
    year_sort[is.na(year_sort)] <- Inf

    # A. 总体
    m <- make_meta(dat, with_subgroup = NULL)
    m_region <- make_meta(dat, with_subgroup = "Region")

    overall <- data.frame(
      Sheet = sheet_name, Comparison = comparison_name,
      Direction = paste0(experimental_label, " minus ", control_label),
      extract_result(m), check.names = FALSE, stringsAsFactors = FALSE
    )
    region <- data.frame(
      Sheet = sheet_name, Comparison = comparison_name,
      extract_region_results(dat, m_region),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    bias <- data.frame(
      Sheet = sheet_name, Comparison = comparison_name,
      extract_publication_bias(m),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    sensitivity <- rbind(
      extract_sensitivity(m, "common"),
      extract_sensitivity(m, "random")
    )

    write_csv_safe(overall, file.path(res_dir, "Table1_overall_results.csv"))
    write_csv_safe(region, file.path(res_dir, "Table2_region_subgroup_results.csv"))
    write_csv_safe(bias, file.path(res_dir, "Table3_publication_bias.csv"))
    write_csv_safe(sensitivity, file.path(res_dir, "Table4_leave_one_out.csv"))

    prefix <- paste(sheet_name, comparison_name, sep = "_")
    plot_analysis(dat, m, m_region, fig_dir, prefix)

    master_overall[[prefix]] <- overall
    master_region[[prefix]] <- region

    # B. Unit 亚组
    dat_unit <- dat[!is.na(dat$Unit) & trimws(as.character(dat$Unit)) != "", ,
                    drop = FALSE]
    if (nrow(dat_unit) >= 2 && length(unique(dat_unit$Unit)) >= 2) {
      m_unit <- make_meta(dat_unit, with_subgroup = "Unit")
      unit <- data.frame(
        Sheet = sheet_name, Comparison = comparison_name,
        extract_subgroup_results(dat_unit, "Unit", m_unit),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      write_csv_safe(unit, file.path(res_dir, "Table5_unit_subgroup_results.csv"))
      year_sort_unit <- suppressWarnings(as.numeric(dat_unit$Year))
      year_sort_unit[is.na(year_sort_unit)] <- Inf
      plot_subgroup(dat_unit, m_unit, "Unit", fig_dir, prefix,
                    experimental_label, control_label, year_sort_unit)
      master_unit[[prefix]] <- unit
    } else {
      message(sheet_name, " / ", comparison_name,
              ": Unit subgroup skipped (fewer than 2 units or 1 level).")
    }

    # C. Sample size 亚组
    if (nrow(dat) >= 2 && length(unique(dat$SampleSizeGroup)) >= 2) {
      m_size <- make_meta(dat, with_subgroup = "SampleSizeGroup")
      size <- data.frame(
        Sheet = sheet_name, Comparison = comparison_name,
        extract_subgroup_results(dat, "SampleSizeGroup", m_size),
        check.names = FALSE, stringsAsFactors = FALSE
      )
      write_csv_safe(size, file.path(res_dir, "Table6_sample_size_subgroup_results.csv"))
      plot_subgroup(dat, m_size, "SampleSizeGroup", fig_dir, prefix,
                    experimental_label, control_label, year_sort)
      master_size[[prefix]] <- size
    } else {
      message(sheet_name, " / ", comparison_name,
              ": Sample size subgroup skipped (only one size level).")
    }

    message("Finished: ", sheet_name, " / ", comparison_name,
            "; datasets = ", nrow(dat))
  }
}

if (length(master_overall) > 0) {
  write_csv_safe(do.call(rbind, master_overall),
                 file.path(output_root, "ALL_overall_results.csv"))
}
if (length(master_region) > 0) {
  write_csv_safe(do.call(rbind, master_region),
                 file.path(output_root, "ALL_region_subgroup_results.csv"))
}
if (length(master_unit) > 0) {
  write_csv_safe(do.call(rbind, master_unit),
                 file.path(output_root, "ALL_unit_subgroup_results.csv"))
}
if (length(master_size) > 0) {
  write_csv_safe(do.call(rbind, master_size),
                 file.path(output_root, "ALL_sample_size_subgroup_results.csv"))
}

message("All analyses completed. Results saved in: ", output_root)
