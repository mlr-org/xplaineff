load_global_runtime_data = function(indir, model_types = c("rf", "toy"), include_mlr3 = FALSE) {
  if (isTRUE(include_mlr3) && !"mlr3_rf" %in% model_types) {
    model_types = c(model_types, "mlr3_rf")
  }

  csv_pattern = "^global_r_runtime_(mlr3_rf|rf|toy)(_(vs_N|vs_D|vs_res))?\\.csv$"
  csv_paths = if (dir.exists(indir)) {
    list.files(indir, pattern = csv_pattern, full.names = TRUE)
  } else {
    character()
  }

  if (!length(csv_paths)) {
    return(data.table())
  }

  read_one_csv = function(path) {
    x = fread(path)
    source_name = basename(path)
    x[, source_file := source_name]
    x[, source_rank := if (grepl("_(vs_N|vs_D|vs_res)\\.csv$", source_name)) 1L else 2L]
    x
  }

  dt = rbindlist(lapply(csv_paths, read_one_csv), use.names = TRUE, fill = TRUE)
  if (!("status" %in% names(dt))) {
    dt[, status := "ok"]
  }
  if (!("sub_experiment" %in% names(dt))) {
    dt[, sub_experiment := "vs_N"]
  }
  dt[is.na(status) | status == "", status := "ok"]
  dt[is.na(sub_experiment) | sub_experiment == "", sub_experiment := "vs_N"]
  dt[package == "gadget", package := "xplaineff"]
  dt[package == "ingredients", package := "DALEX/ingredients"]
  dt = dt[module == "global_r" & model_type %in% model_types]

  integer_key_cols = c("N", "D", "n_grid", "n_intervals", "repetition")
  for (col in intersect(integer_key_cols, names(dt))) {
    dt[, (col) := suppressWarnings(as.integer(fifelse(
      is.na(get(col)) | get(col) == "",
      NA_character_,
      as.character(get(col))
    )))]
  }

  key_cols = c(
    "module", "package", "impl", "method", "model_type", "sub_experiment",
    "N", "D", "n_grid", "n_intervals", "repetition"
  )
  missing_key_cols = setdiff(key_cols, names(dt))
  if (length(missing_key_cols)) {
    stop(sprintf("Raw benchmark CSVs are missing required column(s): %s", paste(missing_key_cols, collapse = ", ")))
  }

  source_report = dt[, .(
    n_sources = uniqueN(source_file),
    source_files = paste(sort(unique(source_file)), collapse = ";")
  ), by = key_cols][n_sources > 1L]
  if (nrow(source_report)) {
    fwrite(source_report, file.path(indir, "global_r_duplicate_sources.csv"))
    message(
      "Found ", nrow(source_report), " duplicated benchmark row key(s); ",
      "keeping per-sub CSV rows before merged rows. See global_r_duplicate_sources.csv."
    )
  }

  setorderv(dt, c("source_rank", "source_file"))
  n_duplicate_rows = sum(duplicated(dt, by = key_cols))
  if (n_duplicate_rows > 0L) {
    dt = unique(dt, by = key_cols)
    message("Dropped ", n_duplicate_rows, " duplicated benchmark row(s) after source prioritization.")
  }

  setorderv(dt, key_cols)
  dt
}
