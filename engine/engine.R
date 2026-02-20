run_voxseq_pipeline <- function(
  uploaded_files_df,
  expt_override = NULL,
  min_calls = 200,
  min_IGI = 0.125,
  min_IBI = 0.225
) {
  stopifnot(is.data.frame(uploaded_files_df))
  stopifnot(all(c("datapath", "name") %in% names(uploaded_files_df)))

  # Create a temp working directory for this run
  run_dir  <- file.path(tempdir(), paste0("voxseq_run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  data_dir <- file.path(run_dir, "data")
  plot_dir <- file.path(run_dir, paste0("Plots_", format(Sys.Date(), "%Y-%m-%d")))
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Copy uploaded Excel files into data_dir
  file.copy(
    from = uploaded_files_df$datapath,
    to   = file.path(data_dir, uploaded_files_df$name),
    overwrite = TRUE
  )

  # Run the pipeline in an isolated environment (so it doesn't pollute Shiny's global env)
  e <- new.env(parent = globalenv())

  # User input variables
  e$data_dir   <- data_dir
  e$plot_dir   <- plot_dir
  e$min_calls  <- min_calls
  e$min_IGI    <- min_IGI
  e$min_IBI    <- min_IBI

  # If you want to override expt (optional)
  if (!is.null(expt_override) && nzchar(expt_override)) {
    e$expt_override <- expt_override
  }

  # Set working directory to engine folder so relative paths work
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(normalizePath(file.path(old_wd, "engine")))

  ## Dependencies
  sys.source("Utils.R", envir = e)
  sys.source("extractUSV.R", envir = e)

  ## Analysis Script
  sys.source("Analysis.R", envir = e)

  # IMPORTANT: Analysis.R should NOT rm(list=ls()) and should not file.choose()
  
  # Return paths for Shiny
  list(
    run_dir = run_dir,
    data_dir = data_dir,
    plot_dir = plot_dir,
    plot_files = list.files(plot_dir, full.names = TRUE)
  )
}

zip_plots <- function(plot_dir, zip_path) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(plot_dir)

  files <- list.files(".", recursive = TRUE)
  if (length(files) == 0) stop("No plot files found to zip.")

  utils::zip(zipfile = zip_path, files = files)
  zip_path
}
