library(shiny)

source("engine/engine.R")

ui <- fluidPage(
   titlePanel("VoxSeq GUI"),

  sidebarLayout(
  
     sidebarPanel(
      h4("Welcome to the VoxSeq GUI!"),
      p("Upload your USV Excel files and set the parameters for analysis. Click 'Run analysis' to process the data and generate plots. You can preview some of the plots to the right and download all of them as a ZIP file."),

      helpText("Note: The analysis may take some time depending on the number and size of the uploaded files. Please be patient!"),
     
      fileInput(
        inputId = "usv_files",
        label   = "Upload USV Excel files (.xlsx)",
        multiple = TRUE,
        accept = c(".xls", ".xlsx"),
        width = "100%"
      ),
      fluidRow(
        column(
          12,
          numericInput("min_calls", "Minimum calls", value = 200, min = 0, step = 1),
          numericInput("min_IGI",   "Min IGI (seconds)", value = 0.125, min = 0, step = 0.001),
          numericInput("min_IBI",   "Min IBI (seconds)", value = 0.225, min = 0, step = 0.001),
          actionButton("run_btn", "Run analysis", class = "btn-primary"),
          hr(),
          verbatimTextOutput("status"),
          downloadButton("download_zip", "Download plots (ZIP)")
        )
      )
     ),
     mainPanel(
      h4("Plot preview:"),
      uiOutput("plot_gallery"),
      hr(),
      h4("Files you selected:"),
      tableOutput("file_table")  
     ),
    position = "left",
    fluid = TRUE
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    result = NULL,
    zip_path = NULL,
    resource_prefix = NULL
  )

  output$file_table <- renderTable({
    req(input$usv_files)
    data.frame(
      filename = input$usv_files$name,
      size_kb  = round(input$usv_files$size / 1024, 1),
      stringsAsFactors = FALSE
    )
  })

  output$status <- renderPrint({
    if (is.null(rv$result)) {
      cat("Waiting for analysis to run.")
    } else {
      cat("Done.\n")
      cat("Plot folder:\n", rv$result$plot_dir, "\n\n")
      cat("Number of plot files:\n", length(rv$result$plot_files), "\n")
    }
  })

  observeEvent(input$run_btn, {
    req(input$usv_files)

    rv$result <- NULL
    rv$zip_path <- NULL
    rv$resource_prefix <- NULL

    withProgress(message = "Running pipeline...", value = 0, {
      incProgress(0.1, detail = "Preparing run...")

      res <- run_voxseq_pipeline(
        uploaded_files_df = input$usv_files,
        min_calls = input$min_calls,
        min_IGI = input$min_IGI,
        min_IBI = input$min_IBI
      )

      incProgress(0.8, detail = "Preparing plot preview + ZIP...")

      # Make plot directory available to the browser
      prefix <- paste0("plots_", as.integer(Sys.time()))
      addResourcePath(prefix, res$plot_dir)
      rv$resource_prefix <- prefix

      # Create zip now (look into creating on download)
      zip_path <- file.path(res$run_dir, "plots.zip")
      zip_plots(res$plot_dir, zip_path)

      rv$result <- res
      rv$zip_path <- zip_path

      incProgress(1, detail = "Done.")
    })
  })

  output$download_zip <- downloadHandler(
    filename = function() {
      paste0("VoxSeq_Plots_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(rv$zip_path)
      file.copy(rv$zip_path, file, overwrite = TRUE)
    }
  )

  output$plot_gallery <- renderUI({
    req(rv$result, rv$resource_prefix)
    plots <- rv$result$plot_files
    if (length(plots) == 0) return(tags$div("No plots found."))

    # Only show PNGs in preview
    pngs <- plots[grepl("\\.png$", plots, ignore.case = TRUE)]
    if (length(pngs) == 0) return(tags$div("No PNG plots found to preview."))

    show <- head(pngs)

    tags$div(
      lapply(show, function(p) {
        rel <- basename(p)
        src <- paste0("/", rv$resource_prefix, "/", rel)
        tags$div(
          style = "display:inline-block; margin:10px; vertical-align:top;",
          tags$div(style="font-size:12px; width:240px; word-wrap:break-word;", rel),
          tags$img(src = src, style = "width:240px; border:1px solid #ddd;")
        )
      })
    )
  })
}

shinyApp(ui, server)