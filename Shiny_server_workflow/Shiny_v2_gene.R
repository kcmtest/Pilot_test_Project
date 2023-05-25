library(shiny)
library(Seurat)
library(BiocParallel)
library(reticulate)
library(Matrix)
library(dplyr)
library(DT)
options(shiny.maxRequestSize = 4e+9)

use_condaenv("Pilot_python_3_10")
sc <- import("scanpy")

# Define UI
ui <- fluidPage(
  titlePanel("Seurat Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("h5adFile", "Upload h5ad file"),
      actionButton("run_button", "Run Analysis"),
      selectInput("gene_input", "Select Gene", choices = NULL),
      downloadButton("export_button_tsv", "Export Data (TSV)"),
      downloadButton("export_button_md", "Export Data (Markdown)")
    ),
    mainPanel(
      plotOutput("feature_plot"),
      DT::dataTableOutput("fetch_data_table")
    )
  )
)

# Define server
server <- function(input, output, session) {
  # Reactive value to store the Seurat object
  seurat_object <- reactiveVal(NULL)
  
  # Function to perform the analysis
  performAnalysis <- function(h5ad_file) {
    adata <- sc$read_h5ad(h5ad_file)
    
    counts <- t(adata$layers["counts"])
    colnames(counts) <- adata$obs_names$to_list()
    rownames(counts) <- adata$var_names$to_list()
    counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
    
    data <- t(adata$layers["winsorized"])
    colnames(data) <- adata$obs_names$to_list()
    rownames(data) <- adata$var_names$to_list()
    data <- Matrix::Matrix(as.matrix(data), sparse = TRUE)
    
    seurat <- CreateSeuratObject(counts)
    seurat <- SetAssayData(seurat, "data", data)
    seurat <- AddMetaData(seurat, adata$obs)
    
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
    seurat <- NormalizeData(seurat)
    seurat <- FindVariableFeatures(seurat, nfeatures = 500)
    seurat <- ScaleData(seurat)
    seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))
    seurat <- RunPCA(seurat, npcs = 50)
    seurat <- RunUMAP(seurat, dims = 1:20)
    seurat <- FindNeighbors(seurat, dims = 1:20)
    seurat <- FindClusters(seurat, resolution = 1)
    cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
    top_gene <- cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
    
    return(seurat)
  }
  
  # Perform analysis when the Run Analysis button is clicked
  observeEvent(input$run_button, {
    req(input$h5adFile)
    h5ad_file <- input$h5adFile$datapath
    seurat <- performAnalysis(h5ad_file)
    seurat_object(seurat)
    
    # Update gene choices for feature plot
    gene_choices <- rownames(seurat@assays$RNA@counts)
    updateSelectInput(session, "gene_input", choices = gene_choices)
  })
  
  # Render the feature plot
  output$feature_plot <- renderPlot({
    seurat <- seurat_object()
    if (!is.null(seurat)) {
      gene <- input$gene_input
      FeaturePlot(seurat, features = gene, label = TRUE)
    }
  })
  
  # Function to fetch data based on specified gene input
  fetchData <- function(seurat_object, gene_input) {
    seurat <- seurat_object()
    if (!is.null(seurat)) {
      FetchData(object = seurat, vars = c("seurat_clusters","louvain",gene_input))
    }
  }
  
  # Render the data table
  output$fetch_data_table <- DT::renderDataTable({
    seurat <- seurat_object()
    gene_input <- input$gene_input
    if (!is.null(seurat) && !is.null(gene_input)) {
      data_table <- fetchData(seurat_object, gene_input)
      DT::datatable(data_table)
    }
  })
  
  # Export the data table as TSV file
  output$export_button_tsv <- downloadHandler(
    filename = function() {
      gene_input <- input$gene_input
      paste0("data_table_", gene_input, "_", format(Sys.time(), "%Y%m%d%H%M%S"), ".tsv")
    },
    content = function(file) {
      seurat <- seurat_object()
      gene_input <- input$gene_input
      if (!is.null(seurat) && !is.null(gene_input)) {
        data_table <- fetchData(seurat_object, gene_input)
        write.table(data_table, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      }
    }
  )
  
  # Export the data table as Markdown file
  output$export_button_md <- downloadHandler(
    filename = function() {
      gene_input <- input$gene_input
      paste0("data_table_", gene_input, "_", format(Sys.time(), "%Y%m%d%H%M%S"), ".md")
    },
    content = function(file) {
      seurat <- seurat_object()
      gene_input <- input$gene_input
      if (!is.null(seurat) && !is.null(gene_input)) {
        data_table <- fetchData(seurat_object, gene_input)
        writeLines(c(
          paste("# Data Table"),
          paste("louvain | seurat_clusters |", gene_input, "|"),
          paste("--- | --- | --- |"),
          paste(apply(data_table, 1, paste, collapse = " | "))
        ), file)
      }
    }
  )
  
  # End the Shiny session and shut down the server when the browser is closed
  session$onSessionEnded(stopApp)
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
