library(shiny)
library(Seurat)
library(BiocParallel)
library(reticulate)
library(Matrix)
library(dplyr)
library(DT)
options(shiny.maxRequestSize = 4e+9)

# set up conda env for scanpy function
use_condaenv("Pilot_python_3_10")
sc <- import("scanpy")

# UI 
ui <- fluidPage(
  titlePanel("Seurat Data Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("h5adFile", "Upload h5ad file"),
      actionButton("run_button", "Run Analysis"),
      selectInput("gene_input", "Select Gene", choices = NULL),
      downloadButton("export_button", "Export Data")
    ),
    mainPanel(
      plotOutput("feature_plot"),
      tableOutput("fetch_data_table")
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive value to store the Seurat object
  seurat_object <- reactiveVal(NULL)
  
  # Analysis function 
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
  
  # Perform analysis when the Run Analysis button clickeed
  observeEvent(input$run_button, {
    req(input$h5adFile)
    h5ad_file <- input$h5adFile$datapath
    seurat <- performAnalysis(h5ad_file)
    seurat_object(seurat)
    
    # Update gene choices for feature plot
    gene_choices <- rownames(seurat@assays$RNA@counts)
    updateSelectInput(session, "gene_input", choices = gene_choices)
  })
  
  # Feature plot rendering
  output$feature_plot <- renderPlot({
    seurat <- seurat_object()
    if (!is.null(seurat)) {
      gene <- input$gene_input
      FeaturePlot(seurat, features = gene, label = TRUE,)
    }
  })
  
  # Function to fetch data based on specified gene input
  fetchData <- function(seurat_object, gene_input) {
    seurat <- seurat_object()
    if (!is.null(seurat)) {
      FetchData(object = seurat, vars = c("seurat_clusters", gene_input))
    }
  }
  

  # Export the data table as TSV file
  output$export_button <- downloadHandler(
    filename = function() {
      gene_input <- input$gene_input
      paste0("data_table_", gene_input, "_", format(Sys.time(), "%Y%m%d%H%M%S"), ".tsv")
    },
    content = function(file) {
      seurat <- seurat_object()
      gene_input <- input$gene_input
      if (!is.null(seurat) && !is.null(gene_input)) {
        data_table <- fetchData(seurat_object, gene_input)
        #colnames(data_table)[1] <- "Gene"  # Rename the first column to "Gene"
        write.table(data_table, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      }
    }
  )
  
  
  # End the Shiny session and shut down the server when the browser is closed
 
  # End the session 
  endSession <- function() {
    stopApp()  # This will end the Shiny session and shut down the server
  }
  
  # end Session browser window closed
  session$onSessionEnded(endSession)
}

  
  
# Run the Shiny app
shinyApp(ui = ui, server = server)
