library(Seurat)
library(BiocParallel)
library(reticulate)
library(Matrix)
library(shiny)
library(dplyr)
library(ggplot2)

# Set the maximum file upload size to 3 GB
options(shiny.maxRequestSize = 3e+9)

# Define the UI
ui <- fluidPage(
  titlePanel("Single-cell RNA-seq Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("h5ad_file", "Upload H5AD File", accept = ".h5ad"),
      actionButton("run_button", "Run Analysis"),
      selectInput("gene_input", "Select Gene for Feature Plot:", choices = NULL),
      selectInput("cluster_input", "Select Cluster for Export:", choices = NULL),
      downloadButton("export_cluster_button", "Export Cluster"),
      downloadButton("export_top_genes_button", "Export Chosen Cluster Genes")
    ),
    mainPanel(
      plotOutput("feature_plot"),
      plotOutput("cluster_plot")
    )
  )
)

# Define the server
server <- function(input, output, session) {
  
  # Function to perform the analysis
  performAnalysis <- function(h5ad_file) {
    sc <- import("scanpy")
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
    top_gene_data <- cl_markers %>% group_by(cluster) %>% top_n(n = 40, wt = avg_log2FC)
    
    return(list(seurat = seurat, top_gene_data = top_gene_data))
  }
  
  # Reactive values to store the Seurat object and top marker genes
  seurat_object <- reactiveVal(NULL)
  top_gene <- reactiveVal(NULL)
  
  # Perform analysis when the Run Analysis button is clicked
  observeEvent(input$run_button, {
    req(input$h5ad_file)
    h5ad_file <- input$h5ad_file$datapath
    analysis_result <- performAnalysis(h5ad_file)
    seurat_object(analysis_result$seurat)
    top_gene(analysis_result$top_gene_data)
    
    # Update gene choices for feature plot
    gene_choices <- rownames(seurat_object()@assays$RNA@counts)
    updateSelectInput(session, "gene_input", choices = gene_choices)
    
    # Update cluster choices for export
    cluster_choices <- unique(seurat_object()$seurat_clusters)
    updateSelectInput(session, "cluster_input", choices = cluster_choices)
  })
  
  # Render the feature plot
  output$feature_plot <- renderPlot({
    seurat <- seurat_object()
    top_gene_data <- top_gene()
    
    if (!is.null(seurat) && !is.null(top_gene_data)) {
      gene <- input$gene_input
      
      # Filter top gene data for the selected gene
      gene_data <- top_gene_data %>% filter(gene == input$gene_input)
      
      # Generate feature plot
      FeaturePlot(seurat, features = gene, label = TRUE) +
        geom_text(data = gene_data, aes(label = gene))
    }
  })
  
  # Render the cluster plot
  output$cluster_plot <- renderPlot({
    seurat <- seurat_object()
    
    if (!is.null(seurat)) {
      DimPlot(seurat, group.by = "seurat_clusters")
    }
  })
  
  # Export selected cluster as TSV file
  output$export_cluster_button <- downloadHandler(
    filename = function() {
      paste("cluster_", input$cluster_input, ".tsv", sep = "")
    },
    content = function(file) {
      seurat <- seurat_object()
      if (!is.null(seurat)) {
        selected_cluster <- input$cluster_input
        selected_cells <- seurat$seurat_clusters == selected_cluster 
        selected_data <- seurat@meta.data[selected_cells, ]
        write.table(selected_data, file, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  )
  
  # Export top marker genes as TSV file
  # Export selected cluster with top genes as TSV file
  # Export top marker genes as TSV file
  output$export_top_genes_button <- downloadHandler(
    filename = function() {
      paste("cluster_", input$cluster_input, "_genes.tsv", sep = "")
    },
    content = function(file) {
      seurat <- seurat_object()
      top_gene_data <- top_gene()
      
      if (!is.null(seurat) && !is.null(top_gene_data)) {
        selected_cluster <- input$cluster_input
        top_genes <- top_gene_data %>% filter(cluster == selected_cluster) %>% select(cluster,gene,p_val_adj,avg_log2FC)
        write.table(top_genes, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      }
    }
  )
  
  
  
  
  
  # Function to end the session
  endSession <- function() {
    stopApp()  # This will end the Shiny session and shut down the server
  }
  
  # Call the endSession function when the session ends (browser window closed)
  session$onSessionEnded(endSession)
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
