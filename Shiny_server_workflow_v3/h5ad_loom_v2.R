library(Seurat)
library(BiocParallel)
library(reticulate)
library(Matrix)
library(shiny)
library(dplyr)
library(ggplot2)
library(loomR)

# Set the maximum file upload size to 3 GB
options(shiny.maxRequestSize = 3e+9)

ui <- fluidPage(
  titlePanel("Single Cell H5ad and loom file analysis"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("fileType", "Choose file type:",
                   choices = c("H5AD", "Loom")),
      conditionalPanel(
        condition = "input.fileType == 'H5AD'",
        fileInput("h5adFile", "Choose an H5AD file", accept = c(".h5ad"))
      ),
      conditionalPanel(
        condition = "input.fileType == 'Loom'",
        fileInput("loomFile", "Choose a Loom file", accept = c(".loom"))
      ),
      actionButton("run_button", "Perform Analysis"),
      selectInput("cluster_input", "Select Cluster", choices = NULL),
      numericInput("fold_change_threshold", "Fold Change Threshold:", value = 0, step = 0.1),
      downloadButton("exportClusters", "Export Clusters"),
      downloadButton("export_top_genes_button", "Export Chosen Cluster Genes")
      
    ),
    mainPanel(
      plotOutput("umapPlot")
    )
  )
)

server <- function(input, output, session) {
  # Initialize cluster choices
  cluster_choices <- reactiveVal(NULL)
  
  # Function to perform H5AD file analysis
  performH5ADAnalysis <- function(h5ad_file) {
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
    cl_markers <- FindAllMarkers(seurat, only.pos = TRUE)
    top_gene_data <- cl_markers 
    
    cluster_choices(levels(seurat$louvain))
    
    return(list(seurat = seurat, top_gene_data = top_gene_data))
  }
  
  # Function to perform Loom file analysis
  performLoomAnalysis <- function(loom_file) {
    l6.immune <- Connect(filename = loom_file, mode = "r")
    vel <- as.Seurat(l6.immune)
    vel <- NormalizeData(vel, verbose = FALSE)
    vel <- ScaleData(vel, verbose = FALSE)
    vel <- RunPCA(vel, features = rownames(vel), verbose = FALSE)
    vel <- FindNeighbors(vel, dims = 1:20, verbose = FALSE)
    vel <- FindClusters(vel, verbose = FALSE)
    vel <- RunUMAP(vel, dims = 1:20, verbose = FALSE)
    cl_markers <- FindAllMarkers(vel, only.pos = TRUE)
    top_gene_data <- cl_markers
    
    cluster_choices(levels(vel$seurat_clusters))
    
    return(list(vel_reduction = vel, top_gene_data = top_gene_data))
  }
  
  # Reactive function to perform the analysis based on the selected file type
  analysisResult <- reactive({
    if (input$fileType == "H5AD") {
      req(input$h5adFile)
      performH5ADAnalysis(input$h5adFile$datapath)
    } else {
      req(input$loomFile)
      performLoomAnalysis(input$loomFile$datapath)
    }
  })
  
  # Perform analysis when the "Perform Analysis" button is clicked
  observeEvent(input$run_button, {
    analysisResult()
  })
  
  # Update cluster choices when analysis result is available
  observe({
    result <- analysisResult()
    if (!is.null(result$seurat)) {
      cluster_choices(levels(result$seurat$seurat_clusters))
    } else if (!is.null(result$vel_reduction)) {
      cluster_choices(levels(result$vel_reduction$seurat_clusters))
    }
  })
  
  # Generate UMAP plot based on the analysis result
  output$umapPlot <- renderPlot({
    result <- analysisResult()
    seurat <- result$seurat
    vel <- result$vel_reduction
    
    if (!is.null(seurat)) {
      DimPlot(seurat, reduction = "umap",group.by = "louvain")
    } else if (!is.null(vel)) {
      DimPlot(vel, reduction = "umap")
    }
  })
  
  # Update the cluster selection input based on available cluster choices
  observeEvent(cluster_choices(), {
    updateSelectInput(session, "cluster_input", choices = cluster_choices())
  })
  
  # Export cluster information as a CSV file
  output$exportClusters <- downloadHandler(
    filename = function() {
      fileType <- input$fileType
      if (fileType == "H5AD") {
        paste("H5ad_cluster_", input$cluster_input, ".tsv", sep = "")
      } else if (fileType == "Loom") {
        paste("loom_clusters_", input$cluster_input, ".tsv", sep = "")
      } else {
        NULL
      }
    },
    content = function(file) {
      fileType <- input$fileType
      result <- analysisResult()
      if (fileType == "H5AD") {
        seurat <- result$seurat
        if (!is.null(seurat)) {
          selected_cluster <- input$cluster_input
          selected_cells <- seurat$louvain == selected_cluster 
          selected_data <- seurat@meta.data[selected_cells, ]
          write.table(selected_data, file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
      } else if (fileType == "Loom") {
        vel <- result$vel_reduction
        if (!is.null(vel)) {
          selected_cluster <- input$cluster_input
          selected_cells <- vel$seurat_clusters == selected_cluster
          selected_data <- vel@meta.data[selected_cells, c("orig.ident", "nCount_RNA", "nFeature_RNA", "TotalUMIs", "RNA_snn_res.0.8", "seurat_clusters")]
          write.table(selected_data, file, sep = "\t", quote = FALSE, row.names = FALSE)
        }
      }
    }
  )
  
  
  
  # Export top marker genes as TSV file
  # Export top marker genes as TSV file
  output$export_top_genes_button <- downloadHandler(
    filename = function() {
      paste("cluster_", input$cluster_input, "_genes.tsv", sep = "")
    },
    content = function(file) {
      if (input$fileType == "H5AD") {
        seurat <- analysisResult()$seurat
        top_gene_data <- analysisResult()$top_gene_data
        
        if (!is.null(seurat) && !is.null(top_gene_data)) {
          selected_cluster <- input$cluster_input
          fold_change_threshold <- input$fold_change_threshold
          top_genes <- top_gene_data %>%
            filter(cluster == selected_cluster, avg_log2FC >= fold_change_threshold) %>%
            select(cluster, gene, p_val_adj, avg_log2FC)
          write.table(top_genes, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        }
      } else {
        vel <- analysisResult()$vel_reduction
        top_gene_data <- analysisResult()$top_gene_data
        
        if (!is.null(vel) && !is.null(top_gene_data)) {
          selected_cluster <- input$cluster_input
          fold_change_threshold <- input$fold_change_threshold
          top_genes <- top_gene_data %>%
            filter(cluster == selected_cluster, avg_log2FC >= fold_change_threshold) %>%
            select(cluster, gene, p_val_adj, avg_log2FC)
          write.table(top_genes, file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        }
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

shinyApp(ui = ui, server = server)
