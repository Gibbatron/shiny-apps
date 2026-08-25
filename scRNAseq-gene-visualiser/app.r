#Load Libraries
library(Seurat)
library(dplyr)
library(shiny)
library(scCustomize)  # For FeaturePlot_scCustom
library(SeuratExtend)  # For VlnPlot2

# Define UI for the application
ui <- fluidPage(

  #format title panel
  titlePanel(
    div(
      h2("Gene Expression Plotter For scRNAseq Seurat Data"),  # Main title
      h5("A tool that takes a processed Seurat object and allows the user to plot the expression of specific genes via a Feature plot and Dot Plot (scCustomize)"),#subtitle
      h4("Author: Dr. Alex Gibbs"),  # Author name
      h5("Last updated: 25/08/2025")  # Manually inserted date
    )
  ),

  #format side panel
  sidebarLayout(
    sidebarPanel(
      tags$style(HTML("
      .section-divider {
        border: none;
        height: 3px;
        background-color: #999;
        margin: 15px 0;
      }
    ")),

      #file upload box
      fileInput("seurat_obj", "1) Upload Seurat object (.rds).", accept = ".rds"),
      helpText("Please upload your Seurat object in .rds format. This file should contain the processed data for analysis in the form of a Seurat object.  Please wait for 'upload complete' before inputting your gene list.",
      style = "margin-top: -15px;"),

  hr(class = "section-divider"),

      #text input box
      textInput("gene_name", "2) Enter your gene for plotting", placeholder = "Gene name"),  # gene input
      helpText("This gene will be visualised via Feature plot and Dot plot. Please ensure the gene name is spelled correctly."),

  hr(class = "section-divider"),

      h5("4) Select your species", style = "font-weight: bold;"), #header line

      #check box
      checkboxInput("human", "Human", value = FALSE),  # Checkbox for Human
      helpText("This will convert your genes to human-specific symbols (all capitals). Example: Gata2 -> GATA2.",
      style = "margin-top: -15px;"),

      #check box
      checkboxInput("mouse", "Mouse", value = FALSE),  # Checkbox for Mouse
      helpText("This will convert your genes to mouse-specific symbols (capitalised and lowercase). Example: PIK3CA -> Pik3ca.",
      style = "margin-top: -15px;"),

  hr(class = "section-divider"),

      #submit button
      actionButton("submit", "Submit"),
    ),

    #format the main panel
    mainPanel(
      verbatimTextOutput("output"),

      plotOutput("feature_plot"), #output for feature plot
      radioButtons("feature_format", "Download Feature Plot as:", choices = c("PNG", "SVG", "PDF"), inline = TRUE),
      downloadButton("downloadFeaturePlot", "Download Feature Plot"),

      tags$hr(), # horizontal divider
      
      plotOutput("dot_plot"),  #output for Dot Plot
      radioButtons("dot_format", "Download Dot Plot as:", choices = c("PNG", "SVG", "PDF"), inline = TRUE),
      downloadButton("downloadDotPlot", "Download Dot Plot")
    )
  )
)


# Define server logic
server <- function(input, output, session) {
  
  # Increase file size limit to 20GB
  options(shiny.maxRequestSize = 20 * 1024^3)
  
  # Reactive expression to store the Seurat object
  seurat_obj <- reactiveVal(NULL)
  
  # Reactive value to store the converted gene for download
  converted_gene <- reactiveVal(NULL)
  
  # Reactive expressions for the plots
  feature_plot_reactive <- reactive({
    req(seurat_obj(), converted_gene())
    seurat_data <- seurat_obj()
    
    # Check if the gene exists in the object
    if (!converted_gene() %in% rownames(seurat_data)) {
      return(NULL)
    }
    
    FeaturePlot_scCustom(seurat_data, features = converted_gene())
  })

dot_plot_reactive <- reactive({
  req(seurat_obj(), converted_gene())
  seurat_data <- seurat_obj()

  if (!converted_gene() %in% rownames(seurat_data)) {
    return(NULL)
  }

  DotPlot_scCustom(
    seurat_data,
    features = converted_gene(),
    assay = "RNA"
  )
})

  # Update variable choices when the Seurat object is uploaded
observeEvent(input$seurat_obj, {
  req(input$seurat_obj)
  seurat_data <- readRDS(input$seurat_obj$datapath)
  seurat_obj(seurat_data)
  
  metadata_cols <- colnames(seurat_data@meta.data)
  updateSelectInput(session, "split_variable", choices = c("", metadata_cols))
})
  
  # Handle submit button click
  observeEvent(input$submit, {
    req(seurat_obj())  # Ensure the Seurat object is uploaded
    
    withProgress(message = 'Processing..', {
      incProgress(0.1)
    
      tryCatch({
        gene_name <- input$gene_name

        if (nchar(gene_name) == 0) {
          output$output <- renderPrint("No genes provided.")
        } else {
          # Convert gene based on species selection
          if (input$human) {
            gene_name <- toupper(gene_name)  # Human: all capitals
          } else if (input$mouse) {
            gene_name <- tools::toTitleCase(tolower(gene_name))  # Mouse: Capitalized
          }
          
          # Store the converted gene for use in plots/download
          converted_gene(gene_name)
          
          incProgress(0.2)
        }
      }, error = function(e) {
        output$output <- renderPrint({ paste("Error: ", e$message) })
      })
    })
  })

  # Render the plots
  output$feature_plot <- renderPlot({
    req(feature_plot_reactive())
    feature_plot_reactive()
  })
  
  output$dot_plot <- renderPlot({
    req(dot_plot_reactive())
    dot_plot_reactive()
  })

  # Download handler for feature plot
  output$downloadFeaturePlot <- downloadHandler(
    filename = function() {
      paste0("Feature_Plot_", Sys.Date(), ".", tolower(input$feature_format))
    },
    content = function(file) {
      format <- input$feature_format
      if (format == "PNG") {
        png(file, width = 1200, height = 1000, res = 150)
        print(feature_plot_reactive())
        dev.off()
      } else if (format == "PDF") {
        pdf(file, width = 10, height = 8)
        print(feature_plot_reactive())
        dev.off()
      } else if (format == "SVG") {
        svg(file, width = 10, height = 8)
        print(feature_plot_reactive())
        dev.off()
      }
    }
  )

  # Download handler for dot plot
  output$downloadDotPlot <- downloadHandler(
    filename = function() {
      paste0("Dot_Plot_", Sys.Date(), ".", tolower(input$dot_format))
    },
    content = function(file) {
      format <- input$dot_format
      if (format == "PNG") {
        png(file, width = 1200, height = 1000, res = 150)
        print(dot_plot_reactive())
        dev.off()
      } else if (format == "PDF") {
        pdf(file, width = 10, height = 8)
        print(dot_plot_reactive())
        dev.off()
      } else if (format == "SVG") {
        svg(file, width = 10, height = 8)
        print(dot_plot_reactive())
        dev.off()
      }
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)
