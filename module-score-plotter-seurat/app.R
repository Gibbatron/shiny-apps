#Install packages
#install.packages('Seurat') #from CRAN
#install.packages("dplyr") #from CRAN
#install.packages("shiny") #from CRAN
#install.packages("scCustomize") #from CRAN
#install SeuratExtend from Github
#if (!requireNamespace("remotes", quietly = TRUE)) {
#    install.packages("remotes")
#}
#remotes::install_github("huayc09/SeuratExtend")

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
      h2("Module Score Feature Plotter & Violin Plotter For scRNAseq Seurat Data"),  # Main title
      h5("A tool that takes a processed Seurat object and list of genes as input, performs AddModuleScore function, and visualises the module score via a Feature plot (scCustomize) and Violin plot (SeuratExtend)."),#subtitle
      h4("Author: Dr. Alex Gibbs"),  # Author name
      h5("Last updated: 07/04/2025")  # Manually inserted date
    )
  ),

  #format side panel
  sidebarLayout(
    sidebarPanel(
      #file upload box
      fileInput("seurat_obj", "1) Upload Seurat object (.rds).", accept = ".rds"),
      helpText("Please upload your Seurat object in .rds format. This file should contain the processed data for analysis in the form of a Seurat object.  Please wait for 'upload complete' before inputting your gene list."),

      #text input box
      textInput("gene_title", "2) Enter a title for the gene list", placeholder = "Gene list title"),  # Title input
      helpText("This title will be used as the name of the module score for plotting"),

      #bigger text input box for gene list input
      textAreaInput("gene_list", "3) Paste your gene list here", rows = 5, placeholder = "Enter one gene per line"),
      helpText("Enter one gene per line. Please use gene symbols (i.e. Gata2 or GATA2). Gene symbols will be converted to the correct species format after selecting the appropriate species in the next step. Note: Users can upload human-specific lists and use the tool to convert them to mouse-specific and vice versa"),

      h4("4) Select your species"), #header line

      #check box
      checkboxInput("human", "Human", value = FALSE),  # Checkbox for Human
      helpText("This will convert your genes to human-specific symbols (all capitals). Example: Gata2 -> GATA2."),

      #check box
      checkboxInput("mouse", "Mouse", value = FALSE),  # Checkbox for Mouse
      helpText("This will convert your genes to mouse-specific symbols (capitalised and lowercase). Example: PIK3CA -> Pik3ca."),

      #dropdown box
      selectInput("vln_variable", "5) Select Variable for Violin Plot", choices = NULL),  # Dropdown for variable selection
      helpText("Choose a metadata variable for grouping cells in the violin plot."),

      #dropdown  box
      selectInput("split_variable", "6) Select Variable to Split Violin Plot", choices = NULL, selected = NULL),  # Dropdown for split variable
      helpText("Optional: Choose a second variable to split the violin plot by group."),

      #submit and download buttons
      actionButton("submit", "Submit"),
      downloadButton("download_genes", "Download Gene List"),# Download button
      helpText("If your gene list was converted and want to download it, click 'Download Gene List'")
    ),

    #format the main panel
    mainPanel(
      verbatimTextOutput("output"),

      plotOutput("feature_plot"), #output for feature plot
      radioButtons("feature_format", "Download FeaturePlot as:", choices = c("PNG", "SVG", "PDF"), inline = TRUE),
      downloadButton("downloadFeaturePlot", "Download FeaturePlot"),

      tags$hr(), # horizontal divider
      
      plotOutput("vln_plot"),  #output for VlnPlot
      radioButtons("vln_format", "Download Violin Plot as:", choices = c("PNG", "SVG", "PDF"), inline = TRUE),
      downloadButton("downloadVlnPlot", "Download Violin Plot")
    )
  )
)


# Define server logic
server <- function(input, output, session) {
  
  # Increase file size limit to 10GB
  options(shiny.maxRequestSize = 10 * 1024^3)
  
  # Reactive expression to store the Seurat object
  seurat_obj <- reactiveVal(NULL)
  
  # Reactive value to store the converted gene list for download
  converted_gene_list <- reactiveVal(NULL)
  
  # Reactive value to store the module score name
  module_score_name <- reactiveVal(NULL)
  
  # Reactive expressions for the plots
  feature_plot_reactive <- reactive({
    req(seurat_obj(), module_score_name())
    seurat_data <- seurat_obj()
    full_feature_name <- paste0(module_score_name(), "1")
    
    # Check if the module score column exists in the metadata
    if(!full_feature_name %in% colnames(seurat_data@meta.data)) {
      return(NULL)
    }
    
    FeaturePlot_scCustom(seurat_data, features = full_feature_name)
  })

  vln_plot_reactive <- reactive({
    req(seurat_obj(), input$vln_variable, module_score_name())
    seurat_data <- seurat_obj()
    full_feature_name <- paste0(module_score_name(), "1")
    
    # Check if the module score column exists in the metadata
    if(!full_feature_name %in% colnames(seurat_data@meta.data)) {
      return(NULL)
    }
    
    VlnPlot2(seurat_data,
             features = full_feature_name,
             assay = "RNA",
             group.by = input$vln_variable,
             split.by = if(nchar(input$split_variable) > 0) input$split_variable else NULL,
             stat.method = "wilcox.test")
  })

  # Update variable choices when the Seurat object is uploaded
  observeEvent(input$seurat_obj, {
    req(input$seurat_obj)
    seurat_data <- readRDS(input$seurat_obj$datapath)
    seurat_obj(seurat_data)  # Store the loaded Seurat object
    
    # Get metadata columns and update the select input choices
    metadata_cols <- colnames(seurat_data@meta.data)
    updateSelectInput(session, "vln_variable", choices = c("", metadata_cols)) #update the variable choices
    updateSelectInput(session, "split_variable", choices = c("", metadata_cols))  #update split variable choices
  })
  
  # Handle submit button click
  observeEvent(input$submit, {
    req(seurat_obj())  # Ensure the Seurat object is uploaded
    
    # Add progress bar
    withProgress(message = 'Processing..', {
      incProgress(0.1)
    
      tryCatch({ #used for error handling
        # Get the Seurat object
        seurat_data <- seurat_obj()
        
        # Get gene list from text area
        gene_list <- unlist(strsplit(input$gene_list, "\n"))
        
        # Ensure that the gene list is not empty
        if (length(gene_list) > 0) {
          # Clean the gene list: trim whitespace and remove empty strings
          gene_list <- trimws(gene_list)
          gene_list <- gene_list[nchar(gene_list) > 0]
          
          # Convert genes based on species selection
          if (input$human) {
            gene_list <- toupper(gene_list)  # Convert to all capitals for Human
          } else if (input$mouse) {
            gene_list <- tools::toTitleCase(tolower(gene_list))  # Capitalize first letter for Mouse
          }
          
          # Store the converted gene list for download
          converted_gene_list(gene_list)
          
          # Use the title provided by the user, defaulting to "ModuleScore" if empty
          module_title <- ifelse(nchar(input$gene_title) > 0, input$gene_title, "ModuleScore")
          
          # Store the module score name in reactive value
          module_score_name(module_title)
          
          # Perform the AddModuleScore
          seurat_data <- AddModuleScore(object = seurat_data, features = list(gene_list), name = module_title)
          
          # Update the Seurat object reactive value with the new object containing module scores
          seurat_obj(seurat_data)
          
          incProgress(0.4)  # Progress after calculating module score
          
          # Render output with the number of found genes
          output$output <- renderPrint({
            total_genes_input <- length(gene_list)
            genes_found <- sum(gene_list %in% rownames(seurat_data))
            list(
              total_genes_input = total_genes_input,
              genes_found = genes_found,
              module_score_column = paste0(module_title, "1")  # Display the exact column name
            )
          })
          
          incProgress(0.2)  # Final progress increment
        } else {
          output$output <- renderPrint("No genes provided.")
        }
      }, error = function(e) {
        # Handle errors and update output
        output$output <- renderPrint({ paste("Error: ", e$message) })
      })
    })
  })
  
  # Download handler for the converted gene list
  output$download_genes <- downloadHandler(
    filename = function() {
      paste("converted_gene_list_", Sys.Date(), ".txt", sep = "")
    },
    content = function(file) {
      req(converted_gene_list())  # Ensure the gene list is available
      writeLines(converted_gene_list(), file)  # Write the gene list to the file
    }
  )

  # Render the plots
  output$feature_plot <- renderPlot({
    req(feature_plot_reactive())
    feature_plot_reactive()
  })
  
  output$vln_plot <- renderPlot({
    req(vln_plot_reactive())
    vln_plot_reactive()
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

  # Download handler for violin plot
  output$downloadVlnPlot <- downloadHandler(
    filename = function() {
      paste0("Violin_Plot_", Sys.Date(), ".", tolower(input$vln_format))
    },
    content = function(file) {
      format <- input$vln_format
      if (format == "PNG") {
        png(file, width = 1200, height = 1000, res = 150)
        print(vln_plot_reactive())
        dev.off()
      } else if (format == "PDF") {
        pdf(file, width = 10, height = 8)
        print(vln_plot_reactive())
        dev.off()
      } else if (format == "SVG") {
        svg(file, width = 10, height = 8)
        print(vln_plot_reactive())
        dev.off()
      }
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)