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
      h5("Last updated: 04/10/2024")  # Manually inserted date
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
      plotOutput("vln_plot")  #output for VlnPlot2
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
  
  # Update variable choices when the Seurat object is uploaded
  observeEvent(input$seurat_obj, {
    req(input$seurat_obj)
    seurat_data <- readRDS(input$seurat_obj$datapath)
    seurat_obj(seurat_data)  # Store the loaded Seurat object
    
    # Get metadata columns and update the select input choices
    metadata_cols <- colnames(seurat_data@meta.data)
    updateSelectInput(session, "vln_variable", choices = metadata_cols) #update the variable choices
    updateSelectInput(session, "split_variable", choices = metadata_cols)  #update split variable choices
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
          module_score_name <- ifelse(nchar(input$gene_title) > 0, input$gene_title, "ModuleScore")
          
          # Perform the FindModuleScore
          seurat_data <- AddModuleScore(object = seurat_data, features = list(gene_list), name = module_score_name)
          incProgress(0.4)  # Progress after calculating module score
          
          # Render output with the number of found genes
          output$output <- renderPrint({
            total_genes_input <- length(gene_list)
            genes_found <- sum(gene_list %in% rownames(seurat_data))
            list(
              total_genes_input = total_genes_input,
              genes_found = genes_found,
              module_score_column = paste0(module_score_name, "1")  # Update to reflect user title
            )
          })
          
          # Plot the module score using FeaturePlot_scCustom
          output$feature_plot <- renderPlot({
            FeaturePlot_scCustom(seurat_data, features = paste0(module_score_name, "1"))  # Use user-defined name
          })
          
          # Plot violin plot using VlnPlot2 with the selected variable
          output$vln_plot <- renderPlot({
            req(input$vln_variable)  # Ensure a variable is selected
            # Split by additional variable if provided
            VlnPlot2(seurat_data, features = paste0(module_score_name, "1"), 
                     assay = "RNA", 
                     group.by = input$vln_variable, 
                     split.by = input$split_variable,  # Additional splitting variable
                     stat.method = "wilcox.test")  # Statistical testing
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
}

# Run the app
shinyApp(ui = ui, server = server)
