This app takes an RDS object as input and allows the user to visualise gene expression via Feature plots and Dot plots
I created this app to allow a user to quickly take a saved Seurat object (.rds) and a gene of interest visualise via scCustomise Feature plot and dot plot functions.

**Prerequisites:**
- The Seurat object must contain processed data up to UMAP/tSNE stage.
- You will need to install the following packages. Simply paste the following into your RStudio console:

```
install.packages('Seurat') #from CRAN
install.packages("dplyr") #from CRAN
install.packages("shiny") #from CRAN
install.packages("scCustomize") #from CRAN

#install SeuratExtend from Github
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

remotes::install_github("huayc09/SeuratExtend")
```

If Mac users want to save the plots as SVG, you will need to ensure you have XQuartz installed on your system.

Assuming you have brew installed on the command line, open a Terminal session:

```
brew install --cask xquartz
```

Then close the Terminal session and restart your RStudio if you have it open.

**Launching the app**

To launch the app, paste the following into the RStudio console:

```
runGitHub(repo = 'Gibbatron/shiny-apps', subdir = 'scRNAseq-gene-visualiser')
```

**Notes:**
- Please allow some time for the app to load your Seurat object, this can take anywhere from 5 seconds to 5 minutes, depending on the size of the object and the speed of your PC.
- Please allow time for the app to parse through your gene list and perform any relevant conversions. This can also take some time depending on the size of your list.
- Once you have selected your Seurat object to upload, please wait for the loading bar to say 'upload completed' before moving onto section 2.
- After the successful loading of the Seurat object, please allow some time for the app to process the metadata and load it into the dropdown boxes in sections 5 and 6.
- If you have used a human-formatted gene list and need to convert it to mouse-specific format, you can download the converted list via the 'Download Gene List' button after clicking 'Submit'
- Depending on the sixe of your Seurat object and gene list, please allow some time for the feature and violin plots to load after the progress bar reaches 100%
- To download the figures, simply right click on the figure and save.
- The AddModuleScore() function by default adds a '1' to the end of the module name (e.g. module-list -> module-list1). Please bear this in mind if you have a modulename that ends with a '1', after processing, it will now end in '11'.
