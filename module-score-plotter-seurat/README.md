I created this app to allow a user to quickly take a saved Seurat object (.rds) and a list of genes and perform Seurats AddModuleScore function, then plot the results in a featureplot and violin plot.

The user wanted to create a module score for various gene lists and assess the module scores expression across their dataset. They also wanted to see the module scores expression between their two sample groups via a Violin Plot.

The scCustomize packages version of the Featureplot (FeaturePlot_scCustom) produces a more informative plot in the means of 'greying' out cells that do not have an expression value for the plotted feature.
The SeuratExtend packages version if the Violin Plot (VlnPlot2) produces nicer looking plots, with optional statistical significance testing between groups.


Prerequisites:
- The Seurat object must contain processed data up to UMAP/tSNE stage.


Notes:
- Please allow some time for the app to load your Seurat object, this can take anywhere from 5 seconds to 5 minutes, depending on the size of the object and the speed of your PC.
- Please allow time for the app to parse through your gene list and perform any relevant conversions. This can also take some time depending on the size of your list.
- Once you have selected your Seurat object to upload, please wait for the loading bar to say 'upload completed' before moving onto section 2.
- After the successful loading of the Seurat object, please allow some time for the app to process the metadata and load it into the dropdown boxes in sections 5 and 6.
- If you have used a human-formatted gene list and need to convert it to mouse-specific format, you can download the converted list via the 'Download Gene List' button after clicking 'Submit'
- Depending on the sixe of your Seurat object and gene list, please allow some time for the feature and violin plots to load after the progress bar reaches 100%
- To download the figures, simply right click on the figure and save.
- The AddModuleScore() function by default adds a '1' to the end of the module name (e.g. module-list -> module-list1). Please bear this in mind if you have a modulename that ends with a '1', after processing, it will now end in '11'.
