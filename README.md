# PICtR

This is the computational workflow to analyse physically interacting cells (PICs) in flow cytometry data using R.  

## Installation
For installation please use the following instructions (installation < 1 minute):
```R
remotes::install_github("agSHaas/PICtR") # Currently the repository is private. It will be made public upon manuscript acceptance
remotes::install_local("/path/to/PICtR_0.1.0.tar.gz", dependencies = T) # After downloading from https://zenodo.org/records/10637096?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjI0ODdiYzc2LWJmYzUtNGE5NS1iMzYzLWY2ZmU5Mzg4MzVmMCIsImRhdGEiOnt9LCJyYW5kb20iOiJlMDdiY2I0MDQzMGE3ZWE3NThiNzc2NGIwMGMzOWQ3MiJ9.rGbkSQ0fvx7ElMbD9HjXLtzen5qcfPonIpTXR1zdrdFQc5yBw5cthsn-yYiHXSWWQbhipHxlm7m8eLgDRXfsGg 
``` 
BPCells can be installed from GitHub using 
```R
remotes::install_github("bnprks/BPCells")
```
and cytoMEM and ComplexHeatmap can be installed from Bioconductor using
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("cytoMEM", "ComplexHeatmap"))
```
## Workflow  
  
1. In brief, the `sketch_wrapper` function is used to
   a) classify cellular events according to their FSC-area to FSC-height ratio
   b) select a representative subset using Seurat's sketch functionality, and
   c) run the standard Seurat workflow.
   It expects a data frame with the dimensionality cells x flow cytometry parameters.  
2. After selecting an appropriate clustering solution in the subset, the results can be predicted for the full dataset with `predict_data`, using Linear Discriminant Analysis as implemented in the MASS package. The overall cellular landscape can be annotated as usual, and clusters representing interacting cells can be selected using `ratio_cluster_plot` and `select_dbt`.
3. Interacting cell clusters can be annotated using cell type exclusive markers. 
4. `plot_functions` contains a selection of helpful wrapper functions for plots.

Please refer to the Method section of the manuscript for a comprehensive description of the analysis pipeline. 

## Demo Data

The package comprises a demo data set of 10% randomly sampled cells from the LCMV infected mice data (subset of data shown in Figure 5 and 6). Data was analysed using high parametric spectral flow cytometry (n=26 markers). Data can be loaded with `data(demo_lcmv)`.

## How to use PICtR

Load PICtR package and demo data:
```R
library(PICtR)
data(demo_lcmv)
```

Run sketch wrapper to identify cells with a high FSC.A/FSC.H ratio, potentially representing interacting cells, and run the standard Seurat workflow. To keep computational time manageable, data will be sketched (more details see manuscript or [Seurat](https://satijalab.org/seurat/articles/seurat5_sketch_analysis)).
```R
demo_obj <- sketch_wrapper(channel = demo_lcmv,
                       meta_data = demo_lcmv,
                       FSC.A = "FSC.A", 
                       FSC.H = "FSC.H", 
                       working_dir=getwd(),
                       resolution = c(0.8, 1),
                       n_sketch_cells=5000,
                       obj_name="demo_data",
                       verbose=T)
```

Generate different plots to evaluate the analysis:
```R
plot_list <- wrapper_for_plots(demo_obj, 
                       feature_plot=TRUE,
                       cluster_plot = TRUE,
                       meta_list = list("ratio_anno"))
```  

Project the remaining cells to obtain cluster information across all cells using either a Seurat object or data frame as query data:
```R
# Seurat object
demo_obj <- predict_data(obj = demo_obj, 
                      data_query = demo_obj, 
                      ref_clusters = "sketch_snn_res.0.8",
                      pred_name = "clusters_predicted_obj",
                      assay_ref = "sketch", 
                      assay_query = "FACS")
                      
# Data frame
demo_obj <- predict_data(obj = demo_obj, 
                      data_query = demo_obj@meta.data[is.na(demo_obj@meta.data$sketch_snn_res.0.8),], 
                      ref_clusters = "sketch_snn_res.0.8",
                      pred_name = "clusters_predicted_df",
                      assay_ref = "sketch", 
                      assay_query = "FACS")
```

By using the cluster label information of sketched and projected cells and combining it with the distribution of ratio high and ratio low cells, interacting cell clusters can be determined: 

```R
ratio_cluster_plot(demo_obj, clusters = "clusters_predicted_obj")
```
Next, the aim is to select clusters with a high number of interacting cells:
```R
demo_obj <- select_dbt(demo_obj, clusters = "clusters_predicted_obj", quantile = 0.95)
```
The selected cluster numbers are stored in `demo_obj@misc$doublet_clusters_q0.95`. The cells within these clusters represent physically interacting cells and can be used for downstream analysis.  
The demo Seurat object is included in the package data to exemplify the expected output of the pipeline. The demo pipeline has a runtime of <30 seconds on a standard desktop computer.

## Dependencies

The core of this framework relies on [Seurat version 5](https://github.com/satijalab/seurat), [BPCells](https://github.com/bnprks/BPCells), [MASS](https://cran.r-project.org/web/packages/MASS/index.html), and the [tidyverse](https://www.tidyverse.org/). The full list of dependencies are listed in the DESCRIPTION file. The package has been tested on the indicated dependency versions. Please make sure that R is operating in version 4.3 or higher. 

## Version 

Current version 0.1.0-alpha
