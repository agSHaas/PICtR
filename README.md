# PICtR

This is the computational workflow to analyse physically interacting cells (PICs) in flow cytometry data using R.  
  
## Workflow  
  
1. In brief, the `sketch_wrapper` function is used to classify cellular events according to their FSC-area to FSC-height ratio, select a representative subset using Seurat's sketch functionality, and run the standard Seurat workflow. It expects a data frame with the dimensionality cells x flow cytometry parameters.  
2. After selecting an appropriate clustering solution in the subset, the results can be predicted for the full dataset with `predict_data`, using Linear Discriminant Analysis as implemented in the MASS package. The overall cellular landscape can be annotated as usual, and clusters representing interacting cells can be selected using `ratio_cluster_plot` and `select_dbt`.
3. Interacting cell clusters can be annotated using cell type exclusive markers. 
4. `plot_functions` contains a selection of helpful wrapper functions for plots.

Please refer to the Method section for a comprehensive description of the analysis pipeline. 
  
## Dependencies 
The core of this framework relies on [Seurat](https://github.com/satijalab/seurat), [BPCells](https://github.com/bnprks/BPCells), [MASS](https://cran.r-project.org/web/packages/MASS/index.html), and the [tidyverse](https://www.tidyverse.org/). Please make sure that Seurat version 4.9.9.9058 or higher is installed on your system.
