[![DOI](https://zenodo.org/badge/671548322.svg)](https://doi.org/10.5281/zenodo.18201450) [![PICtR status badge](https://agshaas.r-universe.dev/PICtR/badges/version)](https://agshaas.r-universe.dev/PICtR)

# PICtR

This is the computational workflow to analyse physically interacting cells (PICs) in flow cytometry data using R.  

## Installation

First, please install [BPCells](https://github.com/bnprks/BPCells). Next, PICtR can be installed using:

```R
install.packages("PICtR", repos = c("https://agshaas.r-universe.dev", "https://cloud.r-project.org"))
```

or directly from GitHub using:

```R
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("agSHaas/PICtR")
```

## Workflow  
  
Please refer to the [Vignette](https://agshaas.github.io/PICtR/articles/PICtR.html) and the Method section of our manuscript for a comprehensive description of the analysis pipeline. 

## Demo Data

The package comprises a demo data set from LCMV infected mice. Data was analysed using high parametric spectral flow cytometry and can be loaded with `data(demo_lcmv)`. See `?demo_lcmv` for details. 


## Dependencies

The core of this framework relies on [Seurat version 5](https://github.com/satijalab/seurat) and [BPCells](https://github.com/bnprks/BPCells). Please make sure that R is operating in version 4.3 or higher.

Some packages are not necessarily required, but expand the functionality of PICtR. Consider installing:

- [flowMeans](https://www.bioconductor.org/packages/release/bioc/html/flowMeans.html)  
- [dbscan](https://github.com/mhahsler/dbscan)  
- [Spectre](https://github.com/ImmuneDynamics/Spectre)
- [ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)  
- [cytoMEM](https://www.bioconductor.org/packages/release/bioc/html/cytoMEM.html)  
- [autothresholdr](https://github.com/rorynolan/autothresholdr)

## Citation

If you use PICtR in your work, please cite

> Vonficht, Jopp-Saile, Yousefian, Flore et al. Ultra-high-scale cytometry-based cellular interaction mapping. _Nature Methods_ (2025). https://doi.org/10.1038/s41592-025-02744-w


## Version 

Current version 1.0.1
