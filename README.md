Repository "first_repo" for AES and ACH's 20.440 term project. 
Project title: Identification and characterization of stem cell niche supporting endothelial cells in 
murine tissues through single cell RNA sequencing
Created 3/28/2020 by AES, last modified 3/30/2020 by AES


Overview:
This repository contains code and data used to generate figures showing UMAP clustering and cluster gene expression. The workflow is as follows, broadly following Sureat "Guided tutorial --- 2,700 PBMCs" (https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html) :
- read in data csv file, adjust data row and column names (https://github.com/satijalab/seurat/issues/1710)
- create SeuratObject
- normalize data using Sureat NormalizeData() (https://satijalab.org/seurat/v3.0/mca.html)
- FindVariableFeatures() (
- Regress out mitochondrial expression using PercetageFeatureSet and ScaleData() ()
- do prinicpal component analysis Sureat  RunPCA() ()
- Create Sureat ElbowPlot() to visualize relevant PCAs ()
- Create dimension heat map with Sureat DimHeatMap() 
- cluster using Sureat FindNeighbors(), Sureat FindClusters()
- visualize clusters using Sureat RunUMAP() and plot, save as pdf to figures folder
- plot gene expression of select genes (Wnt2, Rspo3, Vwf, Jag1) in clusters with Sureat FeaturePlot()
For more information about using Suerat please visit (https://satijalab.org/seurat/) or 
Butler et al., Nature Biotechnology 2018 "Integrating single-cell transcriptomic data across different conditions, technologies, and species".

Data:
Our data set is from Kalucka et al's 2020 Cell publication "Single-Cell Transcriptome Atlas
of Murine Endothelial Cells" DOI: https://doi.org/10.1016/j.cell.2020.01.015
Data is available at https://www.vibcancer.be/software-tools/ec-atlas
For this first analyses only the liver tissue data was downloaded, and is available in this repository as liver_data.csv or through the above website. 
This data set was generated through ensymatic digestion and/or perfusion of murine tissues and FACS isolation of endothelial cells from the resulting cell suspension. 6-8 individual mice (of same genetic background) were pooled for each tissue data set. Libraries were prepared using Chromium Single Cell 3' Reagent Kits, libraries were sequenced on the Illumina HiSeq4000 and then mapped to the mouse genome using CellRanger. In silico, non EC cells were removed. 

Folder Structure:
This repository contains 3 subfolders: Data, which contains the raw data used for this analysis, Code, which contains the scripts to produce the figures, and Figures. 
The Code folder contains two scripts - Liver_Script.R, which was used to produce the figures in the Figures folder, and Small_intestine_Script, another script for use on the small intetsine data set created by ACH. 
The Figures folder contains liver_clustering.pdf, the UMAP clustering of liver endothelial cells, and cluster_expression.pdf, a UMAP clustering with Wnt2, Rspo3, Vwf, and Jag1 gene expressions overlaid. 

Installation:
We reccomend running the .R scripts through Rstudio. Input and output files as well as present working directories should be modified for the users's file structure. 

To run these scripts the user should install R v3.4 or greater (https://repo.miserver.it.umich.edu/cran/), Rstudio (https://rstudio.com/products/rstudio/), and Suerat (instructions here: https://satijalab.org/seurat/install.html , can be executed from RStudio command prompt)








