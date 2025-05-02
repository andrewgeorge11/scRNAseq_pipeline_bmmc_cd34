# Single-Cell RNA-seq Analysis of BMMC and CD34 Samples

- Project Description
This project applies a complete single-cell RNA-seq analysis pipeline to compare bone marrow mononuclear cells (BMMC) and CD34+ hematopoietic stem/progenitor cells using Seurat and other R-based tools. It includes preprocessing, quality control, batch correction, cell type annotation, differential expression analysis, and pathway enrichment.

_ Tools Used
- R (Seurat, DoubletFinder, SingleR, celldex, enrichR)
- ggplot2
- clusterProfiler
- Enrichr
- UMAP, PCA

_ Repository Structure
- script : Complete R pipeline script
- README.md : Project overview
- report: Final summarized report in PDF

_ Key Features
- Doublet removal with DoubletFinder
- Cell type annotation using both SingleR and manual markers
- Batch integration and correction
- Volcano plots and violin plots for DE analysis
- GO enrichment analysis for differentially expressed genes

_ Sample Comparisons
- BMMC vs CD34
- B Cells vs T Cells
- T Cells vs Monocytes
- Monocytes in BMMC vs CD34

