Comparative Transcriptomic Analysis of Breast Cancer Using R

About the Project

I conducted this project to explore gene expression differences between breast tumour and normal breast tissue using publicly available transcriptomic data.

I used the GSE3193 dataset from the NCBI Gene Expression Omnibus (GEO) database and performed the analysis using R and RStudio.

Through this project, I explored differential gene expression, identified significant probes, and visualised gene expression patterns.

Dataset Information

* Database: NCBI GEO
* GEO Accession: GSE3193
* Platform: GPL180
* Sample type: Breast tumour and normal breast tissue
* Samples analysed: 42
* Normal samples: 3
* Tumour samples: 39

I selected samples using the CR1 reference channel for the primary analysis to maintain reference consistency.

Tools and Technologies Used

* R programming
* RStudio
* GEOquery
* limma
* ggplot2
* pheatmap
* dplyr

Work Performed

I completed the following steps during the analysis:

1. Retrieved the GSE3193 dataset from GEO using the GEOquery package.
2. Extracted the expression matrix and sample metadata.
3. Identified and classified the normal and tumour samples.
4. Examined missing values in the expression matrix.
5. Filtered probes based on the availability of expression measurements.
6. Checked the reference channels and selected CR1 samples for the primary comparison.
7. Created a design matrix and contrast using the limma package.
8. Performed differential expression analysis.
9. Applied Benjamini–Hochberg correction to adjust p-values.
10. Annotated significant probes with their corresponding gene symbols.
11. Generated volcano, heatmap and PCA plots.

Results

My differential expression analysis identified 871 significant probes with an adjusted p-value below 0.05.

Some of the top annotated gene symbols included:

* JUN: Associated with cell growth and regulation of gene activity.
* TFDP1: Involved in cell-cycle regulation.
* PLAU: Associated with tissue remodelling.
* PFKFB3: Involved in cellular glucose metabolism.
* FABP4: Associated with fatty-acid transport and lipid metabolism.

These genes provide interesting starting points for further investigation into breast cancer-associated molecular changes.

Visualisations

I generated three visualisations to understand the results:

1. Volcano Plot: To visualise differential expression and statistical significance.
2. Heatmap: To examine expression patterns across the top 30 significant probes.
3. PCA Plot: To explore the overall variation and grouping of samples.

Limitations

* The primary comparison included only three normal samples and 39 tumour samples.
* The findings are exploratory and require validation using independent datasets.
* The dataset contains two-colour microarray log2 ratios rather than direct absolute expression measurements.
* Significant probes may map to the same gene.

Conclusion

Through this project, I gained practical experience in transcriptomic data
 analysis using R.

I learned how to retrieve public gene expression datasets, preprocess expression data, perform differential expression analysis, annotate genes and generate visualisations.

This project strengthened my understanding of computational approaches used in bioinformatics and cancer research.

## Author

**Diana Goldy**

B.Tech Bioinformatics Student  
Amity University Rajasthan